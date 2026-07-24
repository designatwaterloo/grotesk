"""
Waterloo Grotesk variable font builder (phase 2).

The statics pipeline (waterloo-statics.py) uses skia stroke+union per weight,
which produces incompatible point structures across weights -- fine for
statics, fatal for interpolation. This script builds a variable font whose
masters are ALL point-compatible:

1. For every named weight x optical size, generate a lighter-SKELETON
   instance UFO from the (already waterloo-transformed) Inter designspace,
   at the same remapped location the statics use.
2. Bleed each master analytically, preserving point structure exactly:
     - offset every point along its outward normal by R (miter-capped),
       thickening strokes like ink;
     - fillet every corner with a small cubic arc of radius R -- corner
       DECISIONS come from the default master and are applied identically
       everywhere, so all masters keep identical point structures.
   R per master matches the statics formula (graduated at the light end).
3. Assemble a new designspace: 9 wght x 2 opsz masters at their NAMED user
   locations (identity axis map -- the nonlinearity now lives in master
   placement), then fontmake -o variable.

The result interpolates anywhere on (wght, opsz) and agrees with the static
family at every named weight.

Usage:
  python misc/tools/waterloo-var.py roman
  python misc/tools/waterloo-var.py italic
"""
import sys, os, math, argparse, subprocess
import ufoLib2
from fontTools.designspaceLib import (
    DesignSpaceDocument, SourceDescriptor, InstanceDescriptor, AxisDescriptor)
from fontTools.ttLib import TTFont
from fontTools.varLib import instancer
from fontmake.instantiator import Instantiator

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
UFODIR = os.path.join(ROOT, "build", "ufo")
BLEEDDIR = os.path.join(ROOT, "build", "ufo-bleed")
VARDIR = os.path.join(ROOT, "build", "fonts", "var")

WEIGHTS = [
  (100, "Thin"), (200, "ExtraLight"), (300, "Light"), (400, "Regular"),
  (500, "Medium"), (600, "SemiBold"), (700, "Bold"), (800, "ExtraBold"),
  (900, "Black"),
]
OPSZES = [14, 32]
MIN_RADIUS = 16.0
KAPPA = 0.5523  # circular-arc cubic constant
CORNER_SMOOTH_DEG = 168.0  # tangent turns sharper than this get filleted

FAMILY = "Waterloo Grotesk"


# ---------------------------------------------------------------- measurement
# (mirrors waterloo-statics.py; user-space via the compiled InterVariable)

def measure_stem(varfont_path, wght, opsz):
  font = TTFont(varfont_path, lazy=True)
  instancer.instantiateVariableFont(
    font, {"wght": wght, "opsz": opsz}, inplace=True, updateFontNames=False)
  glyf = font["glyf"]
  g = glyf["I"]
  g.recalcBounds(glyf)
  return g.xMax - g.xMin


def build_stem_curve(varfont_path, opsz):
  return [(w, measure_stem(varfont_path, w, opsz)) for w in range(100, 901, 50)]


def invert_stem(curve, stem_target):
  if stem_target <= curve[0][1]:
    return curve[0][0]
  for (w0, s0), (w1, s1) in zip(curve, curve[1:]):
    if s0 <= stem_target <= s1:
      t = 0 if s1 == s0 else (stem_target - s0) / (s1 - s0)
      return w0 + t * (w1 - w0)
  return curve[-1][0]


def stem_at(curve, wght):
  if wght <= curve[0][0]:
    return curve[0][1]
  for (w0, s0), (w1, s1) in zip(curve, curve[1:]):
    if w0 <= wght <= w1:
      t = (wght - w0) / (w1 - w0)
      return s0 + t * (s1 - s0)
  return curve[-1][1]


def user_to_design(axis, user):
  """Apply designspace axis map (piecewise linear)."""
  m = axis.map
  if not m:
    return user
  if user <= m[0][0]:
    return m[0][1]
  for (u0, d0), (u1, d1) in zip(m, m[1:]):
    if u0 <= user <= u1:
      t = 0 if u1 == u0 else (user - u0) / (u1 - u0)
      return d0 + t * (d1 - d0)
  return m[-1][1]


# ---------------------------------------------------------------- bleed core
#
# Contours are handled as flat point lists: (x, y, segmentType, smooth).
# segmentType None = offcurve. All decisions that alter structure are made
# once (on the default master) and replayed identically on every master.

def _anchors(pts):
  return [i for i, p in enumerate(pts) if p[2] is not None]


def _norm(dx, dy):
  d = math.hypot(dx, dy)
  return (dx / d, dy / d) if d > 1e-9 else (0.0, 0.0)


def _tangents(pts, i):
  """Unit tangent arriving at and leaving anchor i (via control polygon)."""
  n = len(pts)
  j = (i - 1) % n
  while pts[j][2] is None and False:
    pass
  # incoming: from previous point (control polygon)
  pin = pts[(i - 1) % n]
  t_in = _norm(pts[i][0] - pin[0], pts[i][1] - pin[1])
  pout = pts[(i + 1) % n]
  t_out = _norm(pout[0] - pts[i][0], pout[1] - pts[i][1])
  return t_in, t_out


def contour_sign(pts):
  """+1 if offsetting along left normal grows signed area, else -1 helper:
  returns signed area of the control polygon."""
  area = 0.0
  n = len(pts)
  for i in range(n):
    x0, y0 = pts[i][0], pts[i][1]
    x1, y1 = pts[(i + 1) % n][0], pts[(i + 1) % n][1]
    area += x0 * y1 - x1 * y0
  return area / 2.0


def corner_indices(pts):
  """Anchor indices whose tangents turn sharply (non-smooth corners)."""
  out = []
  cos_thresh = math.cos(math.radians(180.0 - CORNER_SMOOTH_DEG))
  for i in _anchors(pts):
    if pts[i][3]:  # marked smooth
      continue
    t_in, t_out = _tangents(pts, i)
    if t_in == (0.0, 0.0) or t_out == (0.0, 0.0):
      continue
    dot = t_in[0] * t_out[0] + t_in[1] * t_out[1]
    if dot < math.cos(math.radians(180.0 - CORNER_SMOOTH_DEG)) or dot < 0.98:
      # any turn of more than ~12 degrees counts as a corner
      out.append(i)
  return out


def offset_contour(pts, radius, ink_side):
  """Offset all points away from the ink. The SAME left/right rule applies
  to every contour: because holes are wound opposite to outers, a uniform
  rule grows outers and shrinks holes automatically. ink_side (+1 left of
  travel / -1 right of travel) is calibrated once per font from the 'o'."""
  n = len(pts)
  # normal per anchor via bisector; offcurves inherit neighbouring anchor's
  # displacement (first offcurve after an anchor -> that anchor; last before
  # the next anchor -> next anchor; middle handled by nearest).
  side = ink_side
  disp = [None] * n
  for i in _anchors(pts):
    t_in, t_out = _tangents(pts, i)
    bx, by = t_in[0] + t_out[0], t_in[1] + t_out[1]
    bl = math.hypot(bx, by)
    if bl < 1e-6:  # 180-degree spike; fall back to incoming normal
      nx, ny = -t_in[1], t_in[0]
      k = 1.0
    else:
      bx, by = bx / bl, by / bl
      nx, ny = -by, bx
      # miter factor 1/cos(theta/2), capped at 1.6 (90-degree corners need
      # 1.414; anything spikier gets filleted anyway and a tighter cap
      # avoids self-intersecting offsets at reflex corners)
      cos_half = math.sqrt(max(0.0, (1 + (t_in[0]*t_out[0] + t_in[1]*t_out[1])) / 2))
      k = 1.0 / max(cos_half, 0.625)
    disp[i] = (side * nx * radius * k, side * ny * radius * k)
  anchors = _anchors(pts)
  for ai, i in enumerate(anchors):
    nxt = anchors[(ai + 1) % len(anchors)]
    # offcurves strictly between i and nxt
    j = (i + 1) % n
    betw = []
    while j != nxt:
      if pts[j][2] is None:
        betw.append(j)
      j = (j + 1) % n
    for bi, j in enumerate(betw):
      src = i if bi < (len(betw) + 1) // 2 else nxt
      disp[j] = disp[src]
  out = []
  for p, d in zip(pts, disp):
    dx, dy = d if d else (0.0, 0.0)
    out.append((p[0] + dx, p[1] + dy, p[2], p[3]))
  return out


def fillet_corners(pts, corners, radius):
  """Replace each corner anchor with A + cubic(B): identical structure for
  identical `corners` lists, regardless of geometry."""
  if not corners:
    return pts
  n = len(pts)
  cset = set(corners)
  out = []
  for i, p in enumerate(pts):
    if i not in cset:
      out.append(p)
      continue
    t_in, t_out = _tangents(pts, i)
    # distance to previous/next point limits the cut
    pin = pts[(i - 1) % n]
    pout = pts[(i + 1) % n]
    lin = math.hypot(p[0] - pin[0], p[1] - pin[1])
    lout = math.hypot(pout[0] - p[0], pout[1] - p[1])
    d = min(radius, lin * 0.45, lout * 0.45)
    ax, ay = p[0] - t_in[0] * d, p[1] - t_in[1] * d
    bx, by = p[0] + t_out[0] * d, p[1] + t_out[1] * d
    c1 = (ax + t_in[0] * d * KAPPA, ay + t_in[1] * d * KAPPA, None, False)
    c2 = (bx - t_out[0] * d * KAPPA, by - t_out[1] * d * KAPPA, None, False)
    out.append((ax, ay, p[2], True))
    out.append(c1)
    out.append(c2)
    out.append((bx, by, "curve", True))
  return out


def glyph_points(glyph):
  """ufoLib2 glyph -> list of contours as point lists."""
  res = []
  for c in glyph.contours:
    res.append([(pt.x, pt.y, pt.type, pt.smooth) for pt in c.points])
  return res


def set_glyph_points(glyph, contours):
  from ufoLib2.objects import Contour
  from ufoLib2.objects.point import Point
  glyph.contours.clear()
  for cpts in contours:
    contour = Contour()
    for (x, y, t, s) in cpts:
      contour.points.append(Point(x=round(x), y=round(y), type=t, smooth=s))
    glyph.contours.append(contour)


def ink_side_for_font(ufo):
  """Determine which normal side grows ink, from the 'o' outer contour."""
  o = ufo["o"]
  cts = glyph_points(o)
  outer = max(cts, key=lambda c: abs(contour_sign(c)))
  # try left side (+1): offset and see if |area| grows
  test = offset_contour(outer, 10.0, +1)
  return +1 if abs(contour_sign(test)) > abs(contour_sign(outer)) else -1


def bleed_ufo(ufo, radius, decisions, ink_side):
  """Apply offset+fillet to every contour glyph. `decisions` maps glyph name
  -> list of corner-index lists (one per contour), computed on the default
  master. Returns count of bled glyphs."""
  count = 0
  for glyph in ufo:
    name = glyph.name
    if name not in decisions or not glyph.contours:
      continue
    contours = glyph_points(glyph)
    corner_lists = decisions[name]
    if len(corner_lists) != len(contours):
      continue  # structure mismatch; skip defensively
    newc = []
    for cpts, corners in zip(contours, corner_lists):
      c1 = offset_contour(cpts, radius, ink_side)
      c2 = fillet_corners(c1, corners, radius)
      newc.append(c2)
    set_glyph_points(glyph, newc)
    count += 1
  return count


def compute_decisions(ufo):
  """Corner decisions from the default master; also records which glyphs
  have contours at all."""
  decisions = {}
  for glyph in ufo:
    if not glyph.contours:
      continue
    contours = glyph_points(glyph)
    decisions[glyph.name] = [corner_indices(c) for c in contours]
  return decisions


def structure_signature(ufo, names):
  sig = {}
  for name in names:
    g = ufo[name]
    sig[name] = tuple(
      tuple(pt.type for pt in c.points) for c in g.contours
    )
  return sig


# ---------------------------------------------------------------- pipeline

def main(argv):
  ap = argparse.ArgumentParser()
  ap.add_argument("style", choices=["roman", "italic"])
  ap.add_argument("--fast", action="store_true",
                  help="only wght 100/400/900 masters (debug)")
  args = ap.parse_args(argv[1:])
  italic = args.style == "italic"

  ds_path = os.path.join(
    UFODIR, "Inter-Italic.var.designspace" if italic else "Inter-Roman.var.designspace")
  var_ttf = os.path.join(
    VARDIR, "InterVariable-Italic.ttf" if italic else "InterVariable.ttf")

  print("measuring stem curves...")
  curves = {opsz: build_stem_curve(var_ttf, opsz) for opsz in OPSZES}
  stems_text = dict(curves[14])
  r_full = (stems_text[700] - stems_text[400]) / 2.0
  print("full bleed radius: %.1f" % r_full)

  ds = DesignSpaceDocument.fromfile(ds_path)
  wght_axis = next(a for a in ds.axes if a.tag == "wght")
  ds.loadSourceFonts(ufoLib2.Font.open)
  inst = Instantiator.from_designspace(ds, round_geometry=True)

  os.makedirs(BLEEDDIR, exist_ok=True)
  feat_link = os.path.join(BLEEDDIR, "features")
  if not os.path.islink(feat_link):
    os.symlink(os.path.join("..", "..", "src", "features"), feat_link)

  weights = [w for w in WEIGHTS if not args.fast or w[0] in (100, 400, 900)]

  # 1. generate skeleton instance UFOs
  masters = []  # (path, user_wght, opsz, radius)
  for opsz in OPSZES:
    curve = curves[opsz]
    stem_min = curve[0][1]
    for wval, wname in weights:
      stem_named = stem_at(curve, wval)
      r = max(MIN_RADIUS, min(r_full, (stem_named - stem_min) / 2.0))
      src_user = invert_stem(curve, stem_named - 2 * r)
      src_design = user_to_design(wght_axis, src_user)
      idesc = InstanceDescriptor()
      idesc.familyName = FAMILY
      idesc.styleName = "%s-%d" % (wname, opsz)
      idesc.location = {"Weight": src_design, "Optical size": float(opsz)}
      ufo = inst.generate_instance(idesc)
      suffix = "Italic" if italic else ""
      path = os.path.join(
        BLEEDDIR, "WG-%s%s-%d.ufo" % (wname, suffix, opsz))
      ufo.save(path, overwrite=True)
      masters.append((path, wval, opsz, r))
      print("skeleton %s (user %.1f, design %.1f) R=%.1f" %
            (os.path.basename(path), src_user, src_design, r))

  # 2. bleed
  default_path = next(p for (p, w, o, r) in masters if w == 400 and o == 14)
  default_ufo = ufoLib2.Font.open(default_path)
  decisions = compute_decisions(default_ufo)
  ink_side = ink_side_for_font(default_ufo)
  print("ink side: %+d; %d contour glyphs" % (ink_side, len(decisions)))
  sig_ref = None
  for (path, wval, opsz, r) in masters:
    ufo = ufoLib2.Font.open(path)
    n = bleed_ufo(ufo, r, decisions, ink_side)
    ufo.save(path, overwrite=True)
    sig = structure_signature(ufo, sorted(decisions.keys()))
    if sig_ref is None:
      sig_ref = sig
    else:
      bad = [k for k in sig_ref if sig[k] != sig_ref[k]]
      assert not bad, "structure mismatch in %s: %r" % (path, bad[:10])
    print("bled %s: %d glyphs" % (os.path.basename(path), n))

  # 3. new designspace + fontmake
  out_ds = DesignSpaceDocument()
  a_opsz = AxisDescriptor()
  a_opsz.tag, a_opsz.name = "opsz", "Optical size"
  a_opsz.minimum, a_opsz.default, a_opsz.maximum = 14, 14, 32
  out_ds.addAxis(a_opsz)
  a_wght = AxisDescriptor()
  a_wght.tag, a_wght.name = "wght", "Weight"
  a_wght.minimum, a_wght.default, a_wght.maximum = (
    weights[0][0], 400, weights[-1][0])
  out_ds.addAxis(a_wght)
  for (path, wval, opsz, r) in masters:
    s = SourceDescriptor()
    s.path = path
    s.familyName = FAMILY
    s.styleName = os.path.basename(path)[3:-4]
    s.location = {"Weight": float(wval), "Optical size": float(opsz)}
    if wval == 400 and opsz == 14:
      s.copyInfo = s.copyGroups = s.copyLib = s.copyFeatures = True
    out_ds.addSource(s)
  for wval, wname in weights:
    i = InstanceDescriptor()
    i.familyName = FAMILY
    i.styleName = wname + (" Italic" if italic else "")
    i.location = {"Weight": float(wval), "Optical size": 14.0}
    out_ds.addInstance(i)
  ds_out_path = os.path.join(
    BLEEDDIR, "WaterlooGrotesk-%s.designspace" % ("Italic" if italic else "Roman"))
  out_ds.write(ds_out_path)

  out_ttf = os.path.join(
    VARDIR, "WaterlooGroteskVariable%s.ttf" % ("-Italic" if italic else ""))
  print("fontmake variable...")
  subprocess.run(
    [os.path.join(ROOT, "build", "venv", "bin", "fontmake"),
     "-o", "variable", "-m", ds_out_path, "--output-path", out_ttf,
     "--verbose", "WARNING", "--no-autohint", "--production-names"],
    check=True, cwd=ROOT)
  print("wrote", out_ttf)

  # 4. names
  font = TTFont(out_ttf, lazy=False)
  name = font["name"]
  sub = "Italic" if italic else "Regular"
  ps = "WaterlooGroteskVar" + ("-Italic" if italic else "")
  for nid, val in ((1, FAMILY), (2, sub), (3, "%s Variable %s" % (FAMILY, sub)),
                   (4, "%s %s" % (FAMILY, sub)), (6, ps),
                   (16, FAMILY), (17, sub), (25, ps.replace("-", ""))):
    name.setName(val, nid, 3, 1, 0x409)
    name.setName(val, nid, 1, 0, 0)
  font.save(out_ttf)
  w = TTFont(out_ttf, lazy=False)
  w.flavor = "woff2"
  w.save(out_ttf.replace(".ttf", ".woff2"))
  print("done: %s (+woff2)" % out_ttf)


if __name__ == "__main__":
  main(sys.argv)
