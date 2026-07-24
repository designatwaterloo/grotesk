"""
Waterloo Grotesk static builder.

Takes the built InterVariable TTFs and produces the Waterloo Grotesk static
family with the inkbleed treatment:

  For every named weight, instead of using Inter's own instance we
  instantiate a LIGHTER skeleton and bulk it back up with a stroke
  (round caps/joins) unioned with the fill -- the same recipe used in
  Figma for the Design Waterloo wordmark (e.g. "Bold" = Regular skeleton
  + stroke). Corners round off, joins goop up, counters tighten, but the
  apparent stem weight matches the weight name.

The stroke radius R defaults to (stem(700) - stem(400)) / 2 measured on the
"I" of the text-optical-size roman -- i.e. exactly "Regular + stroke = Bold".

Usage:
  python waterloo-statics.py build/fonts/var/InterVariable.ttf \
      build/fonts/var/InterVariable-Italic.ttf -o build/fonts/waterloo
"""
import sys, os, argparse
from fontTools.ttLib import TTFont
from fontTools.varLib import instancer
from fontTools.pens.ttGlyphPen import TTGlyphPen
import pathops

WEIGHTS = [
  (100, "Thin"),
  (200, "ExtraLight"),
  (300, "Light"),
  (400, "Regular"),
  (500, "Medium"),
  (600, "SemiBold"),
  (700, "Bold"),
  (800, "ExtraBold"),
  (900, "Black"),
]

FAMILY = "Waterloo Grotesk"
VENDOR_SUFFIX = ""  # kept for possible future "Beta" etc.

# opsz axis values in InterVariable
OPSZ_TEXT = 14
OPSZ_DISPLAY = 32


def measure_stem(varfont_path, wght, opsz):
  """Stem width of "I" at a given location, in font units."""
  font = TTFont(varfont_path, lazy=True)
  inst = instancer.instantiateVariableFont(
    font, {"wght": wght, "opsz": opsz}, inplace=True, updateFontNames=False
  )
  glyf = inst["glyf"]
  g = glyf["I"]
  if g.numberOfContours <= 0:
    raise ValueError("I is not a simple glyph")
  g.recalcBounds(glyf)
  return g.xMax - g.xMin


def build_stem_curve(varfont_path, opsz, samples=None):
  """Sampled (wght -> stem) piecewise-linear curve."""
  if samples is None:
    samples = list(range(100, 901, 50))
  return [(w, measure_stem(varfont_path, w, opsz)) for w in samples]


def invert_stem(curve, stem_target):
  """wght coordinate whose stem is stem_target (clamped to curve range)."""
  if stem_target <= curve[0][1]:
    return curve[0][0]
  for (w0, s0), (w1, s1) in zip(curve, curve[1:]):
    if s0 <= stem_target <= s1:
      if s1 == s0:
        return w0
      t = (stem_target - s0) / (s1 - s0)
      return w0 + t * (w1 - w0)
  return curve[-1][0]


def sk_path_from_ttglyph(g, glyf_table):
  path = pathops.Path()
  g.draw(path.getPen(glyphSet=None), glyf_table)
  return path


def inkbleed_glyf(font, radius):
  """Stroke+union every simple glyph. Composites follow their bases."""
  glyf = font["glyf"]
  hmtx = font["hmtx"]
  modified = 0
  for name in font.getGlyphOrder():
    g = glyf[name]
    if g.numberOfContours <= 0:  # empty or composite
      continue
    src = sk_path_from_ttglyph(g, glyf)
    stroked = pathops.Path(src)
    stroked.stroke(
      radius * 2.0,
      pathops.LineCap.ROUND_CAP,
      pathops.LineJoin.ROUND_JOIN,
      4.0,
    )
    # skia stroking emits conics, which op() cannot wind-fix
    stroked.convertConicsToQuads(0.25)
    try:
      result = pathops.op(
        src, stroked, pathops.PathOp.UNION, fix_winding=True
      )
    except pathops.PathOpsError:
      result = src  # keep original on pathological outlines
      print("  warn: union failed for %s; left unbled" % name)
    result.convertConicsToQuads(0.25)
    pen = TTGlyphPen(None)
    result.draw(pen)
    new_glyph = pen.glyph()
    glyf[name] = new_glyph
    new_glyph.recalcBounds(glyf)
    adv, _lsb = hmtx[name]
    hmtx[name] = (adv, new_glyph.xMin)
    modified += 1
  return modified


def set_name(font, nameID, value):
  font["name"].setName(value, nameID, 3, 1, 0x409)
  font["name"].setName(value, nameID, 1, 0, 0)


def rename(font, family, weight_name, weight_value, italic):
  style = weight_name + (" Italic" if italic else "")
  if weight_name == "Regular" and italic:
    style = "Italic"
  # legacy 4-style-linked names
  if weight_name in ("Regular", "Bold"):
    legacy_family = family
    legacy_style = style
  else:
    legacy_family = "%s %s" % (family, weight_name)
    legacy_style = "Italic" if italic else "Regular"
  ps = "%s-%s" % (family.replace(" ", ""), style.replace(" ", ""))
  set_name(font, 1, legacy_family)
  set_name(font, 2, legacy_style)
  set_name(font, 3, "%s %s" % (family, style))
  set_name(font, 4, "%s %s" % (family, style))
  set_name(font, 6, ps)
  set_name(font, 16, family)
  set_name(font, 17, style)
  os2 = font["OS/2"]
  os2.usWeightClass = weight_value
  # fsSelection: bit0 italic, bit5 bold, bit6 regular
  fs = os2.fsSelection & ~0b1100001
  head_mac = font["head"].macStyle & ~0b11
  is_bold = weight_name == "Bold"
  if is_bold:
    fs |= 0b0100000
    head_mac |= 0b01
  if italic:
    fs |= 0b0000001
    head_mac |= 0b10
  if not is_bold and not italic:
    fs |= 0b1000000
  os2.fsSelection = fs
  font["head"].macStyle = head_mac
  return ps


# the lightest weights cannot drop to a lighter skeleton (nothing below
# wght 100 exists), so the bleed radius ramps down there instead --
# each weight gets the largest bleed that still lands on its named stem
# width, floored so even Thin keeps the rounded-corner identity.
MIN_RADIUS = 16.0


def build_static(varfont_path, out_dir, family, opsz, italic, curve, radius):
  built = []
  stem_min = curve[0][1]
  for weight_value, weight_name in WEIGHTS:
    # apparent stem for this weight per stock Inter
    stem_named = dict(curve).get(weight_value)
    if stem_named is None:
      stem_named = invert_stem_value(curve, weight_value)
    r = max(MIN_RADIUS, min(radius, (stem_named - stem_min) / 2.0))
    src_stem = stem_named - 2 * r
    src_wght = invert_stem(curve, src_stem)
    font = TTFont(varfont_path, lazy=False)
    instancer.instantiateVariableFont(
      font, {"wght": src_wght, "opsz": opsz}, inplace=True, updateFontNames=False
    )
    n = inkbleed_glyf(font, r)
    ps = rename(font, family, weight_name, weight_value, italic)
    out_ttf = os.path.join(out_dir, "%s.ttf" % ps)
    font.save(out_ttf)
    woff = TTFont(out_ttf, lazy=False)
    woff.flavor = "woff2"
    out_woff = os.path.join(out_dir, "%s.woff2" % ps)
    woff.save(out_woff)
    built.append(out_ttf)
    print(
      "%s: wght %d <- skeleton %.1f (stem %d -> %d), %d glyphs bled"
      % (ps, weight_value, src_wght, stem_named, int(src_stem), n)
      + "  R=%.1f" % r
    )
  return built


def invert_stem_value(curve, wght):
  """stem at arbitrary wght by linear interpolation of the sample curve."""
  if wght <= curve[0][0]:
    return curve[0][1]
  for (w0, s0), (w1, s1) in zip(curve, curve[1:]):
    if w0 <= wght <= w1:
      t = (wght - w0) / (w1 - w0)
      return s0 + t * (s1 - s0)
  return curve[-1][1]


def main(argv):
  ap = argparse.ArgumentParser()
  ap.add_argument("roman", help="InterVariable.ttf")
  ap.add_argument("italic", nargs="?", help="InterVariable-Italic.ttf")
  ap.add_argument("-o", "--out", default="build/fonts/waterloo")
  ap.add_argument(
    "--bleed",
    type=float,
    default=None,
    help="stroke radius in font units (default: (stem700-stem400)/2)",
  )
  args = ap.parse_args(argv[1:])
  os.makedirs(args.out, exist_ok=True)

  print("measuring stem curve (text roman)...")
  curve_text = build_stem_curve(args.roman, OPSZ_TEXT)
  curve_disp = build_stem_curve(args.roman, OPSZ_DISPLAY)
  radius = args.bleed
  if radius is None:
    stems = dict(curve_text)
    radius = (stems[700] - stems[400]) / 2.0
  print("inkbleed radius: %.1f units (UPM %d)" % (radius, 2048))

  jobs = [(args.roman, FAMILY, OPSZ_TEXT, False, curve_text),
          (args.roman, FAMILY + " Display", OPSZ_DISPLAY, False, curve_disp)]
  if args.italic:
    curve_text_i = build_stem_curve(args.italic, OPSZ_TEXT)
    curve_disp_i = build_stem_curve(args.italic, OPSZ_DISPLAY)
    jobs += [(args.italic, FAMILY, OPSZ_TEXT, True, curve_text_i),
             (args.italic, FAMILY + " Display", OPSZ_DISPLAY, True, curve_disp_i)]

  for path, family, opsz, italic, curve in jobs:
    print("== %s%s (opsz %g) ==" % (family, " Italic" if italic else "", opsz))
    build_static(path, args.out, family, opsz, italic, curve, radius)


if __name__ == "__main__":
  main(sys.argv)
