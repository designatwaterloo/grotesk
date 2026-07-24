"""
Waterloo Grotesk source transform.

Runs on the master UFOs of a designspace, after postprocess-designspace.py:

1. Square punctuation by default: swaps outline geometry between every glyph
   and its .ss07 ("square punctuation") alternate, for pairs where the
   alternate carries real contours. Composite glyphs (e.g. Adieresis.ss07)
   follow automatically through their components, so ss07 becomes a
   "round punctuation" restore-toggle.

2. Goop ligature: builds o_o.dlig in every master as two overlapping "o"
   components -- the Waterloo infinity pair. Referenced by src/features
   dlig.fea (any "oo", opt-in) and calt.fea (auto, only inside the word
   "waterloo"/"Waterloo").

Kept at the UFO stage (not in the .glyphspackage sources) so that upstream
Inter merges stay trivial.
"""
import sys, argparse, math
import defcon
from fontTools.designspaceLib import DesignSpaceDocument

# Goop geometry (matched to the Design Waterloo wordmark SVG):
# the two o's keep their exact normal advance -- tracking is untouched --
# and a smooth ink neck joins the facing bulges at the waist.
GOOP_NAME = "o_o.dlig"
GOOP_ATTACH_FRAC = 0.85  # attachment height as fraction of ink height at flange x
GOOP_NECK_RING = 1.2     # waist half-height = this * ring thickness ...
GOOP_NECK_MAX = 0.6      # ... clamped to this fraction of the o half-height


def glyph_geometry(g):
  contours = [[(p.x, p.y, p.segmentType, p.smooth) for p in c] for c in g]
  components = [(c.baseGlyph, c.transformation) for c in g.components]
  anchors = [dict(a) for a in g.anchors]
  return (contours, components, anchors, g.width)


def set_glyph_geometry(g, geo):
  contours, components, anchors, width = geo
  g.clearContours()
  g.clearComponents()
  g.clearAnchors()
  pen = g.getPointPen()
  for contour in contours:
    pen.beginPath()
    for x, y, segmentType, smooth in contour:
      pen.addPoint((x, y), segmentType=segmentType, smooth=smooth)
    pen.endPath()
  for baseGlyph, transformation in components:
    pen.addComponent(baseGlyph, transformation)
  for a in anchors:
    g.appendAnchor(dict(a))
  g.width = width


SWAP_MARKER = "com.designwaterloo.ss07swapped"


def swap_square_defaults(ufo):
  if ufo.lib.get(SWAP_MARKER):
    return []  # already swapped; never double-apply
  ufo.lib[SWAP_MARKER] = True
  swapped = []
  for name in list(ufo.keys()):
    if not name.endswith(".ss07"):
      continue
    base_name = name[: -len(".ss07")]
    if base_name not in ufo:
      continue
    alt = ufo[name]
    if len(alt) == 0:  # pure composite; follows its components
      continue
    base = ufo[base_name]
    geo_base = glyph_geometry(base)
    geo_alt = glyph_geometry(alt)
    set_glyph_geometry(base, geo_alt)
    set_glyph_geometry(alt, geo_base)
    swapped.append(base_name)
  return swapped


def _flatten_contour(contour, steps=24):
  """Contour -> polyline point list (approximating curves)."""
  segs = contour.segments
  pts = []
  if not segs:
    return pts
  # find start point: last on-curve of last segment
  prev = None
  for seg in segs:
    on = [p for p in seg if p.segmentType is not None]
    if not on:
      return pts  # all-offcurve (TrueType-style); not expected in Inter UFOs
    prev = (on[-1].x, on[-1].y)
  for seg in segs:
    on = [p for p in seg if p.segmentType is not None][-1]
    offs = [(p.x, p.y) for p in seg if p.segmentType is None]
    p3 = (on.x, on.y)
    if len(offs) == 2:  # cubic
      p0, (x1, y1), (x2, y2) = prev, offs[0], offs[1]
      for i in range(1, steps + 1):
        t = i / steps
        mt = 1 - t
        x = mt**3 * p0[0] + 3 * mt**2 * t * x1 + 3 * mt * t**2 * x2 + t**3 * p3[0]
        y = mt**3 * p0[1] + 3 * mt**2 * t * y1 + 3 * mt * t**2 * y2 + t**3 * p3[1]
        pts.append((x, y))
    elif len(offs) == 1:  # quadratic
      p0, (x1, y1) = prev, offs[0]
      for i in range(1, steps + 1):
        t = i / steps
        mt = 1 - t
        x = mt**2 * p0[0] + 2 * mt * t * x1 + t**2 * p3[0]
        y = mt**2 * p0[1] + 2 * mt * t * y1 + t**2 * p3[1]
        pts.append((x, y))
    else:  # line
      pts.append(p3)
    prev = p3
  return pts


def measure_o(ufo):
  """Returns (xMin, xMax, yMin, yMax, ring_thickness) of the 'o' glyph."""
  o = ufo["o"]
  polylines = [_flatten_contour(c) for c in o]
  all_pts = [p for poly in polylines for p in poly]
  xmin = min(p[0] for p in all_pts)
  xmax = max(p[0] for p in all_pts)
  ymin = min(p[1] for p in all_pts)
  ymax = max(p[1] for p in all_pts)
  ymid = (ymin + ymax) / 2
  # x positions where contours cross ymid
  crossings = []
  for poly in polylines:
    n = len(poly)
    for i in range(n):
      (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % n]
      # half-open test so a crossing at a shared vertex counts once
      if (y0 <= ymid < y1) or (y1 <= ymid < y0):
        t = (ymid - y0) / (y1 - y0)
        crossings.append(x0 + t * (x1 - x0))
  crossings.sort()
  if len(crossings) >= 4:
    thickness = crossings[1] - crossings[0]  # left ring wall
  else:
    thickness = (xmax - xmin) * 0.18  # fallback guess
  return xmin, xmax, ymin, ymax, thickness


def _outer_winding_sign(o):
  polylines = [_flatten_contour(c) for c in o]
  outer = max(polylines, key=lambda poly: abs(_shoelace(poly)))
  return 1 if _shoelace(outer) >= 0 else -1


def _shoelace(poly):
  area = 0.0
  n = len(poly)
  for i in range(n):
    x0, y0 = poly[i]
    x1, y1 = poly[(i + 1) % n]
    area += x0 * y1 - x1 * y0
  return area / 2.0


def add_goop_ligature(ufo):
  """o_o.dlig: two o components at the o's EXACT normal advance (tracking is
  untouched), plus an ink-neck contour joining the facing bulges at the
  waist -- per the Design Waterloo wordmark. The neck scales with the ring
  thickness, so it goes wiry at Thin and chunky at Black."""
  if GOOP_NAME in ufo:
    del ufo[GOOP_NAME]
  o = ufo["o"]
  xmin, xmax, ymin, ymax, t = measure_o(ufo)
  dx = o.width
  cx = (xmin + xmax) / 2.0
  cy = (ymin + ymax) / 2.0
  a_out = (xmax - xmin) / 2.0
  b_out = (ymax - ymin) / 2.0
  h = min(GOOP_NECK_RING * t, GOOP_NECK_MAX * b_out)  # waist half-height
  # flange x: mid ring wall, but always beyond the counter so the bridge's
  # side edges never cross into it (counter max half-width = a_out - t)
  fx = max(a_out - 0.5 * t, a_out - t + 6.0)
  xl = cx + fx          # flange x on left o (facing right)
  xr = dx + cx - fx     # flange x on right o (facing left)
  # ink height at flange x (outer ellipse), attachment inside it
  y_outer = b_out * math.sqrt(max(0.0, 1.0 - (fx / a_out) ** 2))
  y_att = max(GOOP_ATTACH_FRAC * y_outer, min(h * 1.15, 0.95 * y_outer))
  # cubic through midpoint at waist: yc = (8h - 2*y_att)/6
  yc = (8.0 * h - 2.0 * y_att) / 6.0
  cxl = xl + 0.45 * (xr - xl)
  cxr = xr - 0.45 * (xr - xl)
  # counter-clockwise construction (positive shoelace), starting top-left
  bridge_ccw = [
    (xl, cy + y_att, "line", False),
    (cxl, cy + yc, None, False),
    (cxr, cy + yc, None, False),
    (xr, cy + y_att, "curve", False),
    (xr, cy - y_att, "line", False),
    (cxr, cy - yc, None, False),
    (cxl, cy - yc, None, False),
    (xl, cy - y_att, "curve", False),
  ]
  # the bridge must wind the SAME way as the o's outer contour, or the
  # overlap cancels under nonzero fill; flip construction if needed
  bridge_sign = 1 if _shoelace([(p[0], p[1]) for p in bridge_ccw]) >= 0 else -1
  if bridge_sign != _outer_winding_sign(o):
    bridge_ccw = [
      (xl, cy - y_att, "line", False),
      (cxl, cy - yc, None, False),
      (cxr, cy - yc, None, False),
      (xr, cy - y_att, "curve", False),
      (xr, cy + y_att, "line", False),
      (cxr, cy + yc, None, False),
      (cxl, cy + yc, None, False),
      (xl, cy + y_att, "curve", False),
    ]
  g = ufo.newGlyph(GOOP_NAME)
  g.width = 2 * o.width
  pen = g.getPointPen()
  pen.addComponent("o", (1, 0, 0, 1, 0, 0))
  pen.addComponent("o", (1, 0, 0, 1, dx, 0))
  pen.beginPath()
  for (x, y, seg, smooth) in bridge_ccw:
    pen.addPoint((round(x), round(y)), segmentType=seg, smooth=smooth)
  pen.endPath()
  order = ufo.lib.get("public.glyphOrder")
  if order is not None and GOOP_NAME not in order:
    order.append(GOOP_NAME)
    ufo.lib["public.glyphOrder"] = order
  return dx, t


def main(argv):
  ap = argparse.ArgumentParser()
  ap.add_argument("designspace", help="designspace whose source UFOs to transform")
  args = ap.parse_args(argv[1:])
  ds = DesignSpaceDocument.fromfile(args.designspace)
  seen = set()
  for source in ds.sources:
    if source.path in seen:
      continue
    seen.add(source.path)
    ufo = defcon.Font(source.path)
    swapped = swap_square_defaults(ufo)
    dx, thickness = add_goop_ligature(ufo)
    ufo.save()
    print(
      "waterloo-transform: %s: %d square-punct swaps; %s dx=%d (ring=%.0f)"
      % (source.filename, len(swapped), GOOP_NAME, dx, thickness)
    )


if __name__ == "__main__":
  main(sys.argv)
