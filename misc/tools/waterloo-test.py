"""
Waterloo Grotesk sanity tests + visual proofs.

Shaping tests (uharfbuzz):
  - "Waterloo" goops (o_o.dlig) with default features
  - "book"/"loop" do NOT goop by default
  - any "oo" goops with dlig on
  - "?!" becomes interrobang with dlig on
  - period is square by default; ss07 restores the round dot
  - ss04 (disambiguation) exists but is off by default

Proof sheets: SVG renders of sample lines for eyeballing.

Usage: python misc/tools/waterloo-test.py [font.ttf ...]
"""
import sys, os
import uharfbuzz as hb
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_FONTS = [os.path.join(ROOT, "build", "fonts", "waterloo", f)
                 for f in ("WaterlooGrotesk-Regular.ttf", "WaterlooGrotesk-Bold.ttf")]

FAILS = []


def shape(path, text, features=None):
  with open(path, "rb") as f:
    data = f.read()
  face = hb.Face(data)
  font = hb.Font(face)
  buf = hb.Buffer()
  buf.add_str(text)
  buf.guess_segment_properties()
  hb.shape(font, buf, features or {})
  ttf = TTFont(path, lazy=True)
  order = ttf.getGlyphOrder()
  return [order[i.codepoint] for i in buf.glyph_infos]


def check(desc, ok):
  print(("  PASS  " if ok else "  FAIL  ") + desc)
  if not ok:
    FAILS.append(desc)


def run_checks(path):
  print(os.path.basename(path))
  g = shape(path, "Waterloo")
  check("Waterloo goops by default (calt): %r" % g, "o_o.dlig" in g)
  g = shape(path, "waterloo")
  check("lowercase waterloo goops too", "o_o.dlig" in g)
  g = shape(path, "book loop good")
  check("book/loop/good do NOT goop by default", "o_o.dlig" not in g)
  g = shape(path, "book", {"dlig": True})
  check("dlig goops any oo", "o_o.dlig" in g)
  g = shape(path, "Waterloo", {"calt": False})
  check("calt off disables the goop", "o_o.dlig" not in g)
  g = shape(path, "?!", {"dlig": True})
  check("?! interrobang with dlig", "interrobang" in g or "uni203D" in g)

  # square punctuation default: a square(ish) dot fills ~90%+ of its bbox,
  # a round dot ~78% -- robust even after the bleed pass rounds corners.
  ttf = TTFont(path, lazy=True)
  glyphset = ttf.getGlyphSet()

  def fill_ratio(name):
    from fontTools.pens.areaPen import AreaPen
    from fontTools.pens.boundsPen import BoundsPen
    ap = AreaPen(glyphset)
    bp = BoundsPen(glyphset)
    glyphset[name].draw(ap)
    glyphset[name].draw(bp)
    (x0, y0, x1, y1) = bp.bounds
    return abs(ap.value) / ((x1 - x0) * (y1 - y0))

  rb, ra = fill_ratio("period"), fill_ratio("period.ss07")
  check("period is square by default (fill %.2f)" % rb, rb > 0.85)
  check("period.ss07 is round (fill %.2f, restore toggle)" % ra, ra < 0.85)
  g = shape(path, ".", {"ss07": True})
  check("ss07 remaps period -> period.ss07", g == ["period.ss07"])
  gsub_feats = set()
  if "GSUB" in ttf:
    for fr in ttf["GSUB"].table.FeatureList.FeatureRecord:
      gsub_feats.add(fr.FeatureTag)
  check("ss04 present but opt-in", "ss04" in gsub_feats)
  check("calt present", "calt" in gsub_feats)


def proof(path, out_svg, lines=None, size=110):
  """Render sample lines to an SVG proof sheet (no positioning niceties)."""
  if lines is None:
    lines = [
      ("Waterloo Grotesk", {}),
      ("Design Waterloo 47", {}),
      ("book loop good.", {}),
      ("book loop good.", {"dlig": True}),
      ("«Hello.» (period?!)", {}),
      ("0123456789 ?!", {}),
    ]
  ttf = TTFont(path, lazy=True)
  upm = ttf["head"].unitsPerEm
  scale = size / upm
  glyphset = ttf.getGlyphSet()
  order = ttf.getGlyphOrder()
  with open(path, "rb") as f:
    face = hb.Face(f.read())
  hbfont = hb.Font(face)
  rows = []
  width = 100
  y = size * 1.25
  for text, feats in lines:
    buf = hb.Buffer()
    buf.add_str(text)
    buf.guess_segment_properties()
    hb.shape(hbfont, buf, feats or {})
    x = 40
    parts = []
    for info, pos in zip(buf.glyph_infos, buf.glyph_positions):
      name = order[info.codepoint]
      pen = SVGPathPen(glyphset)
      glyphset[name].draw(pen)
      d = pen.getCommands()
      if d:
        parts.append(
          '<path transform="translate(%.1f %.1f) scale(%.5f -%.5f)" d="%s"/>'
          % (x + pos.x_offset * scale, y - pos.y_offset * scale, scale, scale, d)
        )
      x += pos.x_advance * scale
    width = max(width, x + 40)
    label = text.replace("&", "&amp;").replace("<", "&lt;")
    ftxt = ",".join(k for k, v in (feats or {}).items() if v) or "default"
    rows.append(
      '<text x="40" y="%.1f" font-size="12" fill="#888" font-family="monospace">%s  [%s]</text>'
      % (y - size * 1.02, label, ftxt)
    )
    rows.append("<g>%s</g>" % "".join(parts))
    y += size * 1.55
  svg = (
    '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
    '<rect width="100%%" height="100%%" fill="white"/>%s</svg>'
    % (int(width), int(y), int(width), int(y), "".join(rows))
  )
  with open(out_svg, "w") as f:
    f.write(svg)
  print("proof: %s" % out_svg)


def main(argv):
  fonts = argv[1:] or DEFAULT_FONTS
  for p in fonts:
    run_checks(p)
    out = os.path.join(
      os.path.dirname(p), "proof-%s.svg" % os.path.basename(p).replace(".ttf", "")
    )
    proof(p, out)
  if FAILS:
    print("\n%d FAILURES" % len(FAILS))
    sys.exit(1)
  print("\nall checks passed")


if __name__ == "__main__":
  main(sys.argv)
