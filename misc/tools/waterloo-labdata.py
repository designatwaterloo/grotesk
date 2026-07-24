"""
Generates lab/ data for the Waterloo Grotesk tester site:

- copies the built woff2 fonts into lab/fonts/
- writes lab/data/glyphs-<Weight>.json for the text roman weights, with an
  SVG path for EVERY glyph (encoded or not) for the debug glyph grid
- writes woff2 versions of the stock-geometry InterVariable for the
  before/after compare mode

Usage: python misc/tools/waterloo-labdata.py
"""
import os, json, shutil, glob
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
FONTDIR = os.path.join(ROOT, "build", "fonts", "waterloo")
VARDIR = os.path.join(ROOT, "build", "fonts", "var")
LAB = os.path.join(ROOT, "lab")

WEIGHTS = ["Thin", "ExtraLight", "Light", "Regular", "Medium", "SemiBold",
           "Bold", "ExtraBold", "Black"]


def glyph_data(ttf_path):
  font = TTFont(ttf_path, lazy=True)
  upm = font["head"].unitsPerEm
  hhea = font["hhea"]
  cmap = font.getBestCmap()
  rev = {}
  for cp, name in cmap.items():
    rev.setdefault(name, cp)
  glyphset = font.getGlyphSet()
  out = []
  for name in font.getGlyphOrder():
    g = glyphset[name]
    pen = SVGPathPen(glyphset, ntos=lambda f: ("%.1f" % f).rstrip("0").rstrip("."))
    try:
      g.draw(pen)
      d = pen.getCommands()
    except Exception:
      d = ""
    out.append({
      "n": name,
      "u": rev.get(name),
      "a": g.width,
      "d": d,
    })
  return {
    "upm": upm,
    "ascender": hhea.ascender,
    "descender": hhea.descender,
    "glyphs": out,
  }


def main():
  os.makedirs(os.path.join(LAB, "fonts"), exist_ok=True)
  os.makedirs(os.path.join(LAB, "data"), exist_ok=True)

  for f in glob.glob(os.path.join(FONTDIR, "*.woff2")):
    shutil.copy(f, os.path.join(LAB, "fonts", os.path.basename(f)))
  print("copied %d woff2 files" % len(glob.glob(os.path.join(LAB, "fonts", "*.woff2"))))

  # Waterloo variable fonts (phase 2), if built
  for f in glob.glob(os.path.join(VARDIR, "WaterlooGroteskVariable*.woff2")):
    shutil.copy(f, os.path.join(LAB, "fonts", os.path.basename(f)))

  # stock InterVariable for compare mode
  for name in ("InterVariable.ttf", "InterVariable-Italic.ttf"):
    src = os.path.join(VARDIR, name)
    if os.path.exists(src):
      f = TTFont(src, lazy=False)
      f.flavor = "woff2"
      dst = os.path.join(LAB, "fonts", name.replace(".ttf", ".woff2"))
      f.save(dst)
      print("wrote", os.path.basename(dst))

  for w in WEIGHTS:
    ttf = os.path.join(FONTDIR, "WaterlooGrotesk-%s.ttf" % w)
    if not os.path.exists(ttf):
      print("missing", ttf)
      continue
    data = glyph_data(ttf)
    out = os.path.join(LAB, "data", "glyphs-%s.json" % w)
    with open(out, "w") as fp:
      json.dump(data, fp, separators=(",", ":"))
    print("wrote %s (%d glyphs)" % (os.path.basename(out), len(data["glyphs"])))


if __name__ == "__main__":
  main()
