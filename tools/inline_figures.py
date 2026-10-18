"""Write a copy of manuscript.qmd with each figure and its legend placed right
after the first paragraph that cites it (for reading and co-author review;
the submission file keeps the figures at the end).

Usage (from the project root):
    python3 -I tools/inline_figures.py
    quarto render manuscript/manuscript_figures_inline.qmd --to docx
"""
import pathlib
import re

MS = pathlib.Path("manuscript")
qmd = (MS / "manuscript.qmd").read_text("utf-8")

# Expand the {{< include >}} shortcodes
text = re.sub(r"\{\{< include (\S+) >\}\}", lambda m: (MS / m.group(1)).read_text("utf-8"), qmd)

# Collect figure images and legends ("**Figure X.** ..."), then remove them
paras = re.split(r"\n\s*\n", text)
images, legends, kept = [], {}, []
for p in paras:
    s = p.strip()
    img = re.fullmatch(r"!\[\]\((\S+)\)", s)
    leg = re.match(r"\*\*Figure (\w+)\.\*\*", s)
    if img:
        images.append(img.group(1))
    elif leg:
        legends[leg.group(1)] = s
    elif re.fullmatch(r"#+ Figure (legends|B\d legend)", s) or s.startswith("Figure B2: legend in Box 2"):
        continue
    else:
        kept.append(p)

fig_ids = {"fig1": "1", "fig2": "2", "fig3": "3", "fig4": "4", "fig5": "5",
           "figB1": "B1", "box_design": "B2"}
image_of = {}
for path in images:
    key = next(k for k in fig_ids if pathlib.Path(path).name.startswith(k))
    image_of[fig_ids[key]] = path
assert set(image_of) == set(legends), (sorted(image_of), sorted(legends))

# Insert each figure after the first paragraph (outside the YAML header and
# HTML comments) that cites it
def cites(p, k):
    return re.search(rf"Figs?\.\s[^)]*?(?<![\w.]){re.escape(k)}(?![0-9])", p) is not None

body_start = next(i for i, p in enumerate(kept) if p.lstrip().startswith("#") and not p.lstrip().startswith("---"))
out = list(kept)
for k in sorted(legends, key=lambda x: (x.startswith("B"), x)):
    i = next(i for i, p in enumerate(out) if i > body_start and not p.lstrip().startswith("<!--")
             and not p.lstrip().startswith("!") and cites(p, k))
    out.insert(i + 1, f"![]({image_of[k]}){{width=16cm}}\n\n{legends[k]}")

(MS / "manuscript_figures_inline.qmd").write_text("\n\n".join(out) + "\n", "utf-8")
print("figures placed:", ", ".join(sorted(legends)))
