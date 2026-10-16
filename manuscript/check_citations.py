"""Cross-check author-year citations in the manuscript sections against
references.json. Usage: python3 manuscript/check_citations.py"""
import json, re, pathlib

here = pathlib.Path(__file__).parent
refs = json.loads((here / "references.json").read_text())
R = set()
for d in refs:
    a = d["author"][0]
    R.add((a.get("family", a.get("literal", "")).split()[-1], str(d["issued"]["date-parts"][0][0])))

body = ""
for f in ["abstract", "introduction", "methods_results", "discussion", "box2_design"]:
    body += re.sub(r"<!--.*?-->", "", (here / f"{f}.md").read_text(), flags=re.S)

name = r"([A-Z][A-Za-zÀ-ÿ'\-]+)(?: et al\.| and [A-Z][A-Za-zÀ-ÿ'\-]+)?"
C = {(m[1], m[2]) for m in re.finditer(name + r",? \(?((?:19|20)\d{2})", body)}
C |= {(m[1], m[3]) for m in re.finditer(name + r" \(?(\d{4}), (\d{4})", body)}

print("Cited but missing from references.json:", sorted(C - R) or "none")
print("In references.json but not cited:", sorted(R - C) or "none")
