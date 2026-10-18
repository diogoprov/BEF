"""Replace personal names in .xlsx files (document properties and cell-comment
authors) by neutral labels, in place. Cell contents are not touched.

Usage: python3 -I tools/anonymize_xlsx.py file.xlsx [...]
"""
import re
import shutil
import sys
import tempfile
import zipfile

LABELS = [
    (r"Amy L\. Downing", "Data originator 1"),
    (r"Amy Downing", "Data originator 1"),
    (r"[Mm]att [Bb]runs", "Data originator 2"),
    (r"Diogo( B\.)? (Borges )?Provete", ""),
    (r"Ohio Wesleyan University", ""),
]

for path in sys.argv[1:]:
    with zipfile.ZipFile(path) as zin, tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
        with zipfile.ZipFile(tmp, "w") as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename.endswith((".xml", ".vml", ".rels")):
                    text = data.decode("utf-8")
                    for pat, rep in LABELS:
                        text = re.sub(pat, rep, text)
                    data = text.encode("utf-8")
                zout.writestr(item, data, compress_type=item.compress_type)
    shutil.move(tmp.name, path)
    print("anonymized", path)
