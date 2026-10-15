"""Extract initial (stocked) composition from Table 1 of the Downing & Leibold
(2002, Nature) Supplementary Information.

Table legend: '-' = not stocked; '0' = stocked but went extinct; a value =
stocked and survived; italic value = not stocked, invaded by final census.

Usage (from the project root; the .doc is first converted with macOS textutil):
    textutil -convert html -output /tmp/si.html refs/Downing_Leibold_2002_Nature_SI.doc
    python3 -I analysis/00a_initial_composition.py /tmp/si.html data/derived/si_table1_long.csv
"""
import csv
import re
import sys
from html.parser import HTMLParser

SPECIES = {
    "A": "Potamogeton_zosteriformis",  # 'Potamageton compressus' in the SI
    "B": "Utricularia_vulgaris",
    "C": "Myriophyllum_verticullatum",
    "D": "Elodea_nuttallii",           # 'Elodea occidentale' in the SI
    "E": "Vallisneria_americana",
    "F": "Ceratophyllum_demersum",
    "G": "Potamogeton_crispus",
    "H": "Crangonyx_richmondensis",
    "I": "Lithobates_catesbeianus",    # Rana catesbeiana
    "J": "Lithobates_clamitans",       # Rana clamitans
    "K": "Physa_gyrina",
    "L": "Hyalella_azteca",
    "M": "Helisoma_trivolvis",
    "N": "Trichocorixa_sp",
    "O": "Neoplea_striola",
    "P": "Ambrysus_sp",
    "Q": "Gyrinus_sp",
    "R": "Acilius_sp",
    "S": "Belostoma_flumireum",
    "T": "Notonecta_undulata",
    "U": "Notonecta_sp",
}
FUNCTIONAL_GROUP = {**{k: "macrophyte" for k in "ABCDEFG"},
                    **{k: "grazer" for k in "HIJKLMN"},
                    **{k: "predator" for k in "OPQRSTU"}}


class RowParser(HTMLParser):
    """Collect paragraphs as lists of (cell_text, italic) split at tab spans."""

    def __init__(self):
        super().__init__()
        self.paras, self.cells = [], None
        self.italic = 0
        self.buf, self.buf_italic = "", False

    def _flush(self):
        if self.cells is not None:
            self.cells.append((self.buf.strip(), self.buf_italic))
        self.buf, self.buf_italic = "", False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "p":
            self.cells, self.buf, self.buf_italic = [], "", False
        elif tag == "i":
            self.italic += 1
        elif tag == "span" and a.get("class") == "Apple-tab-span":
            self._flush()

    def handle_endtag(self, tag):
        if tag == "i":
            self.italic -= 1
        elif tag == "p" and self.cells is not None:
            self._flush()
            self.paras.append(self.cells)
            self.cells = None

    def handle_data(self, data):
        if self.cells is None:
            return
        # whitespace also separates cells (some are space-separated); italics
        # are tracked per token, not per cell
        for ch in data:
            if ch.isspace():
                if self.buf:
                    self._flush()
            else:
                self.buf += ch
                self.buf_italic = self.buf_italic or self.italic > 0


def main(html_path, out_path):
    p = RowParser()
    p.feed(open(html_path, encoding="utf-8").read())

    rows, block = [], 0
    for cells in p.paras:
        first = cells[0][0] if cells else ""
        if re.fullmatch(r"\d+", first) or (first == "" and len(cells) > 1 and cells[1][0] == "1"):
            continue
        letter = first.strip()
        if letter in SPECIES:
            values = [c for c in cells[1:] if c[0] != ""]
            if letter == "A":
                block += 1
            comps = {1: range(1, 8), 2: range(8, 15), 3: range(15, 22)}[block]
            if len(values) > 7:  # SD line merged into the same paragraph
                print(f"note: row {letter}, block {block}: dropping trailing SDs {values[7:]}")
                values = values[:7]
            if len(values) != 7:
                sys.exit(f"Row {letter}, block {block}: expected 7 cells, got {len(values)}: {values}")
            for comp, (txt, ital) in zip(comps, values):
                if txt == "-":
                    status = "not_stocked"
                elif ital:
                    status = "invader"
                elif re.fullmatch(r"0", txt):
                    status = "stocked_extinct"
                else:
                    status = "stocked_survived"
                rows.append({"comp": comp, "code": letter, "species": SPECIES[letter],
                             "group": FUNCTIONAL_GROUP[letter], "cell": txt, "status": status})

    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows ({len({r['comp'] for r in rows})} compositions)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
