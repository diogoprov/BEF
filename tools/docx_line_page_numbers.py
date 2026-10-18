"""Add continuous line numbers, page numbers (centred footer) and A4 page with 2.5 cm
margins to docx files rendered by Quarto/pandoc (pandoc drops the section properties of
the reference document, so they are injected after rendering).

Usage (from the repo root, after quarto render):
    python3 -I tools/docx_line_page_numbers.py manuscript/manuscript.docx manuscript/title_page.docx
"""
import os
import re
import shutil
import sys
import tempfile
import zipfile

FOOTER = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>1</w:t></w:r>'
    '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>'
)
FOOTER_RID = "rIdPageNumberFooter"
FOOTER_REF = f'<w:footerReference w:type="default" r:id="{FOOTER_RID}"/>'
SECT_LAYOUT = (
    '<w:pgSz w:w="11906" w:h="16838"/>'
    '<w:pgMar w:top="1418" w:right="1418" w:bottom="1418" w:left="1418" '
    'w:header="709" w:footer="709" w:gutter="0"/>'
    '<w:lnNumType w:countBy="1" w:distance="283" w:restart="continuous"/>'
)


def process(path):
    with zipfile.ZipFile(path) as z:
        files = {n: z.read(n) for n in z.namelist()}

    doc = files["word/document.xml"].decode("utf-8")
    if "<w:lnNumType" in doc:
        print(f"{path}: already processed")
        return
    if 'xmlns:r="' not in doc[:3000]:
        doc = doc.replace("<w:document ", '<w:document xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" ', 1)
    # last (body-level) section properties; schema order: footerReference first,
    # then footnotePr (written by pandoc), then pgSz, pgMar, lnNumType
    idx = doc.rfind("<w:sectPr")
    end = doc.find(">", idx) + 1
    close = doc.find("</w:sectPr>", end)
    doc = doc[:end] + FOOTER_REF + doc[end:close] + SECT_LAYOUT + doc[close:]
    files["word/document.xml"] = doc.encode("utf-8")
    files["word/footer_pagenum.xml"] = FOOTER.encode("utf-8")

    rels = files["word/_rels/document.xml.rels"].decode("utf-8")
    rels = rels.replace(
        "</Relationships>",
        f'<Relationship Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
        f'Id="{FOOTER_RID}" Target="footer_pagenum.xml"/></Relationships>')
    files["word/_rels/document.xml.rels"] = rels.encode("utf-8")

    ct = files["[Content_Types].xml"].decode("utf-8")
    ct = ct.replace(
        "</Types>",
        '<Override PartName="/word/footer_pagenum.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>')
    files["[Content_Types].xml"] = ct.encode("utf-8")

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".docx").name
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for n, b in files.items():
            z.writestr(n, b)
    shutil.move(tmp, path)
    os.chmod(path, 0o644)
    print(f"{path}: line and page numbers added")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        process(p)
