#!/bin/zsh
# Render the submission files to docx (manuscript/), add line and page numbers
# to all except the cover letter.
#   manuscript.docx                   anonymized manuscript, figures at the end
#   manuscript_figures_inline.docx    same, figures after their first citation
#   appendix_s1.docx                  supporting information (anonymized)
#   title_page.docx                   authors, acknowledgements, funding
#   cover_letter_oikos.docx           cover letter and suggested reviewers
set -euo pipefail
cd "${0:A:h}/.."
python3 manuscript/check_citations.py
python3 -I tools/inline_figures.py
for f in manuscript.qmd manuscript_figures_inline.qmd appendix_s1.qmd title_page.qmd cover_letter_oikos.md; do
  quarto render manuscript/$f --to docx --quiet
done
python3 -I tools/docx_line_page_numbers.py manuscript/{manuscript,manuscript_figures_inline,appendix_s1,title_page}.docx
