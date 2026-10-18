#!/bin/zsh
# Build the anonymized Zenodo deposit (double-anonymous review) in zenodo/:
#   zenodo/zenodo_deposit.zip  (files to upload)
#   zenodo/README.html         (paste into the Zenodo "Description" field)
# Copies only the files that the pipeline reads or writes, renames the legacy
# data file, removes personal names from the spreadsheet metadata and comment
# authors, and renders the analysis report without author names.
# Usage: zsh tools/build_zenodo_deposit.sh   (from the project root)
set -euo pipefail
cd "${0:A:h}/.."
OUT=zenodo; STAGE=$OUT/deposit
rm -rf $OUT; mkdir -p $STAGE

copy() { for f in "$@"; do mkdir -p "$STAGE/${f:h}"; cp -p "$f" "$STAGE/$f"; done; }

copy run_all.R renv.lock R/brms_helpers.R R/phylo_metrics.R R/t99_data.R
copy analysis/0*.R analysis/00a_initial_composition.py analysis/1*.R analysis/2*.R analysis/report.qmd
copy "data/raw/nature2001/Tanks 98 (Nature 2001 paper).xlsx" data/raw/nature2001/treat.txt \
     "data/raw/ecology2005/Tanks 99 Functional groups Ecology.xlsx" data/raw/ecology2005/Ecology2005.txt \
     data/raw/phylogeny/phylog_dated_2002.txt data/raw/phylogeny/phylog_dated_2002_21spp_4b4da4d.txt \
     data/raw/phylogeny/phylog_dated_2005.txt data/raw/phylogeny/timetree_calibrations.csv
copy data/derived/*(.) output/tables/**/*(.) output/draws/*(.)
copy figures/manuscript/*(.) figures/box_design_simulation.* figures/ppc/*(.) figures/ppc_detail/*(.) figures/prior_predictive/*(.)
cp docs/zenodo/README.md docs/zenodo/LICENSE $STAGE/
touch $STAGE/.here

# Legacy data file: neutral name
cp -p data/raw/nature2001/data_diogo.csv $STAGE/data/raw/nature2001/data_legacy_2014.csv
sed -i '' 's/data_diogo\.csv/data_legacy_2014.csv/g' $STAGE/analysis/00_build_data.R $STAGE/analysis/report.qmd
sed -i '' "s/from Amy Downing's original spreadsheet/from the original spreadsheet/" $STAGE/analysis/00_build_data.R
sed -i '' "s/from Amy Downing's original spreadsheets/from the original spreadsheets/" $STAGE/analysis/report.qmd
sed -i '' '/^author: /d' $STAGE/analysis/report.qmd

# Step 00a needs the publisher's SI file, which is not redistributed
python3 -I - $STAGE/run_all.R <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
a = 'run("analysis/00a_initial_composition.py")   # initial composition from the Nature SI\n'
b = ('if (file.exists(here("refs", "Downing_Leibold_2002_Nature_SI.doc"))) {\n'
     '  run("analysis/00a_initial_composition.py")   # initial composition from the Nature SI\n'
     '} else message("Nature SI not found in refs/: using data/derived/si_table1_long.csv")\n')
assert s.count(a) == 1; open(p, "w").write(s.replace(a, b))
PY

# Spreadsheets: personal names in document properties and comment authors
python3 -I tools/anonymize_xlsx.py $STAGE/data/raw/*/*.xlsx

# Report without author names (rendered from a temporary copy next to the original)
tmp=analysis/report_anonymous.qmd
sed -e '/^author: /d' -e "s/from Amy Downing's original spreadsheets/from the original spreadsheets/" \
    -e 's/data_diogo\.csv/data_legacy_2014.csv/g' analysis/report.qmd > $tmp
quarto render $tmp --quiet
mkdir -p $STAGE/output/reports; mv analysis/report_anonymous.html $STAGE/output/reports/report.html
rm -rf $tmp analysis/report_anonymous_files analysis/report_anonymous_cache

# Check: no personal names or local paths left (text files; embedded base64
# images in the report are skipped; xlsx parts are checked after unzipping)
python3 -I - $STAGE <<'PY'
import pathlib, re, sys, zipfile
names = re.compile(r"(?i)provete|diogo|\bamy\b|mathew|matt bruns|ohio wesleyan|ufms|provetelab")
paths = re.compile(r"/Users/")
bad = []
for f in pathlib.Path(sys.argv[1]).rglob("*"):
    if not f.is_file() or f.suffix in {".png", ".rds", ".pdf"}:
        continue
    if f.suffix == ".xlsx":
        with zipfile.ZipFile(f) as z:
            texts = [z.read(n).decode("utf-8", "ignore") for n in z.namelist() if n.endswith((".xml", ".vml"))]
    else:
        texts = [re.sub(r"data:[^\"')\s]+", "", f.read_text("utf-8", "ignore"))]
    for t in texts:
        for m in list(names.finditer(t)) + list(paths.finditer(t)):
            bad.append(f"{f}: ...{t[max(0, m.start() - 40):m.end() + 40]}...")
if bad:
    print("\n".join(bad[:20])); sys.exit("Identifying strings found")
print("anonymity check passed")
PY

(cd $STAGE && zip -q -r -X ../zenodo_deposit.zip . -x ".DS_Store")
quarto pandoc docs/zenodo/README.md -f gfm -t html -o $OUT/README.html
ls -la $OUT; du -sh $STAGE
