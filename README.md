# Phylogenetic diversity and ecosystem functioning in freshwater mesocosms

Research compendium for a PhD thesis chapter (D. B. Provete) on the effects of
phylogenetic diversity (PD, MPD, MNTD) on ecosystem functioning (productivity,
respiration, decomposition) in multitrophic pond mesocosms, using data from
the experiment by A. L. Downing (Downing & Leibold 2002, *Nature*).

## Experimental design

Derived from the original spreadsheet
(`data/raw/nature2001/Tanks 98 (Nature 2001 paper).xlsx`):

- 84 tanks = 21 species compositions × 4 replicate tanks.
- Richness treatment among compositions: `div` = 1, 9, 15 codes compositions
  with 3, 9, and 15 initial species (`INIT_SP`); 7 compositions per level.
- Acid pulse disturbance (`dist` = 0/1), crossed with composition (2 tanks of
  each composition per level), applied after sampling period 2.
- Six sampling periods per tank (`time`); 504 observations.
- Responses: productivity and respiration (diurnal oxygen), decomposition (`dc`).

Phylogenetic metrics are computed for the **initial** (stocked) composition,
parsed from Table 1 of the Downing & Leibold (2002) Supplementary Information,
and for the **final** composition (`treat.txt`), which equals the union of the
species surviving in the four tanks of each composition. Each composition has
a unique value, so these metrics are composition-level predictors. The second
experiment (T99; Downing 2005, Ecology) is in `data/raw/ecology2005/`.

## Workflow

```r
source("run_all.R")   # full pipeline (~3-4 h from scratch) and report
```

`run_all.R` runs, in order: `00a`/`00` (T98 data and metrics), `01`-`03`
(T98 models), `10`-`13` (T99 data, models, influence analysis and the
exploratory species-identity Box), `14` (prior predictive and prior
sensitivity checks), and renders `analysis/report.qmd` to a self-contained
HTML report (BARG-style; not versioned).

| Script | Output |
|---|---|
| `analysis/00_build_data.R` | `data/derived/bef_long.csv`, `metrics_composition.csv`, `metrics_tank.csv`; `output/tables/pd_checks.csv` |
| `analysis/01_models.R` | `output/tables/*.csv`, `figures/ppc/*.png`; fitted models in `output/models/` (git-ignored) |

Software: R 4.6.1, brms 2.23.0 with CmdStan 2.39.0, picante, ape.

## Data corrections relative to the 2014 analyses

1. Faith's PD was recomputed with the re-dated pseudochronogram
   (`phylog_dated_2002.txt`, commit `c47b468`). The 2014 values came from the
   earlier tree with *Utricularia* and *Vallisneria* pruned via `drop.tip()`,
   which underestimated PD by 21.1 Myr in nine compositions.
2. PSV of composition 6 was 0.620 (copied from composition 5); the correct
   value is 0.762.
3. `dist` for tank 5, period 6 was 6 in `data_diogo.csv`; the original
   spreadsheet has 1.
4. The legacy scripts assigned MPD with `rep(mpd, each = 12)`, but each
   composition has 24 rows (4 tanks × 6 periods). R recycled the vector
   silently, so MPD and MNTD were wrong in 480 of 504 rows in all legacy
   MPD/MNTD analyses (including manuscript Tables 2, 4 and 6). PD analyses used
   the `PD` column of the data file and are not affected. The new pipeline
   joins metrics by composition.

## Layout

| Folder | Content |
|---|---|
| `data/raw/` | Original data and phylogenies. **Read-only** |
| `data/derived/` | Data produced by `analysis/00_build_data.R` |
| `R/` | Functions (`phylo_metrics.R`) |
| `analysis/` | Analysis scripts; `analysis/legacy/` holds the 2013–2019 scripts |
| `output/` | Tables and logs; fitted models are git-ignored |
| `manuscript/` | Manuscript versions, tables, and the thesis chapter |
| `figures/` | Figures; `figures/legacy/` holds old outputs |
| `docs/` | Notes on the statistical approach, SAS code, experimental design |
| `refs/` | Third-party PDFs (git-ignored; local only) |

## History

This repository merges `diogoprov/BEF` (code and data) with
`diogoprov/Capitulo-4---B-EF` (manuscript and SAS files; archived on GitHub).
