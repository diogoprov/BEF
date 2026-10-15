# Phylogenetic diversity and ecosystem functioning in freshwater mesocosms

Research compendium for a PhD thesis chapter (D. B. Provete) on the effects of
phylogenetic diversity (PD, MPD, MNTD) on ecosystem functioning (productivity,
respiration, decomposition) in multitrophic pond mesocosms, using data from the
experiment by A. L. Downing.

## Design (from `data/raw/nature2001/data_diogo.csv`)

- 84 tanks; 21 species compositions with 4 tank replicates each
- 3 species-richness levels (`div` = 1, 9, 15) × 2 levels of a second factor (`treat` 1–6)
- 6 repeated sampling occasions (`time`) per tank (504 rows)
- Each composition has a unique PD/MPD value, so phylogenetic metrics are
  composition-level predictors

## Layout

| Folder | Content |
|---|---|
| `data/raw/` | Original data and phylogenies. **Read-only**: never edit in place |
| `data/derived/` | Data produced by scripts (to be created) |
| `R/` | Reusable functions |
| `analysis/` | Analysis scripts; `analysis/legacy/` holds the 2013–2019 scripts |
| `manuscript/` | Manuscript versions, tables, and the thesis chapter |
| `figures/` | Figures; `figures/legacy/` holds old outputs |
| `docs/` | Notes on the statistical approach, SAS code, and experimental design |
| `refs/` | Third-party PDFs (git-ignored; local only) |

## History

This repository merges `diogoprov/BEF` (code and data) with
`diogoprov/Capitulo-4---B-EF` (manuscript and SAS files; now archived).
