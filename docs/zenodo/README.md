# README

## Data and code for: Whole-community phylogenetic diversity restates trophic structure in experimental pond food webs

**Note for reviewers:** *Oikos* operates under double-anonymous peer review. Author names, affiliations, funding agencies and acknowledgements have been removed from this deposit and will be added to the metadata after acceptance. In the two original spreadsheets, the names of the people who wrote cell comments were replaced by neutral labels; the comments themselves are unchanged.

## GENERAL INFORMATION

### 1. Paper citation

==to be fulfilled upon acceptance== (DOI pending).

### 2. Brief abstract

We reanalysed two pond mesocosm experiments in which the richness and composition of rooted macrophytes, periphyton grazers and invertebrate predators were manipulated simultaneously (T98: 21 species compositions in 84 tanks, Downing and Leibold 2002; T99: 26 compositions in 104 tanks, Downing 2005), and in which productivity, respiration and decomposition were measured repeatedly and half of the tanks received an acid pulse. We computed Faith's phylogenetic diversity (PD), mean pairwise distance (MPD) and the variance of pairwise distances (VPD) of the stocked and realized compositions on pseudochronograms dated with TimeTree 5 (Kumar et al. 2022), and estimated their effects with Bayesian mixed models (brms) that accounted for the design, the repeated measures and the acid pulse. We also fitted a model of phylogenetically structured species effects and simulated alternative experimental designs. This deposit contains the original spreadsheets of both experiments, the phylogenies and node-age calibrations, the derived data needed to refit every model, the R scripts of the whole pipeline, the posterior draws used in the figures, all output tables, the figures and an HTML report of the analyses.

### 3. Originators

==to be fulfilled upon acceptance==

### 4. Contact information

==to be fulfilled upon acceptance==

### 5. Date of data collection

No new field data were collected. The experiments were run outdoors in 1998 (T98) and 1999 (T99); field and laboratory methods are described in Downing and Leibold (2002, 2010) and Downing (2005). Node ages were obtained from TimeTree 5 (timetree.org) on 16 October 2026. Derived data, models and figures were produced between 15 and 18 October 2026.

### 6. Geographic location(s)

Kellogg Biological Station, Hickory Corners, Michigan, USA (outdoor 300-L mesocosms).

### 7. Funding

==to be fulfilled upon acceptance==

## ACCESS INFORMATION

**Licence:** data under Creative Commons Attribution 4.0 International (CC-BY 4.0); code under the MIT licence (file `LICENSE`).

**Data derived from other sources** (please cite the original sources when reusing the derived files):

- Ecosystem functioning, design and final biomass of both experiments: original spreadsheets of Downing and Leibold (2002, 2010) and Downing (2005).
- Stocked composition of T98: Table 1 of the Supplementary Information of Downing and Leibold (2002). **The Supplementary Information file is not redistributed.** Its table is included as `data/derived/si_table1_long.csv`; step 00a, which parses it, runs only if the file `Downing_Leibold_2002_Nature_SI.doc` (downloaded from the publisher) is saved in `refs/`.
- Stocked composition of T99: Appendix A of Downing (2005), transcribed in `data/raw/ecology2005/Ecology2005.txt`.
- Node ages: TimeTree 5 (Kumar et al. 2022), median pairwise divergence times (`data/raw/phylogeny/timetree_calibrations.csv`).

## DATA FILES AND VARIABLES

The folder structure mirrors the R project, so that the scripts run with paths built by `here::here()` (the empty file `.here` marks the project root). Tables are comma-separated UTF-8 CSV files, except where noted. Species names follow current nomenclature and are written as tip labels (`Genus_species`); some differ from the names used in the original publications (e.g., *Lithobates catesbeianus* for *Rana catesbeiana*, *Elodea nuttallii* for *E. occidentale*, *Potamogeton zosteriformis* for *P. compressus*). Phylogenetic distances are in million years (PD, MPD and MNTD in Myr, VPD in Myr²; PSV is unitless).

### 1. `data/raw/` (original data, read-only)

- `nature2001/Tanks 98 (Nature 2001 paper).xlsx`: original T98 spreadsheet. Sheets `ecosystem functioning T98 wide` and `... long` (one row per tank × sampling period: `tank`, `treat`, `comp` (composition, 1–21), `div` (richness level: 1, 9 and 15 code 3, 9 and 15 stocked species), `dist` (acid pulse: 1, pulsed; 0, control), `time` (sampling period, 1–6), `dc` (decomposition), `ph` (phytoplankton, chlorophyll *a*), `prbio` (periphyton biomass), `prod` (productivity, diurnal oxygen), `resp` (respiration, diurnal oxygen) and other variables not used here), `species biomass T98` (final biomass per tank) and `codes in spreadsheets` (definitions of the codes, written by the originators).
- `nature2001/treat.txt`: tab-delimited species × composition matrix (presence/absence) of the realized (final) T98 compositions.
- `nature2001/data_legacy_2014.csv`: data set used in earlier, unpublished analyses (2014); read only by step 00 to document the corrections listed in the paper.
- `ecology2005/Tanks 99 Functional groups Ecology.xlsx`: original T99 spreadsheet. Sheets `FG expt design T99` (tank, composition, acid pulse, total richness `spdv`, number of functional groups `fgdv`, presence of macrophytes `mcr`, grazers `grz` and predators `prd`, and stocked species codes), `species biomasses T99`, `ecosystem function wide/long T99` (productivity `prd_corr`, respiration `rsp_corr`, decomposition `dec_A1`, by tank × sampling period, with `days` since the first sampling), `ecosystem biomass wide/long T99` and `spreadsheet codes`.
- `ecology2005/Ecology2005.txt`: tab-delimited species × composition matrix of the stocked T99 compositions, transcribed from Appendix A of Downing (2005); used to check the species codes of the design sheet.
- `phylogeny/phylog_dated_2002.txt`, `phylog_dated_2002_21spp_4b4da4d.txt`, `phylog_dated_2005.txt`: Newick pseudochronograms of 2013 (T98 with 19 and 21 species; T99 with 24 species), dated in Phylocom; used as topologies by step 00b.
- `phylogeny/timetree_calibrations.csv`: node calibrations. `node` (label), `tip_a`, `tip_b` (two tips whose most recent common ancestor is the calibrated node), `age_mya` (TimeTree median age, Myr), `ci_low`, `ci_high` (TimeTree confidence interval), `n_studies`, `taxon_a_queried`, `taxon_b_queried` (names queried in TimeTree), `note`.

### 2. `data/derived/` (written by the scripts; included so that any step can be rerun alone)

- `phylo_timetree_t98_19spp.tre`, `phylo_timetree_t98_21spp.tre`, `phylo_timetree_t99.tre`: Newick trees re-dated with TimeTree ages (step 00b).
- `si_table1_long.csv`: T98 stocked and final composition (step 00a). `comp`, `code`, `species`, `group` (macrophyte, grazer, predator), `cell` (symbol in the original table), `status` (`stocked_survived`, `stocked_extinct`, `invader`, `not_stocked`).
- `comm_initial.csv`, `t99_comm_initial.csv`: composition × species matrices (0/1) of the stocked compositions of T98 and T99.
- `metrics_initial.csv`, `metrics_composition.csv`, `metrics_tank.csv` (T98) and `t99_metrics_initial.csv`, `t99_metrics_final.csv` (T99): phylogenetic metrics of the stocked (`initial`), realized (`composition`, `final`) and tank-level compositions. `comp` or `tank`, `SR` (species richness), `PD` (Faith's PD, including the root), `MPD`, `MNTD`, `VPD`, `PSV`.
- `bef_long.csv`: T98 analysis data, one row per tank × sampling period (504 rows). Columns of the original sheet (see above), the metrics of the realized (`*_comp`), stocked (`*_init`) and tank-level (`*_tank`) compositions, `post` (1 from period 3, after the pulse) and `dist_post` (`dist × post`).
- `t99_long.csv`: T99 analysis data, one row per tank × sampling period (1,404 rows; composition 0, without stocked species, is excluded from the analyses). `tank`, `time` (1–13), `days`, `prod`, `resp`, `dc`, `comp`, `dist`, `spdv`, `fgdv`, `fgcomp` (functional groups present: M, G, MG, GP, MGP), `mcr`, `grz`, `prd`, metrics of the stocked (`*_init`) and realized (`*_final`) compositions, `cell` (design cell, `fgcomp_spdv`), `post` (1 from period 5), `early` (periods 5–6, resistance phase) and `late` (periods 7–13, recovery phase).

### 3. `output/tables/`

Model estimates, diagnostics (R-hat, bulk and tail effective sample sizes, divergent transitions, tree-depth hits) and checks, named by step: `pd_checks`, `predictor_correlations` (step 00); `phylo_redating_*` (00b: old vs new node ages and metric correlations); `fixef_*`, `pd_slopes`, `ar1_check`, `diagnostics` (01); `step2_*` (02: slopes within richness levels, variance within strata, among-composition variance explained); `step3_*` (03: acid pulse); `t99_*` (10–13 and 16: within-cell variance, branch-length sensitivity, turnover, slopes, design cells vs metrics, pulse, leave-one-composition-out influence, species-effects model, variance-structure sensitivity); `within_trophic_*` (17); `design_sim_*` (18–19: design simulation); `prior_predictive_summary`, `prior_sensitivity` (14); `ppc_grouped_pvalues` (15). Columns `q2.5`/`q97.5` and `q0.83`/`q99.17` are the limits of the 95% and 98.33% credible intervals, and `p_gt0` is the posterior probability that a parameter is positive.

### 4. `output/draws/`

Posterior draws of the slopes (`slope_draws_t98.rds`, `slope_draws_t99.rds`) and of the acid-pulse effects (`pulse_draws_t98.rds`, `pulse_draws_t99.rds`) used in Figs 3 and 5 (R `.rds` files; written by steps 02, 03 and 11).

### 5. `figures/`

`manuscript/` (Figs 1–5 and B1, PDF and 600-dpi PNG), `box_design_simulation.*` (Fig. B2), `ppc/` and `ppc_detail/` (posterior predictive checks of every model) and `prior_predictive/` (prior predictive checks).

### 6. `output/reports/report.html`

Self-contained HTML report of all analyses, written following the Bayesian Analysis Reporting Guidelines (priors, convergence, posterior predictive checks, prior sensitivity and estimates), rendered with Quarto from `analysis/report.qmd`.

## CODE SCRIPTS AND WORKFLOW

Shared functions are in `R/` and the steps in `analysis/`; `run_all.R` runs every step in order from the project root (`Rscript run_all.R`). Each step reads from `data/` and writes to `data/derived/`, `output/` or `figures/`; seeds are set for every stochastic step. Fitted models (`output/models/*.rds`) are not included because of their size; they are recreated by the scripts (about 3–4 h on four or more cores; later runs refit only models whose data or code changed).

| Step | Script | Input | Output |
|---|---|---|---|
| 00a | `00a_initial_composition.py` | Supplementary Information of Downing and Leibold (2002), converted to HTML with macOS `textutil` (not redistributed; optional) | `si_table1_long.csv` |
| 00b | `00b_redate_phylogenies.R` | 2013 trees, `timetree_calibrations.csv` | `phylo_timetree_*.tre`, `phylo_redating_*.csv` |
| 00 | `00_build_data.R` | T98 spreadsheet, `treat.txt`, `si_table1_long.csv`, trees | `bef_long.csv`, T98 metrics |
| 01 | `01_models.R` | 00 | metrics of the realized composition, AR(1) check, multivariate model |
| 02 | `02_initial_composition.R` | 00 | slopes within richness levels (T98) |
| 03 | `03_disturbance.R` | 00 | acid-pulse effects and their modulation (T98) |
| 10 | `10_build_data_t99.R` | T99 spreadsheet, tree | `t99_long.csv`, T99 metrics |
| 11 | `11_models_t99.R` | 10 | slopes within design cells, design cells vs metrics, acid pulse (T99) |
| 12 | `12_influence_t99.R` | 11 | leave-one-composition-out refits (about 30 min) |
| 13 | `13_species_identity_t99.R` | 10 | phylogenetically structured species effects (Box 1) |
| 17 | `17_within_trophic_check.R` | 00, 10 | MPD within functional groups (descriptive) |
| 18 | `18_design_simulation.R` | T99 tree | design simulation (Box 2; set `FULL_SIM=1` for 1,000 replicates) |
| 19 | `19_design_simulation_brms.R` | 18 | Bayesian check of the simulation |
| 20 | `20_figure_design_box.R` | 18 | Fig. B2 |
| 21 | `21_manuscript_figures.R` | trees, `output/tables/`, `output/draws/` | Figs 1–5 and B1 |
| 14 | `14_prior_checks.R` | 02, 11 | prior predictive checks and prior sensitivity |
| 15 | `15_ppc_detail.R` | 02, 11 | grouped posterior predictive checks |
| 16 | `16_variance_sensitivity_t99.R` | 11 | alternative variance structures (T99) |

Stan models use four chains of 2,000 iterations (1,000 warm-up; 4,000 for models refitted after failing the convergence criteria) run in parallel.

## SOFTWARE VERSIONS

R 4.6.1 (R Core Team 2026); CmdStan 2.39.0 through cmdstanr 0.9.0.9002 (Gabry et al. 2025); Quarto 1.8.26; Python 3 (standard library only; step 00a). Package versions are recorded in `renv.lock` (phylocomr and ggtree are missing from it). Main packages: brms 2.23.0 (Bürkner 2017), posterior 1.7.0, loo 2.10.1, priorsense 1.4.0 (Kallioinen et al. 2024), bayesplot 1.16.0, tidybayes 3.0.7, ggdist 3.3.3, ape 5.8-1 (Paradis and Schliep 2019), picante 1.8.4 (Kembel et al. 2010), phytools 2.5-2 (Revell 2024), phylocomr 0.3.4 (Ooms et al. 2023; runs bladj of Phylocom, Webb et al. 2008), ggtree 4.2.0, tidyverse 2.0.0 (dplyr 1.2.1, ggplot2 4.0.3), readxl 1.5.0.1, patchwork 1.3.2, ragg 1.5.2, withr 3.0.3, here 1.0.2.

## REFERENCES

- Bürkner, P.-C. 2017. brms: an R package for Bayesian multilevel models using Stan. – J. Stat. Softw. 80: 1–28. https://doi.org/10.18637/jss.v080.i01
- Downing, A. L. 2005. Relative effects of species composition and richness on ecosystem properties in ponds. – Ecology 86: 701–715. https://doi.org/10.1890/03-0239
- Downing, A. L. and Leibold, M. A. 2002. Ecosystem consequences of species richness and composition in pond food webs. – Nature 416: 837–841. https://doi.org/10.1038/416837a
- Downing, A. L. and Leibold, M. A. 2010. Species richness facilitates ecosystem resilience in aquatic food webs. – Freshwater Biol. 55: 2123–2137. https://doi.org/10.1111/j.1365-2427.2010.02472.x
- Gabry, J., Češnovar, R., Johnson, A. and Bronder, S. 2025. cmdstanr: R interface to "CmdStan". – R package version 0.9.0.9002. https://github.com/stan-dev/cmdstanr
- Kallioinen, N., Paananen, T., Bürkner, P.-C. and Vehtari, A. 2024. Detecting and diagnosing prior and likelihood sensitivity with power-scaling. – Stat. Comput. 34: 57. https://doi.org/10.1007/s11222-023-10366-5
- Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H., Ackerly, D. D., Blomberg, S. P. and Webb, C. O. 2010. Picante: R tools for integrating phylogenies and ecology. – Bioinformatics 26: 1463–1464. https://doi.org/10.1093/bioinformatics/btq166
- Kumar, S., Suleski, M., Craig, J. M., Kasprowicz, A. E., Sanderford, M., Li, M., Stecher, G. and Hedges, S. B. 2022. TimeTree 5: an expanded resource for species divergence times. – Mol. Biol. Evol. 39: msac174. https://doi.org/10.1093/molbev/msac174
- Ooms, J., Chamberlain, S. and Sanchez Reyes, L. 2023. phylocomr: interface to "Phylocom". – R package version 0.3.4. https://doi.org/10.32614/CRAN.package.phylocomr
- Paradis, E. and Schliep, K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. – Bioinformatics 35: 526–528. https://doi.org/10.1093/bioinformatics/bty633
- R Core Team 2026. R: a language and environment for statistical computing. – R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
- Revell, L. J. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). – PeerJ 12: e16505. https://doi.org/10.7717/peerj.16505
- Webb, C. O., Ackerly, D. D. and Kembel, S. W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. – Bioinformatics 24: 2098–2100. https://doi.org/10.1093/bioinformatics/btn358
