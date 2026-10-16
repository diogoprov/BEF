# Build analysis data from Amy Downing's original spreadsheet and recompute
# phylogenetic diversity metrics.
#
# Inputs (read-only): data/raw/nature2001/, data/raw/phylogeny/
# Requires: data/derived/si_table1_long.csv (analysis/00a_initial_composition.py,
#           Table 1 of the Downing & Leibold 2002 Supplementary Information)
# Outputs: data/derived/bef_long.csv, data/derived/metrics_composition.csv,
#          data/derived/metrics_initial.csv, data/derived/metrics_tank.csv,
#          data/derived/comm_initial.csv, output/tables/pd_checks.csv

library(tidyverse)
library(readxl)
library(ape)
library(here)

source(here("R", "phylo_metrics.R"))
dir.create(here("data", "derived"), showWarnings = FALSE)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)

xlsx <- here("data", "raw", "nature2001", "Tanks 98 (Nature 2001 paper).xlsx")

# ---- Ecosystem functioning (tank x sampling period) ----
ef <- read_excel(xlsx, sheet = "ecosystem functioning T98 long") |>
  select(tank, treat, comp, div, dist, time, dc, ph, prbio, prod, prpro, resp)

stopifnot(
  nrow(ef) == 504,
  n_distinct(ef$tank) == 84,
  all(count(ef, tank)$n == 6),
  all(count(distinct(ef, comp, tank), comp)$n == 4),
  # each tank has a single treatment/disturbance code
  all(count(distinct(ef, tank, treat, dist), tank)$n == 1)
)

# ---- Phylogeny: pseudochronogram without the two macrophytes that went
# extinct, re-dated with TimeTree ages (analysis/00b_redate_phylogenies.R) ----
tree <- read.tree(here("data", "derived", "phylo_timetree_t98_19spp.tre"))
stopifnot(is.ultrametric(tree, tol = 1e-4), Ntip(tree) == 19)

# ---- Composition-level community (final composition, D&L 2002) ----
comm_comp <- read.table(here("data", "raw", "nature2001", "treat.txt"),
                        header = TRUE, row.names = 1) |>
  as.matrix() |>
  t()
colnames(comm_comp) <- fix_species_names(colnames(comm_comp))

# ---- Tank-level community from final biomass ----
biomass <- read_excel(xlsx, sheet = "species biomass T98")
comm_tank <- biomass_presence(biomass)

# Validate species mapping: union of the 4 tanks == composition matrix
union_comp <- map(1:21, \(k) as.integer(colSums(comm_tank[biomass$comp == k, , drop = FALSE]) > 0)) |>
  do.call(what = rbind)
colnames(union_comp) <- colnames(comm_tank)
stopifnot(all(union_comp[, colnames(comm_comp)] == comm_comp))

# ---- Metrics ----
met_comp <- phylo_metrics(comm_comp, tree) |>
  mutate(comp = as.integer(str_remove(site, "^T")), .keep = "unused", .before = 1)

met_tank <- phylo_metrics(comm_tank, tree) |>
  mutate(tank = as.integer(site), .keep = "unused", .before = 1)

# ---- Initial (stocked) composition: SI Table 1 of Downing & Leibold 2002 ----
si <- read_csv(here("data", "derived", "si_table1_long.csv"), show_col_types = FALSE)
comm_init <- si |>
  mutate(stocked = as.integer(status %in% c("stocked_survived", "stocked_extinct"))) |>
  select(comp, species, stocked) |>
  pivot_wider(names_from = species, values_from = stocked) |>
  arrange(comp) |>
  column_to_rownames("comp") |>
  as.matrix()

# Design check: 1, 3 or 5 species per functional group (3, 9, 15 in total)
stopifnot(all(rowSums(comm_init) == rep(c(3, 9, 15), each = 7)))
# Final composition in the SI (survivors + invaders) == treat.txt
final_si <- si |>
  mutate(present = as.integer(status %in% c("stocked_survived", "invader"))) |>
  select(comp, species, present) |>
  pivot_wider(names_from = species, values_from = present) |>
  arrange(comp) |>
  column_to_rownames("comp") |>
  as.matrix()
stopifnot(all(final_si[, colnames(comm_comp)] == comm_comp),
          all(final_si[, setdiff(colnames(final_si), colnames(comm_comp))] == 0))

# The initial pool includes Utricularia and Vallisneria (extinct in all tanks),
# so it requires the 21-species pseudochronogram
tree21 <- read.tree(here("data", "derived", "phylo_timetree_t98_21spp.tre")) |>
  fix_tip_labels()
stopifnot(is.ultrametric(tree21, tol = 1e-4), setequal(tree21$tip.label, colnames(comm_init)))

met_init <- phylo_metrics(comm_init, tree21) |>
  mutate(comp = as.integer(site), .keep = "unused", .before = 1)
stopifnot(all(met_init$SR == rep(c(3, 9, 15), each = 7)))
write_csv(as_tibble(comm_init, rownames = "comp"), here("data", "derived", "comm_initial.csv"))
write_csv(met_init, here("data", "derived", "metrics_initial.csv"))

# Sensitivity to branch lengths: Grafen (topology only) and unit branch lengths
tree_grafen <- compute.brlen(tree, method = "Grafen")
tree_unit <- compute.brlen(tree, 1)
sens <- map(list(bladj = tree, grafen = tree_grafen, unit = tree_unit),
            \(tr) phylo_metrics(comm_comp, tr))
sens_cor <- map_dfr(c("PD", "MPD", "MNTD", "VPD", "PSV"), \(v) tibble(
  metric = v,
  r_grafen = cor(sens$bladj[[v]], sens$grafen[[v]]),
  r_unit = cor(sens$bladj[[v]], sens$unit[[v]])
))
print(sens_cor)

# ---- Compare with values in the legacy data_diogo.csv ----
legacy <- read.csv(here("data", "raw", "nature2001", "data_diogo.csv")) |>
  distinct(comp, PD_legacy = PD, PSV_legacy = PSV)
pd_checks <- met_comp |>
  left_join(legacy, by = "comp") |>
  mutate(dPD = PD - PD_legacy, dPSV = PSV - PSV_legacy)
write_csv(pd_checks, here("output", "tables", "pd_checks.csv"))
print(pd_checks, n = 21)

# ---- Assemble long data ----
bef <- ef |>
  left_join(met_comp |> rename_with(\(x) paste0(x, "_comp"), -comp), by = "comp") |>
  left_join(met_init |> rename_with(\(x) paste0(x, "_init"), -comp), by = "comp") |>
  left_join(met_tank |> rename_with(\(x) paste0(x, "_tank"), -tank), by = "tank") |>
  mutate(
    # acid pulse applied after sampling period 2
    post = as.integer(time >= 3),
    dist_post = dist * post
  ) |>
  arrange(comp, tank, time)

write_csv(bef, here("data", "derived", "bef_long.csv"))
write_csv(met_comp, here("data", "derived", "metrics_composition.csv"))
write_csv(met_tank, here("data", "derived", "metrics_tank.csv"))

# Correlations among composition-level predictors (21 compositions)
cor_tab <- met_init |>
  rename_with(\(x) paste0(x, "_init"), -comp) |>
  left_join(met_comp |> rename_with(\(x) paste0(x, "_final"), -comp), by = "comp") |>
  left_join(distinct(ef, comp, div), by = "comp") |>
  select(div, PD_init, MPD_init, MNTD_init, VPD_init, SR_final, PD_final, MPD_final, MNTD_final, VPD_final) |>
  cor(use = "pairwise.complete.obs") |>
  round(2)
print(cor_tab)
write_csv(as_tibble(cor_tab, rownames = "var"), here("output", "tables", "predictor_correlations.csv"))
