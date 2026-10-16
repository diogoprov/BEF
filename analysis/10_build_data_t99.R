# Build analysis data for the 1999 functional-group experiment (T99;
# Downing 2005, Ecology) and compute phylogenetic diversity metrics.
#
# Design: 27 compositions x 4 tanks (2 with acid pulse); composition 0 is a
# no-species control (excluded from the phylogenetic analyses). Species
# richness 6 or 18; functional-group (FG) composition: M, G, MG, GP, MGP.
# Thirteen sampling periods with irregular spacing (days 0-72); decomposition
# measured in nine of them.
#
# Inputs (read-only): data/raw/ecology2005/, data/derived/phylo_timetree_t99.tre
# Outputs: data/derived/t99_long.csv, data/derived/t99_metrics_initial.csv,
#          data/derived/t99_metrics_final.csv, output/tables/t99_*.csv

library(tidyverse)
library(readxl)
library(ape)
library(here)

source(here("R", "phylo_metrics.R"))

xlsx <- here("data", "raw", "ecology2005", "Tanks 99 Functional groups Ecology.xlsx")

# Species codes in "FG expt design T99" -> tip labels. Resolved by matching
# presence patterns across the 26 compositions against Ecology2005.txt (24 of
# 24 codes matched one species in >= 25 of 26 compositions). The two
# disagreements are errors in Ecology2005.txt: composition 5 lacks
# L. clamitans (it lists only 5 of 6 species) and composition 16 lists
# Buenoa instead of Dineutus; final biomass confirms the design sheet
# (L. clamitans in 2 tanks of comp. 5; Dineutus in 3 tanks, Buenoa in none,
# of comp. 16).
design_codes <- c(
  cer = "Ceratophyllum_demersum", elo = "Elodea_nuttallii", mir = "Myriophyllum_verticullatum",
  po1 = "Potamogeton_crispus", po2 = "Potamogeton_natans", sag = "Sagittaria_rigida",
  utr = "Utricularia_vulgaris", val = "Vallisneria_americana",
  bull = "Lithobates_catesbeianus", gre = "Lithobates_clamitans", gam = "Crangonyx_richmondensis",
  hya = "Hyalella_azteca", hel = "Helisoma_trivolvis", phy = "Physa_gyrina",
  cor = "Hesperocorixa_sp", co2 = "Trichocorixa_sp",
  nau = "Ambrysus_sp", bel = "Belostoma_flumineum", no1 = "Notonecta_undulata",
  no2 = "Buenoa_sp", ple = "Neoplea_striola", dit = "Coptotomus_sp",
  whi = "Gyrinus_sp", hyd = "Dineutus_sp"
)

# Final biomass columns -> tip labels (invaders not in the phylogeny, e.g.
# damselflies, leeches, mites, salamanders, are excluded from PD metrics)
biomass_codes <- c(
  cer_g = "Ceratophyllum_demersum", elo_g = "Elodea_nuttallii", mir_g = "Myriophyllum_verticullatum",
  pocriso_g = "Potamogeton_crispus", ponat_g = "Potamogeton_natans", sag_g = "Sagittaria_rigida",
  utr_g = "Utricularia_vulgaris", val_g = "Vallisneria_americana",
  bul_g = "Lithobates_catesbeianus", gre_g = "Lithobates_clamitans", gam_g = "Crangonyx_richmondensis",
  hya_g = "Hyalella_azteca", hel_g = "Helisoma_trivolvis", phy_g = "Physa_gyrina",
  cobhesp_g = "Hesperocorixa_sp", cosmtri_g = "Trichocorixa_sp",
  amb_nau_g = "Ambrysus_sp", bel_g = "Belostoma_flumineum", nund_g = "Notonecta_undulata",
  buen_g = "Buenoa_sp", ple_g = "Neoplea_striola", dit_g = "Coptotomus_sp",
  whb_gyr_g = "Gyrinus_sp", whs_din_g = "Dineutus_sp"
)

tree <- read.tree(here("data", "derived", "phylo_timetree_t99.tre")) |> fix_tip_labels()
stopifnot(is.ultrametric(tree, tol = 1e-4), setequal(tree$tip.label, design_codes),
          setequal(tree$tip.label, biomass_codes))

# ---- Design and initial composition ----
design <- read_excel(xlsx, sheet = "FG expt design T99") |>
  filter(!is.na(tank)) |>
  mutate(fgcomp = case_when(
    mcr == 0 & grz == 0 & prd == 0 ~ "none",
    TRUE ~ paste0(if_else(mcr == 1, "M", ""), if_else(grz == 1, "G", ""), if_else(prd == 1, "P", ""))
  ))
stopifnot(nrow(design) == 108, all(count(design, comp)$n == 4),
          all(count(distinct(design, tank, dist), tank)$n == 1))

init_long <- design |>
  distinct(comp, .keep_all = TRUE) |>
  select(comp, matches("^(mac|graz|pred)\\.\\.\\.")) |>
  pivot_longer(-comp, values_to = "code") |>
  filter(code != "0") |>
  mutate(species = design_codes[code])
stopifnot(!anyNA(init_long$species))

comm_init <- init_long |>
  filter(comp > 0) |>
  mutate(v = 1L) |>
  distinct(comp, species, v) |>
  pivot_wider(names_from = species, values_from = v, values_fill = 0L) |>
  arrange(comp) |>
  column_to_rownames("comp") |>
  as.matrix()
cells <- distinct(design, comp, spdv, fgdv, fgcomp) |> filter(comp > 0) |> arrange(comp)
stopifnot(all(rowSums(comm_init) == cells$spdv))

# ---- Final composition (union of species with final biomass > 0 over tanks) ----
biomass <- read_excel(xlsx, sheet = "species biomasses T99") |> filter(!is.na(tank))
pres_tank <- as.matrix(biomass[, names(biomass_codes)])
pres_tank[is.na(pres_tank)] <- 0
pres_tank <- (pres_tank > 0) * 1L
colnames(pres_tank) <- unname(biomass_codes)
comm_final <- map(cells$comp, \(k) as.integer(colSums(pres_tank[biomass$comp == k, , drop = FALSE]) > 0)) |>
  do.call(what = rbind)
dimnames(comm_final) <- list(cells$comp, colnames(pres_tank))

# Stocked species that persisted vs invaders (within the phylogeny)
sp <- colnames(comm_final)
turnover <- tibble(
  comp = cells$comp,
  stocked = rowSums(comm_init[, sp]),
  persisted = rowSums(comm_init[, sp] * comm_final),
  invaders = rowSums((1 - comm_init[, sp]) * comm_final)
)
write_csv(turnover, here("output", "tables", "t99_turnover.csv"))

# ---- Metrics ----
met_init <- phylo_metrics(comm_init, tree) |>
  mutate(comp = as.integer(site), .keep = "unused", .before = 1)
met_final <- phylo_metrics(comm_final[rowSums(comm_final) > 0, , drop = FALSE], tree) |>
  mutate(comp = as.integer(site), .keep = "unused", .before = 1)

# Sensitivity to branch lengths (initial composition)
sens <- map(list(bladj = tree, grafen = compute.brlen(tree, method = "Grafen"),
                 unit = compute.brlen(tree, 1)),
            \(tr) phylo_metrics(comm_init, tr))
sens_cor <- map_dfr(c("PD", "MPD", "MNTD", "VPD"), \(v) tibble(
  metric = v,
  r_grafen = cor(sens$bladj[[v]], sens$grafen[[v]], use = "complete.obs"),
  r_unit = cor(sens$bladj[[v]], sens$unit[[v]], use = "complete.obs")
))
write_csv(sens_cor, here("output", "tables", "t99_branch_length_sensitivity.csv"))
print(sens_cor)

# Share of each metric's variance within design cells (FG composition x richness)
within_share <- cells |>
  left_join(met_init |> select(comp, PD, MPD, MNTD, VPD) |> rename_with(\(x) paste0(x, "_init"), -comp), by = "comp") |>
  left_join(met_final |> select(comp, PD, MPD, MNTD, VPD) |> rename_with(\(x) paste0(x, "_final"), -comp), by = "comp") |>
  mutate(cell = paste(fgcomp, spdv, sep = "_")) |>
  pivot_longer(matches("_(init|final)$"), names_to = "metric") |>
  group_by(metric) |>
  summarise(r2_cell = summary(lm(value ~ cell))$r.squared,
            r2_spdv = summary(lm(value ~ factor(spdv)))$r.squared, .groups = "drop") |>
  mutate(within_cell_share = 1 - r2_cell)
write_csv(within_share, here("output", "tables", "t99_within_cell_variance.csv"))
print(within_share)

# ---- Ecosystem functioning (tank x sampling period) ----
ef <- read_excel(xlsx, sheet = "ecosystem function long T99") |>
  filter(!is.na(tank)) |>
  transmute(tank, time, days, prod = prd_corr, resp = rsp_corr, dc = dec_A1)
stopifnot(nrow(ef) == 108 * 13)

t99 <- ef |>
  left_join(design |> select(tank, comp, dist, spdv, fgdv, fgcomp, mcr, grz, prd), by = "tank") |>
  left_join(met_init |> rename_with(\(x) paste0(x, "_init"), -comp), by = "comp") |>
  left_join(met_final |> rename_with(\(x) paste0(x, "_final"), -comp), by = "comp") |>
  mutate(
    cell = paste(fgcomp, spdv, sep = "_"),
    # acid pulse between periods 4 (day 19) and 5 (day 22): inferred from
    # treatment means, which diverge only from period 5 onwards
    post = as.integer(time >= 5),
    early = as.integer(time %in% 5:6),
    late = as.integer(time >= 7)
  ) |>
  arrange(comp, tank, time)

write_csv(t99, here("data", "derived", "t99_long.csv"))
write_csv(met_init, here("data", "derived", "t99_metrics_initial.csv"))
write_csv(met_final, here("data", "derived", "t99_metrics_final.csv"))
write_csv(as_tibble(comm_init, rownames = "comp"), here("data", "derived", "t99_comm_initial.csv"))

print(turnover, n = Inf)
