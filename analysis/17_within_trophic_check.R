# Within-trophic-group phylogenetic diversity (T98 and T99), descriptive only.
#
# Whole-community metrics are nearly determined by richness and functional-group
# composition because plant-animal divergences dominate the distances. This
# script asks whether MPD computed *within* each trophic group (macrophytes,
# grazers, predators; stocked composition) varies among compositions within
# design cells, and which phylogenetic splits generate that variation.
#
# Inputs (read-only): data/derived/comm_initial.csv, si_table1_long.csv,
#   t99_comm_initial.csv; data/raw/phylogeny/; T99 design sheet.
# Outputs: output/tables/within_trophic_{group_summary,composition,lineages}.csv

suppressPackageStartupMessages({
  library(tidyverse); library(ape); library(readxl); library(here)
})
source(here("R", "phylo_metrics.R"))

# MPD and PD (root excluded) of each trophic group in each composition
group_metrics <- function(comm, tree, groups) {
  D <- cophenetic(tree)
  map_dfr(rownames(comm), \(k) {
    sp <- colnames(comm)[comm[k, ] > 0]
    map_dfr(split(sp, groups[sp]), \(s) tibble(
      n   = length(s),
      mpd = if (length(s) > 1) mean(D[s, s][lower.tri(D[s, s])]) else NA_real_,
      pd  = if (length(s) > 1) sum(keep.tip(tree, s)$edge.length) else NA_real_
    ), .id = "group") |> mutate(comp = k)
  })
}

# Fraction of the variance of x that lies within levels of cell
within_share <- function(x, cell) {
  ok <- !is.na(x); x <- x[ok]; cell <- cell[ok]
  sum((x - ave(x, cell))^2) / sum((x - mean(x))^2)
}

# Counts of the taxa that sit on the deep splits within each group
lineages <- function(comm) {
  has <- \(sp) rowSums(comm[, intersect(sp, colnames(comm)), drop = FALSE])
  tibble(
    comp = rownames(comm),
    n_frog    = has(c("Lithobates_catesbeianus", "Lithobates_clamitans")),
    n_snail   = has(c("Physa_gyrina", "Helisoma_trivolvis")),
    n_beetle  = has(c("Gyrinus_sp", "Acilius_sp", "Coptotomus_sp", "Dineutus_sp")),
    n_monocot = has(c("Potamogeton_zosteriformis", "Potamogeton_crispus", "Potamogeton_natans",
                      "Elodea_nuttallii", "Vallisneria_americana", "Sagittaria_rigida"))
  )
}

# ---- T98 ----
tree98 <- read.tree(here("data", "derived", "phylo_timetree_t98_21spp.tre")) |>
  fix_tip_labels()
grp98 <- read_csv(here("data", "derived", "si_table1_long.csv"), show_col_types = FALSE) |>
  distinct(species, group) |> deframe()
comm98 <- read_csv(here("data", "derived", "comm_initial.csv"), show_col_types = FALSE) |>
  column_to_rownames("comp") |> as.matrix()
stopifnot(setequal(colnames(comm98), tree98$tip.label), all(colnames(comm98) %in% names(grp98)))
cell98 <- tibble(comp = rownames(comm98), cell = paste0("S", rowSums(comm98)))

g98 <- group_metrics(comm98, tree98, grp98) |> left_join(cell98, by = "comp") |>
  mutate(experiment = "T98")

# ---- T99 ----
tree99 <- read.tree(here("data", "derived", "phylo_timetree_t99.tre")) |> fix_tip_labels()
codes99 <- c(
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
design99 <- suppressMessages(read_excel(
  here("data", "raw", "ecology2005", "Tanks 99 Functional groups Ecology.xlsx"),
  sheet = "FG expt design T99")) |>
  filter(!is.na(tank), comp > 0)
grp99 <- design99 |>
  distinct(comp, .keep_all = TRUE) |>
  select(matches("^(mac|graz|pred)\\.\\.\\.")) |>
  pivot_longer(everything()) |>
  filter(value != "0") |>
  mutate(group = recode(sub("\\.\\.\\..*", "", name),
                        mac = "macrophyte", graz = "grazer", pred = "predator"),
         species = codes99[value]) |>
  distinct(species, group) |> deframe()
stopifnot(length(grp99) == 24, setequal(names(grp99), tree99$tip.label))
comm99 <- read_csv(here("data", "derived", "t99_comm_initial.csv"), show_col_types = FALSE) |>
  column_to_rownames("comp") |> as.matrix()
cell99 <- design99 |>
  distinct(comp, mcr, grz, prd, spdv) |>
  transmute(comp = as.character(comp),
            cell = paste0(if_else(mcr == 1, "M", ""), if_else(grz == 1, "G", ""),
                          if_else(prd == 1, "P", ""), "_", spdv))

g99 <- group_metrics(comm99, tree99, grp99) |> left_join(cell99, by = "comp") |>
  mutate(experiment = "T99")

# ---- Summaries ----
g <- bind_rows(g98, g99)

group_summary <- g |>
  filter(n > 1) |>
  group_by(experiment, group, n) |>
  summarise(n_comp = n(), n_distinct = n_distinct(round(mpd, 1)),
            mpd_mean = mean(mpd), mpd_sd = sd(mpd), mpd_min = min(mpd), mpd_max = max(mpd),
            cv = sd(mpd) / mean(mpd), .groups = "drop")

# Mean within-group MPD per composition: community-level metric that excludes
# cross-trophic distances
composition <- g |>
  group_by(experiment, comp, cell) |>
  summarise(mpd_within_groups = mean(mpd, na.rm = TRUE), .groups = "drop") |>
  mutate(mpd_within_groups = na_if(mpd_within_groups, NaN)) |>
  left_join(bind_rows(mutate(lineages(comm98), experiment = "T98"),
                      mutate(lineages(comm99), experiment = "T99")),
            by = c("experiment", "comp"))

share <- composition |>
  group_by(experiment) |>
  summarise(within_cell_share = within_share(mpd_within_groups, cell))

print(group_summary, n = Inf, width = Inf)
print(share)
print(composition |> group_by(experiment, cell) |>
        summarise(n = n(), mean = mean(mpd_within_groups), sd = sd(mpd_within_groups),
                  .groups = "drop"))

write_csv(group_summary, here("output", "tables", "within_trophic_group_summary.csv"))
write_csv(composition, here("output", "tables", "within_trophic_composition.csv"))
write_csv(share, here("output", "tables", "within_trophic_share.csv"))
