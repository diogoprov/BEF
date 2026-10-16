# Re-date the three pseudochronograms with TimeTree node ages (bladj).
#
# The 2013 pseudochronograms were dated in Phylocom (bladj) with very few
# calibrated nodes, so bladj placed shallow undated nodes halfway between deep
# calibrations and the tips (e.g., Lithobates catesbeianus vs L. clamitans at
# 391 Myr, Physa vs Helisoma at 293 Myr). Here the topologies follow TimeTree
# (see below) and every node that TimeTree 5 dates is calibrated with its median pairwise
# age (timetree.org/api/pairwise, queried 2026-10-16; ages, CIs and number of
# studies in data/raw/phylogeny/timetree_calibrations.csv). Remaining nodes are
# interpolated by bladj, run through the R package phylocomr.
#
# Topology: the full TimeTree pairwise age matrix of the 27 species (queried
# 2026-10-16) agrees with the 2013 trees except in two places, both changed:
#   - Nepomorpha: TimeTree groups Belostomatidae with Corixidae (116.6 Myr)
#     and Naucoridae with Notonectidae (196.8 Myr), whereas the 2013 trees had
#     Belostoma sister to all other bugs. TimeTree has no data for Neoplea
#     (Pleidae), which is kept as sister to Notonectidae, as in the 2013 trees.
#   - T99: the 2013 tree grouped Sagittaria with Vallisneria, making
#     Hydrocharitaceae (Elodea + Vallisneria) non-monophyletic; Sagittaria
#     (Alismataceae) is sister to Hydrocharitaceae.
# No TimeTree data within Lithobates s.l., Notonectidae, Gyrinidae, Corixidae
# or for Neoplea; these nodes are interpolated, except that the Lithobates
# pair is fixed at the age of Lithobates s.l. vs Rana (an upper bound) instead
# of being interpolated from Bilateria.
#
# Inputs (read-only): data/raw/phylogeny/phylog_dated_{2002,2002_21spp_4b4da4d,2005}.txt,
#   data/raw/phylogeny/timetree_calibrations.csv
# Outputs: data/derived/phylo_timetree_{t98_19spp,t98_21spp,t99}.tre,
#   output/tables/phylo_redating_{node_ages,metric_correlations}.csv

suppressPackageStartupMessages({
  library(tidyverse); library(ape); library(phylocomr); library(here)
})
source(here("R", "phylo_metrics.R"))

cal <- read_csv(here("data", "raw", "phylogeny", "timetree_calibrations.csv"),
                show_col_types = FALSE)

old <- list(
  t98_19spp = read.tree(here("data", "raw", "phylogeny", "phylog_dated_2002.txt")),
  t98_21spp = read.tree(here("data", "raw", "phylogeny", "phylog_dated_2002_21spp_4b4da4d.txt")),
  t99       = read.tree(here("data", "raw", "phylogeny", "phylog_dated_2005.txt"))
) |> map(fix_tip_labels)

# Topology only, rearranged to the TimeTree topology
topology <- function(tr) {
  tr$edge.length <- NULL; tr$node.label <- NULL; tr$root.edge <- NULL
  nwk <- sub("((Sagittaria_rigida,Vallisneria_americana),Elodea_nuttallii)",
             "(Sagittaria_rigida,(Vallisneria_americana,Elodea_nuttallii))",
             write.tree(tr), fixed = TRUE)
  tr <- read.tree(text = nwk)
  # Nepomorpha: ((Belostoma, Corixidae), (Ambrysus, (Neoplea, Notonectidae)))
  sub_nwk <- \(tips) if (length(tips) == 1) tips else
    sub(";$", "", write.tree(extract.clade(tr, getMRCA(tr, tips))))
  corix <- intersect(c("Trichocorixa_sp", "Hesperocorixa_sp"), tr$tip.label)
  old_clade <- sub_nwk(c("Belostoma_flumineum", "Ambrysus_sp"))
  new_clade <- sprintf("((Belostoma_flumineum,%s),%s)", sub_nwk(corix),
                       sub_nwk(c("Ambrysus_sp", "Notonecta_undulata")))
  stopifnot(grepl(old_clade, nwk, fixed = TRUE))
  out <- read.tree(text = sub(old_clade, new_clade, nwk, fixed = TRUE))
  stopifnot(setequal(out$tip.label, tr$tip.label),
            is.monophyletic(out, c("Belostoma_flumineum", corix)))
  out
}

redate <- function(tr) {
  nt <- Ntip(tr)
  tr$node.label <- paste0("n", seq_len(tr$Nnode))
  use <- cal |>
    filter(tip_a %in% tr$tip.label, tip_b %in% tr$tip.label) |>
    mutate(nd = map2_int(tip_a, tip_b, \(a, b) getMRCA(tr, c(a, b)))) |>
    distinct(nd, .keep_all = TRUE)
  stopifnot(!anyDuplicated(use$node), (nt + 1) %in% use$nd)
  ages <- tibble(a = tr$node.label[use$nd - nt], b = use$age_mya)
  out <- read.tree(text = ph_bladj(ages = as.data.frame(ages), phylo = write.tree(tr)))
  out$root.edge <- NULL
  # Phylocom lowercases labels
  out$tip.label <- tr$tip.label[match(out$tip.label, tolower(tr$tip.label))]
  stopifnot(setequal(out$tip.label, tr$tip.label), is.ultrametric(out, tol = 1e-4))
  # Calibrated ages are reproduced exactly
  bt <- branching.times(out)
  chk <- map2_dbl(use$tip_a, use$tip_b, \(a, b) unname(bt[getMRCA(out, c(a, b)) - Ntip(out)]))
  stopifnot(all(abs(chk - use$age_mya) < 1e-3))
  out
}

new <- map(old, \(tr) redate(topology(tr)))
iwalk(new, \(tr, nm) write.tree(tr, here("data", "derived", paste0("phylo_timetree_", nm, ".tre"))))

# ---- Old vs new node ages (clades defined by their tips) ----
node_table <- function(tr) {
  nt <- Ntip(tr); bt <- branching.times(tr)
  map_dfr((nt + 1):(nt + tr$Nnode), \(nd) tibble(
    clade = paste(sort(extract.clade(tr, nd)$tip.label), collapse = ";"),
    age = unname(bt[nd - nt])))
}
node_ages <- imap_dfr(new, \(tr, nm) {
  full_join(node_table(old[[nm]]) |> rename(age_2013 = age),
            node_table(tr) |> rename(age_timetree = age), by = "clade") |>
    mutate(tree = nm, n_tips = str_count(clade, ";") + 1,
           calibrated = round(age_timetree, 3) %in% round(cal$age_mya, 3))
}) |> relocate(tree, n_tips)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
write_csv(node_ages, here("output", "tables", "phylo_redating_node_ages.csv"))

# ---- Old vs new metrics for the compositions analysed ----
# (needs the community matrices written by 00_build_data.R and
# 10_build_data_t99.R; skipped in a run from scratch)
if (!all(file.exists(here("data", "derived", c("comm_initial.csv", "t99_comm_initial.csv"))))) {
  message("Community matrices not found; skipping the metric comparison.")
  quit(save = "no")
}
comm98_init <- read_csv(here("data", "derived", "comm_initial.csv"), show_col_types = FALSE) |>
  column_to_rownames("comp") |> as.matrix()
colnames(comm98_init) <- fix_species_names(colnames(comm98_init))
comm99_init <- read_csv(here("data", "derived", "t99_comm_initial.csv"), show_col_types = FALSE) |>
  column_to_rownames("comp") |> as.matrix()
colnames(comm99_init) <- fix_species_names(colnames(comm99_init))

compare <- function(comm, tr_old, tr_new, label) {
  m_old <- phylo_metrics(comm, tr_old); m_new <- phylo_metrics(comm, tr_new)
  map_dfr(c("PD", "MPD", "MNTD", "VPD", "PSV"), \(v) tibble(
    set = label, metric = v,
    r = cor(m_old[[v]], m_new[[v]], use = "complete.obs"),
    rho = cor(m_old[[v]], m_new[[v]], use = "complete.obs", method = "spearman")))
}
metric_cor <- bind_rows(
  compare(comm98_init, old$t98_21spp, new$t98_21spp, "T98 stocked (21 spp)"),
  compare(comm99_init[rowSums(comm99_init) > 0, ], old$t99, new$t99, "T99 stocked (24 spp)")
)
write_csv(metric_cor, here("output", "tables", "phylo_redating_metric_correlations.csv"))
print(metric_cor, n = Inf)
print(node_ages |> filter(abs(age_timetree - age_2013) > 25 | is.na(age_2013)) |>
        mutate(clade = str_trunc(clade, 60)), n = Inf, width = 200)
