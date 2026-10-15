# Phylogenetic diversity metrics for the Downing & Leibold (2002) mesocosms

# Biomass columns in "species biomass T98" mapped to tip labels of the
# 2002 phylogeny. Validated in analysis/00_build_data.R: the union of tank-level
# presences per composition reproduces data/raw/nature2001/treat.txt exactly.
biomass_to_species <- c(
  po1_dw_g     = "Potamogeton_zosteriformis",
  elo_dw_g     = "Elodea_nuttallii",
  pcrisp_dw_g  = "Potamogeton_crispus",
  cer_dw_g     = "Ceratophyllum_demersum",
  mir_dw_g     = "Myriophyllum_verticullatum",
  bul_dw       = "Lithobates_catesbeianus",
  gre_dw       = "Lithobates_clamitans",
  gam_dw       = "Crangonyx_richmondensis",
  hya_dw       = "Hyalella_azteca",
  hel_dw_tot_g = "Helisoma_trivolvis",
  phy_dw       = "Physa_gyrina",
  cor_dw       = "Trichocorixa_sp",
  nau_dw       = "Ambrysus_sp",
  ple_dw       = "Neoplea_striola",
  not2sm_dw    = "Notonecta_undulata",
  not1big_dw   = "Notonecta_sp",
  bel_dw       = "Belostoma_flumireum",
  whi_dw       = "Gyrinus_sp",
  dit_dw       = "Acilius_sp"
)

# Presence/absence matrix (sites x species) from final dry-weight biomass
biomass_presence <- function(biomass, site_col = "tank") {
  m <- as.matrix(biomass[, names(biomass_to_species)])
  m[is.na(m)] <- 0
  m <- (m > 0) * 1
  colnames(m) <- unname(biomass_to_species)
  rownames(m) <- biomass[[site_col]]
  m
}

# Spelling used in the 2014 trees -> current tip labels
fix_tip_labels <- function(tree) {
  tree$tip.label <- dplyr::recode(tree$tip.label,
    Potamageton_zosteriformis = "Potamogeton_zosteriformis",
    Potamageton_crispus = "Potamogeton_crispus",
    Hyallela_azteca = "Hyalella_azteca",
    Helisoma_trivolis = "Helisoma_trivolvis"
  )
  tree
}

# Variance of pairwise phylogenetic distances (regularity dimension;
# Tucker et al. 2017). NA for sites with fewer than two species.
vpd <- function(comm, dist_mat) {
  apply(comm, 1, \(x) {
    sp <- names(x)[x > 0]
    if (length(sp) < 2) return(NA_real_)
    d <- dist_mat[sp, sp]
    var(d[lower.tri(d)])
  })
}

# Faith's PD (root included), MPD, MNTD, VPD, PSV and richness for each site.
# Unweighted (presence/absence) metrics; MPD/MNTD/VPD are NA for
# single-species sites.
phylo_metrics <- function(comm, tree) {
  mm <- picante::match.phylo.comm(tree, comm)
  dist_mat <- ape::cophenetic.phylo(mm$phy)
  pd <- picante::pd(mm$comm, mm$phy, include.root = TRUE)
  tibble::tibble(
    site = rownames(mm$comm),
    SR   = pd$SR,
    PD   = pd$PD,
    MPD  = picante::mpd(mm$comm, dist_mat),
    MNTD = picante::mntd(mm$comm, dist_mat),
    VPD  = vpd(mm$comm, dist_mat),
    PSV  = picante::psv(mm$comm, mm$phy)$PSVs
  )
}
