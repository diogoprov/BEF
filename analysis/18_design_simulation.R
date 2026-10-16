# Box (design): which experimental designs can separate a relatedness effect
# from trophic structure and species identity?
#
# Generative model (composition-mean scale; response standardized), following
# the parametrization of analysis/13_species_identity_t99.R:
#   y_c = X_c b + theta * z(MPD_wg,c) + u_phy,c + u_id,c + e_c
#   u_phy ~ MVN(0, s_phy^2 K_phy), K_phy = P C P' (C: within-group phylogenetic
#     correlation, block-diagonal across trophic groups)
#   u_id  ~ MVN(0, s_id^2 K_id),   K_id  = P P'
#   e_c   ~ N(0, s_comp^2 + s_tank^2 / n_tank)   (non-additive + tank noise)
# Kernels are projected orthogonally to the fixed effects and scaled to unit
# mean diagonal (as in script 13), so s_phy^2 and s_id^2 are the same
# between-composition variances in every design. MPD_wg = mean of the
# within-group MPDs. theta is the effect of 1 SD of MPD_wg among *random*
# compositions from the design's candidate set (what the species pool offers);
# selecting compositions along MPD increases the realized spread and power.
#
# Designs (n_comp compositions x n_tank tanks):
#   D1  Downing-like: T99 pools, 3 groups, k = 2 or 6 species per group (half
#       the compositions each), random draws; richness as fixed effect.
#   D2  Gradient: T99 pools, k = 3 per group, compositions spread evenly along
#       MPD_wg.
#   D3  Gradient + balanced incidence: as D2, chosen to keep the MPD_wg spread
#       while minimizing max |cor(species incidence, MPD_wg)|.
#   D4b Focal group = T99 macrophytes (8 spp., k = 3), grazers and predators
#       fixed in all tanks; selection as in D3. Separates "other groups fixed"
#       from "large pool".
#   D4  Focal group = 20 species on a simulated pure-birth tree (crown age of
#       the T99 macrophytes), k = 4, other groups fixed; selection as in D3;
#       10 simulated trees, cycled over replicates.
# With only 56 distinct macrophyte triplets, D4b repeats compositions when
# n_comp > 56; repeats share species effects but not the non-additive term
# (slightly optimistic for D4b at n_comp = 96).
#
# Estimation per replicate:
#   naive: OLS on composition means (= LMM with iid composition effect when
#          tanks per composition are balanced)
#   full:  REML with K_phy, K_id and iid terms (as in script 13), Wald test
#
# FULL_SIM=1 runs the full grid; otherwise a small calibration run.
# Outputs: output/tables/design_sim_summary[_small].csv,
#          output/models/design_sim_reps[_small].rds (git-ignored),
#          output/tables/design_sim_designs.csv (realized design properties)

suppressPackageStartupMessages({library(tidyverse); library(ape); library(here); library(parallel)})
source(here("R", "phylo_metrics.R"))

FULL   <- identical(Sys.getenv("FULL_SIM"), "1")
N_REP  <- if (FULL) 1000 else 50
N_COMPS <- if (FULL) c(24, 48, 96) else c(24, 96)
N_TANK <- 4
N_CAND <- 5000
N_TREES <- 10
N_CORES <- max(1, detectCores() - 1)
RNGkind("L'Ecuyer-CMRG")
set.seed(20261015)

# ---- Variance components (scenario values; T98/T99 sd_comp ~ 0.15-0.30) ----
vc_scen <- tibble(
  scenario = c("phy+id", "phy only", "id only"),
  s_phy = c(0.15, 0.21, 0), s_id = c(0.15, 0, 0.21)
)
S_COMP <- 0.10    # non-additive composition variation
S_TANK <- 0.30    # tank-level residual of the time-averaged response
THETAS <- c(0, 0.1, 0.2)

# ---- Species pools (T99) ----
tree99 <- read.tree(here("data", "derived", "phylo_timetree_t99.tre")) |> fix_tip_labels()
pools <- list(
  macrophyte = c("Ceratophyllum_demersum", "Elodea_nuttallii", "Myriophyllum_verticullatum",
                 "Potamogeton_crispus", "Potamogeton_natans", "Sagittaria_rigida",
                 "Utricularia_vulgaris", "Vallisneria_americana"),
  grazer = c("Lithobates_catesbeianus", "Lithobates_clamitans", "Crangonyx_richmondensis",
             "Hyalella_azteca", "Helisoma_trivolvis", "Physa_gyrina", "Hesperocorixa_sp",
             "Trichocorixa_sp"),
  predator = c("Ambrysus_sp", "Belostoma_flumireum", "Notonecta_undulata", "Buenoa_sp",
               "Neoplea_striola", "Coptotomus_sp", "Gyrinus_sp", "Dineutus_sp")
)
stopifnot(setequal(unlist(pools), tree99$tip.label))

group_mat <- \(t) list(D = cophenetic(t), C = vcv(t, corr = TRUE))
mats99 <- map(pools, \(s) group_mat(keep.tip(tree99, s)))

# Simulated focal pools for D4 (crown age = T99 macrophyte crown age)
mac_crown <- max(mats99$macrophyte$D) / 2
trees_sim <- map(seq_len(N_TREES), \(i) {
  t <- phytools::pbtree(n = 20, scale = mac_crown)
  t$tip.label <- paste0("focal_", seq_len(20))
  t
})

# ---- Composition builders ----
mpd_of <- function(D, s) mean(D[s, s][lower.tri(D[s, s])])

# A composition is a named list of species vectors, one per group
draw_comp <- function(mats, k) imap(mats, \(m, g) sample(rownames(m$D), k[[g]]))

comp_stats <- function(comps, mats) {
  sp_all <- unlist(map(mats, \(m) rownames(m$D)), use.names = FALSE)
  P <- t(vapply(comps, \(cp) { x <- as.numeric(sp_all %in% unlist(cp)); x / sum(x) },
                numeric(length(sp_all))))
  colnames(P) <- sp_all
  C <- as.matrix(Matrix::bdiag(map(mats, "C"))); dimnames(C) <- list(sp_all, sp_all)
  mpd <- vapply(comps, \(cp) mean(unlist(imap(cp, \(s, g)
    if (length(s) > 1) mpd_of(mats[[g]]$D, s) else NA_real_)), na.rm = TRUE), numeric(1))
  list(P = P, mpd = mpd, C = C, inc = (P > 0) * 1)
}

add_kernels <- function(st, X) {
  M <- diag(nrow(X)) - X %*% solve(crossprod(X)) %*% t(X)
  scale1 <- \(K) { K <- M %*% K %*% M; K / mean(diag(K)) }
  c(st, list(X = X, K_phy = scale1(st$P %*% st$C %*% t(st$P)),
             K_id = scale1(st$P %*% t(st$P))))
}

max_cor_inc <- function(inc, mpd) {
  I <- inc[, apply(inc, 2, var) > 0, drop = FALSE]
  max(abs(cor(I, mpd)))
}

# Spread compositions evenly along MPD_wg
pick_gradient <- function(mpd, n) {
  idx <- integer(0)
  for (t in seq(min(mpd), max(mpd), length.out = n)) {
    idx <- c(idx, setdiff(order(abs(mpd - t)), idx)[1])
  }
  idx
}

# Swap search: keep MPD spread, minimize the strongest incidence-MPD correlation
pick_balanced <- function(mpd, inc, n, ref_sd, iter = 4000, lambda = 2) {
  score <- \(i) sd(mpd[i]) / ref_sd - lambda * max_cor_inc(inc[i, , drop = FALSE], mpd[i])
  cur <- pick_gradient(mpd, n); s_cur <- score(cur)
  for (it in seq_len(iter)) {
    prop <- cur; prop[sample(n, 1)] <- sample(setdiff(seq_along(mpd), cur), 1)
    s_prop <- score(prop)
    if (s_prop > s_cur) { cur <- prop; s_cur <- s_prop }
  }
  cur
}

# Distinct candidate compositions (species sorted within groups)
candidates <- function(draw, mats, n = N_CAND) {
  comps <- replicate(n, map(draw(), sort), simplify = FALSE)
  comps <- comps[!duplicated(map_chr(comps, \(cp) paste(unlist(cp), collapse = "|")))]
  c(comp_stats(comps, mats), list(comps = comps))
}

# Selection when n may exceed the number of distinct candidates: use every
# distinct composition, then add repeats chosen by the same rule
pick_n <- function(pick, mpd, n) {
  if (n <= length(mpd)) return(pick(n))
  c(seq_along(mpd), pick(n - length(mpd)))
}

select_design <- function(cand, mats, idx) {
  st <- comp_stats(cand$comps[idx], mats)
  c(add_kernels(st, matrix(1, length(idx), 1)), list(ref_sd = sd(cand$mpd)))
}

# ---- Fixed designs (D2, D3, D4b, D4) for each n_comp ----
k3 <- list(macrophyte = 3, grazer = 3, predator = 3)
cand3 <- candidates(\() draw_comp(mats99, k3), mats99)

fixed_rest <- list(grazer = sample(pools$grazer, 3), predator = sample(pools$predator, 3))
mats4b <- mats99
cand4b <- candidates(\() c(list(macrophyte = sample(pools$macrophyte, 3)), fixed_rest), mats4b)

mats4 <- map(trees_sim, \(t) c(list(focal = group_mat(t)), mats99[c("grazer", "predator")]))
cand4 <- map(mats4, \(m) candidates(\() c(list(focal = sample(rownames(m$focal$D), 4)),
                                         fixed_rest), m))

designs_fixed <- map(set_names(N_COMPS), \(n) list(
  D2  = list(select_design(cand3, mats99, pick_gradient(cand3$mpd, n))),
  D3  = list(select_design(cand3, mats99, pick_balanced(cand3$mpd, cand3$inc, n, sd(cand3$mpd)))),
  D4b = list(select_design(cand4b, mats4b, pick_n(
    \(m) pick_balanced(cand4b$mpd, cand4b$inc, m, sd(cand4b$mpd)), cand4b$mpd, n))),
  D4  = map2(cand4, mats4, \(c4, m4) select_design(c4, m4,
                                                   pick_balanced(c4$mpd, c4$inc, n, sd(c4$mpd))))
))

# D1 is redrawn in every replicate (random assembly, as in the experiments);
# theta is scaled by the SD of MPD_wg among random compositions within levels
ref_sd1 <- local({
  kk <- \(k) list(macrophyte = k, grazer = k, predator = k)
  sqrt((var(candidates(\() draw_comp(mats99, kk(2)), mats99, 2000)$mpd) +
        var(candidates(\() draw_comp(mats99, kk(6)), mats99, 2000)$mpd)) / 2)
})
draw_d1 <- function(n) {
  ks <- rep(c(2, 6), each = n / 2)
  comps <- map(ks, \(k) draw_comp(mats99, list(macrophyte = k, grazer = k, predator = k)))
  c(add_kernels(comp_stats(comps, mats99), model.matrix(~ factor(ks))), list(ref_sd = ref_sd1))
}

design_props <- imap_dfr(designs_fixed, \(dl, n) imap_dfr(dl, \(ds, nm) map_dfr(ds, \(d) tibble(
  n_comp = as.integer(n), design = nm,
  realized_sd = sd(d$mpd) / d$ref_sd,
  max_cor_inc = max_cor_inc(d$inc, d$mpd),
  n_distinct = nrow(unique(d$P)),
  kernel_alignment = sum(d$K_phy * d$K_id) / sqrt(sum(d$K_phy^2) * sum(d$K_id^2))
))))
# D1 is random: summarize 200 draws per n_comp (own seed, so the simulation's
# random stream is unchanged)
design_props <- bind_rows(design_props, withr::with_seed(1, map_dfr(N_COMPS, \(n) map_dfr(1:200, \(i) {
  d <- draw_d1(n)
  tibble(n_comp = n, design = "D1", realized_sd = sd(d$mpd) / d$ref_sd,
         max_cor_inc = max_cor_inc(d$inc, d$mpd), n_distinct = nrow(unique(d$P)),
         kernel_alignment = sum(d$K_phy * d$K_id) / sqrt(sum(d$K_phy^2) * sum(d$K_id^2)))
}))))
write_csv(design_props, here("output", "tables", "design_sim_designs.csv"))
dir.create(here("output", "models"), showWarnings = FALSE)
saveRDS(designs_fixed, here("output", "models", "design_sim_designs.rds"))
if (identical(Sys.getenv("DESIGNS_ONLY"), "1")) quit(save = "no")
print(design_props |> group_by(n_comp, design) |> summarise(across(everything(), mean), .groups = "drop"))

# ---- Estimation ----
reml_fit <- function(y, Xf, Ks) {
  n <- length(y)
  build_V <- \(s2) Reduce(`+`, map2(s2[-length(s2)], Ks, `*`)) + diag(s2[length(s2)], n)
  nll <- function(lp) {
    R <- tryCatch(chol(build_V(exp(lp))), error = \(e) NULL)
    if (is.null(R)) return(1e10)
    Vi <- chol2inv(R)
    ViX <- Vi %*% Xf
    XtViX <- crossprod(Xf, ViX)
    b <- solve(XtViX, crossprod(ViX, y))
    r <- y - Xf %*% b
    0.5 * (2 * sum(log(diag(R))) + determinant(XtViX)$modulus + sum(r * (Vi %*% r)))
  }
  o <- optim(rep(log(0.02), length(Ks) + 1), nll, method = "Nelder-Mead",
             control = list(maxit = 2000))
  s2 <- exp(o$par)
  Vi <- chol2inv(chol(build_V(s2)))
  cov_b <- solve(crossprod(Xf, Vi %*% Xf))
  b <- cov_b %*% crossprod(Xf, Vi %*% y)
  list(est = b[ncol(Xf)], se = sqrt(cov_b[ncol(Xf), ncol(Xf)]), s2 = s2)
}

rmvn <- function(K) {
  e <- eigen(K, symmetric = TRUE)
  drop(e$vectors %*% (sqrt(pmax(e$values, 0)) * rnorm(nrow(K))))
}

one_rep <- function(des, theta, s_phy, s_id) {
  n <- nrow(des$P)
  z <- (des$mpd - mean(des$mpd)) / des$ref_sd
  y <- theta * z + s_phy * rmvn(des$K_phy) + s_id * rmvn(des$K_id) +
    rnorm(n, 0, sqrt(S_COMP^2 + S_TANK^2 / N_TANK))
  Xf <- cbind(des$X, z)
  ols <- summary(lm(y ~ 0 + Xf))$coefficients[ncol(Xf), ]
  full <- reml_fit(y, Xf, list(des$K_phy, des$K_id))
  c(naive_est = ols[[1]], naive_p = ols[[4]],
    full_est = full$est, full_p = 2 * pnorm(-abs(full$est / full$se)),
    s2_phy = full$s2[1], s2_id = full$s2[2])
}

grid <- expand_grid(n_comp = N_COMPS, design = c("D1", "D2", "D3", "D4b", "D4"),
                    vc_scen, theta = THETAS)

run_cell <- function(i) {
  g <- grid[i, ]
  res <- t(vapply(seq_len(N_REP), \(r) {
    des <- if (g$design == "D1") draw_d1(g$n_comp) else {
      ds <- designs_fixed[[as.character(g$n_comp)]][[g$design]]
      ds[[(r - 1) %% length(ds) + 1]]
    }
    one_rep(des, g$theta, g$s_phy, g$s_id)
  }, numeric(6)))
  bind_cols(g[rep(1, N_REP), ], as_tibble(res), rep = seq_len(N_REP))
}

t0 <- Sys.time()
sims <- bind_rows(mclapply(seq_len(nrow(grid)), run_cell, mc.cores = N_CORES,
                           mc.preschedule = FALSE))
elapsed <- difftime(Sys.time(), t0, units = "mins")

mc_se <- \(p) sqrt(p * (1 - p) / N_REP)
summary_tab <- sims |>
  group_by(n_comp, design, scenario, s_phy, s_id, theta) |>
  summarise(
    reject_naive = mean(naive_p < 0.05), reject_full = mean(full_p < 0.05),
    bias_full = mean(full_est) - first(theta), sd_full = sd(full_est),
    rmse_full = sqrt(mean((full_est - first(theta))^2)),
    phylo_signal_med = median(s2_phy / (s2_phy + s2_id), na.rm = TRUE),
    phylo_signal_q10 = quantile(s2_phy / (s2_phy + s2_id), 0.1, na.rm = TRUE),
    phylo_signal_q90 = quantile(s2_phy / (s2_phy + s2_id), 0.9, na.rm = TRUE),
    .groups = "drop") |>
  mutate(mcse_full = mc_se(reject_full))

print(summary_tab, n = Inf, width = Inf)
message("Replicates per cell: ", N_REP, "; cores: ", N_CORES, "; elapsed: ", round(elapsed, 1), " min")

suffix <- if (FULL) "" else "_small"
write_csv(summary_tab, here("output", "tables", paste0("design_sim_summary", suffix, ".csv")))
dir.create(here("output", "models"), showWarnings = FALSE)
saveRDS(sims, here("output", "models", paste0("design_sim_reps", suffix, ".rds")))
