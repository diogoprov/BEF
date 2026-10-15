# Exploratory (Box): additive species-identity effects with and without
# phylogenetic structure, T99.
#
# If species i has effect beta_i and composition c contains species with equal
# proportions p_ci = 1/S_c (initial, stocked composition), the composition
# effect is u = P beta. With beta ~ MVN(0, s_phy^2 C + s_id^2 I), C the
# phylogenetic correlation matrix, u has covariance
#   s_phy^2 P C P' + s_id^2 P P',
# i.e. two composition-level random effects with known covariance kernels
# (as in an "animal model"), plus an iid composition term for non-additive
# (interaction) effects. This is the identity part of Diversity-Interaction
# models (Kirwan et al. 2009) with a phylogenetic prior on species effects
# (cf. Connolly et al. 2011 for phylogenetically structured interactions).
#
# Fixed effects: richness, sampling period, disturbance (no design cells), so
# the species terms can absorb functional-group composition effects. Kernels
# are projected orthogonally to the richness means (fixed effects) and scaled
# to unit mean diagonal. Identifiability is checked by comparing posterior and
# prior-only fits.

library(tidyverse)
library(brms)
library(posterior)
library(ape)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))
source(here("R", "phylo_metrics.R"))
source(here("R", "t99_data.R"))

d <- prep_t99()

# ---- Kernels ----
comm <- read_csv(here("data", "derived", "t99_comm_initial.csv"), show_col_types = FALSE) |>
  column_to_rownames("comp") |>
  as.matrix()
P <- comm / rowSums(comm)
tree <- read.tree(here("data", "raw", "phylogeny", "phylog_dated_2005.txt")) |> fix_tip_labels()
C <- vcv(tree, corr = TRUE)[colnames(P), colnames(P)]

rich_c <- distinct(d, comp, spdv) |> arrange(as.integer(as.character(comp)))
stopifnot(all(rownames(P) == as.character(rich_c$comp)))
X <- model.matrix(~ factor(rich_c$spdv))
M <- diag(nrow(P)) - X %*% solve(crossprod(X)) %*% t(X)   # projection orthogonal to richness means

make_kernel <- function(K) {
  K <- M %*% K %*% M
  K <- K / mean(diag(K))
  K <- K + diag(1e-4, nrow(K))     # jitter: projected kernels are rank-deficient
  dimnames(K) <- list(rownames(P), rownames(P))
  K
}
K_phy <- make_kernel(P %*% C %*% t(P))
K_id  <- make_kernel(P %*% t(P))

alignment <- sum(K_phy * K_id) / sqrt(sum(K_phy^2) * sum(K_id^2))
message("Kernel alignment (phy vs id): ", round(alignment, 3))

d <- d |> mutate(comp_phy = comp, comp_id = comp)

rhs_sp_for <- \(resp) paste0("0 + rich + ftime + dist + dist:post",
                             " + (1 | gr(comp_phy, cov = K_phy)) + (1 | gr(comp_id, cov = K_id))",
                             random_for(resp))

fit_sp <- function(resp, prior_only = FALSE) {
  r <- responses99[[resp]]
  sig <- sigma99[[resp]]
  rhs_sp <- rhs_sp_for(resp)
  pri <- if (grepl("ar(", rhs_sp, fixed = TRUE)) priors_uni else priors_uni[priors_uni$class != "ar", ]
  if (!is.null(sig)) pri <- c(pri[pri$class != "sigma", ], priors_sigma)
  if (families[[r]]$family == "student") pri <- c(pri, prior_nu)
  brm(
    if (is.null(sig)) bf(paste(r, "~", rhs_sp)) else bf(as.formula(paste(r, "~", rhs_sp)), as.formula(paste("sigma ~", sig))),
    data = d, data2 = list(K_phy = K_phy, K_id = K_id), family = families[[r]], prior = pri,
    sample_prior = if (prior_only) "only" else "no",
    iter = 2 * ITER, chains = 4, seed = SEED, control = list(adapt_delta = 0.99),
    file = here("output", "models", paste0("t99sp_", resp, if (prior_only) "_prior")),
    file_refit = "on_change"
  )
}

fits <- expand_grid(resp = names(responses99), prior_only = c(FALSE, TRUE)) |>
  mutate(fit = map2(resp, prior_only, fit_sp),
         name = paste0("t99sp_", resp, if_else(prior_only, "_prior", "")))

diag <- pmap_dfr(list(fits$fit[!fits$prior_only], fits$name[!fits$prior_only]), diagnose)
print(diag)
write_csv(diag, here("output", "tables", "t99_species_identity_diagnostics.csv"))
walk2(fits$fit[!fits$prior_only], fits$name[!fits$prior_only], save_ppc)

# ---- Variance components and phylogenetic share ----
vc <- pmap_dfr(list(fits$fit, fits$resp, fits$prior_only), \(f, r, po) {
  dr <- as_draws_df(f)
  s_phy <- dr$sd_comp_phy__Intercept; s_id <- dr$sd_comp_id__Intercept; s_c <- dr$sd_comp__Intercept
  tot <- s_phy^2 + s_id^2 + s_c^2
  list(sd_phy = s_phy, sd_id = s_id, sd_comp_iid = s_c,
       share_phy = s_phy^2 / tot, share_id = s_id^2 / tot, share_nonadditive = s_c^2 / tot,
       phylo_signal = s_phy^2 / (s_phy^2 + s_id^2)) |>
    imap_dfr(\(v, k) tibble(quantity = k, median = median(v),
                            q2.5 = quantile(v, 0.025), q97.5 = quantile(v, 0.975))) |>
    mutate(response = r, fit = if_else(po, "prior", "posterior"), .before = 1)
})
write_csv(vc, here("output", "tables", "t99_species_identity_variance.csv"))
print(vc |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf)

# Posterior-to-prior contraction of the phylogenetic share (1 = fully informed)
contraction <- vc |>
  filter(quantity %in% c("phylo_signal", "share_nonadditive")) |>
  mutate(width = q97.5 - q2.5) |>
  select(response, fit, quantity, width) |>
  pivot_wider(names_from = fit, values_from = width) |>
  mutate(contraction = 1 - posterior / prior)
write_csv(contraction, here("output", "tables", "t99_species_identity_contraction.csv"))
print(contraction)
