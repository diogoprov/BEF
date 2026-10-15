# Steps 1-3 replicated in the 1999 functional-group experiment (T99).
#
# (a) Step 2 analogue: PD, MPD, VPD (initial composition; MPD also final)
#     standardized within design cells (FG composition x species richness):
#     y ~ 0 + cell + zw_X + ftime + dist + dist:post + (1|comp) + (1|tank) + AR(1)
# (b) Trophic-structure test: can a single phylogenetic metric replace the
#     design cells? Among-composition variance (sd_comp^2) is compared between
#     a baseline with richness only, models adding a globally standardized
#     metric, and the model with design cells.
# (c) Step 3 analogue: resistance (periods 5-6, right after the acid pulse) and
#     recovery phase (periods 7-13) by species richness, and their modulation
#     by zw_X.
#
# Responses are standardized (prod and resp contain negative values, so log is
# not an option). prod and resp: Gaussian with residual SD varying by sampling
# period (sigma ~ ftime; their SD varies 3-5 fold over time, which a first
# Student-t fit absorbed as nu ~ 2). dc: Student-t with constant sigma.
# AR(1) is defined on the sampling index; spacing is irregular (1-10 days).
# Sampling period uses sum-to-zero contrasts.
# Composition 0 (no species) is excluded.

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))

source(here("R", "t99_data.R"))
metrics <- metrics99
d <- prep_t99()

# ---- (a) within-cell slopes ----
rhs_a <- \(x, resp) paste0("0 + cell + ", x, " + ftime + dist + dist:post", random_for(resp))
fits_a <- expand_grid(resp = names(responses99), x = paste0("zw_", metrics)) |>
  mutate(name = paste("t99a", resp, x, "sig", sep = "_"),
         fit = pmap(list(resp, x, name), \(r, x, nm)
           fit_uni(responses99[[r]], data = d, name = nm, rhs_override = rhs_a(x, r), sigma_rhs = sigma99[[r]])))

slopes_a <- pmap_dfr(list(fits_a$fit, fits_a$resp, fits_a$x), \(f, r, x)
  bind_cols(tibble(response = r, predictor = x),
            summarise_draws_ci(as_draws_df(f)[[paste0("b_", x)]])))

# ---- (b) can a phylogenetic metric replace the design cells? ----
rhs_b <- list(
  rich_only = "0 + rich + ftime + dist + dist:post",
  rich_PD   = "0 + rich + zg_PD_init + ftime + dist + dist:post",
  rich_MPD  = "0 + rich + zg_MPD_init + ftime + dist + dist:post",
  rich_VPD  = "0 + rich + zg_VPD_init + ftime + dist + dist:post",
  cells     = "0 + cell + ftime + dist + dist:post"
)
fits_b <- expand_grid(resp = names(responses99), model = names(rhs_b)) |>
  mutate(name = paste("t99b", resp, model, "sig", sep = "_"),
         fit = pmap(list(resp, model, name), \(r, m, nm)
           fit_uni(responses99[[r]], data = d, name = nm, rhs_override = paste0(rhs_b[[m]], random_for(r)),
                   sigma_rhs = sigma99[[r]])))

sd_comp <- \(f) as_draws_df(f)$sd_comp__Intercept
set.seed(SEED)
r2_b <- fits_b |>
  group_by(resp) |>
  group_modify(\(g, k) {
    s0 <- sd_comp(g$fit[[which(g$model == "rich_only")]])
    map_dfr(seq_len(nrow(g)), \(i) {
      s <- sd_comp(g$fit[[i]])
      v <- 1 - s^2 / sample(s0^2, length(s), replace = TRUE)
      tibble(model = g$model[i], sd_comp = median(s),
             r2_comp = 1 - mean(s^2) / mean(s0^2),
             q2.5 = quantile(v, 0.025), q97.5 = quantile(v, 0.975))
    })
  }) |>
  ungroup()

slopes_b <- fits_b |>
  filter(str_detect(model, "_(PD|MPD|VPD)$")) |>
  pmap_dfr(\(resp, model, name, fit) {
    x <- paste0("zg_", str_remove(model, "rich_"), "_init")
    bind_cols(tibble(response = resp, model = model),
              summarise_draws_ci(as_draws_df(fit)[[paste0("b_", x)]]))
  })

# ---- (c) disturbance: resistance and recovery, modulation by zw_X ----
fits_c <- expand_grid(resp = names(responses99), x = paste0("zw_", metrics)) |>
  mutate(name = paste("t99c", resp, x, "sig", sep = "_"),
         fit = pmap(list(resp, x, name), \(r, x, nm)
           fit_uni(responses99[[r]], data = d, name = nm, rhs_override = rhs_c(x, r), sigma_rhs = sigma99[[r]])))

effects_c <- pmap_dfr(list(fits_c$fit, fits_c$resp, fits_c$x), \(f, r, x) {
  dr <- as_draws_df(f)
  nms <- names(dr)
  rich <- map_dfr(c("S6", "S18"), \(lv) {
    e <- dr[[find_col(nms, c(paste0("rich", lv), "dist", "early"))]]
    l <- dr[[find_col(nms, c(paste0("rich", lv), "dist", "late"))]]
    bind_rows(
      bind_cols(tibble(term = "resistance (dist x early)", rich = lv), summarise_draws_ci(e)),
      bind_cols(tibble(term = "recovery phase (dist x late)", rich = lv), summarise_draws_ci(l)),
      bind_cols(tibble(term = "recovery (late - early)", rich = lv), summarise_draws_ci(l - e))
    )
  })
  phy <- bind_rows(
    bind_cols(tibble(term = "x modulates resistance", rich = "all"),
              summarise_draws_ci(dr[[find_col(nms, c("dist", "early", x))]])),
    bind_cols(tibble(term = "x modulates recovery phase", rich = "all"),
              summarise_draws_ci(dr[[find_col(nms, c("dist", "late", x))]]))
  )
  bind_rows(rich, phy) |> mutate(response = r, predictor = x, .before = 1)
})

# ---- Diagnostics and outputs ----
all_fits <- bind_rows(select(fits_a, name, fit), select(fits_b, name, fit), select(fits_c, name, fit))
diag <- pmap_dfr(list(all_fits$fit, all_fits$name), diagnose)
walk2(all_fits$fit, all_fits$name, save_ppc)

write_csv(diag, here("output", "tables", "t99_diagnostics.csv"))
write_csv(slopes_a, here("output", "tables", "t99_step2_slopes.csv"))
write_csv(r2_b, here("output", "tables", "t99_cells_vs_phylo_r2.csv"))
write_csv(slopes_b, here("output", "tables", "t99_cells_vs_phylo_slopes.csv"))
write_csv(effects_c, here("output", "tables", "t99_step3_disturbance.csv"))

r2 <- \(df) mutate(df, across(where(is.numeric), \(x) round(x, 2)))
print(diag, n = Inf)
print(r2(slopes_a), n = Inf, width = Inf)
print(r2(r2_b), n = Inf, width = Inf)
print(r2(slopes_b), n = Inf, width = Inf)
print(r2(effects_c), n = Inf, width = Inf)
