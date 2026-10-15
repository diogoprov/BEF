# Step 3: does phylogenetic diversity modulate the response to the acid pulse?
#
# Following Downing & Leibold (2010, Freshwater Biology), the disturbance
# response is split into resistance (period 3, the first sampling after the
# pulse applied after period 2) and the recovery phase (periods 4-6).
# Disturbance effects are estimated per richness level (replicating the
# richness x resilience test) and their modulation by within-richness
# standardized phylogenetic metrics (zw_X):
#
#   y ~ 0 + div + div:zw_X + ftime + div:dist + div:dist:early + div:dist:late
#       + dist:early:zw_X + dist:late:zw_X + (1|comp) + (1|tank) + AR(1)

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))

# ---- Data ----
bef <- read_csv(here("data", "derived", "bef_long.csv"), show_col_types = FALSE)

comp_preds <- bef |>
  distinct(comp, div, PD_init, MPD_init, VPD_init, MPD_final = MPD_comp) |>
  mutate(across(-c(comp, div), \(x) z_within(x, div), .names = "zw_{.col}"))

d <- bef |>
  left_join(select(comp_preds, comp, starts_with("zw_")), by = "comp") |>
  mutate(early = as.integer(time == 3), late = as.integer(time >= 4)) |>
  prep_data()

predictors <- c("zw_MPD_init", "zw_PD_init", "zw_VPD_init", "zw_MPD_final")

rhs_dist <- function(x) {
  paste0("0 + div + div:", x, " + ftime + div:dist + div:dist:early + div:dist:late",
         " + dist:early:", x, " + dist:late:", x,
         " + (1 | comp) + (1 | tank) + ar(time = time, gr = tank, p = 1)")
}

fits <- expand_grid(resp = names(responses), x = predictors) |>
  mutate(name = paste0(paste("dist", resp, x, sep = "_"), if_else(resp %in% c("prod", "resp"), "_sig", "")),
         fit = pmap(list(resp, x, name), \(r, x, nm)
           fit_uni(responses[[r]], data = d, name = nm, rhs_override = rhs_dist(x),
                   sigma_rhs = sigma98[[responses[[r]]]])))

diag <- pmap_dfr(list(fits$fit, fits$name), diagnose)
print(diag, n = Inf)
write_csv(diag, here("output", "tables", "step3_diagnostics.csv"))
walk2(fits$fit, fits$name, save_ppc)

# ---- Effects ----
# Columns are matched by their component terms, whatever order brms uses
find_col <- function(nms, parts) {
  hit <- nms[map_lgl(str_split(str_remove(nms, "^b_"), ":"), \(p) setequal(p, parts))]
  stopifnot(length(hit) == 1)
  hit
}

effects <- pmap_dfr(list(fits$fit, fits$resp, fits$x), \(f, r, x) {
  dr <- as_draws_df(f)
  nms <- names(dr)
  rich <- map_dfr(c("S3", "S9", "S15"), \(lv) {
    e <- dr[[find_col(nms, c(paste0("div", lv), "dist", "early"))]]
    l <- dr[[find_col(nms, c(paste0("div", lv), "dist", "late"))]]
    bind_rows(
      bind_cols(tibble(term = "resistance (dist x early)", div = lv), summarise_draws_ci(e)),
      bind_cols(tibble(term = "recovery phase (dist x late)", div = lv), summarise_draws_ci(l)),
      bind_cols(tibble(term = "recovery (late - early)", div = lv), summarise_draws_ci(l - e))
    )
  })
  phy <- bind_rows(
    bind_cols(tibble(term = "x modulates resistance", div = "all"),
              summarise_draws_ci(dr[[find_col(nms, c("dist", "early", x))]])),
    bind_cols(tibble(term = "x modulates recovery phase", div = "all"),
              summarise_draws_ci(dr[[find_col(nms, c("dist", "late", x))]]))
  )
  bind_rows(rich, phy) |> mutate(response = r, predictor = x, .before = 1)
})
write_csv(effects, here("output", "tables", "step3_disturbance_effects.csv"))
print(effects |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf, width = Inf)
