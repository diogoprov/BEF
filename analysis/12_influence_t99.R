# Influence of single compositions on the T99 result "MPD modulates the
# recovery phase after the acid pulse" (analysis/11_models_t99.R, part c).
#
# Each of the 4 models (prod, resp x MPD_init, MPD_final) is refitted 26
# times, each time leaving one composition out (predictor values are kept as
# standardized in the full data). The compiled model is reused via update().
# Summaries are cached per fit in output/models/influence/ so the script can
# resume.

library(tidyverse)
library(brms)
library(posterior)
library(here)

options(brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))
source(here("R", "t99_data.R"))

d <- prep_t99()
cache_dir <- here("output", "models", "influence")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

targets <- expand_grid(resp = c("prod", "resp"), x = c("zw_MPD_init", "zw_MPD_final"))
comps <- levels(d$comp)

summ <- function(fit, x) {
  dr <- as_draws_df(fit)
  col <- find_col(names(dr), c("dist", "late", x))
  v <- dr[[col]]
  bind_cols(summarise_draws_ci(v),
            tibble(rhat = posterior::rhat(extract_variable_matrix(fit, col)),
                   ess_bulk = posterior::ess_bulk(extract_variable_matrix(fit, col))))
}

influence <- pmap_dfr(targets, \(resp, x) {
  full <- fit_uni(responses99[[resp]], data = d, name = paste("t99c", resp, x, "sig", sep = "_"),
                  rhs_override = rhs_c(x), sigma_rhs = sigma99[[resp]])
  # Fits loaded from file carry no compiled CmdStan executable: compile once
  # with the first leave-one-out fit and reuse it as template
  drop_fit <- function(template, k, recompile) {
    f_cache <- file.path(cache_dir, paste0(resp, "_", x, "_drop", k, ".rds"))
    if (file.exists(f_cache) && !recompile) return(readRDS(f_cache))
    fit_k <- update(template, newdata = droplevels(filter(d, comp != k)), recompile = recompile,
                    chains = 4, cores = 4, iter = ITER, seed = SEED, refresh = 0)
    out <- summ(fit_k, x) |> mutate(dropped = k)
    saveRDS(out, f_cache)
    if (recompile) attr(out, "fit") <- fit_k
    out
  }
  first <- drop_fit(full, comps[1], recompile = TRUE)
  template <- attr(first, "fit")
  attr(first, "fit") <- NULL
  res <- parallel::mclapply(comps[-1], \(k) drop_fit(template, k, recompile = FALSE), mc.cores = 4)
  bad <- !map_lgl(res, is.data.frame)
  if (any(bad)) stop("Leave-one-out fits failed: ", paste(comps[-1][bad], collapse = ", "), "\n", res[bad][[1]])
  res <- c(list(first), res)
  bind_rows(summ(full, x) |> mutate(dropped = "none"), bind_rows(res)) |>
    mutate(response = resp, predictor = x, .before = 1)
})

write_csv(influence, here("output", "tables", "t99_influence_mpd_recovery.csv"))

influence |>
  group_by(response, predictor) |>
  summarise(full = mean[dropped == "none"],
            min_mean = min(mean[dropped != "none"]), max_mean = max(mean[dropped != "none"]),
            most_influential = dropped[dropped != "none"][which.max(abs(mean[dropped != "none"] - full))],
            n_ci95_excl0 = sum(q97.5[dropped != "none"] < 0 | q2.5[dropped != "none"] > 0),
            n_ci98_excl0 = sum(q99.17[dropped != "none"] < 0 | q0.83[dropped != "none"] > 0),
            max_rhat = max(rhat), .groups = "drop") |>
  print(width = Inf)
