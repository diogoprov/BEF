# Data preparation and model components shared by the T99 scripts
# (analysis/11_models_t99.R, 12_influence_t99.R). Requires R/brms_helpers.R.

responses99 <- c(prod = "y99_prod", resp = "y99_resp", dc = "y99_dc")
families <- c(families, list(y99_prod = gaussian(), y99_resp = gaussian(), y99_dc = student()))
sigma99 <- list(prod = "ftime", resp = "ftime", dc = NULL)
metrics99 <- c("PD_init", "MPD_init", "VPD_init", "MPD_final")

prep_t99 <- function(path = here::here("data", "derived", "t99_long.csv")) {
  t99 <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::filter(comp > 0)

  comp_preds <- t99 |>
    dplyr::distinct(comp, cell, dplyr::across(dplyr::all_of(metrics99))) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(metrics99), \(x) z_within(x, cell), .names = "zw_{.col}"),
                  dplyr::across(dplyr::all_of(metrics99), z, .names = "zg_{.col}"))

  d <- t99 |>
    dplyr::left_join(dplyr::select(comp_preds, comp, dplyr::starts_with("zw_"), dplyr::starts_with("zg_")),
                     by = "comp") |>
    dplyr::mutate(
      y99_prod = z(prod), y99_resp = z(resp), y99_dc = z(dc),
      cell = factor(cell), rich = factor(spdv, levels = c(6, 18), labels = c("S6", "S18")),
      ftime = factor(time), comp = factor(comp), tank = factor(tank)
    )
  # Sum-to-zero contrasts for sampling period: with treatment contrasts the
  # cell means and period effects traded off a common level (R-hat up to 1.019)
  contrasts(d$ftime) <- contr.sum(nlevels(d$ftime))
  d
}

random_ar <- " + (1 | comp) + (1 | tank) + ar(time = time, gr = tank, p = 1)"

rhs_c <- \(x) paste0("0 + cell + ", x, " + ftime + rich:dist + rich:dist:early + rich:dist:late",
                     " + dist:early:", x, " + dist:late:", x, random_ar)

# Match a coefficient by its component terms, whatever order brms uses
find_col <- function(nms, parts) {
  hit <- nms[purrr::map_lgl(stringr::str_split(stringr::str_remove(nms, "^b_"), ":"), \(p) setequal(p, parts))]
  stopifnot(length(hit) == 1)
  hit
}
