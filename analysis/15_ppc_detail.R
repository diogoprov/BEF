# Detailed posterior predictive checks for the primary models (initial-MPD
# models of T98 and T99): density on the observed range, mean and SD by
# sampling period and by design stratum, and LOO-PIT calibration.

library(tidyverse)
library(brms)
library(bayesplot)
library(patchwork)
library(here)

options(mc.cores = 4, brms.backend = "cmdstanr")
source(here("R", "brms_helpers.R"))
source(here("R", "t99_data.R"))
dir.create(here("figures", "ppc_detail"), showWarnings = FALSE, recursive = TRUE)

bef <- read_csv(here("data", "derived", "bef_long.csv"), show_col_types = FALSE)
d98 <- bef |>
  left_join(bef |> distinct(comp, div, MPD_init) |>
              mutate(zw_MPD_init = z_within(MPD_init, div)) |> select(comp, zw_MPD_init),
            by = "comp") |>
  prep_data()
d99 <- prep_t99()
rhs_a99 <- \(resp) paste0("0 + cell + zw_MPD_init + ftime + dist + dist:post", random_for(resp))

models <- bind_rows(
  tibble(exp = "T98", resp = names(responses), y = unname(responses)),
  tibble(exp = "T99", resp = names(responses99), y = unname(responses99))
)

ppc_detail <- function(exp, resp, y) {
  if (exp == "T98") {
    f <- fit_uni(y, "zw_MPD_init", data = d98, sigma_rhs = sigma98[[y]])
    strata <- "div"
  } else {
    f <- fit_uni(y, data = d99, name = paste0("t99a_", resp, "_zw_MPD_init_sig"),
                 rhs_override = rhs_a99(resp), sigma_rhs = sigma99[[resp]])
    strata <- "cell"
  }
  dat <- f$data
  yobs <- dat[[y]]
  set.seed(SEED)
  ids <- sample(ndraws(f), 400)
  yrep <- posterior_predict(f, draw_ids = ids)
  lims <- range(yobs) + c(-1, 1)

  p1 <- ppc_dens_overlay(yobs, yrep[1:100, ]) + coord_cartesian(xlim = lims) + ggtitle("Density")
  p2 <- ppc_stat_grouped(yobs, yrep, group = dat$ftime, stat = "mean") + ggtitle("Mean by period")
  p3 <- ppc_stat_grouped(yobs, yrep, group = dat$ftime, stat = "sd") + ggtitle("SD by period")
  p4 <- ppc_stat_grouped(yobs, yrep, group = dat[[strata]], stat = "sd") + ggtitle(paste("SD by", strata))
  # LOO-PIT only without AR(1): with AR terms the pointwise log-likelihood is
  # conditional on the previous observation and PSIS-LOO is not valid
  has_ar <- "ar[1]" %in% variables(f)
  p5 <- if (!has_ar) {
    lw <- weights(loo::psis(-log_lik(f, draw_ids = ids), r_eff = NA))
    ppc_loo_pit_overlay(yobs, yrep, lw = lw) + ggtitle("LOO-PIT")
  } else {
    ggplot() + annotate("text", x = 0, y = 0, size = 3.5,
                        label = "LOO-PIT not shown:\nPSIS-LOO is not valid\nfor AR(1) models") + theme_void()
  }

  n_per <- nlevels(factor(dat$ftime))
  h_per <- ceiling(n_per / 4)
  p <- (p1 | p5) / p2 / p3 / p4 +
    plot_layout(heights = c(1.4, h_per, h_per, 1.2)) +
    plot_annotation(title = paste(exp, resp, "— initial-MPD model")) &
    theme(legend.position = "none", axis.text.x = element_text(size = 6))
  ggsave(here("figures", "ppc_detail", paste0(exp, "_", resp, ".png")), p,
         width = 10, height = 6 + 2.6 * h_per, dpi = 110)

  # Numerical summary: observed statistic vs posterior predictive interval
  stat_check <- function(g, fun) {
    map_dfr(levels(factor(g)), \(lv) {
      i <- which(g == lv)
      s_rep <- apply(yrep[, i, drop = FALSE], 1, fun)
      tibble(level = lv, obs = fun(yobs[i]), p = mean(s_rep >= fun(yobs[i])))
    })
  }
  bind_rows(
    stat_check(dat$ftime, mean) |> mutate(stat = "mean", group = "period"),
    stat_check(dat$ftime, sd) |> mutate(stat = "sd", group = "period"),
    stat_check(dat[[strata]], sd) |> mutate(stat = "sd", group = strata)
  ) |> mutate(exp = exp, response = resp, .before = 1)
}

ppp <- pmap_dfr(models, ppc_detail)
write_csv(ppp, here("output", "tables", "ppc_grouped_pvalues.csv"))
ppp |> filter(p < 0.05 | p > 0.95) |> print(n = Inf)
