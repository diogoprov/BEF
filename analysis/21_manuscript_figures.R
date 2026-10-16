# Publication figures for the manuscript (Figs 1-5 and B1). Reads the trees
# in data/derived/, the tables in output/tables/ and the draws in
# output/draws/, so it can be rerun without refitting any model. Figure B2
# is drawn by 20_figure_design_box.R.
#
# Palettes checked for colour-vision deficiency (minimum pairwise CIEDE2000
# under deutan/protan/tritan simulation, colorspace + farver): functional
# groups 31.0, stocked vs realized 30.3, T98 vs T99 29.9. Groups and
# experiments are also told apart by shape, so identity never relies on colour.
#
# Inputs (read-only): data/derived/phylo_timetree_{t98_21spp,t99}.tre,
#   data/derived/si_table1_long.csv, output/tables/*.csv,
#   output/draws/{pulse,slope}_draws_{t98,t99}.rds
# Outputs: figures/manuscript/{fig1..fig5,figB1}_*.{pdf,png}
#   (180 mm wide = double column; pdf vector, png 600 dpi)

suppressPackageStartupMessages({
  library(tidyverse); library(ape); library(ggtree); library(patchwork); library(here)
})
source(here("R", "phylo_metrics.R"))

tab <- \(f) read_csv(here("output", "tables", paste0(f, ".csv")), show_col_types = FALSE)
out_dir <- here("figures", "manuscript")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, height, width = 180) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height,
         units = "mm", device = cairo_pdf, bg = "white")
  ragg::agg_png(file.path(out_dir, paste0(name, ".png")), width = width, height = height,
                units = "mm", res = 600, background = "white")
  print(p); invisible(dev.off())
}

# ---- Shared style (same inks as 20_figure_design_box.R) ----
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"
theme_ms <- theme_minimal(base_size = 9) +
  theme(text = element_text(colour = ink),
        axis.text = element_text(colour = ink2),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.3),
        panel.spacing = unit(4, "mm"),
        strip.text = element_text(colour = ink, face = "bold", hjust = 0),
        plot.tag = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", colour = NA))

fg_col <- c(Macrophytes = "#1baf7a", Grazers = "#6250d6", Predators = "#eda100")
fg_shp <- c(Macrophytes = 21, Grazers = 22, Predators = 24)
comp_col <- c(Stocked = "#3d3c39", Realized = "#eda100")
exp_col <- c(T98 = "#6250d6", T99 = "#e34948")
exp_shp <- c(T98 = 21, T99 = 24)

resp_lab <- c(prod = "Productivity", resp = "Respiration", dc = "Decomposition")
lab_response <- \(x) factor(resp_lab[x], resp_lab)
# zw_MPD_init -> metric "MPD", composition "Stocked"
split_predictor <- \(df) df |>
  mutate(predictor = str_remove(predictor, "^z[wg]_")) |>
  separate(predictor, c("metric", "composition"), sep = "_") |>
  mutate(metric = factor(metric, c("PD", "MPD", "VPD")),
         composition = factor(c(init = "Stocked", final = "Realized")[composition],
                              c("Stocked", "Realized")))

# Point + 95% (thick) and 98.33% (thin) credible intervals
# (fill = NULL when fill and shape are mapped to a variable)
ci_layers <- function(dodge = 0, fill = ink) {
  pos <- position_dodge(width = dodge)
  pt <- if (is.null(fill)) geom_point(size = 1.9, stroke = 0.5, colour = "white", position = pos) else
    geom_point(size = 1.9, stroke = 0.5, colour = "white", shape = 21, fill = fill, position = pos)
  list(geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3, linetype = "22", show.legend = FALSE),
       geom_linerange(aes(xmin = q0.83, xmax = q99.17), linewidth = 0.35, position = pos),
       geom_linerange(aes(xmin = q2.5, xmax = q97.5), linewidth = 1, position = pos),
       pt)
}
effect_lab <- "Effect (SD of response per SD of metric)"

# ---- Fig. 1: TimeTree-dated pseudochronograms ----
groups <- read_csv(here("data", "derived", "si_table1_long.csv"), show_col_types = FALSE) |>
  distinct(species, group)
plants <- "Potamogeton|Myriophyllum|Elodea|Ceratophyllum|Vallisneria|Utricularia|Sagittaria"
grazers <- "Lithobates|Crangonyx|Hyalella|Helisoma|Physa|corixa"
fg_of <- \(x) case_when(str_detect(x, plants) ~ "Macrophytes", str_detect(x, grazers) ~ "Grazers",
                        TRUE ~ "Predators")
stopifnot(all(fg_of(groups$species) ==
                c(macrophyte = "Macrophytes", grazer = "Grazers", predator = "Predators")[groups$group]))

tree_panel <- function(file, title) {
  tr <- read.tree(here("data", "derived", file)) |> fix_tip_labels()
  age <- max(node.depth.edgelength(tr))
  tips <- tibble(label = tr$tip.label, group = fg_of(label),
                 name = str_replace(label, "_sp$", " sp.") |> str_replace_all("_", " "))
  ggtree(tr, linewidth = 0.35, colour = ink2) %<+% tips +
    geom_tippoint(aes(fill = group, shape = group), size = 1.8, stroke = 0.3, colour = "white") +
    geom_tiplab(aes(label = name), fontface = "italic", size = 2.3, offset = age * 0.02, colour = ink) +
    scale_fill_manual(values = fg_col, breaks = names(fg_col), name = NULL) +
    scale_shape_manual(values = fg_shp, breaks = names(fg_col), name = NULL) +
    scale_x_continuous(breaks = age - seq(0, 1500, 250), labels = seq(0, 1500, 250),
                       expand = expansion(mult = c(0.02, 0.55))) +
    coord_cartesian(clip = "off") +
    theme_tree2() +
    theme(text = element_text(colour = ink, size = 8), axis.text.x = element_text(colour = ink2, size = 7),
          axis.line.x = element_line(colour = ink2, linewidth = 0.3),
          axis.title.x = element_text(colour = ink, size = 8, margin = margin(t = 4)),
          plot.title = element_text(face = "bold", size = 9), legend.position = "bottom") +
    # after theme_tree2(), which drops the axis title
    labs(title = title, x = "Million years ago")
}
fig1 <- (tree_panel("phylo_timetree_t98_21spp.tre", "a  T98 (21 species)") |
           tree_panel("phylo_timetree_t99.tre", "b  T99 (24 species)")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(fig1, "fig1_phylogenies", height = 140)

# ---- Fig. 2: variance of each metric within design strata ----
fig2_dat <- bind_rows(
  tab("step2_within_div_variance") |> transmute(exp = "T98", metric, share = within_div_share),
  tab("t99_within_cell_variance") |> transmute(exp = "T99", metric, share = within_cell_share)
) |>
  filter(!str_detect(metric, "MNTD")) |>
  rename(predictor = metric) |> split_predictor() |>
  mutate(exp = factor(exp, c("T98", "T99"),
                      c("T98: within richness levels", "T99: within design cells")))
fig2 <- ggplot(fig2_dat, aes(share, fct_rev(metric), fill = composition)) +
  geom_col(position = position_dodge2(width = 0.8, padding = 0.15, reverse = TRUE), width = 0.7) +
  geom_text(aes(label = scales::percent(share, accuracy = 1), group = composition),
            position = position_dodge2(width = 0.7, reverse = TRUE), hjust = -0.15, size = 2.4,
            colour = ink2) +
  facet_wrap(~exp) +
  scale_fill_manual(values = comp_col, name = "Composition") +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Variance among compositions within design strata", y = NULL) +
  theme_ms
save_fig(fig2, "fig2_within_strata_variance", height = 70)

# ---- Fig. 3: within-strata slopes, T98 and T99 ----
# Full posteriors from draws saved by 02_initial_composition.R (T98: average
# of the within-richness slopes) and 11_models_t99.R (T99: within-cell slope).
# Each slab is scaled to its own maximum, so the narrow T99 posteriors do not
# flatten the T98 ones.
fig3_dat <- bind_rows(
  readRDS(here("output", "draws", "slope_draws_t98.rds")) |> mutate(exp = "T98"),
  readRDS(here("output", "draws", "slope_draws_t99.rds")) |> mutate(exp = "T99")
) |>
  split_predictor() |>
  mutate(response = lab_response(response), exp = factor(exp, c("T99", "T98")),
         row = fct_rev(interaction(metric, composition, sep = ", ", lex.order = TRUE)))
fig3 <- ggplot(fig3_dat, aes(value, row, fill = exp)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3, linetype = "22") +
  ggdist::stat_halfeye(point_interval = ggdist::mean_qi, .width = c(0.95, 1 - 0.05 / 6),
                       normalize = "groups", position = position_dodge(width = 0.55),
                       height = 0.55, slab_alpha = 0.8, colour = ink, point_size = 1.1,
                       interval_size_range = c(0.35, 0.9)) +
  facet_wrap(~response) +
  scale_fill_manual(values = exp_col, breaks = c("T98", "T99"), name = NULL) +
  labs(x = effect_lab, y = NULL) +
  theme_ms
save_fig(fig3, "fig3_slopes", height = 110)

# ---- Fig. 4: design cells vs single metrics (T99) ----
# Point estimates of the variance ratio can be slightly negative (a metric
# that adds nothing); intervals are truncated at -10%.
fig4_dat <- tab("t99_cells_vs_phylo_r2") |>
  filter(model != "rich_only") |>
  mutate(model = factor(model, c("rich_PD", "rich_MPD", "rich_VPD", "cells"),
                        c("PD", "MPD", "VPD", "Design cells")),
         response = lab_response(resp))
fig4 <- ggplot(fig4_dat, aes(r2_comp, fct_rev(model))) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3) +
  geom_linerange(aes(xmin = pmax(q2.5, -0.1), xmax = q97.5), linewidth = 1, colour = "#b9b8b3") +
  geom_point(aes(fill = model == "Design cells"), shape = 21, size = 2.2, stroke = 0.5, colour = "white") +
  scale_fill_manual(values = c(`TRUE` = ink, `FALSE` = exp_col[["T99"]]), guide = "none") +
  facet_wrap(~response) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(-0.1, 1)) +
  labs(x = "Among-composition variance explained (relative to richness only)", y = NULL) +
  theme_ms
save_fig(fig4, "fig4_design_vs_metrics_t99", height = 55)

# ---- Fig. 5: acid pulse by richness level, T98 (a) and T99 (b) ----
# Full posteriors (half-eye: density, mean, 95% thick and 98.33% thin
# intervals) from draws saved by 03_disturbance.R and 11_models_t99.R.
# Fill follows stocked richness on a multi-hue sequential ramp (viridis,
# begin = 0.15, end = 0.9: yellow-green -> teal -> purple), monotonic in
# lightness, so it also reads in grey scale; richness is on the y axis too.
rich_fill <- c(S3 = "#BBDF27", S6 = "#49C16D", S9 = "#1F968B", S15 = "#30698E", S18 = "#463480")
pulse_panel <- function(file, levels, tag_title) {
  readRDS(here("output", "draws", file)) |>
    pivot_longer(c(resistance, recovery), names_to = "phase", values_to = "value") |>
    mutate(phase = factor(phase, c("recovery", "resistance"), c("Recovery", "Resistance")),
           response = lab_response(response),
           level = factor(level, rev(levels))) |>
    ggplot(aes(value, level, fill = level)) +
    geom_vline(xintercept = 0, colour = ink2, linewidth = 0.3, linetype = "22") +
    ggdist::stat_halfeye(point_interval = ggdist::mean_qi, .width = c(0.95, 1 - 0.05 / 6),
                         slab_alpha = 0.85, height = 0.9, colour = ink, point_size = 1.3,
                         interval_size_range = c(0.4, 1.1)) +
    facet_grid(phase ~ response) +
    scale_fill_manual(values = rich_fill, guide = "none") +
    scale_y_discrete(labels = \(x) str_remove(x, "^S")) +
    labs(title = tag_title, x = "Effect of the pulse (SD of response)", y = "Stocked richness") +
    theme_ms + theme(plot.title = element_text(face = "bold", size = 9),
                     strip.text.y = element_text(angle = 0))
}
fig5 <- pulse_panel("pulse_draws_t98.rds", c("S3", "S9", "S15"), "a  T98") /
  pulse_panel("pulse_draws_t99.rds", c("S6", "S18"), "b  T99")
save_fig(fig5, "fig5_pulse", height = 160)

# ---- Fig. B1: prior vs posterior of the variance shares (T99) ----
figB1_dat <- tab("t99_species_identity_variance") |>
  filter(quantity %in% c("phylo_signal", "share_nonadditive")) |>
  mutate(quantity = factor(quantity, c("phylo_signal", "share_nonadditive"),
                           c("Phylogenetic share of species effects", "Non-additive share")),
         fit = factor(str_to_sentence(fit), c("Prior", "Posterior")),
         response = lab_response(response))
figB1 <- ggplot(figB1_dat, aes(median, fct_rev(response), colour = fit, fill = fit, shape = fit)) +
  geom_linerange(aes(xmin = q2.5, xmax = q97.5), linewidth = 1, position = position_dodge(0.5)) +
  geom_point(size = 2, stroke = 0.5, colour = "white", position = position_dodge(0.5)) +
  facet_wrap(~quantity) +
  scale_colour_manual(values = c(Prior = "#3d3c39", Posterior = "#eda100"), name = NULL) +
  scale_fill_manual(values = c(Prior = "#3d3c39", Posterior = "#eda100"), name = NULL) +
  scale_shape_manual(values = c(Prior = 21, Posterior = 24), name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  guides(colour = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE),
         shape = guide_legend(reverse = TRUE)) +
  labs(x = "Share of variance (median and 95% interval)", y = NULL) +
  theme_ms
save_fig(figB1, "figB1_species_effects", height = 60)

message("Figures written to ", out_dir)
