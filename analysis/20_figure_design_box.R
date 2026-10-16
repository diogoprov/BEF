# Figure for the design Box (analysis/18_design_simulation.R output).
# Palette checked for colour-vision deficiency (worst pair OKLab Delta E:
# normal 20.8, deutan 6.8, protan 9.1, tritan 17.2); designs also differ in
# point shape and are labelled directly, so identity never relies on colour.
#
# Outputs: figures/box_design_simulation.{png,pdf}

suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(here)})

sim <- read_csv(here("output", "tables", "design_sim_summary.csv"), show_col_types = FALSE)
props <- read_csv(here("output", "tables", "design_sim_designs.csv"), show_col_types = FALSE)

design_levels <- c("D1", "D2", "D3", "D4b", "D4")
design_labels <- c(
  D1  = "D1 Random assembly",
  D2  = "D2 MPD gradient",
  D3  = "D3 Gradient + balanced incidence",
  D4b = "D4b Focal group, 8-species pool",
  D4  = "D4 Focal group, 20-species pool"
)
pal <- c(D1 = "#3d3c39", D2 = "#eda100", D3 = "#e34948", D4b = "#1baf7a", D4 = "#6250d6")
shp <- c(D1 = 21, D2 = 22, D3 = 23, D4b = 24, D4 = 25)

ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"
theme_box <- theme_minimal(base_size = 9) +
  theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 9.5),
    plot.subtitle = element_text(colour = ink2, size = 8),
    strip.text = element_text(colour = ink2, face = "bold", hjust = 0),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA)
  )

sim <- sim |> mutate(design = factor(design, design_levels))
x_scale <- scale_x_continuous(breaks = c(24, 48, 96), trans = "log2",
                              expand = expansion(mult = c(0.08, 0.08)))
common <- list(scale_colour_manual(values = pal, labels = design_labels, name = NULL, drop = FALSE),
               scale_fill_manual(values = pal, labels = design_labels, name = NULL, drop = FALSE),
               scale_shape_manual(values = shp, labels = design_labels, name = NULL, drop = FALSE),
               guides(alpha = "none"))

# (a) Who drives the relatedness gradient?
pa_dat <- props |>
  filter(design != "D1", n_comp %in% c(24, 96)) |>
  group_by(n_comp, design) |>
  summarise(max_cor = mean(max_cor_inc), .groups = "drop") |>
  mutate(design = factor(design, design_levels))
d1_cor <- props |> filter(design == "D1") |> pull(max_cor_inc) |> mean()
pa <- ggplot(pa_dat,
             aes(x = max_cor, y = fct_rev(design), colour = design, fill = design, shape = design)) +
  geom_vline(xintercept = d1_cor, colour = pal["D1"], linetype = "22", linewidth = 0.4) +
  annotate("text", x = d1_cor, y = 4.45, label = "random\nassembly", hjust = -0.08, vjust = 1,
           size = 2.4, colour = ink2, lineheight = 0.9) +
  geom_line(aes(group = design), linewidth = 0.6, alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(alpha = factor(n_comp)), size = 2.6, stroke = 0.6, colour = "white",
             show.legend = FALSE) +
  scale_alpha_manual(values = c(`24` = 0.45, `96` = 1)) +
  scale_y_discrete(labels = \(x) x) +
  scale_x_continuous(limits = c(0.2, 0.8), breaks = seq(0.2, 0.8, 0.2)) +
  common +
  labs(title = "a  One species can carry the gradient",
       subtitle = "Strongest correlation between a species' presence and MPD\n(light: 24, dark: 96 compositions)",
       x = "max |r| (species incidence, MPD)", y = NULL) +
  theme_box

# (b) Power at theta = 0.1 SD, phy+id scenario, kernel model
pb_dat <- sim |> filter(scenario == "phy+id", theta == 0.1)
pb <- ggplot(pb_dat, aes(n_comp, reject_full, colour = design, fill = design, shape = design)) +
  geom_hline(yintercept = 0.8, colour = ink2, linetype = "22", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4, stroke = 0.6, colour = "white") +
  geom_text(data = filter(pb_dat, n_comp == 24),
            aes(label = design), colour = ink, nudge_x = -0.32, size = 2.5, fontface = "bold") +
  x_scale + scale_y_continuous(limits = c(0.3, 1), labels = scales::percent) +
  common +
  labs(title = "b  Selecting compositions buys power",
       subtitle = "Detection rate for an effect of 0.1 SD per SD of MPD",
       x = "Number of compositions (4 tanks each)", y = "Power") +
  theme_box

# (c) Type I error: independent compositions vs species-level kernels
pc_dat <- sim |>
  filter(scenario == "phy+id", theta == 0) |>
  pivot_longer(c(reject_naive, reject_full), names_to = "model", values_to = "rate") |>
  mutate(model = factor(model, c("reject_naive", "reject_full"),
                        c("Compositions treated as independent", "Shared species modelled (kernels)")))
pc <- ggplot(pc_dat, aes(n_comp, rate, colour = design, fill = design, shape = design)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.05, fill = "#f1f0ec") +
  geom_hline(yintercept = 0.05, colour = ink2, linewidth = 0.3) +
  geom_line(linewidth = 0.7, show.legend = FALSE) +
  geom_point(size = 2.2, stroke = 0.6, colour = "white", show.legend = FALSE) +
  facet_wrap(~model) +
  x_scale + scale_y_continuous(limits = c(0, 0.55), labels = scales::percent) +
  common +
  labs(title = "c  Shared species inflate false positives",
       subtitle = "Rate of 'significant' MPD effects when the true effect is zero (nominal 5%)",
       x = "Number of compositions", y = "False-positive rate") +
  theme_box

# (d) Can phylogenetically structured species effects be told apart from
# idiosyncratic identity effects? 10-90% range of the estimated share
pd_dat <- sim |>
  filter(theta == 0, scenario %in% c("phy only", "id only")) |>
  mutate(scenario = factor(scenario, c("phy only", "id only"),
                           c("Truth: phylogenetically structured (share = 1)",
                             "Truth: idiosyncratic identity (share = 0)")),
         truth = if_else(str_detect(scenario, "= 1"), 1, 0),
         xpos = log2(n_comp) + (as.integer(design) - 3) * 0.14)
pd <- ggplot(pd_dat, aes(x = xpos, colour = design, fill = design, shape = design)) +
  geom_hline(aes(yintercept = truth), colour = ink2, linewidth = 0.3, linetype = "22") +
  geom_linerange(aes(ymin = phylo_signal_q10, ymax = phylo_signal_q90), linewidth = 0.9,
                 show.legend = FALSE) +
  geom_point(aes(y = phylo_signal_med), size = 2.1, stroke = 0.6, colour = "white", show.legend = FALSE) +
  facet_wrap(~scenario) +
  scale_x_continuous(breaks = log2(c(24, 48, 96)), labels = c(24, 48, 96)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  common +
  labs(title = "d  Phylogenetic signal in species effects is rarely identifiable",
       subtitle = "Estimated phylogenetic share of species-effect variance (median, 10-90% of replicates)",
       x = "Number of compositions", y = "Phylogenetic share") +
  theme_box

fig <- wrap_plots(a = pa, b = pb, c = pc, d = pd, design = "ab\ncc\ndd") +
  plot_layout(guides = "collect") &
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE)) &
  theme(legend.position = "bottom", legend.justification = "center",
        legend.text = element_text(size = 8, margin = margin(r = 10)))

ragg::agg_png(here("figures", "box_design_simulation.png"), width = 180, height = 230,
              units = "mm", res = 300, background = "white")
print(fig); invisible(dev.off())
ggsave(here("figures", "box_design_simulation.pdf"), fig, width = 180, height = 230,
       units = "mm", device = cairo_pdf, bg = "white")
