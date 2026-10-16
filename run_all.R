# Run the full analysis pipeline from raw data. Fitted models are cached in
# output/models/ (git-ignored): a fresh run refits everything (~3-4 h on 4+
# cores); later runs only refit models whose data or code changed.
#
# Requires macOS textutil (step 00a converts the Nature SI .doc to HTML) and
# the SI file in refs/ (not versioned; see README).

library(here)

run <- function(script) {
  message("\n==== ", script, " ====")
  t0 <- Sys.time()
  if (grepl("\\.py$", script)) {
    html <- file.path(tempdir(), "si.html")
    stopifnot(system2("textutil", c("-convert", "html", "-output", shQuote(html),
                                    shQuote(here("refs", "Downing_Leibold_2002_Nature_SI.doc")))) == 0)
    stopifnot(system2("python3", c("-I", shQuote(here(script)), shQuote(html),
                                   shQuote(here("data", "derived", "si_table1_long.csv")))) == 0)
  } else {
    stopifnot(system2("Rscript", shQuote(here(script))) == 0)
  }
  message("done in ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
}

# T98 (Downing & Leibold 2002)
run("analysis/00a_initial_composition.py")   # initial composition from the Nature SI
run("analysis/00_build_data.R")              # derived data and phylogenetic metrics
run("analysis/01_models.R")                  # final-composition metrics; AR(1) check; multivariate model
run("analysis/02_initial_composition.R")     # step 2: metrics fixed a priori, within richness
run("analysis/03_disturbance.R")             # step 3: resistance and recovery

# T99 (Downing 2005)
run("analysis/10_build_data_t99.R")
run("analysis/11_models_t99.R")              # steps 2-3 and cells vs phylogenetic metrics
run("analysis/12_influence_t99.R")           # leave-one-composition-out (effects excluding zero)
run("analysis/13_species_identity_t99.R")    # Box: phylogenetically structured species effects
run("analysis/17_within_trophic_check.R")    # MPD within trophic groups (descriptive, T98 and T99)
Sys.setenv(FULL_SIM = "1")
run("analysis/18_design_simulation.R")       # Box: design simulation (~7 min on 17 cores)
run("analysis/19_design_simulation_brms.R")  # Bayesian check of the simulation
run("analysis/20_figure_design_box.R")       # Box figure

# BARG checks
run("analysis/14_prior_checks.R")            # prior predictive and prior sensitivity
run("analysis/15_ppc_detail.R")              # grouped posterior predictive checks
run("analysis/16_variance_sensitivity_t99.R") # alternative variance structures (T99)

# Report
stopifnot(system2("quarto", c("render", shQuote(here("analysis", "report.qmd")))) == 0)
