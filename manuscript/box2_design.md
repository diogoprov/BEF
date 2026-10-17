# Box 2. What would an experiment need to separate relatedness from trophic structure and species identity?

<!-- Numbers from output/tables/design_sim_summary.csv and design_sim_brms.csv
     (analysis/18-20). Scenario shown: phylogenetic + identity species effects
     (s_phy = s_id = 0.15), effect size 0.1 SD of the response per SD of MPD. -->

Phylogenetic diversity has been proposed as a predictor of ecosystem
functioning from single trophic levels to complex food webs (Srivastava et
al. 2012). We simulated experiments under the same model used in Box 1: each
species has an effect on the response, these effects covary with
phylogenetic relatedness (variance σ²_phy) or vary idiosyncratically (σ²_id),
as in phylogenetic mixed models of community data (Ives and Helmus 2011),
and mean pairwise distance within trophic groups (MPD) may add an effect of
its own, analogous to the community phylogenetic diversity term added to
Diversity–Interactions models (Kirwan et al. 2009, Connolly et al. 2011). Variance
components were of the magnitude estimated in the two experiments. We
compared five designs of equal size: random assembly at two richness levels,
as in the original experiments (D1); compositions spread along MPD with the
species pools of the 1999 experiment (D2); the same gradient with the
incidence of each species balanced along it (D3); and a single focal group
varying while the rest of the food web was held constant, drawn from a pool
of 8 (D4b) or 20 species (D4).

With eight species per group, a gradient in MPD is largely a gradient in the
presence of a few deeply divergent taxa (Fig. B2a). In D2, MPD correlated at
0.61 with the presence of *Buenoa* sp. and at about 0.52 with that of
*Physa gyrina*, *Coptotomus* sp. and *Lithobates catesbeianus*. Balancing incidence
reduced this correlation, and selecting compositions along MPD more than
doubled power relative to random assembly (Fig. B2b; 0.97 vs 0.37 with 24
compositions). Selecting compositions along relatedness at fixed richness
has been done within single trophic levels (Cadotte 2013, Venail and Vives
2013), and balancing species incidence across treatments follows the logic
of designs that separate species identity from richness (Bell et al. 2009).

Treating compositions as independent units inflated false positives, and the
inflation increased with the number of compositions because compositions
share species (Fig. B2c; 10–25% with 24 compositions, 16–52% with 96, for a
nominal 5%), as also reported for naive regressions across communities
that overlap in species and phylogeny (Tsang and Cadotte 2025). Modelling
shared species with covariance kernels reduced the
rate to 4–8% with 96 compositions, and to 5–7% with 24 compositions when
variance-component uncertainty was propagated in Bayesian fits (D2 and D4).

Phylogenetically structured species effects could not be distinguished from
idiosyncratic identity effects with 24 compositions in any design, nor
with 48 in most designs (Fig. B2d); estimates began to separate only with
96 compositions, and not even then when a single group of eight species
varied (D4b). Testing whether evolutionary history adds to species identity
therefore requires large species pools within a trophic level, compositions
chosen along relatedness with balanced species incidence, and analyses that
model shared species explicitly.


## Figure B2 legend

**Figure B2.** Simulated performance of five experimental designs for
detecting an effect of within-trophic-group phylogenetic relatedness (MPD) on
ecosystem functioning. D1, random assembly with two or six species per
trophic group, as in the original experiments; D2, 24–96 compositions spread
evenly along MPD, with three species per group drawn from the pools of the
1999 experiment; D3, as D2, with compositions chosen to minimize the
correlation between the incidence of any species and MPD; D4b, a single focal
group (the eight macrophytes) varying while grazers and predators are held
constant; D4, as D4b, with a focal pool of 20 species on simulated
pure-birth trees. Each composition was replicated in four tanks.
(a) Strongest correlation between the incidence of a single species and MPD
across the selected compositions (light symbols, 24 compositions; dark
symbols, 96); the dashed line marks the mean for random assembly.
(b) Power to detect an effect of 0.1 SD of the response per SD of MPD, with
species effects modelled through covariance kernels; the dotted line marks
80%. (c) False-positive rate when the true effect is zero, treating
compositions as independent (left) or modelling shared species (right;
restricted maximum likelihood with Wald tests); the solid line marks the
nominal 5%. (d) Estimated phylogenetic share of the variance in species
effects, σ²_phy / (σ²_phy + σ²_id), when the truth is 1 (left) or 0 (right);
points are medians and bars span the 10th–90th percentiles of 1000
replicates. Variance components: σ_phy = σ_id = 0.15 in (a)–(c); σ_phy = 0.21
or σ_id = 0.21 in (d); non-additive composition SD 0.10, tank SD 0.30.
