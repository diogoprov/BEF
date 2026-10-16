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
0.72 with the presence of *Lithobates catesbeianus*. Balancing incidence
reduced this correlation, and selecting compositions along MPD more than
doubled power relative to random assembly (Fig. B2b; 0.86 vs 0.33 with 24
compositions). Selecting compositions along relatedness at fixed richness
has been done within single trophic levels (Cadotte 2013, Venail and Vives
2013), and balancing species incidence across treatments follows the logic
of designs that separate species identity from richness (Bell et al. 2009).

Treating compositions as independent units inflated false positives, and the
inflation increased with the number of compositions because compositions
share species (Fig. B2c; 15–22% with 24 compositions, 26–52% with 96, for a
nominal 5%), as also reported for naive regressions across communities
that overlap in species and phylogeny (Tsang and Cadotte 2025). Modelling
shared species with covariance kernels reduced the
rate to 6–8% with 96 compositions, and to 3–5% with 24 compositions when
variance-component uncertainty was propagated in Bayesian fits (D2 and D4).

Phylogenetically structured species effects could not be distinguished from
idiosyncratic identity effects in almost any design (Fig. B2d); separation
began only with a 20-species focal pool and 96 compositions. Testing whether evolutionary history adds to species identity
therefore requires large species pools within a trophic level, compositions
chosen along relatedness with balanced species incidence, and analyses that
model shared species explicitly.

## References

Bell, T., Lilley, A. K., Hector, A., Schmid, B., King, L. and Newman, J. A.
2009. A linear model method for biodiversity–ecosystem functioning
experiments. – Am. Nat. 174: 836–849. https://doi.org/10.1086/647931

Cadotte, M. W. 2013. Experimental evidence that evolutionarily diverse
assemblages result in higher productivity. – Proc. Natl Acad. Sci. USA 110:
8996–9000. https://doi.org/10.1073/pnas.1301685110

Connolly, J., Cadotte, M. W., Brophy, C., Dooley, Á., Finn, J., Kirwan, L.,
Roscher, C. and Weigelt, A. 2011. Phylogenetically diverse grasslands are
associated with pairwise interspecific processes that increase biomass. –
Ecology 92: 1385–1392. https://doi.org/10.1890/10-2270.1

Ives, A. R. and Helmus, M. R. 2011. Generalized linear mixed models for
phylogenetic analyses of community structure. – Ecol. Monogr. 81: 511–525.
https://doi.org/10.1890/10-1264.1

Kirwan, L., Connolly, J., Finn, J. A., Brophy, C., Lüscher, A., Nyfeler, D.
and Sebastià, M. T. 2009. Diversity–interaction modeling: estimating
contributions of species identities and interactions to ecosystem function.
– Ecology 90: 2032–2038. https://doi.org/10.1890/08-1684.1

Srivastava, D. S., Cadotte, M. W., MacDonald, A. A. M., Marushia, R. G. and
Mirotchnick, N. 2012. Phylogenetic diversity and the functioning of
ecosystems. – Ecol. Lett. 15: 637–648.
https://doi.org/10.1111/j.1461-0248.2012.01795.x

Tsang, T. P. N. and Cadotte, M. W. 2025. Species overlap and phylogenetic
relatedness result in community statistical non-independence (and what to
do about it). – Ecol. Lett. 28: e70267. https://doi.org/10.1111/ele.70267

Venail, P. A. and Vives, M. J. 2013. Phylogenetic distance and species
richness interactively affect the productivity of bacterial communities. –
Ecology 94: 2529–2536. https://doi.org/10.1890/12-2002.1

<!-- Verified 2026-10-16 via Scholar Gateway, Consensus and Undermind
     (metadata; full text read for Cadotte 2013, Connolly et al. 2011,
     Tsang & Cadotte 2025). Initials "A. A. M." for MacDonald in Srivastava
     et al. 2012 not confirmed by the tools: check before submission. -->

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
