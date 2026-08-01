# Colocalisation analysis: state-of-the-art approaches

Survey of methods used by ImageJ/Fiji and other established colocalisation
tools, done to scope the design of this module. 2026-08-01.

## Two main families, plus a newer third

1. **Pixel/intensity-based** — treat each channel as an intensity image and
   correlate intensities across all pixels/voxels, with no segmentation step.
2. **Object-based** — segment each channel into discrete objects first, then
   compare objects (overlap fraction, centroid/nearest-neighbour distance).
3. **Point-pattern / spatial-statistics based** — treat detected objects or
   localisations as point sets and test whether their cross-channel spatial
   arrangement departs from complete spatial randomness (CSR). This is the
   newest family and the one closest to code AnalyzER already has.

Which family is appropriate depends on the biology: diffuse/overlapping
compartments suit pixel-based metrics; discrete puncta or organelles suit
object-based or point-pattern metrics.

## Pixel-based methods

- **Pearson's Correlation Coefficient (PCC)** — linear correlation of the two
  channels' intensities, −1..+1. Reference metric in essentially every tool
  (Coloc2, JACoP, CellProfiler, Imaris, Huygens). Insensitive to absolute
  intensity/gain differences, but blind to spatial structure and easily
  distorted by background, noise, and channel bleed-through — background
  subtraction is a prerequisite, not optional.
- **Manders' Coefficients (M1, M2)** — fraction of channel A's signal
  overlapping thresholded channel B, and vice versa; asymmetric and
  biologically interpretable ("38% of A's signal overlaps B"). The most
  widely reported metric after PCC, but the result is only as good as the
  threshold used to define "signal" in each channel.
- **Costes' automatic thresholding** — removes the subjectivity of manual
  Manders' thresholds by finding, per channel, the highest intensity value
  below which PCC over the remaining pixels is ≈0 (i.e. that population is
  indistinguishable from uncorrelated background). Standard in Coloc2,
  CellProfiler, Imaris; should be the default rather than user-dragged
  sliders.
- **Costes' randomisation significance test** — block-scrambles one channel
  and recomputes PCC/Manders many times (Coloc2 defaults to ~100 iterations)
  to build a null distribution, then reports where the real value falls
  against it. A bare coefficient without this kind of null-model comparison
  is close to meaningless — this is the biggest gap between "quick colour
  overlay by eye" and a defensible quantitative result.
- **Li's Intensity Correlation Quotient (ICQ)** — sign-test based nonparametric
  measure; distinguishes segregated staining, random co-occurrence, and
  dependent (correlated) staining, complementing PCC.
- **van Steensel's Cross-Correlation Function (CCF)** — recomputes PCC while
  progressively shifting one channel pixel-by-pixel; a symmetric peak at
  zero shift indicates genuine spatial association, a flat or off-centre
  profile indicates coincidental overlap. Good for catching cases where PCC
  alone is misleading. A 2024 extension applies the same shift-correlation
  idea directly to super-resolution/SMLM point data rather than pixel grids.
- **Spearman's rank correlation** — nonparametric alternative to PCC when the
  intensity relationship isn't linear.

## Object-based methods

Segment each channel into discrete objects (thresholding, watershed, or
increasingly deep-learning segmentation such as Cellpose/StarDist — directly
applicable here since AnalyzER's organelle strand already uses Cellpose),
then compare across channels by:
- fractional area/volume overlap between matched objects, or
- centroid-to-centroid / boundary-to-boundary nearest-neighbour distance.

Established tools: **DiAna** (ImageJ, 3D object-based colocalisation +
distance analysis), **ComDet** (spot detection + colocalisation, built for
punctate structures), **EzColocalization** and **ComDet** (both handle >2
channels), CellProfiler's object-based `MeasureColocalization`/`RelateObjects`
modules, Imaris's spot/surface colocalisation (distance-transform based).

Advantage: gives a spatially interpretable answer ("N% of mitochondria
contact ER") and tolerates differing label stoichiometry between channels
much better than pixel correlation. Disadvantage: result quality is capped by
segmentation quality, and it's a poor fit for diffuse, non-punctate signal.

## Point-pattern / distance-based methods (newest family)

Rather than binary overlap, treat each channel's detected objects/localisations
as a point pattern and ask whether the *cross-channel* spatial relationship
(nearest-neighbour distance distribution, cross pair-correlation, cross-Ripley's
K/L) differs from what's expected under CSR or from a matched random-label
permutation null. This generalises cleanly to single-molecule
localisation microscopy, where images are point lists rather than intensity
matrices and classic pixel-based methods don't apply directly.

**SPACE** (Spatial Pattern Analysis using Closest Events, 2024) is a recent
example, reporting better sensitivity than pixel/object metrics for punctate
data by working directly on nearest-event distances rather than forcing a
threshold/overlap decision.

This family is conceptually identical to machinery already in this codebase:
`analyzerRipleyL.m` (NetworkCommon, cross-type Ripley's L clustering) and the
ray-cast nearest-distance functions in OrganelleDistances_sandbox
(`organelleD2OrganelleCompute.m` etc.) already compute cross-object nearest
distances with path-obstruction awareness. A colocalisation point-pattern
metric would mostly be a null-model/significance layer on top of distance
distributions these repos can already produce, not a new algorithm from
scratch.

## Practical requirements common to all state-of-the-art tools

- **Registration first.** Sub-pixel chromatic misalignment between channels
  destroys pixel-based metrics before any statistics are computed — channel
  registration/chromatic-aberration correction must happen upstream of any
  correlation step.
- **Background handling.** Consistent background subtraction (rolling-ball or
  otherwise) before correlating; raw-camera-noise floors bias PCC toward
  spurious positive correlation.
- **Prefer automatic, reproducible thresholding** (Costes) over manual
  threshold-dragging, which is a common source of irreproducible published
  colocalisation numbers.
- **Always pair a coefficient with a significance test** (Costes
  randomisation for pixel-based; label/point permutation for object- or
  point-pattern-based) — this is the single most common gap between casual
  and rigorous colocalisation analysis in the literature (Bolte &
  Cordelières 2006; Dunn et al. 2011 are the standard citable guides for this).
- **3D-native where possible.** Many older tools are 2D-slice-only; both the
  object-based and point-pattern families increasingly operate natively in
  3D/4D, which matches how this project's other analysis modules
  (OrganelleDistances_sandbox) are already built.

## Suggested starting scope for this module

Given AnalyzER already has organelle segmentation (Cellpose/threshold-based)
and both a ray-cast distance pipeline and Ripley's-L clustering elsewhere in
the toolset, the lowest-duplication path is:

1. **Pixel-based core**: PCC + Manders' M1/M2 with Costes automatic
   thresholding and Costes randomisation significance — this is the
   universally-expected "Coloc2-equivalent" baseline every reviewer will look
   for.
2. **Object-based extension**: reuse existing segmented-object output
   (wherever the calling strand already has masks/labels) to compute
   fractional overlap and nearest-neighbour object-to-object distance,
   sharing code with OrganelleDistances_sandbox's distance primitives rather
   than reimplementing them.
3. **Point-pattern/CCF layer**: van Steensel CCF and/or a cross-Ripley's-L
   significance layer as a phase-2 addition, since it overlaps conceptually
   with `analyzerRipleyL.m` already in NetworkCommon_sandbox.

## Sources

- [Coloc 2 — ImageJ](https://imagej.net/plugins/coloc-2)
- [Colocalization Analysis — ImageJ](https://imagej.net/imaging/colocalization-analysis)
- [DesignNotes · fiji/Colocalisation_Analysis Wiki](https://github.com/fiji/Colocalisation_Analysis/wiki/DesignNotes)
- [Development of a Novel Automated Workflow in Fiji ImageJ for Batch Analysis... (Manders Coefficient)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11986846/)
- [Colocalization by cross-correlation, a new method suited for super-resolution microscopy](https://pmc.ncbi.nlm.nih.gov/articles/PMC10837882/)
- [Statistical analysis of molecule colocalization in bioimaging — Lagache et al., Cytometry A 2015](https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22629)
- [SVI — Colocalization Coefficients in Brief](https://svi.nl/ColocalizationCoefficientsInBrief)
- [Intensity versus Object Based Colocalization — Oxford Instruments (Imaris)](https://imaris.oxinst.com/learning/view/article/intensity-versus-object-based-coloc)
- [A versatile toolbox for semi-automatic cell-by-cell object-based colocalization analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7643144/)
- [CellProfiler Colocalization Tutorial](https://cellprofiler-website.s3.amazonaws.com/outreach/content/Example_Colocalization_Tutorial.pdf)
- [DiAna, an ImageJ tool for object-based 3D co-localization and distance analysis](https://www.researchgate.net/publication/310817317_DiAna_an_ImageJ_tool_for_object-based_3D_co-localization_and_distance_analysis)
- [Spatial Pattern Analysis using Closest Events (SPACE) — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11057815/)
- [Analysis of conditional colocalization relationships and hierarchies in three-color microscopy images — JCB](https://rupress.org/jcb/article/221/7/e202106129/213216/Analysis-of-conditional-colocalization)
