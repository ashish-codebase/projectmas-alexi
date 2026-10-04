# Paper Summary — Hain & Anderson (2017), *GRL* 44, 9723–9733

> **Source PDF:** `Geophysical Research Letters - 2017 - Hain - Estimating morning change in land surface temperature from MODIS day night (1).pdf`
> **Extracted text:** `_hain2017_text.txt`
> **DOI:** 10.1002/2017GL074952 · Received 24 Jul 2017 · Accepted 10 Aug 2017 · Published 2017 (US Government work, public domain)
> **Authors:** C. R. Hain (NASA MSFC, Earth Science Office) and M. C. Anderson (USDA-ARS Hydrology Remote Sensing Laboratory)
> **Support:** NASA Applied Sciences Water Resources Program, grant NNX12AK90G

---

## 1. What the paper is about

- **Problem.** The Atmosphere–Land Exchange Inverse (**ALEXI**) surface energy-balance model requires `LST_ALEXI` — the *mid-morning rise* in land surface temperature over a window from ~1.5 h after sunrise to ~1.5 h before local noon. That quantity is only directly observable from **geostationary** sensors (GOES, ≤1 h cadence), which cannot serve the globe:
  - 4+ geostationary data streams needed to mosaic 60°N–60°S;
  - each sensor has different (and changing) thermal retrieval algorithms, plus spatial/temporal/spectral differences;
  - oblique view angles prevent accurate LST beyond ±60° latitude.
- **Core idea.** Polar-orbiting sensors (MODIS / AVHRR / VIIRS) give only **two** LST looks per day (one night, one day). A physically constrained diurnal-cycle fit (e.g. Gottsche–Olesen DTC, Holmes 2013, Inamdar 2008) needs ≥4 observations/day. Instead, **data-mine** the day/night pair to predict `LST_ALEXI` directly from a **single** sensor.
- **Why single-sensor beats multi-sensor.** ~25% higher probability of a clear day–night acquisition pair from one sensor (Aqua MODIS) than of 4 clear retrievals in a day from two sensors (Terra + Aqua); less operational risk from sensor failure; far less data volume and processing.

## 2. Methodology

- **Algorithm.** Rule-based regression-tree data mining — **RuleQuest Cubist**, deliberately capped at **20 rules per model** to curb overfitting.
- **Training target.** GOES-observed `LST_ALEXI` over 2007–2014, treated as ground truth.
- **Predictor variables.** Ten candidates tested; **five retained** because they showed the strongest relationship to the training data:
  1. MODIS day–night LST difference,
  2. MODIS daytime LST,
  3. MODIS nighttime LST,
  4. leaf area index (LAI),
  5. topographic variability (1 km GTOPO DEM).
  - Dropped as weak: land cover classification (UMD global land cover), emissivity (ASTER Global Emissivity Dataset), surface albedo (MOD43C), ALEXI window length, mean annual precipitation.
- **Model stratification.** Separate models per regime of **maximum vegetation fraction** `fcmax` (0–25, 25–50, 50–75, 75–100 %) × **mean annual precipitation** (0–500, 500–1000, 1000–5000 mm/yr) → **12 sub-models**.
  - `fc` derived by inverting LAI with a canopy gap-fraction method (Norman 1995); `fcmax` = maximum sampled over 2003–2015 per pixel.
  - Rationale: `fcmax` and `P` inject climate/land-surface information, letting the tree find better relationships within subsets.
- **Data pipeline.**
  - *GOES:* 11 µm brightness temperature inverted to LST (Price 1983 algorithm), single-channel, GOES-East 75°W / GOES-West 105°W; atmospheric temperature + water-vapour profiles from the NOAA **CFSR reanalysis**; cloud-contaminated retrievals flagged and rejected via visible reflectance plus a **bispectral composite threshold** on the 3.7 µm and 11 µm Tb (Jedlovec 2008); 30-min temporal resolution; native 4 km at nadir → **aggregated to 8 km** to include high-latitude Canada; then temporally interpolated per pixel over the ALEXI window.
  - *MODIS:* twice daily (1:30 AM/PM local solar) LST from standard **Aqua** products, generalized **split-window atmospheric compensation** (Wan 2002/2004/2007); extracted from `MYD11C1` L3 Global 0.05° Climate-Modeling-Grid LST; resampled to 8 km; product QC flags applied; pixels with satellite view angle > 60° masked.
  - All training/predicator layers remapped to a North America domain at **8 km**.
- **Quality control.** Product-level QA plus extra filtering: extreme outliers removed; physically impossible values removed (`LST_ALEXI` < 0 or day–night LST difference < 0 → symptom of undetected cloud). Analysis restricted to **snow-free** pixels via a dynamic mask from the 16 km NOAA **IMS** snow-cover product (Ramsey 1998).
- **Validation protocol.** **10-fold cross-validation** over the 2007–2014 period; predictions generated on a daily time step; **n = 226,419,494** pixel-day observations compared.
- **Second-order test (utility, not just error).** ALEXI run over **2003–2015** at 8 km over CONUS twice — once with GOES-observed, once with MODIS-predicted `LST_ALEXI` — producing **Evaporative Stress Index (ESI)** time series (Anderson 2011a convention: standardized anomalies of the clear-sky ratio `fRET = ET / ETref`, with `ETref` = FAO-96 Penman–Monteith reference ET for grass, which normalizes out non-moisture drivers such as seasonal radiation load). Compared against an independent reference: **NLDAS-2 Noah LSM** soil moisture, layer-weighted over the first three layers (0–10, 10–40, 40–100 cm) to a 0–100 cm estimate, converted to normalized σ-anomalies over the same period of record, resampled to the 8 km grid. Agreement quantified as spatial anomaly correlation *r* per warm-season month (May–Sep) per year.

## 3. Findings

- **Agreement is high.**
  - Coefficient of determination **R² = 0.96** overall.
  - Seasonal R²: 0.91 (Jan–Mar), 0.94 (Apr–Jun), 0.94 (Jul–Sep), 0.92 (Oct–Dec) — no seasonal collapse in skill.
- **Error statistics.**
  - Domain-averaged **bias ≈ 0** (±0.01 K); most regions within ±0.5 K.
  - **Normalized RMSE ≈ 9%** (absolute RMSE = 1.84 K) domain-averaged; seasonal subsets 8–15%; spatially 5–20% (conclusions quote 8–20%).
  - Absolute RMSE is smaller in cool seasons, larger in warm seasons — it tracks mean `LST_ALEXI`.
- **Physical magnitudes behave as expected.** `LST_ALEXI` ≈ **3–9 K** over wet/highly vegetated regions vs **15–27 K** over semiarid–arid western US and Mexico — consistent with the "drier surfaces heat faster under rising morning solar load" signal ALEXI exploits.
- **Downstream equivalence.** Warm-season spatial anomaly correlation with Noah soil moisture: **0.63 (GOES)** vs **0.61 (MODIS-predicted)** averaged 2003–2015 — **no statistically significant difference**. Per-year MODIS is slightly lower in 2004–2008, 2010, 2011, 2013–2014; slightly higher in 2003, 2009, 2012, 2015 — all within noise.
- **Conclusion drawn by the authors.** Single-sensor `LST_ALEXI` has accuracy comparable to current geostationary products and can complement (not replace) diurnal-fit methods that require >2 observations/day. They recommend a **hybrid**: exploit multiple daily observations where available, and fill spatiotemporal gaps with the data-mining approach where only two retrievals are possible per day.

## 4. Limitations

### 4.1 Stated by the authors

- Cold-season discrimination between a cold land surface and low cloud degrades **both** GOES and MODIS, especially for **nighttime** MODIS observations.
- Worst nRMSE regions are attributed to **inadequate snow-cover screening over Canada** in the cold season and **inadequate cloud screening over Canada and South America**.
- Biases of 1–2 K in the western US coincide with **complex terrain** — plausibly surface-emissivity error plus differences in **satellite view-angle geometry** between GOES and MODIS. Largest absolute biases occur where `LST_ALEXI` is already high, so *relative* bias is not much worse than elsewhere.
- GOES (single-channel) and MODIS (split-window) use different atmospheric-correction algorithms; the authors argue this is tolerable because the target is a **time differential**, but this is asserted, not demonstrated.

### 4.2 Methodological — worth flagging when reusing this work

- **Empirical, not physical.** No energy-budget or diurnal-curve constraint; skill is interpolation inside the 2007–2014 training regimes. Extrapolation to unseen climate/land-surface regimes carries unbounded risk.
- **Ground truth is GOES.** Reported accuracy is *relative to GOES*, not absolute truth; shared-mode errors and GOES bias are invisible to the metric.
- **Static predictors.** `fcmax` is a *maximum over 2003–2015*; LAI/albedo/GTOPO are climatologies. Intra-season phenology, crop growth stage, planting, irrigation, and disturbance are not represented.
- **No self-detection of bad predictions.** With only two looks/day there is no internal redundancy to flag a day where cloud contamination slipped past QC.
- **Coverage ≠ accuracy.** Snow masking collapses high-latitude coverage in cold seasons; cloud gaps leave holes. The "global" claim is an **enabling** statement — validation is over North America (plus partial Canada/South America), not a demonstrated global product.
- **Coarse model** (20 rules, 12 regimes) trades fidelity for robustness; the actual rule sets live in **Supporting Information S1**, so the method is not reproducible from the letter alone.
- **Reference-data circularity.** ESI-vs-Noah is a model-vs-model comparison; Noah LSM has its own biases, and the 0.63/0.61 equivalence partly reflects a shared climate signal rather than proving equal LST quality.
- **Venue constraints.** *Geophysical Research Letters* is a short-form letter: no per-pixel uncertainty field, no cost/throughput analysis for operational deployment, limited methodological detail.

## 5. Relevance to this project

- This paper specifies the **ALEXI input problem** that `projectmas-alexi-master` works on: a defensible way to feed ALEXI from MODIS day/night pairs when GOES-class temporal sampling is unavailable.
- Its caveats map directly onto failure modes to gate in a product pipeline: **snow cover**, **low cloud / fog** (especially at night), **complex terrain / emissivity**, and **static vegetation fraction**.
- Key transferable design points: stratify by `fcmax` and precipitation regime; keep the model deliberately low-order; validate with k-fold cross-validation against a geostationary reference; and evaluate downstream (ESI vs an independent soil-moisture reference), not just pixel error.

## 6. Key numbers (quick reference)

| Quantity | Value |
| --- | --- |
| ALEXI sampling window | 1.5 h after sunrise → 1.5 h before local noon |
| MODIS acquisition times | 1:30 AM / PM local solar (Aqua) |
| Spatial resolution (all layers) | 8 km |
| Training / validation period | 2007–2014 (10-fold cross-validation) |
| Observations compared | 226,419,494 |
| Sub-models | 12 (`fcmax` × precipitation regimes), ≤20 rules each |
| R² (overall / seasonal) | 0.96 / 0.91–0.94 |
| Bias (domain average) | ±0.01 K |
| RMSE (absolute / normalized) | 1.84 K / ~9% (spatially 5–20%) |
| `LST_ALEXI` range, wet–vegetated | 3–9 K |
| `LST_ALEXI` range, semiarid–arid | 15–27 K |
| ESI↔Noah correlation, GOES / MODIS | 0.63 / 0.61 (no significant difference) |
| ESI / Noah comparison period | 2003–2015, warm season May–Sep |

---

## Appendix A — RuleQuest Cubist (the data-mining engine used above)

> Sources: RuleQuest documentation (`cubist-unix.html`, `cubist-info.html`, `cr210.html`, `cubist-previous.html`, `download.html`), the R `Cubist` vignette (topepo/CRAN), and the tidymodels `rules` Cubist-vs-RuleFit comparison.

### A.1 What it is

- **Cubist** is a RuleQuest Research tool from Ross Quinlan's group — the numeric-prediction sibling of See5/C5.0 (which predicts categories). Runs on Windows 8/10/11 and Linux; a **single-threaded Linux C source is released under GNU GPL (Release 2.07 GPL Edition)**; the demo build is capped at 200 cases.
- **Built for scale:** designed for "hundreds of thousands to millions of records and tens to thousands of numeric or nominal fields", using up to 8 cores — which is why it could ingest the paper's **226 million** pixel-day samples.
- **Model form:** a **rulebook** — a set of `if <conditions> then <multivariate linear model>` rules. A case matching a rule's conditions is scored by that rule's linear model. Chosen for *intelligibility*: stronger than plain multivariate linear regression, far more readable than a neural network.

### A.2 Algorithm stages

- **Base:** an extension of Quinlan's **M5 model tree** (Quinlan 1992), plus nearest-neighbour corrections from **Quinlan 1993a**.
- **Grow a model tree:** terminal leaves hold **linear regression models** built only from the predictors used in the splits above them; a linear model is also saved at **every intermediate node**.
- **Smoothed prediction:** the leaf model's value is blended recursively with the parent models up to the root, so a prediction is a linear combination of *all* models along the path. (Consequence: the variable-usage statistics printed by `summary()` do not correspond to the printed terminal models.)
- **Tree → rules:** each root-to-leaf path becomes a rule, then is **simplified** (Quinlan 1987): subsumed conditions removed (`A<10 & A<7 → A<7`), and conditions that do not improve performance are **pruned**.
- **Feature selection** runs inside every linear model, so a rule's model may use fewer variables than its conditions mention; the coefficients attached to a rule are blends of the ancestor models.
- **Optional add-ons** (part of the tool, not evidenced in the paper):
  - **Composite models** — find the *n* most similar training cases and combine their known values, their model-predicted values, and the model's prediction for the new case.
  - **Committee models** — a boosting-like ensemble of model trees; later trees fit residuals via adjusted pseudo-outcomes, and the final prediction is a **simple average** (no stage weights).
  - **Unbiased rules** — force each rule's mean prediction to equal its training mean, at slightly higher mean absolute error; recommended when the target has a preponderant single value.

### A.3 Paper settings ↔ Cubist options

| Paper statement | Cubist option | Effect |
| --- | --- | --- |
| "Cubist was limited to developing 20 rules for each model" | `-r 20` | Max rules in a model (**default 500**; 100 in older releases). The docs present `-r` as the *simplicity–accuracy trade-off* control, so 20 rules × 12 regimes is a deliberate, very coarse, interpretable model. |
| "a f-fold (10 folds used) cross-validation technique" | `-X 10` | Cubist's **built-in** cross-validation: cases are split into *f* blocks matched for size and target distribution; a model is built per fold and scored on the held-out block, so every case is tested exactly once. **The paper's validation protocol is the tool's own option, not an independent one.** |
| (unstated, but active by default) | `-e 5` | **Extrapolation clamp.** Each rule records the min/max target of its training cases; a computed value outside that range is snapped to the nearer bound, extended by 5% of the range. Predictions are therefore effectively **bounded by the 2007–2014 training envelope** — a real constraint on the paper's claim of predicting unseen years. |
| (would be needed) | names file | `continuous` / ordered-discrete declarations, derived attributes (`x := expr`), `attributes included:/excluded:`, `case weight`, `?` = missing, `N/A` = not applicable. |
| (not used, per the paper) | `-i` / `-a` / `-n 1–9`, `-C 5`, `-u`, `-S`, `-I` | composite models, nearest neighbours, committees, unbiased rules, sampling. |

### A.4 Version history (relevant to a 2017 paper)

- **2.03** — model complexity controlled *only* by a rule-count limit (default 100); the "minimum case cover" parameter dropped; per-case weighting added.
- **2.02** — faster composite models; multi-core / hyper-thread support.
- **2.07** — faster composite models on large data, higher-precision neighbour distances, multithreaded public model-reading code; this is the **GPL edition** still shipped today.
- **2.09** — targeted at millions of cases: one 30 M-case job ran **30× faster** (22 min vs >11 h) with 20% less memory; revised **Aug 2016**.
- **2.10** — improved composite models and nearest-neighbour search; documentation banner dated **Jun 2019**.
- ⇒ The 2017 paper almost certainly ran **2.08/2.09**. The exact build and the rulebooks are **not stated**; the rule sets live only in **Supporting Information S1**, so the method is **not reproducible** from the article alone. (SI Figure S2 shows the same method at 1 km MODIS; `LST_ALEXI` fields are archived at NASA MSFC, FTP on request from the author.)

### A.5 Using it today

- **R:** `Cubist` on CRAN (topepo) — a port of the GPL 2.07 C code, with `caret`/`tidymodels` integration.
- **Python:** `cubist` on PyPI (pjaselin) — wrapper around the Cubist **v2.07** C code, Python ≥3.10, GPL-3.0, scikit-learn-style plotting.
- **C:** RuleQuest's **public C source** reads/interprets `.model` files (ASCII since release 1.08, cross-platform) for embedding models in your own systems — the natural target for the `csharp/` conversion work.
- **Caveats:** no support is offered for GPL releases; commercial use of the Windows/Linux binaries requires the licence agreement; the R port omits predictor binning, non-numeric/non-categorical attribute types, and the "let the C code decide about instances" option.

### A.6 Modern equivalents of the Cubist approach

| Goal | Nearest equivalent | Notes |
| --- | --- | --- |
| Drop-in, same algorithm | **M5P** — Weka `weka.classifiers.trees.M5P`, R `RWeka::M5P` / mlr3 `regr.m5p`, Python `m5py` (`smarie/python-m5p`, scikit-learn API) | Wang & Witten (1997) "rational reconstruction" of Quinlan 1992: model trees with linear leaves. Same family, no committees/instances. |
| Rules as *features* in one global model | **RuleFit** (Friedman & Popescu 2008) — R `xrf` (XGBoost builds rules, `glmnet` fits a sparse linear model), tidymodels `rules` | Opposite composition: one linear model over raw features **plus** rule-derived binary features, rather than one model per rule. |
| Interpretable additive model | **GAM** — R `mgcv`, Python `pygam` | Per-feature spline/piecewise-linear terms combined linearly; a natural fit for this problem (day–night difference, LAI, topography as smooth terms). |
| GAM with interactions, scales to many features | **CAT / TAM** (concept-based Taylor additive models, *SIGKDD* 30) | Addresses GAM's parameter count and overfitting on real-world tabular data. |
| Ensemble that is inherently explainable | **shallow-tree ensembles** (arXiv 2410.19098, "Inherently Interpretable Tree Ensemble Learning") | Reframes a boosted-forest as an additive model over shallow trees. |
| Neural analogue of "a linear model per region" | **Mesomorphic NNs / hypernetworks** (NeurIPS 2024) | Trains a network that *generates* a linear model in the original feature space per sample — closest learned-model counterpart to a Cubist rulebook. |
| Plain black-box baseline | XGBoost / LightGBM | Better raw accuracy at scale, but no rulebook and known to under-perform in the tails; model trees do not have that defect because their leaves regress over the full range. |
| LLM / foundation model | **not** an equivalent | No calibrated numeric output, no reproducible rulebook, no guarantee of monotone/physical behaviour; usable at most as a semantic feature-augmentation layer, not as the regressor. |
