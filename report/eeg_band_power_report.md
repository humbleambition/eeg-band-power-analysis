Resting EEG band power in adolescents with schizophrenia-spectrum
disorders
================
Jennifer Hsieh, PhD
October 2026

- [Question](#question)
- [Data](#data)
- [Methods](#methods)
- [Results](#results)
  - [Power spectra](#power-spectra)
  - [Channel-by-band effect sizes](#channel-by-band-effect-sizes)
  - [Region-level mixed models](#region-level-mixed-models)
  - [Alpha peak frequency](#alpha-peak-frequency)
  - [Classification](#classification)
- [Interpretation](#interpretation)
- [Limitations](#limitations)
- [Reference](#reference)
- [Session info](#session-info)

## Question

Is resting, eyes-closed EEG spectral power different in adolescents with
schizophrenia-spectrum disorders compared with healthy adolescents, in
which frequency bands and where on the scalp, and how well do these
differences separate the two groups out of sample?

## Data

The MSU adolescent schizophrenia EEG database (Borisov et al., 2005): 39
healthy schoolboys and 45 boys with schizophrenia-spectrum diagnoses
(ICD-10 F20, F21, F25), aged about 11 to 14, unmedicated at the time of
recording. Each record is 60 s of wakeful, eyes-closed EEG from 16 scalp
electrodes (10-20 system, linked-ear reference), sampled at 128 Hz and
supplied as artefact-free segments. Source:
<http://brain.bio.msu.ru/eeg_schizophrenia.htm>.

## Methods

**Cleaning and QC.** Each channel was mean-centred. Recordings were
checked for missing values and flat channels (none). The 60 s were split
into 29 overlapping 4-s segments, and a segment was rejected if any
channel’s peak-to-peak amplitude exceeded that channel’s within-subject
median + 5 MAD. Median rejected segments per subject: 0 (max 4 of 29),
so no subject met the pre-specified exclusion rule (\> 25% rejected). 7
channel recordings had outlying overall variance (\|robust z\| \> 3.5,
from 2 subjects); they were kept, because log band power absorbs gain
differences and the data were supplied as clean.

**Spectra and features.** Welch PSD per channel (4-s Hann windows, 50%
overlap, linear detrend, 0.25 Hz resolution). Band power was integrated
over delta (1-4 Hz), theta (4-8), alpha (8-13) and beta (13-30), and
expressed as **log10 absolute power** and **relative power** (share of
1-30 Hz). Individual **alpha peak frequency** was the maximum of the
mean posterior spectrum (O1, O2, P3, P4, Pz) between 7 and 13 Hz,
refined by parabolic interpolation; 2 adolescents whose maximum sat on
the 7 Hz edge (no identifiable peak) were left out of that test.

**Statistics.** For each of 16 channels x 4 bands x 2 measures: Welch
t-test, with a Wilcoxon rank-sum test as a robustness check, and Hedges’
*g* with 95% CI (schizophrenia minus control). P-values were corrected
with Benjamini-Hochberg FDR within each measure (64 tests). For a
region-level view, a linear mixed model per band,
`log10 power ~ group * region + (1 | subject)` (sum contrasts, Type III,
Satterthwaite df; `lmerTest`).

**Classifier.** Elastic-net logistic regression (`glmnet`, alpha = 0.5)
on all 128 band-power features, with lambda chosen by inner 5-fold CV
inside each outer training fold (nested CV; standardisation learned on
the training fold only). Performance: 10 x repeated stratified 5-fold
CV. Baseline: logistic regression on one feature, mean relative alpha
across channels. Significance: 200 label permutations of the full
pipeline.

## Results

### Power spectra

<img src="../figures/fig1_power_spectra.png" width="100%" style="display: block; margin: auto;" />

The healthy group shows the expected posterior alpha peak; in the
schizophrenia group it is lower and shifted to slower frequencies, while
delta and theta power look similar.

### Channel-by-band effect sizes

<img src="../figures/fig2_effect_size_heatmap.png" width="100%" style="display: block; margin: auto;" />

| Measure      | Band  | Channels FDR q\<.05 (t-test) | Channels FDR q\<.05 (Wilcoxon) | Hedges g range |
|:-------------|:------|-----------------------------:|-------------------------------:|:---------------|
| Log absolute | delta |                            0 |                              0 | -0.13 to 0.35  |
| Log absolute | theta |                            0 |                              0 | -0.27 to 0.29  |
| Log absolute | alpha |                           12 |                             14 | -1.02 to -0.15 |
| Log absolute | beta  |                            1 |                              3 | -0.60 to -0.09 |
| Relative     | delta |                           13 |                              7 | 0.37 to 1.05   |
| Relative     | theta |                            9 |                              4 | 0.20 to 0.82   |
| Relative     | alpha |                           15 |                             15 | -1.25 to -0.41 |
| Relative     | beta  |                            0 |                              1 | -0.29 to 0.45  |

- **Alpha is reduced** in schizophrenia: lower absolute alpha at 12 of
  16 channels and lower relative alpha at 15 of 16. The largest effect
  is relative alpha at O1 (g = -1.25, 95% CI -1.72 to -0.78).
- **Absolute delta and theta do not differ** (all \|g\| \< 0.36, none
  significant), so the higher *relative* delta and theta in
  schizophrenia mostly reflects the smaller alpha share, not extra
  slow-wave power. The relative delta/theta effects are less robust to
  the rank-based test.
- **Beta** trends lower in absolute terms (one channel, P3, survives
  FDR).

<img src="../figures/fig3_scalp_effect_maps.png" width="100%" style="display: block; margin: auto;" />

### Region-level mixed models

| Band  | Term         | F     | df      | p       | FDR q (group) |
|:------|:-------------|:------|:--------|:--------|:--------------|
| delta | group        | 0.13  | 1, 82   | 0.719   | 0.766         |
| delta | group:region | 2.99  | 4, 1252 | 0.018   |               |
| theta | group        | 0.09  | 1, 82   | 0.766   | 0.766         |
| theta | group:region | 4.94  | 4, 1252 | 5.9e-04 |               |
| alpha | group        | 15.29 | 1, 82   | 1.9e-04 | 7.6e-04       |
| alpha | group:region | 10.04 | 4, 1252 | 5.3e-08 |               |
| beta  | group        | 3.88  | 1, 82   | 0.052   | 0.104         |
| beta  | group:region | 4.99  | 4, 1252 | 5.4e-04 |               |

Absolute alpha power by region (region-averaged log10 power)

| Region    | Power ratio SCZ/CTL | Hedges g \[95% CI\]    | FDR q   |
|:----------|:--------------------|:-----------------------|:--------|
| Frontal   | 0.59                | -0.96 \[-1.42, -0.51\] | 3.9e-04 |
| Temporal  | 0.73                | -0.46 \[-0.89, -0.02\] | 0.100   |
| Central   | 0.62                | -0.76 \[-1.20, -0.31\] | 0.005   |
| Parietal  | 0.54                | -0.73 \[-1.18, -0.29\] | 0.005   |
| Occipital | 0.46                | -0.88 \[-1.33, -0.43\] | 8.9e-04 |

Only alpha shows a main effect of group after FDR. Group x region
interactions are significant in every band, i.e. the group difference is
not uniform over the scalp: for alpha it is largest over occipital and
frontal sites, where the schizophrenia group has roughly 45-60% of
control alpha power.

<img src="../figures/fig4_relative_power_by_region.png" width="100%" style="display: block; margin: auto;" />

### Alpha peak frequency

<img src="../figures/fig5_alpha_peak_frequency.png" width="80%" style="display: block; margin: auto;" />

Mean alpha peak was 9.83 Hz (SD 0.70) in controls and 9.26 Hz (SD 0.86)
in the schizophrenia group: Welch t(79.2) = -3.28, p = 0.002, g = -0.71
\[-1.16, -0.26\].

### Classification

<img src="../figures/fig6_classifier.png" width="100%" style="display: block; margin: auto;" />

| Model                           | AUC mean (range over repeats) | Balanced accuracy | Permutation p |
|:--------------------------------|:------------------------------|:------------------|:--------------|
| Baseline: global relative alpha | 0.76 (0.74-0.77)              | 0.65              |               |
| Elastic net (128 features)      | 0.79 (0.72-0.83)              | 0.71              | 0.005         |

The elastic net reached a cross-validated AUC of 0.79, well above the
permutation null (95th percentile 0.55; none of 200 permutations reached
the observed value), but only slightly better than the single
global-relative-alpha feature (AUC 0.76). The most consistently selected
features were rel_alpha_O1, rel_alpha_O2, rel_alpha_F3,
abs_log10_alpha_F4, rel_delta_O1: alpha again carries most of the
signal.

## Interpretation

Unmedicated adolescents with schizophrenia-spectrum disorders show a
**reduced and slower alpha rhythm** at rest, with otherwise similar
absolute slow-wave power. This sits alongside the original authors’
report, in this same sample, of altered alpha-segment structure and
reduced structural synchrony of alpha activity; here the effect shows up
in plain spectral power as well. A spectral signature like this
separates groups at roughly 70% balanced accuracy, which is useful for
research but nowhere near a diagnostic test.

## Limitations

- Small, single-site, all-male sample from one 2005 database; no age, IQ
  or symptom data per subject, so age (which strongly shapes alpha
  frequency in adolescence) cannot be modelled directly. The authors
  report matched mean ages (12.3 years in both groups).
- 60 s per subject; recordings were pre-selected as artefact-free by the
  data providers.
- The amplitude scale in the files is larger than typical scalp EEG in
  microvolts, so absolute power is reported on a log scale and
  interpreted only relative to the other group; relative power and peak
  frequency are unaffected by gain.
- 50 Hz is barely visible in the spectra, which suggests the data were
  low-pass filtered before sharing; analysis was restricted to 1-30 Hz.
- No licence is stated on the source page; the data are not
  redistributed here.

## Reference

Borisov SV, Kaplan AYa, Gorbachevskaya NL, Kozlova IA (2005). Analysis
of EEG structural synchrony in adolescents with schizophrenic disorders.
*Human Physiology* 31(3):255-261.
<https://doi.org/10.1007/s10747-005-0042-z>

## Session info

    ## R version 4.5.0 (2025-04-11)
    ## Platform: aarch64-unknown-linux-gnu
    ## Running under: Debian GNU/Linux 13 (trixie)
