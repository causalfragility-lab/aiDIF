# aiDIF: Differential Item Functioning for AI-Scored Assessments

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)]()
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

**aiDIF** extends the robust DIF framework from
[robustDIF](https://github.com/phalpin/robustDIF) to address a modern
measurement challenge: **does AI scoring introduce subgroup-dependent item
bias?**

As AI systems increasingly score essays, short answers, and structured
responses in educational and psychological assessments, a critical fairness
question arises: does the AI scoring engine shift item difficulties
differently for different demographic groups — even when no human-scoring DIF
exists?

`aiDIF` provides:

- **Robust DIF analysis** under both human and AI scoring conditions
  (inheriting M-estimation from `robustDIF`)
- **Differential AI Scoring Bias (DASB) test** — a new Wald test for
  item-level scoring shifts that differ across groups
- **AI-effect summaries** classifying items as: *stable*, *masked*,
  *introduced*, or *new_direction* across scoring conditions
- **Anchor weight diagnostics** under potential AI contamination
- **Simulation utilities** for benchmarking DIF methods in AI-scored settings

---

## Installation

```r
# Install robustDIF dependency first (from CRAN or GitHub)
# install.packages("robustDIF")

# Install aiDIF from GitHub (once published)
# devtools::install_github("yourname/aiDIF")

# Or install from local source
devtools::install_local("path/to/aiDIF")
```

---

## Quick Start

```r
library(aiDIF)

# Use the built-in example dataset
# aidif.eg contains paired human/AI scoring parameters for 6 items, 2 groups.
# Item 1 has DIF in human scoring; Item 3 has DASB (AI-induced group-specific bias).
data(aidif.eg)

# Fit the model
mod <- fit_aidif(
  human_mle = aidif.eg$human,
  ai_mle    = aidif.eg$ai
)

# Compact summary
print(mod)

# Full report
summary(mod)

# Visualisations
plot(mod, type = "dif_forest")  # Forest plot: human vs AI DIF estimates
plot(mod, type = "dasb")        # Bar chart of DASB with error bars
plot(mod, type = "weights")     # Anchor weights in each scoring condition
plot(mod, type = "rho")         # Bi-square objective for human scoring
```

---

## Core Concepts

### Differential AI Scoring Bias (DASB)

For item *i* and group *g*, define the **scoring shift**:

```
δ_ig = d_ig^AI − d_ig^Human
```

where *d_ig* is the IRT intercept (difficulty) parameter. The DASB is:

```
DASB_i = δ_i2 − δ_i1
```

Under H₀: DASB_i = 0, a Wald test is conducted using the asymptotic
variance derived from the delta method (assuming independent groups and
scoring conditions):

```
Var(DASB_i) = σ²(d_i1^H) + σ²(d_i2^H) + σ²(d_i1^AI) + σ²(d_i2^AI)
```

A significant result means the AI scoring engine does **not** merely
re-scale all items uniformly — it disadvantages (or advantages) one group at
specific items.

### AI-Effect Classification

`ai_effect_summary()` compares DIF flagging patterns between scoring
conditions:

| Status | Meaning |
|---|---|
| `stable_clean` | Not flagged in either condition |
| `stable_dif` | Flagged in both (same direction) |
| `introduced` | Flagged only under AI scoring |
| `masked` | Flagged only under human scoring |
| `new_direction` | Flagged in both, but bias reverses sign |

---

## From Existing IRT Fits

If you have fitted IRT models in `mirt` or `lavaan`, use
`get_ai_parms()` to extract parameters:

```r
library(mirt)

# Fit multigroup 2PL under human scoring
human_fit <- mirt(human_data, model = 1, itemtype = "2PL",
                  group = "group", SE = TRUE)

# Fit multigroup 2PL under AI scoring
ai_fit <- mirt(ai_data, model = 1, itemtype = "2PL",
               group = "group", SE = TRUE)

# Extract and bundle
dat <- get_ai_parms(human_fit, ai_fit)

# Fit aiDIF model
mod <- fit_aidif(dat$human, dat$ai)
```

---

## Simulation

```r
# Generate synthetic data with known DIF and DASB
dat <- simulate_aidif_data(
  n_items    = 10,
  n_obs      = 500,
  impact     = 0.5,      # 0.5 SD group mean difference
  dif_items  = c(1, 2),  # items with human-scoring DIF
  dif_mag    = 0.5,
  dasb_items = 5,        # item with AI-induced differential bias
  dasb_mag   = 0.4,
  ai_drift   = 0.1       # uniform AI calibration offset
)

mod <- fit_aidif(dat$human, dat$ai)
summary(mod)
```

---

## Example Dataset

`aidif.eg` is a list with elements `$human` and `$ai`, each formatted
identically to the output of `robustDIF::get_model_parms()`:

- **6 items**, 2 groups (reference and focal)
- **Item 1**: DIF in human scoring (intercept bias +0.5 in focal group)
- **Item 3**: DASB — AI scoring adds +0.4 to focal group intercept only
- **Impact**: 0.5 SD (focal group higher on latent trait)

---

## Package Architecture

```
aiDIF/
├── R/
│   ├── read_functions.R    # read_ai_scored(), get_ai_parms()
│   ├── aidif_core.R        # fit_aidif() — main estimation engine
│   ├── scoring_bias.R      # scoring_bias_test(), ai_effect_summary(),
│   │                       #   anchor_weights()
│   ├── simulate.R          # simulate_aidif_data()
│   ├── validate_inputs.R   # Internal validation helpers
│   └── class_functions.R   # print/summary/plot S3 methods
├── data/
│   └── aidif.eg.rda        # Built-in example dataset
├── tests/
│   └── testthat/
│       └── test-aidif.R    # Unit tests
└── DESCRIPTION
```

---

## Relationship to robustDIF

`aiDIF` depends on `robustDIF` and reuses its core infrastructure:

| `robustDIF` | `aiDIF` usage |
|---|---|
| `get_model_parms()` | Called internally by `get_ai_parms()` |
| `rdif()` | Called twice in `fit_aidif()` (human + AI) |
| `dif_test()` | Called for both scoring conditions |
| `delta_test()` | Stored in each `rdif` sub-object |

The novel contribution of `aiDIF` is the **DASB test** and the
AI-effect summary infrastructure.

---

## Citation

If you use `aiDIF` in published research, please cite both this package and
the underlying `robustDIF` methodology:

```
[Author] (2025). aiDIF: Differential Item Functioning for AI-Scored
Assessments. R package version 0.1.0.

Halpin, P., Nickodem, K., & Eagle, J. (2024). robustDIF: Differential Item
Functioning Using Robust Scaling. R package version 0.2.0.
```

---

## License

GPL (>= 3)
