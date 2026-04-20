# aiDIF: Differential Item Functioning for AI-Scored Assessments

[![R-CMD-check](https://github.com/causalfragility-lab/aiDIF/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/aiDIF/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

**aiDIF** addresses a modern measurement fairness challenge: does AI scoring
introduce subgroup-dependent item bias?

As AI systems increasingly score essays, short answers, and structured
responses in educational and psychological assessments, a critical question
arises: does the AI scoring engine shift item difficulties differently for
different demographic groups — even when no human-scoring DIF exists?

aiDIF provides:

- **Robust DIF analysis** under both human and AI scoring conditions,
  using a fully self-contained M-estimation engine (IRLS with bi-square loss)
- **Differential AI Scoring Bias (DASB) test** — a novel Wald test for
  item-level scoring shifts that differ across groups
- **AI-effect summaries** classifying items as: `stable_clean`, `stable_dif`,
  `introduced`, `masked`, or `new_direction` across scoring conditions
- **Anchor weight diagnostics** under potential AI contamination
- **Simulation utilities** for benchmarking DIF methods in AI-scored settings

## Installation

```r
# Install from GitHub
devtools::install_github("causalfragility-lab/aiDIF")

# Or install from local source
devtools::install_local("path/to/aiDIF")
```

## Quick Start

```r
library(aiDIF)

# Generate synthetic data with known DIF and DASB
dat <- simulate_aidif_data(n_items = 6, seed = 1)

# Fit the model
mod <- fit_aidif(
  human_mle = dat$human,
  ai_mle    = dat$ai
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

## Core Concepts

### Differential AI Scoring Bias (DASB)

For item *i* and group *g*, define the scoring shift:

```
delta_ig = d_ig^AI - d_ig^Human
```

where `d_ig` is the IRT intercept (difficulty) parameter. The DASB is:

```
DASB_i = delta_i2 - delta_i1
```

Under H₀: DASB_i = 0, a Wald test is conducted using the asymptotic
variance derived from the delta method (assuming independent groups and
scoring conditions):

```
Var(DASB_i) = Var(d_i1^H) + Var(d_i2^H) + Var(d_i1^AI) + Var(d_i2^AI)
```

A significant result means the AI scoring engine does not merely re-scale
all items uniformly — it disadvantages (or advantages) one group at
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

## From Existing IRT Fits

If you have fitted IRT models in **mirt**, use `read_ai_scored()` to bundle
your parameter estimates into the format `fit_aidif()` expects:

```r
library(mirt)

# Fit multigroup 2PL under human scoring
human_fit <- mirt(human_data, model = 1, itemtype = "2PL",
                  group = "group", SE = TRUE)

# Extract parameters manually and bundle
# (see ?read_ai_scored for the required list structure)
dat <- read_ai_scored(human_mle, ai_mle)

# Fit aiDIF model
mod <- fit_aidif(dat$human, dat$ai)
```

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

## Package Architecture

```
aiDIF/
├── R/
│   ├── read_functions.R    # read_ai_scored()
│   ├── aidif_core.R        # fit_aidif() — main estimation wrapper
│   ├── robust_engine.R     # estimate_robust_scale(), Wald tests, IRLS engine
│   ├── scoring_bias.R      # scoring_bias_test(), ai_effect_summary(),
│   │                       #   anchor_weights()
│   ├── simulate.R          # simulate_aidif_data()
│   ├── validate_inputs.R   # Internal validation helpers
│   └── class_functions.R   # print/summary/plot S3 methods
├── tests/
│   └── testthat/
│       └── test-aidif.R
└── DESCRIPTION
```

## Citation

If you use aiDIF in published research, please cite:

```
Hait, S. (2026). aiDIF: Differential Item Functioning for AI-Scored
Assessments. R package version 0.1.0.
https://github.com/causalfragility-lab/aiDIF
```

## License

GPL (>= 3)