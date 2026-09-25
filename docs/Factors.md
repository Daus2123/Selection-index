

## 1. Executive Summary & Scope

The **Agricultural & Experimental Design Analytical Pipeline** is an R Shiny application designed to automate statistical modeling, hypothesis testing, mean separation, and interaction plotting for multi-factor experiments. The engine handles **2-factor and 3-factor** layouts under two primary experimental designs: **Factorial Designs** and **Split-Plot (and Split-Split-Plot) Designs**.

The pipeline supports both **Completely Randomized Designs (CRD)** and **Randomized Complete Block Designs (RCBD)** as base randomization structures, with **Fisher’s Least Significant Difference (LSD)** and **Tukey’s Honestly Significant Difference (HSD)** for post-hoc pairwise mean comparisons.

---

## 2. Core Functional Requirements

### 2.1 Design & Statistical Engine Matrix

| Feature / Model | 2-Factor Factorial | 3-Factor Factorial | 2-Factor Split-Plot | 3-Factor Split-Split-Plot |
| --- | --- | --- | --- | --- |
| **Base Randomization** | CRD / RCBD | CRD / RCBD | CRD / RCBD (at Whole-Plot) | CRD / RCBD (at Whole-Plot) |
| **Factors Evaluated** | Factor A, Factor B | Factor A, Factor B, Factor C | Whole-Plot (A), Sub-Plot (B) | Whole-Plot (A), Sub-Plot (B), Sub-Sub-Plot (C) |
| **Model Type** | Fixed-Effects ANOVA | Fixed-Effects ANOVA | Linear Mixed Model (LMM) | Linear Mixed Model (LMM) |
| **Primary Computational Engine** | `stats::aov()` / `lm()` | `stats::aov()` / `lm()` | `lmerTest::lmer()` | `lmerTest::lmer()` |
| **Error Partitioning** | Residual error | Residual error | Whole-plot error ($\text{Block} \times A$) + Sub-plot residual | Whole-plot error ($\text{Block} \times A$) + Sub-plot error ($\text{Block} \times A \times B$) + Residual |

---

### 2.2 Mathematical Specifications

#### 1. Factorial Designs

* **2-Factor RCBD Model:**

$$y_{ijk} = \mu + b_k + \alpha_i + \beta_j + (\alpha\beta)_{ij} + \epsilon_{ijk}$$

Where $b_k \sim N(0, \sigma^2_b)$ is the block effect, $\alpha_i$ and $\beta_j$ are main effects, $(\alpha\beta)_{ij}$ is the interaction effect, and $\epsilon_{ijk} \sim N(0, \sigma^2_e)$ is the residual error. Set $b_k = 0$ for CRD.
* **3-Factor RCBD Model:**

$$y_{ijkl} = \mu + b_l + \alpha_i + \beta_j + \gamma_k + (\alpha\beta)_{ij} + (\alpha\gamma)_{ik} + (\beta\gamma)_{jk} + (\alpha\beta\gamma)_{ijk} + \epsilon_{ijkl}$$


#### 2. Split-Plot & Split-Split-Plot Designs

* **2-Factor Split-Plot Model (RCBD at Whole-Plot level):**

$$y_{ijk} = \mu + b_k + \alpha_i + \eta_{ik} + \beta_j + (\alpha\beta)_{ij} + \epsilon_{ijk}$$



Where $\eta_{ik} \sim N(0, \sigma^2_{wp})$ represents the Whole-Plot Error (Block $\times$ Factor A interaction).
* **3-Factor Split-Split-Plot Model:**

$$y_{ijkl} = \mu + b_l + \alpha_i + \eta_{il} + \beta_j + (\alpha\beta)_{ij} + \delta_{ijl} + \gamma_k + (\alpha\gamma)_{ik} + (\beta\gamma)_{jk} + (\alpha\beta\gamma)_{ijk} + \epsilon_{ijkl}$$



Where $\eta_{il}$ is Whole-Plot Error ($\text{Block} \times A$) and $\delta_{ijl}$ is Sub-Plot Error ($\text{Block} \times A \times B$).

---

### 2.3 Post-Hoc Comparison Engine

The post-hoc module computes Estimated Marginal Means (`emmeans`) across main effects and interaction factors, supporting two primary mean separation techniques:

1. **Fisher’s Least Significant Difference (LSD):**
* Computes unadjusted pairwise $p$-values and critical difference threshold:

$$LSD = t_{\alpha/2, df} \times SED$$


* Where $SED = \sqrt{2 \cdot MSE / n}$ is the Standard Error of Differences.


2. **Tukey’s Honestly Significant Difference (HSD):**
* Adjusts for family-wise error rates using the Studentized Range distribution:

$$HSD = q_{\alpha, k, df} \times \sqrt{\frac{MSE}{n}}$$


* **Compact Letter Display (CLD):** Generates group lettering (e.g., `a`, `ab`, `b`) assigned to treatments based on significance thresholds ($\alpha = 0.05$).

---

## 3. Data Input Schema & Input Requirements

| Input Field | Permitted Data Types | Design Applicability | Validation Rules |
| --- | --- | --- | --- |
| `Response Variable` | Continuous Numeric | All | Must contain non-null numeric values |
| `Design Type` | Categorical (`CRD`, `RCBD`) | All | Defines presence/absence of blocking factor |
| `Block / Replication` | Factor / Character / Integer | Mandatory for RCBD | Must have $\ge 2$ unique levels |
| `Factor A` | Factor / Character | All (Factorial A / Whole-Plot) | Minimum 2 levels |
| `Factor B` | Factor / Character | All (Factorial B / Sub-Plot) | Minimum 2 levels |
| `Factor C` | Factor / Character | Optional (3-Factor runs) | Minimum 2 levels when selected |

---

## 4. Visualizations & Outputs

1. **ANOVA / LMM Table:**
* Degree of Freedom ($df$), Sum of Squares ($SS$), Mean Squares ($MS$), $F$-statistic, $p$-value ($P(>F)$).
* For Split-Plot designs: Displays Satterthwaite approximation degrees of freedom and Type III ANOVA tables.


2. **Post-Hoc Tables:**
* Estimated Marginal Means, Standard Errors, Lower/Upper Confidence Limits, and Tukey/LSD Group Lettering.
 Layout Option A: Combined Mean ± SD / SE Header

| Group / Variety | BLUE ± SE | Compact Letter Display |
| --- | --- | --- |
| Variety A | 28.40 ± 0.42 | a |
| Variety B | 26.70 ± 0.38 | ab |
| Variety C | 24.20 ± 0.45 | b |
| **CV (%)** | **7.90%** |  |
| **Critical Difference ($\text{LSD}_{0.05}$, $\text{HSD}_{0.05}$, or $\text{Dunn}_{\text{crit}}$)** | **1.85** |  |
| **Heritability ($h^2$)** | **0.68** |  |

3. **Interaction Plots:**
* 2-Way profile plots (Factor A on x-axis, response on y-axis, grouped by Factor B).
* 3-Way faceted profile plots (Faceted by Factor C).
* Error bars displaying standard errors per factor combination.

# System Component Diagram

```
 +-----------------------------------------------------------------------+
 |                            Shiny UI (app.R)                           |
 |  [ Data Upload & Mapping ]  |  [ Factorial Tab ]  |  [ Split-Plot Tab ]|
 +-----------------------------------------------------------------------+
                                     |
                                     v
 +-----------------------------------------------------------------------+
 |                     Data Validation (utils_data_val.R)                |
 | - Coerce factors | Validate missing values | Check balance            |
 +-----------------------------------------------------------------------+
                                     |
                 +-------------------+-------------------+
                 |                                       |
                 v                                       v
 +-------------------------------+       +-------------------------------+
 | Factorial Module (utils_models) |       | Split-Plot Module(utils_models)|
 | - Fit standard ANOVA          |       | - Construct error strata      |
 | - Handle CRD / RCBD formulas  |       | - Fit LMM via lmerTest        |
 +-------------------------------+       +-------------------------------+
                 |                                       |
                 +-------------------+-------------------+
                                     |
                                     v
 +-----------------------------------------------------------------------+
 |                     Post-Hoc Engine (utils_posthoc.R)                 |
 | - Calculate EMMeans across factors & interactions                     |
 | - Execute LSD / Tukey pairwise tests                                  |
 | - Extract Compact Letter Display (CLD) groupings                      |
 +-----------------------------------------------------------------------+
                                     |
                                     v
 +-----------------------------------------------------------------------+
 |                   Output Generation & Rendering                       |
 | - Formatted ANOVA tables | EMMeans + Lettering | ggplot2 Plots        |
 +-----------------------------------------------------------------------+



#### Layout Option A: Combined Mean ± SD / SE Header

| Group / Variety | BLUE ± SE | Compact Letter Display |
| --- | --- | --- |
| Variety A | 28.40 ± 0.42 | a |
| Variety B | 26.70 ± 0.38 | ab |
| Variety C | 24.20 ± 0.45 | b |
| **CV (%)** | **7.90%** |  |
| **Critical Difference ($\text{LSD}_{0.05}$, $\text{HSD}_{0.05}$, or $\text{Dunn}_{\text{crit}}$)** | **1.85** |  |
| **Heritability ($h^2$)** | **0.68** |  |

