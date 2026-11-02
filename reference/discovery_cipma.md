# Run Data-Driven cIPMA with Causal Discovery

This function performs causal discovery using constraint-based,
score-based, or hybrid algorithms to learn the structural model from
data, then runs cIPMA on the discovered structure. This exploratory
approach is useful when theoretical guidance is limited or when
researchers want to compare theory-driven and data-driven models.

## Usage

``` r
discovery_cipma(
  data,
  measurement_model,
  target_construct,
  scales,
  target_level = 85,
  nca_rep = NULL,
  pi = 1,
  B = 1999L,
  algorithm = "hc",
  blacklist = NULL,
  whitelist = NULL,
  algorithm.args = list(),
  nboot = 1000
)
```

## Arguments

- data:

  The raw data frame containing indicator variables.

- measurement_model:

  A seminr measurement model object (defined with constructs()).

- target_construct:

  String. The name of the dependent variable for cIPMA.

- scales:

  A named list of theoretical min/max for each indicator.

- target_level:

  Numeric. Desired target construct level (0-100). Default is 85.

- nca_rep:

  Deprecated. Number of permutations; use `B`.

- pi:

  Tolerance of the bottleneck, passed to
  [`cipma`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).
  Default 1.

- B:

  Number of permutations of the significance test, passed to
  [`cipma`](https://marekdejauj.github.io/cIPMA/reference/cipma.md).
  Default 1999.

- algorithm:

  Character. Causal discovery algorithm. Options:

  - Score-based: `"hc"` (default), `"tabu"`

  - Constraint-based: `"pc"`, `"gs"`, `"iamb"`, `"fast.iamb"`,
    `"inter.iamb"`, `"mmpc"`, `"hpc"`

  - Hybrid: `"mmhc"`, `"h2pc"`, `"rsmax2"`

- blacklist:

  Optional. A data frame with columns 'from' and 'to' specifying
  forbidden arcs (e.g., based on temporal ordering or theory).

- whitelist:

  Optional. A data frame with columns 'from' and 'to' specifying
  required arcs (e.g., well-established theoretical relationships).

- algorithm.args:

  Optional. A list of additional arguments to pass to the bnlearn
  algorithm (e.g., `list(restart = 10)` for hc, `list(alpha = 0.01)` for
  constraint-based algorithms).

- nboot:

  Integer. Number of bootstrap samples for PLS-SEM. Default 1000. Use
  lower values (e.g., 500) for faster exploratory runs.

## Value

A cipma object with additional elements:

- `learned_graph`: The bnlearn network object

- `discovered_arcs`: Data frame of discovered relationships

- `algorithm`: The algorithm used for discovery

## Details

This function implements the causal discovery approach for PLS-SEM as
discussed by Richter et al. (2020), who argue that combining causal
discovery algorithms with PLS-SEM can be valuable for:

- Exploratory research where theory is underdeveloped

- Model comparison (theory-driven vs. data-driven)

- Identifying potentially missing paths in existing models

- Cross-validation of theoretical assumptions

The function uses proxy scores (standardized mean of indicators) for
structure learning, which is consistent with the composite-based nature
of PLS-SEM. After structure learning, proper PLS-SEM estimation is
performed on the original indicators.

## Available Algorithms

- Score-based (fast):

  - `"hc"`: Hill-Climbing (default) - greedy search, fastest

  - `"tabu"`: Tabu Search - avoids local optima

- Constraint-based:

  - `"pc"`: PC algorithm (stable version)

  - `"gs"`: Grow-Shrink

  - `"iamb"`: Incremental Association Markov Blanket

  - `"fast.iamb"`: Fast IAMB

  - `"inter.iamb"`: Interleaved IAMB

  - `"mmpc"`: Max-Min Parents and Children

  - `"hpc"`: Hybrid Parents and Children

- Hybrid (recommended for robustness):

  - `"mmhc"`: Max-Min Hill-Climbing

  - `"h2pc"`: Hybrid HPC

  - `"rsmax2"`: General 2-phase Restricted Maximization

## Runtime Note

The main computational cost comes from PLS-SEM bootstrapping (Step 4),
not structure learning. Consider reducing `nboot` (e.g., 500) for
initial exploration, then increase for final results.

## Caution

Results from causal discovery should be interpreted as exploratory
hypotheses requiring theoretical justification and replication. The
discovered structure represents statistical associations that may
reflect causal relationships but cannot establish causality from
observational data alone.

## References

Richter, N. F., Schubring, S., Hauff, S., Ringle, C. M., & Sarstedt, M.
(2020). When predictors of outcomes are necessary: Guidelines for the
combined use of PLS-SEM and NCA. Industrial Management & Data Systems,
120(12), 2243-2267.

Scutari, M. (2010). Learning Bayesian Networks with the bnlearn R
Package. Journal of Statistical Software, 35(3), 1-22.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with whitelist (required arcs) and blacklist (forbidden arcs)

# Whitelist: Theory strongly supports these paths
whitelist <- data.frame(
  from = c("Quality", "Satisfaction"),
  to = c("Satisfaction", "Loyalty")
)

# Blacklist: These paths are theoretically impossible
blacklist <- data.frame(
  from = c("Loyalty", "Loyalty"),
  to = c("Quality", "Satisfaction")
)

# Run with hybrid algorithm (recommended)
result <- discovery_cipma(
  data = survey_data,
  measurement_model = mm,
  target_construct = "Loyalty",
  scales = scales,
  algorithm = "mmhc",
  whitelist = whitelist,
  blacklist = blacklist,
  nboot = 500  # Lower for faster exploration
)

# View discovered structure
print(result$discovered_arcs)
} # }
```
