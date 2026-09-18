# TTM

TTM R package for terminal trend model.

## Overview

The TTM package fits a retrospective-time longitudinal model for outcomes measured before death using weighted generalized estimating equations (WGEE). The main idea is to model the longitudinal outcome on a retrospective time scale,

$$
Y = \beta_\mu(t^*) + A \times \beta_A(t^*) + \beta_X^T X + \epsilon,
$$

where $t^*$ is the retrospective time, $A$ is a treatment indicator, $X$ are baseline covariates, and $\epsilon$ is a random error term. The retrospective time is defined as

$$
t^* = OS\_time - t_{pros},
$$

where $OS\_time$ is the observed survival time and $t_{pros}$ is the prospective measurement time. This makes it natural to study how a longitudinal outcome changes as a subject approaches death.

The package supports two model interfaces:

- `TTM()`: the default natural-spline version with automatic internal-knot selection over a candidate list.
- `TTM_linear()`: a linear spline version that uses a simpler basis and can be specified via `knots` and `spline_type`.

The target trajectories are:

- $\beta_\mu(t^*)$: the baseline retrospective trajectory,
- $\beta_A(t^*)$: the treatment effect trajectory over retrospective time.

The weighted GEE step uses inverse-probability-type weights built from a Cox dropout model and censoring information, and uncertainty is estimated using delete-group jackknife replicates.

## Data structure

The input data should be in long format, with one row per subject-visit. Each row should include:

- subject identifier (`id`),
- treatment indicator (`treatment`),
- longitudinal outcome (`outcome`),
- baseline covariates (`covariates`),
- observed survival or censoring time (`os_time`),
- event indicator (`event`),
- prospective follow-up time (`t_pros`).

The event variable is assumed to be coded as:

- 1 = death,
- 2 = dropout,
- other values = neither death nor dropout.

## Example dataset

The package includes a long-format example dataset:

```r
library(TTM)

data("TTM_data", package = "TTM")
```

A typical analysis uses:

- `id` as subject ID,
- `A` as treatment indicator,
- `Y` as the longitudinal outcome,
- `X1`, `X2` as baseline covariates,
- `OS_time` as observed survival/censoring time,
- `event` as the event indicator,
- `t_pros` as the prospective observation time.

## Main fitting function: `TTM()`

The default entry point is `TTM()`. It fits a retrospective-time model using a natural-spline basis for time and selects the number of internal knots from `n_inner_knot_list` using the supplied model-selection criterion.

The original `TTM()` implementation is the package's primary model-fitting function. It begins by cleaning the long-format data and constructing retrospective time as $t^* = OS\_time - t_{pros}$. It then fits a Cox model for dropout, computes the corresponding subject-level weights, and uses a weighted GEE among observations from subjects with observed death events to estimate the baseline trajectory $\beta_\mu(t^*)$ and the treatment effect trajectory $\beta_A(t^*)$. Uncertainty is obtained with delete-group jackknife replicates.

The most important arguments are:

- `data_long`: long-format data frame,
- `id`: subject identifier,
- `treatment`: treatment indicator,
- `outcome`: longitudinal outcome,
- `covariates`: baseline covariates,
- `os_time`: observed survival or censoring time,
- `event`: event indicator coded as 1 = death, 2 = dropout, and other values treated as neither death nor dropout,
- `t_pros`: prospective follow-up time,
- `n_inner_knot_list`: candidate values for the number of internal spline knots,
- `jk_block_size`: delete-group jackknife block size,
- `corstr`: working correlation structure passed to `geepack::geeglm()`,
- `criteria`: criterion used for choosing the final knot count.

```r
fit <- TTM(data_long = TTM_data,
								id = "id",
								treatment = "A",
								outcome = "Y",
								covariates = c("X1", "X2"),
								os_time = "OS_time",
								event = "event",
								t_pros = "t_pros",
								jk_block_size = NULL,
								n_inner_knot_list = range(10),
								corstr = "ar1",
								criteria = "RJC")
```

### Difference between `TTM()` and `TTM_linear()`

The two functions share the same overall objective and weighting framework, but they differ in the spline specification:

- `TTM()` uses the original natural-spline formulation with internal knots chosen from `n_inner_knot_list` and therefore returns the selected knot count in `fit$n_inner_knots`.
- `TTM_linear()` uses a simpler linear spline specification and allows direct control through `knots` and `spline_type`.
- `TTM()` is the default package workflow for flexible retrospective trajectories.
- `TTM_linear()` is a simpler alternative intended for a lower-dimensional or more directly specified basis.

In other words, `TTM()` is the more general original implementation, while `TTM_linear()` is a streamlined version built on the same retrospective-time model structure.

The function returns a list with components such as:

- `fit$dropout_formula`
- `fit$dropout_model`
- `fit$longitudinal_formula`
- `fit$weighted_gee`
- `fit$spline_boundary_knots`
- `fit$spline_inner_knots`
- `fit$n_inner_knots`
- `fit$variable_name_map`
- `fit$weights`
- `fit$data`
- `fit$gee_metrics`

The weighted GEE results are stored in `fit$weighted_gee`, including:

- `fixed_effects`
- `fixed_effects_vcov`
- `fixed_effects_sd`

## Linear-spline variant: `TTM_linear()`

`TTM_linear()` is a lower-dimensional alternative that constructs a spline basis using a linear spline representation. This can be useful when a simpler basis is desired or when you want to specify `knots` directly rather than using the automatic knot-selection procedure in `TTM()`.

```r
fit_linear <- TTM_linear(
	data_long = TTM_data,
	id = "id",
	treatment = "A",
	outcome = "Y",
	covariates = c("X1", "X2"),
	os_time = "OS_time",
	event = "event",
	t_pros = "t_pros",
	knots = 6,
	jk_block_size = NULL,
	corstr = "ar1",
	spline_type = "bs"
)
```

`TTM_linear()` has the same overall structure as `TTM()` but stores additional information such as:

- `fit_linear$type` (set to `"linear"`),
- `fit_linear$knots`,
- `fit_linear$spline_type`.

This is a convenient option when the main interest is a simpler retrospective-time specification while still using the same weighting and jackknife framework.

## Extracting trajectory estimates

Once a model has been fit, the estimated retrospective trajectories can be evaluated at specific time points with `get_beta_mu()` and `get_beta_A()`.

To evaluate the baseline trajectory $\beta_\mu(t^*)$ at selected retrospective times:

```r
get_beta_mu(new_times = 1:10, results = fit, conf.level = 0.95)
```

To evaluate the treatment-effect trajectory $\beta_A(t^*)$ at selected retrospective times:

```r
get_beta_A(new_times = 1:10, results = fit, conf.level = 0.95)
```

The same functions also work for `TTM_linear()` results, because the helper functions reconstruct the appropriate spline basis from the fitted object.

These functions return a list containing the evaluated time grid, estimated trajectory, pointwise standard error, and lower/upper confidence limits.

## Visualization

The package includes plotting functions for the fitted retrospective trajectories.

To plot the estimated baseline trajectory $\beta_\mu(t^*)$ with pointwise confidence intervals:

```r
plot_beta_mu(results = fit, conf.level = 0.95)
```

To plot the estimated treatment-effect trajectory $\beta_A(t^*)$ with pointwise confidence intervals:

```r
plot_beta_A(results = fit, conf.level = 0.95)
```

The same plotting functions can be used with `fit_linear` as well:

```r
plot_beta_mu(results = fit_linear, conf.level = 0.95)
plot_beta_A(results = fit_linear, conf.level = 0.95)
```

## Test slope difference for `TTM_linear()` object

Test the difference in the slope of $\beta_\mu(t^*)$ between and after the specified knots. `new_times` should be a vector of three retrospective times, where the first one is before the knot and the last is after the knot. The middle one is the knot in the `TTM_linear()` function.

```r
test_beta_slope(new_times = c(5, 6, 7), fit = fit_linear, var_name = "X0")
```

Test the difference in the slope of $\beta_A(t^*)$ between and after the specified knots.

```r
test_beta_slope(new_times = c(5, 6, 7), fit = fit_linear, var_name = "A")
```

## Summary

The package provides two complementary retrospective-time model fits:

- `TTM()` for the default natural-spline model with automatic knot selection,
- `TTM_linear()` for a simpler linear-spline formulation with direct control of knots and spline type.

Both approaches estimate treatment and baseline retrospective trajectories using weighted GEE and return trajectory summaries and plots through the helper functions in the package.
