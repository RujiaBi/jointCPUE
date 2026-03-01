---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# jointCPUE <img src="man/figures/logo.png" align="right" height="136" alt="jointCPUE logo" />

> Joint spatiotemporal CPUE standardization in TMB

<!-- badges: start -->
[![R-CMD-check](https://github.com/RujiaBi/jointCPUE/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/RujiaBi/jointCPUE/actions/workflows/R-CMD-check.yml)
<!-- badges: end -->

**jointCPUE** is an R package for spatiotemporal CPUE standardization using
Template Model Builder (TMB) and SPDE-based random fields. I wrote it for
multi-fleet settings where observed CPUE reflects both underlying population
density and fleet-specific catchability differences.

## Overview

The current implementation is built around a positive-CPUE model with:

- time effects
- population-level spatial random fields
- population-level spatiotemporal random fields
- optional fleet-specific catchability differences in:
  - system
  - time
  - space

For aggregated neon flying squid (NFS) data, the proportion of zero catches can
be too low to justify a stable delta-model parameterization. In the workflows
shown here, I therefore treat positive CPUE as the baseline response and retain
population and fleet-specific structure in space and time.

A typical workflow is:

1. prepare a CPUE data set with 0-based time and fleet indices
2. convert lon/lat to projected coordinates with `make_utm()`
3. build a mesh with `make_mesh()`
4. fit a model with `jointCPUE()`
5. extract a standardized index with `get_index()`
6. inspect predicted values and residual patterns

The package vignettes show the full model-comparison workflow, including yearly
and monthly analyses for the left and right sides of `170E`.

## Model principle

At a high level, the package treats observed CPUE as the product of population
biomass density and catchability:

```text
CPUE[f,s,t] = Biomass[f,s,t] * q[f,s,t]
```

where:

- `f` = fleet
- `s` = space
- `t` = time

On the log scale, the model is decomposed into a population-density component
and optional fleet-specific catchability differences:

```text
log(CPUE[f,s,t]) = y[t] + omega[s] + epsilon[s,t] + q_diffs[f,s,t]
```

Equivalently,

```text
log(Density[s,t]) = y[t] + omega[s] + epsilon[s,t]
log(q[f,s,t])     = q_diffs[f,s,t]
```

where:

- `y[t]` is a fixed time effect
- `omega[s]` is a time-constant spatial field
- `epsilon[s,t]` is a spatiotemporal field
- `q_diffs[f,s,t]` captures fleet-specific catchability differences

The first three terms determine the latent population density surface. The
`q_diffs` terms adjust the observation model for fleet-specific catchability,
but do not enter the standardized index itself.

### SPDE-based spatial fields

The spatial fields are represented with the SPDE approach on a triangular mesh.
In practice, this means the package approximates a continuous spatial Gaussian
field with a sparse precision matrix, which makes large spatial and
spatiotemporal models computationally feasible. `omega` is the time-constant
field, and `epsilon` is the time-varying field. The spatiotemporal field can be
fit either as independent time slices (`"iid"`) or as a first-order random walk
through time (`"rw"`).

Using an RW specification for `epsilon` often lowers the CV of the standardized
index because adjacent time steps share information instead of being estimated
independently. This stabilizes the latent density surface when some periods
have sparse observations or weak spatial coverage, and it usually produces a
smoother index than an iid specification.

For the monthly NFS workflows provided here, I do not assume an RW structure by
default. Once the data are filtered and split, observed months are not fully
continuous, and an RW process can force dependence across gaps or smooth over
real month-to-month changes. Monthly models are therefore fit with an iid
spatiotemporal field unless a continuous monthly time axis is retained and
temporal correlation is explicitly intended.

### Standardized CPUE index

After fitting the model, the latent population-density surface is projected onto
an extrapolation grid and area-weighted over space to obtain a standardized
index:

```text
Index[t] = sum_g Area[g] * Density[g,t]
```

This is the quantity returned by `get_index()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("RujiaBi/jointCPUE")
```

```r
library(jointCPUE)
```

## Required data structure

A minimal input data frame must contain:

- `cpue`: positive CPUE values
- `lon`, `lat`: geographic coordinates
- `tid`: 0-based contiguous time index
- `fleetid`: 0-based contiguous fleet index, with `0` as the reference fleet

Example:

```r
data_input <- data.frame(
  cpue = my_data$cpue,
  lon = my_data$lon,
  lat = my_data$lat,
  tid = my_data$tid,
  fleetid = my_data$fleetid
)
```

### Important notes

- `tid` must be `0:(n_t - 1)` with no gaps
- `fleetid` must be `0:(n_f - 1)` with no gaps
- the reference fleet must be coded as `0`
- the current baseline workflow assumes positive CPUE input
- this reflects the current NFS application, where aggregated CPUE data often
  have a relatively low proportion of zero catches

## Coordinate projection

`make_utm()` converts lon/lat to projected UTM coordinates and adds scaled
coordinates used internally for mesh construction and fitting.

It automatically:

- detects longitude convention (`-180..180` or `0..360`)
- chooses a UTM zone when `utm_zone = NULL`
- scales projected coordinates for numerical stability

```r
utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
data_utm <- utm$data_utm
```

The returned object includes:

- `utm_x`, `utm_y`
- `utm_x_scale`, `utm_y_scale`
- `utm_scale`, `utm_epsg`, `utm_zone`

## Build a spatial mesh

The mesh should be built with `utm_x_scale` and `utm_y_scale`.

### Cutoff mesh

```r
mesh <- make_mesh(
  data_utm,
  xy_cols = c("utm_x_scale", "utm_y_scale"),
  type = "cutoff",
  cutoff = 0.1
)
plot(mesh)
```

### K-means mesh

```r
mesh <- make_mesh(
  data_utm,
  xy_cols = c("utm_x_scale", "utm_y_scale"),
  type = "kmeans",
  n_knots = 50
)
plot(mesh)
```

### Tailored mesh

```r
mesh <- make_mesh(
  data_utm,
  xy_cols = c("utm_x_scale", "utm_y_scale"),
  type = "tailored",
  convex = -0.1,
  max.edge = c(0.5, 2),
  offset = c(0.1, 0.5),
  cutoff = 0.05
)
plot(mesh)
```

### Custom mesh

```r
mesh_custom <- fmesher::fm_rcdt_2d_inla(
  as.matrix(data_utm[, c("utm_x_scale", "utm_y_scale")]),
  cutoff = 0.1
)

mesh <- make_mesh(
  data_utm,
  xy_cols = c("utm_x_scale", "utm_y_scale"),
  type = "custom",
  mesh = mesh_custom
)
plot(mesh)
```

## Fit the model

A typical fit looks like this:

```r
fit <- jointCPUE(
  data_utm = data_utm,
  mesh = mesh,
  pop_spatial = "on",
  pop_spatiotemporal = "on",
  pop_spatiotemporal_type = "rw",
  q_diffs_system = "on",
  q_diffs_time = "off",
  q_diffs_spatial = "off",
  ncores = 4
)
```

### Population components

- `pop_spatial = "on"` includes `omega[s]`
- `pop_spatial = "off"` removes `omega[s]` from both the likelihood and the
  standardized index
- `pop_spatiotemporal = "on"` includes `epsilon[s,t]`
- `pop_spatiotemporal = "off"` removes `epsilon[s,t]` from both the likelihood
  and the standardized index
- `pop_spatiotemporal_type = "iid"` fits independent time slices
- `pop_spatiotemporal_type = "rw"` fits a first-order random walk through time

### Fleet-specific catchability components

- `q_diffs_system`: fleet-specific systematic differences
- `q_diffs_time`: fleet-specific temporal deviations
- `q_diffs_spatial`: fleet-specific spatial deviations

These terms are centered relative to the reference fleet and are used to adjust
catchability, not the standardized abundance index.

## Extract the standardized index

```r
index <- get_index(fit)
plot_index(index)
```

If you want calendar labels on the x-axis:

```r
year_lookup <- data.frame(
  Year = c(1997, 1998, 1999, 2001, 2002),
  t_year = 0:4
)

plot_index(
  index,
  time_values = year_lookup$Year,
  time_positions = year_lookup$Year
)
```

You can also compute an approximate marginal AIC:

```r
calc_marginal_aic(fit)
```

## Predicted values and residual diagnostics

Observation-level predicted values can be extracted with:

```r
pred <- get_predicted(
  fit,
  data = data_utm
)
```

This returns fitted values and residuals on both the response and log scales.

To visualize residual diagnostics:

```r
plots <- plot_residuals(
  pred,
  observed_col = "cpue",
  x_col = "utm_x_scale",
  y_col = "utm_y_scale",
  residual = "log",
  bins = 40
)

plots$observed_predicted
plots$spatial_residual
```

The spatial residual plot shows binned mean residuals, which is generally more
useful than plotting all residual points directly when observations overlap in
space.

## Fleet-specific q-differences

If q-difference components are turned on, they can be extracted and plotted:

```r
fleet_labels <- c("Reference", "Fleet 1", "Fleet 2")

plot_q_diffs_system(fit, fleet_labels = fleet_labels)
plot_q_diffs_time(fit, fleet_labels = fleet_labels)
plot_q_diffs_spatial(fit, fleet_labels = fleet_labels)
```

## Contact

For questions, suggestions, or collaboration:

- **Rujia Bi** — <bikayla5@gmail.com>

## Citation

If you use **jointCPUE**, please cite it as software:

> Bi, R. (2026). *jointCPUE: Joint spatiotemporal CPUE standardization using TMB*. R package (v0.1.0). <https://github.com/RujiaBi/jointCPUE>.
