# Final yearly workflow for NFS data west of 170E.
# Run from the package root after installing jointCPUE, or source interactively.

data_file <- "NFS_jointCPUE_Data_split170E.RData"
ncores <- 8L
floor_value <- 1e-5
mesh_cutoff <- 0.5
fleet_labels <- c("CHN", "JPN", "KOR", "RUS", "TWN_hand", "TWN_machine")

library(jointCPUE)

if (!file.exists(data_file)) {
  stop("Cannot find data file: ", data_file, call. = FALSE)
}

load(data_file)

required_objects <- c("dd_left_170", "year_lookup_left_170")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects)) {
  stop(
    "The loaded data file is missing required objects: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

dd <- dd_left_170
year_lookup <- year_lookup_left_170

data_input <- data.frame(
  cpue = dd$CPUE + floor_value,
  lon = dd$X,
  lat = dd$Y,
  month = dd$Month,
  tid = dd$t_year,
  fleetid = dd$Fleet_code
)

utm <- make_utm(data_input, utm_zone = NULL, coord_scale = "auto")
data_utm <- utm$data_utm

mesh <- make_mesh(
  data_utm,
  xy_cols = c("utm_x_scale", "utm_y_scale"),
  type = "cutoff",
  cutoff = mesh_cutoff
)

fit <- jointCPUE(
  data_utm = data_utm,
  mesh = mesh,
  index = "yearly",
  pop_spatial = "on",
  pop_spatiotemporal = "on",
  pop_spatiotemporal_type = "rw",
  q_diffs_system = "on",
  q_diffs_time = "off",
  q_diffs_spatial = "on",
  ncores = ncores
)

convergence <- check_convergence(fit)
marginal_aic <- calc_marginal_aic(fit)

index_yearly <- get_index(fit)
index_plot <- plot_index(
  index_yearly,
  time_values = year_lookup$Year,
  time_positions = year_lookup$Year
)

predicted <- get_predicted(
  fit,
  data = data_utm,
  drop_floor = TRUE,
  floor_value = floor_value
)

residual_plots <- plot_residuals(
  predicted,
  observed_col = "cpue",
  x_col = "utm_x_scale",
  y_col = "utm_y_scale",
  residual = "log",
  bins = 40
)

q_system_plot <- plot_q_diffs_system(
  fit,
  fleet_labels = fleet_labels
)

q_spatial_plot <- plot_q_diffs_spatial(
  fit,
  fleet_labels = fleet_labels
)

result_left_170_yearly <- list(
  side = "left_170",
  fit = fit,
  convergence = convergence,
  marginal_aic = marginal_aic,
  index = index_yearly,
  index_plot = index_plot,
  predicted = predicted,
  residual_plots = residual_plots,
  q_system_plot = q_system_plot,
  q_spatial_plot = q_spatial_plot
)

print(convergence)
print(marginal_aic)
print(index_plot)

