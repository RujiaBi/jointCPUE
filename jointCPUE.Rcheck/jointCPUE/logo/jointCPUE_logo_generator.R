# ============================================================
# jointCPUE hex logo generator
# Output:
#   - man/figures/logo.png
# ============================================================

library(sf)
library(ggplot2)
library(dplyr)
library(scales)

sf_use_s2(FALSE)
set.seed(2026)

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)])
script_path <- gsub("~\\+~", " ", script_path)
script_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
pkg_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
logo_dir <- file.path(pkg_root, "man", "figures")
dir.create(logo_dir, recursive = TRUE, showWarnings = FALSE)

make_hex <- function(r = 1, rot = pi / 2) {
  ang <- rot + (0:5) * pi / 3
  xy <- cbind(x = r * cos(ang), y = r * sin(ang))
  rbind(xy, xy[1, ])
}

crs_use <- 3857
hex_outer <- st_sf(geometry = st_sfc(st_polygon(list(make_hex(1.00))), crs = crs_use))
hex_inner <- st_sf(geometry = st_sfc(st_polygon(list(make_hex(0.962))), crs = crs_use))

# ------------------------------------------------------------
# Sea surface field: stronger contrast, diagonal light lane,
# softened edge darkening, and non-flat texture.
# ------------------------------------------------------------
gx <- seq(-0.98, 0.98, length.out = 320)
gy <- seq(-0.98, 0.98, length.out = 340)
sea_grid <- expand.grid(x = gx, y = gy)

sea_grid$z <- with(sea_grid, {
  0.52 +
    0.20 * (0.55 * x + 0.85 * y) +
    0.09 * sin(2.6 * x + 0.5) +
    0.07 * cos(2.8 * y - 0.4) +
    0.08 * sin(3.4 * (x + 0.6 * y) + 1.2) +
    0.06 * cos(4.5 * (x - 0.8 * y) - 0.7)
})

blobs <- data.frame(
  cx = runif(8, -0.75, 0.75),
  cy = runif(8, -0.75, 0.75),
  amp = runif(8, -0.16, 0.16),
  sx = runif(8, 0.18, 0.42),
  sy = runif(8, 0.18, 0.38)
)

for (k in seq_len(nrow(blobs))) {
  sea_grid$z <- sea_grid$z + with(
    blobs[k, ],
    amp * exp(-((sea_grid$x - cx)^2 / sx^2 + (sea_grid$y - cy)^2 / sy^2))
  )
}

# Diagonal specular lane, like sunlight dragging across the sea surface.
sea_grid$glow <- with(sea_grid, {
  exp(-((y - (0.10 + 0.35 * x))^2) / (2 * 0.13^2))
})

# Vignette to frame the hex without making the edges too heavy.
sea_grid$edge <- with(sea_grid, sqrt(x^2 + y^2))
sea_grid$z <- sea_grid$z + 0.18 * sea_grid$glow - 0.17 * sea_grid$edge^1.55
sea_grid$z <- scales::rescale(sea_grid$z, to = c(0.16, 1.00))
sea_grid$glow <- scales::rescale(sea_grid$glow, to = c(0, 0.24))

sea_pts <- st_as_sf(sea_grid, coords = c("x", "y"), crs = crs_use, remove = FALSE)
inside <- lengths(st_within(sea_pts, hex_inner)) > 0
sea_grid <- sea_grid[inside, , drop = FALSE]

tile_w <- diff(gx[1:2])
tile_h <- diff(gy[1:2])

sea_pal <- c(
  "#123054",
  "#1B4271",
  "#275694",
  "#3A74B6",
  "#5A9FD6",
  "#88C2E1",
  "#B7E4EC"
)

# ------------------------------------------------------------
# Wave lines
# ------------------------------------------------------------
x <- seq(-1.25, 1.25, length.out = 900)
line_levels <- c(0.76, 0.64, 0.54, 0.45, 0.36, 0.28, 0.20, 0.11, 0.03, -0.08, -0.20, -0.35, -0.52, -0.72)

flow_lines <- lapply(seq_along(line_levels), function(i) {
  k <- 2.1 + i * 0.08
  amp <- 0.010 + i * 0.0013
  ph <- i * 0.65
  yy <- line_levels[i] +
    amp * sin(k * x + ph) +
    0.004 * sin((k * 1.7) * x + ph * 0.4)
  st_linestring(cbind(x, yy))
})

flow_sf <- st_sf(
  wave_alpha = seq(0.06, 0.18, length.out = length(flow_lines)),
  wave_size = seq(0.35, 0.85, length.out = length(flow_lines)),
  geometry = st_sfc(flow_lines, crs = crs_use)
)
flow_clip <- suppressWarnings(st_intersection(flow_sf, st_make_valid(hex_inner)))

# ------------------------------------------------------------
# Sparkles concentrated around the diagonal glow
# ------------------------------------------------------------
spark_core <- st_sample(hex_inner, size = 70, type = "random") |> st_sf()
sxy <- st_coordinates(spark_core)
spark_df <- data.frame(
  x = sxy[, 1],
  y = sxy[, 2]
)
spark_df$lane <- exp(-((spark_df$y - (0.10 + 0.35 * spark_df$x))^2) / (2 * 0.18^2))
spark_df$size <- rescale(spark_df$lane, to = c(0.12, 1.10))
spark_df$alpha <- rescale(spark_df$lane, to = c(0.02, 0.18))

# ------------------------------------------------------------
# A few brighter crest streaks near the light lane
# ------------------------------------------------------------
crest_x <- seq(-0.65, 0.75, length.out = 420)
crest_lines <- lapply(c(0.34, 0.10, -0.16), function(level) {
  yy <- level + 0.015 * sin(3.0 * crest_x + level * 8) + 0.010 * sin(7.0 * crest_x - 0.6)
  st_linestring(cbind(crest_x, yy))
})
crest_sf <- st_sf(geometry = st_sfc(crest_lines, crs = crs_use))
crest_clip <- suppressWarnings(st_intersection(crest_sf, st_make_valid(hex_inner)))

# ------------------------------------------------------------
# Build plot
# ------------------------------------------------------------
p <- ggplot() +
  geom_sf(data = hex_outer, fill = "#26384B", color = NA) +
  geom_sf(data = hex_inner, fill = "#295A94", color = NA) +
  geom_tile(
    data = sea_grid,
    aes(x, y, fill = z),
    width = tile_w,
    height = tile_h,
    inherit.aes = FALSE
  ) +
  geom_tile(
    data = sea_grid,
    aes(x, y, alpha = glow),
    width = tile_w,
    height = tile_h,
    fill = "#D3F1F3",
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(colors = sea_pal, guide = "none") +
  geom_sf(
    data = flow_clip,
    aes(alpha = wave_alpha, linewidth = wave_size),
    color = "#E7FBFF",
    lineend = "round",
    show.legend = FALSE
  ) +
  scale_linewidth_identity() +
  geom_sf(
    data = crest_clip,
    color = alpha("#F7FEFF", 0.18),
    linewidth = 1.8,
    lineend = "round"
  ) +
  geom_point(
    data = spark_df,
    aes(x, y, size = size, alpha = alpha),
    color = "white"
  ) +
  scale_size_identity() +
  scale_alpha_identity()

# Text: stronger shadow, slight ocean-glow under the title.
p <- p +
  annotate(
    "text",
    x = 0.010,
    y = -0.012,
    label = "jointCPUE",
    color = alpha("#04284A", 0.70),
    size = 31.8,
    family = "sans",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0.004,
    y = -0.004,
    label = "jointCPUE",
    color = alpha("#8FEAF1", 0.22),
    size = 31.0,
    family = "sans",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 0,
    y = 0,
    label = "jointCPUE",
    color = "white",
    size = 30.2,
    family = "sans",
    fontface = "bold"
  )

p <- p +
  geom_sf(data = hex_inner, fill = NA, color = "#4D6887", linewidth = 2.0) +
  coord_sf(
    xlim = c(-1.04, 1.04),
    ylim = c(-1.04, 1.04),
    expand = FALSE,
    clip = "on"
  ) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(logo_dir, "logo.png"), p, width = 8, height = 8, dpi = 1200, bg = "transparent")

message("Saved: ", file.path(logo_dir, "logo.png"))
