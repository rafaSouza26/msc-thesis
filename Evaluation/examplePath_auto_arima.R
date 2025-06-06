# Install and load necessary packages
if (!requireNamespace("forecast", quietly = TRUE)) {
  install.packages("forecast")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

library(forecast)
library(ggplot2)
library(stringr)

# --- Configuration ---
example_data <- WWWusage
max_search_p <- 5
max_search_q <- 5
max_plot_order <- 7
plot_aspect_ratio <- 1
true_model_coords <- NULL # Example: c(p = 2, q = 6)

# --- 1. Run auto.arima and Capture Trace ---
set.seed(123)
trace_output <- capture.output({
  fit <- forecast::auto.arima(example_data,
                              stepwise = TRUE,
                              trace = TRUE,
                              approximation = FALSE,
                              seasonal = FALSE,
                              max.p = max_search_p,
                              max.q = max_search_q,
                              allowdrift = TRUE,
                              allowmean = TRUE)
})
print("--- Fitted Model ---")
print(fit)

# --- 2. Parse Trace Output to Extract Path ---
path_p <- integer()
path_q <- integer()
path_d <- integer()
model_fit_regex <- "^\\s*ARIMA\\((\\d+),(\\d+),(\\d+)\\)(?:\\(0,0,0\\)\\[0\\])?(?:\\s+with\\s+(?:non-)?zero\\s+mean|\\s+with\\s+drift)?\\s*:\\s*(-?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?|Inf|NA|-?NaN)"
for (line in trace_output) {
  match <- stringr::str_match(line, model_fit_regex)
  if (!is.na(match[1,1])) {
    p_val <- as.integer(match[1,2]); d_val <- as.integer(match[1,3]); q_val <- as.integer(match[1,4])
    path_p <- c(path_p, p_val); path_d <- c(path_d, d_val); path_q <- c(path_q, q_val)
  }
}
if (length(path_p) == 0) stop("No models found in trace.")
path_coords_df <- data.frame(p = path_p, q = path_q, d = path_d, step = 1:length(path_p))

# --- Create data frame for segments (for arrows) ---
if (nrow(path_coords_df) > 1) {
  segments_df <- data.frame(
    p_start = head(path_coords_df$p, -1),
    q_start = head(path_coords_df$q, -1),
    p_end   = tail(path_coords_df$p, -1),
    q_end   = tail(path_coords_df$q, -1)
  )
} else {
  segments_df <- data.frame(p_start=numeric(), q_start=numeric(), p_end=numeric(), q_end=numeric()) # Empty df
}

# --- 3. Prepare Data for Plotting Special Points ---
final_model_p <- fit$arma[1]; final_model_q <- fit$arma[2]; final_model_d <- fit$arma[6]
points_for_plot <- data.frame(p = path_coords_df$p, q = path_coords_df$q, type = "Path Point")
special_points_list <- list()
if (nrow(path_coords_df) > 0) {
  special_points_list[["Start Point"]] <- data.frame(p = path_coords_df$p[1], q = path_coords_df$q[1], type = "Start Point")
}
special_points_list[["End Point"]] <- data.frame(p = final_model_p, q = final_model_q, type = "End Point")
if (!is.null(true_model_coords) && length(true_model_coords) == 2 && !anyNA(true_model_coords)) {
  special_points_list[["True Model"]] <- data.frame(p = true_model_coords['p'], q = true_model_coords['q'], type = "True Model")
}
if (length(special_points_list) > 0) {
  priority_levels = c("True Model", "End Point", "Start Point", "Path Point")
  all_points_temp <- do.call(rbind, special_points_list); all_points_temp <- rbind(all_points_temp, points_for_plot)
  all_points_temp$type <- factor(all_points_temp$type, levels = priority_levels)
  all_points_temp <- all_points_temp[order(all_points_temp$p, all_points_temp$q, all_points_temp$type), ]
  final_plot_points <- all_points_temp[!duplicated(all_points_temp[, c("p", "q")]), ]
} else { final_plot_points <- points_for_plot[!duplicated(points_for_plot[, c("p", "q")]), ] }
legend_order_levels <- c("Start Point", "End Point", "Path Point")
if (!is.null(true_model_coords) && length(true_model_coords) == 2 && !anyNA(true_model_coords)) {
  if ("True Model" %in% final_plot_points$type) legend_order_levels <- c(legend_order_levels, "True Model")
}
final_plot_points$type <- factor(final_plot_points$type, levels = legend_order_levels)

# --- 4. Define Aesthetics ---
color_map <- c("Path Point" = "black", "Start Point" = "green", "End Point" = "blue", "True Model" = "red")
shape_map <- c("Path Point" = 16, "Start Point" = 16, "End Point" = 16, "True Model" = 16)
size_map  <- c("Path Point" = 2.0, "Start Point" = 2.5, "End Point" = 2.5, "True Model" = 3.0)
present_types <- levels(droplevels(final_plot_points$type))
color_map <- color_map[names(color_map) %in% present_types]; shape_map <- shape_map[names(shape_map) %in% present_types]; size_map <- size_map[names(size_map) %in% present_types]

# --- 5. Generate the Plot ---
plot_title <- paste0("Stepwise Search Path to Model (", final_model_p, ",", final_model_d[1], ",", final_model_q, ")")
search_path_plot <- ggplot() +
  # MODIFICATION: Use geom_segment for arrows on each segment
  geom_segment(data = segments_df,
               aes(x = p_start, y = q_start, xend = p_end, yend = q_end),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               color = "black", linewidth = 1.0) +
  # Points (drawn on top of arrow ends)
  geom_point(data = final_plot_points,
             aes(x = p, y = q, color = type, shape = type, size = type)) +
  scale_color_manual(values = color_map, name = NULL, drop = FALSE) +
  scale_shape_manual(values = shape_map, name = NULL, drop = FALSE) +
  scale_size_manual(values = size_map, name = NULL, drop = FALSE) +
  labs(title = plot_title, x = "p order", y = "q order") +
  
  # ----- MODIFIED LINES FOR THE NEW LIMITS -----
scale_x_continuous(breaks = 0:4, limits = c(-0.5, 4.5), expand = c(0,0)) +
  scale_y_continuous(breaks = 0:4, limits = c(-0.5, 4.5), expand = c(0,0)) +
  # -----------------------------------------------

theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "lightgrey"), panel.grid.minor = element_blank(),
    legend.position = "bottom", legend.box = "horizontal",
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.background = element_rect(fill="white", color=NA),
    legend.margin = margin(t = 5, unit = "pt"),
    aspect.ratio = if(!is.null(plot_aspect_ratio)) plot_aspect_ratio else NULL
  ) +
  guides(color = guide_legend(nrow = 1))

# Re-running these two lines is essential
print(search_path_plot)
ggsave(output_filename, plot = search_path_plot, width = 8, height = 8.5, dpi = 300)