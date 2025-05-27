library(tidyverse)
library(ggplot2)

#' Read and validate data
read_and_validate_data <- function(file_path) {
  data <- tryCatch(
    {
      read.csv(file_path, header = TRUE, check.names = FALSE)
    },
    error = function(e) {
      stop(paste("Error reading file:", e$message))
    }
  )
  return(data)
}

#' Transform data to long format
#' @param data Raw data
#' @param group_pattern Grouping pattern for the replicate names
#' @param group_names Group names for the two populations
reshape_to_long_format <- function(data, group_pattern, group_names) {
  data_long <- data %>%
    rename(Replicate = 1) %>%
    pivot_longer(
      cols = -Replicate,
      names_to = "Time",
      values_to = "Value"
    ) %>%
    mutate(
      Time = as.numeric(Time),
      Group = case_when(
        grepl(group_pattern[[1]], Replicate, ignore.case = TRUE) ~ group_names[1],
        grepl(group_pattern[[2]], Replicate, ignore.case = TRUE) ~ group_names[2],
        TRUE ~ "Other"
      )
    ) %>%
    filter(Group != "Other")
  
  return(data_long)
}

#' Calculate summary statistics
#' @param data_long Long format data
calculate_summary_stats <- function(data_long) {
  summary_data <- data_long %>%
    group_by(Group, Time) %>%
    summarise(
      Mean = median(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      SE = SD / sqrt(n()),
      CI_lower = quantile(Value, probs = 0.025, na.rm = TRUE),
      CI_upper = quantile(Value, probs = 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary_data)
}

#' Process a single statistic
process_single_stat <- function(file_path, stat_name, group_pattern, group_names) {
  data <- read_and_validate_data(file_path)
  data_long <- reshape_to_long_format(data, group_pattern, group_names)
  summary_data <- calculate_summary_stats(data_long)
  summary_data$Statistic <- stat_name
  return(summary_data)
}


#' Combine plots of each statistics
plot_combined_stats <- function(combined_data, scenario_name, 
  colors = c("Constant" = "#2bbac1", "Decline" = "#f2635b"),
  crash_time = 36) {
  
  # Plot the data
  p <- ggplot(combined_data, aes(x = Time, y = Mean, color = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 1.5, alpha = 0.8) +
  facet_wrap(~Statistic, scales = "free_y", ncol = 3) +

  labs(
  title = NULL,
  x = "Generations ago (past -> present)",
  y = NULL,
  color = "Population",
  linetype = NULL,
  caption = str_wrap(
    "Figure 1. Summary statistics of different population simulations across time. Each plot compares a constant and a declined population. Dashed lines represent the time of population crash. Error bars represent 95% confidence intervals. Note: The x-axis is reversed to show the most recent generations on the right.", 
    width = 180
  )
  ) +

  scale_color_manual(values = colors) +
  scale_x_reverse() +
  theme_minimal(base_size = 13) +
  theme(
  legend.position = c(0.75, 0.2),
  plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  strip.text = element_text(face = "bold", size = 13),
  strip.placement = "outside",
  panel.spacing = unit(1, "lines"),

  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16),
  plot.caption = element_text(
    size = 14,
    hjust = 0,  # Align to the left
    vjust = 1,  # Align to the bottom
    margin = margin(t = 10, b = 5, r = 5),  # Add margin to the top, bottom and right
    face = "italic"
  ),
  plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),  # Margin: top, right, bottom, left
  plot.background = element_rect(colour = "black", fill = NA, size = 1)  # Add border
  )
  
  # Add vertical line for population crash
  p <- p + 
  geom_vline(aes(xintercept = crash_time, linetype = "Population crash"), 
   color = "darkred", alpha = 0.7) +
  scale_linetype_manual(values = c("Population crash" = "dashed"))

  return(p)
}

#' Main function to process a scenario
process_scenario <- function(scenario_dir, scenario_name, 
                           group_pattern = list(constant = "constant", decline = "decline"),
                           group_names = c("Constant", "Decline")) {
  
  setwd(scenario_dir)
  
  # Define the statistics to process
  stats_list <- list(
    list(file = "Tajimas_D.csv", name = "Tajima's D"),
    list(file = "Density_of_segregating_sites.csv", name = "Density of segregating sites"),
    list(file = "Nucleotide_diversity_pi.csv", name = "Nucleotide diversity (π)"),
    list(file = "Number_of_mutations.csv", name = "Number of mutations"),
    list(file = "Number_of_trees.csv", name = "Number of trees")
  )
  
  # Process each statistic and combine the results
  combined_data <- map_dfr(stats_list, function(stat) {
    process_single_stat(stat$file, stat$name, group_pattern, group_names)
  })
  
  # Plotting
  p <- plot_combined_stats(combined_data, scenario_name)
  
  # Save the plot
  output_dir <- "../combined_plots"
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  output_file <- file.path(output_dir, paste0(scenario_name, "_combined_stats.pdf"))
  ggsave(output_file, plot = p, width = 18, height = 10, device = cairo_pdf)
  
  output_file_png <- file.path(output_dir, paste0(scenario_name, "_combined_stats.png"))
  ggsave(output_file_png, plot = p, width = 18, height = 10, dpi = 300, bg = "white")

  return(p)
}

# ' Define the scenarios to process
scenarios <- list(
  list(dir = "~/Documents/MresProject/outputs/statistics", name = NULL)
)

# ' Process all scenarios
walk(scenarios, function(scenario) {
  process_scenario(scenario$dir, scenario$name)
})
