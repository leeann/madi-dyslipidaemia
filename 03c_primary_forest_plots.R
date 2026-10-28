# Create publication forest plots from Model 1 and Model 2 estimates.

if (!exists("load_analysis_packages", mode = "function")) {
  source(file.path("R", "00b_functions.R"))
}

load_analysis_packages(c("broom", "dplyr", "ggplot2"))
create_output_directories()

primary_models <- readRDS(
  file.path("results", "models", "primary-models.rds")
)

forest_rows <- lapply(names(outcome_labels), function(outcome) {
  dplyr::bind_rows(
    tidy_madi_terms(
      primary_models[[outcome]]$models$model1,
      outcome,
      "Model 1"
    ),
    tidy_madi_terms(
      primary_models[[outcome]]$models$model2,
      outcome,
      "Model 2"
    )
  )
})

forest_data <- dplyr::bind_rows(forest_rows) |>
  dplyr::mutate(
    outcome = factor(outcome, levels = unname(outcome_labels)),
    contrast = dplyr::case_when(
      grepl("Q2", term) ~ "Q2 vs Q1",
      grepl("Q3", term) ~ "Q3 vs Q1",
      grepl("Q4", term) ~ "Q4 vs Q1",
      grepl("Q5", term) ~ "Q5 vs Q1"
    ),
    contrast = factor(
      contrast,
      levels = c("Q5 vs Q1", "Q4 vs Q1", "Q3 vs Q1", "Q2 vs Q1")
    ),
    specification = factor(specification, levels = c("Model 1", "Model 2")),
    estimate_ci = sprintf(
      "%s: %.2f (%.2f, %.2f)",
      specification, estimate, conf.low, conf.high
    )
  )

make_forest_plot <- function(data, logarithmic = TRUE, q4_q5_only = FALSE) {
  if (q4_q5_only) {
    data <- data |>
      dplyr::filter(contrast %in% c("Q4 vs Q1", "Q5 vs Q1")) |>
      droplevels()
    header_y <- 2.48
    figure_title <- "Q4 and Q5 versus Q1 (logarithmic x-axis)"
  } else {
    header_y <- 4.48
    figure_title <- if (logarithmic) {
      "Logarithmic x-axis"
    } else {
      "Linear x-axis"
    }
  }

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = estimate, y = contrast, shape = specification)
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = 2,
      colour = "grey55"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0,
      colour = "black",
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::geom_point(
      size = 2.2,
      colour = "black",
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = 2.02,
        y = contrast,
        label = estimate_ci,
        group = specification
      ),
      position = ggplot2::position_dodge(width = 0.55),
      hjust = 0,
      size = 2.0,
      colour = "black",
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 2.02,
      y = header_y,
      label = "Model: PR (95% CI)",
      hjust = 0,
      fontface = "bold",
      size = 2.2
    ) +
    ggplot2::facet_wrap(
      ~outcome,
      ncol = 1,
      labeller = ggplot2::labeller(
        outcome = ggplot2::label_wrap_gen(width = 30)
      )
    ) +
    ggplot2::coord_cartesian(xlim = c(0.75, 2.0), clip = "off") +
    ggplot2::scale_shape_manual(
      values = c("Model 1" = 17, "Model 2" = 15)
    ) +
    ggplot2::labs(
      title = figure_title,
      x = "Prevalence ratio",
      y = NULL,
      shape = "Model",
      caption = paste0(
        "Q1 is the reference.\n",
        "M1 adjusts for age, sex, and ethnicity.\n",
        "M2 additionally adjusts for marital status, education,\n",
        "occupation, and income."
      )
    ) +
    ggplot2::theme_minimal(base_size = 7.5) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 7.5,
        hjust = 0,
        margin = ggplot2::margin(l = -31)
      ),
      strip.clip = "off",
      plot.title = ggplot2::element_text(size = 9, face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 6),
      legend.title = ggplot2::element_text(size = 6.5),
      legend.text = ggplot2::element_text(size = 6.5),
      axis.text = ggplot2::element_text(size = 6.5),
      axis.title = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(4, 95, 4, 4)
    )

  if (logarithmic) {
    plot <- plot + ggplot2::scale_x_log10(
      breaks = c(0.75, 1.0, 1.25, 1.5, 2.0),
      labels = c("0.75", "1.00", "1.25", "1.50", "2.00")
    )
  } else {
    plot <- plot + ggplot2::scale_x_continuous(
      breaks = c(0.75, 1.0, 1.25, 1.5, 2.0),
      labels = c("0.75", "1.00", "1.25", "1.50", "2.00")
    )
  }

  plot
}

forest_log <- make_forest_plot(forest_data, logarithmic = TRUE)
forest_linear <- make_forest_plot(forest_data, logarithmic = FALSE)
forest_q4_q5 <- make_forest_plot(
  forest_data,
  logarithmic = TRUE,
  q4_q5_only = TRUE
)

save_plot_outputs(forest_log, "forest-primary-log", width = 3.5, height = 10.5)
save_plot_outputs(
  forest_linear,
  "forest-primary-linear",
  width = 3.5,
  height = 10.5
)
save_plot_outputs(
  forest_q4_q5,
  "forest-primary-q4-q5-log",
  width = 3.5,
  height = 8
)

utils::write.csv(
  forest_data,
  file.path("results", "data", "primary-forest-estimates.csv"),
  row.names = FALSE
)

message("Saved primary forest plots and plotting data.")
