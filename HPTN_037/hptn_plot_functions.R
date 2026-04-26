

# Auxiliary function that generate the figures of the results
library(ggplot2)
library(latex2exp)
library(scales)

sa_ie_plot <- function(sa_ie_data, 
                       sa_ie_naive_data,
                       palette,
                       plot_labels,
                       m_vec_alters){
  
  min_y <- round(min(sa_ie_data$ci_low),1)
  max_y <- round(max(sa_ie_data$ci_high),1)
  
  subtitle_text <- paste0("Naive Estimate [95% CI]: ",
                          round(sa_ie_naive_data$ie_rd[1],3),
                          " [",
                          round(sa_ie_naive_data$ci_low[1],3),
                          ", ",
                          round(sa_ie_naive_data$ci_high[1],3),
                          "]")
  ie_plot <- ggplot() + 
    geom_ribbon(data = sa_ie_data,
                aes(
                  x = pi_param,
                  ymin = ci_low,
                  ymax = ci_high,
                  fill = model_type
                ),
                alpha = 0.3 
    ) +
    geom_line(
      data = sa_ie_data,
      aes(
        x = pi_param,
        y = ie_rd,
        color = model_type
      ),
      linewidth = 1.8
    ) +
    # Add the Naive point estimates
    geom_point(
      data = sa_ie_naive_data,
      aes(
        x = pi_param,
        y = ie_rd
        # shape = aug,
        # group = aug
      ),
      color = "black",
      size = 4,
      # position = position_dodge(12)
    ) +
    # Add the Naive error bars
    geom_errorbar(
      data = sa_ie_naive_data,
      aes(
        x = pi_param,
        ymin = ci_low,
        ymax = ci_high
        # group = aug
      ),
      color = "black",
      width = 10,
      linewidth = 0.7,
      # position = position_dodge(12)
    ) +
    # Horizontal line at y=0 (no effect)
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey48"
    ) +
    # Use the custom color palette for lines and fills
    scale_color_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    scale_fill_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    # Customize the shape legend
    # scale_shape_manual(
    #   name = "Naive Model:",
    #   values = c("TRUE" = 17, "FALSE" = 16), # Triangle and Circle
    #   labels = c("TRUE" = "Augmented", "FALSE" = "Unaugmented")
    # ) +
    # Add labels (using latex2exp for the x-axis)
    labs(
      title = "Indirect Effect",
      subtitle = subtitle_text,
      y = "Estimated IE",
      x = TeX("$m^a$ (Expected number of missing alter-ego edges)")
    ) +
    
    # Adjust x-axis scale
    scale_x_continuous(breaks = seq(0, max(m_vec_alters), 50),
                       labels = seq(0, max(m_vec_alters), 50)) +
    # scale_y_continuous(breaks = c(-0.45, -0.35, -0.25, -0.15, -0.05,0,0.05,0.15)) +
    # scale_y_continuous(breaks = seq(-0.5,0.2,0.1)) +
    scale_y_continuous(breaks = seq(min_y,max_y,0.2), 
                       labels = number_format(accuracy = 0.1)) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 13),
      axis.text = element_text(size = 13)
    )
  
  return(ie_plot)
}

sa_de_plot_given_kappa <- function(sa_de_data, 
                                   sa_de_naive_data,
                                   kappa,
                                   palette,
                                   plot_labels,
                                   m_vec_egos){
  subtitle_text <- paste0("Naive Estimate [95% CI]: ",
                          round(sa_de_naive_data$de_rd[1],3),
                          " [",
                          round(sa_de_naive_data$ci_low[1],3),
                          ", ",
                          round(sa_de_naive_data$ci_high[1],3),
                          "]")
  de_plot <- ggplot() + 
    geom_ribbon(data = sa_de_data,
                aes(
                  x = pi_param,
                  ymin = ci_low,
                  ymax = ci_high,
                  fill = spec
                ),
                alpha = 0.3 
    ) +
    geom_line(
      data = sa_de_data,
      aes(
        x = pi_param,
        y = de_rd,
        color = spec
      ),
      linewidth = 1.8
    ) +
    # Add the Naive point estimates
    geom_point(
      data = sa_de_naive_data,
      aes(
        x = pi_param,
        y = de_rd
        # shape = aug,
        # group = aug
      ),
      color = "black",
      size = 4,
      # position = position_dodge(12)
    ) +
    # Add the Naive error bars
    geom_errorbar(
      data = sa_de_naive_data,
      aes(
        x = pi_param,
        ymin = ci_low,
        ymax = ci_high
        # group = aug
      ),
      color = "black",
      width = 10,
      linewidth = 0.7,
      # position = position_dodge(12)
    ) +
    # Horizontal line at y=0 (no effect)
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey48"
    ) +
    # Use the custom color palette for lines and fills
    scale_color_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    scale_fill_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    # Customize the shape legend
    # scale_shape_manual(
    #   name = "Naive Model:",
    #   values = c("TRUE" = 17, "FALSE" = 16), # Triangle and Circle
    #   labels = c("TRUE" = "Augmented", "FALSE" = "Unaugmented")
    # ) +
    # Add labels (using latex2exp for the x-axis)
    labs(
      title = paste0("Direct Effect; kappa = ", kappa),
      subtitle = subtitle_text,
      y = "Estimated DE",
      x = TeX("$m^e$ (Expected number of missing ego-ego edges)")
    ) +
    
    # Adjust x-axis scale
    scale_x_continuous(breaks = seq(0, max(m_vec_egos), 50),
                       labels = seq(0, max(m_vec_egos), 50)) +
    # scale_y_continuous(breaks = c(-0.45, -0.35, -0.25, -0.15, -0.05,0,0.05,0.15)) +
    # scale_y_continuous(breaks = seq(min(sa_de_data$ci_low)-0.01,
                                    # max(sa_de_data$ci_high)+0.01,
                                    # 0.1)) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 13),
      axis.text = element_text(size = 13)
    )
  
  return(de_plot)
}

sa_de_plot_given_m_e <- function(sa_de_data, 
                                 sa_de_naive_data,
                                 m_e,
                                 palette,
                                 plot_labels){
  subtitle_text <- paste0("Naive Estimate [95% CI]: ",
                          round(sa_de_naive_data$de_rd[1],3),
                          " [",
                          round(sa_de_naive_data$ci_low[1],3),
                          ", ",
                          round(sa_de_naive_data$ci_high[1],3),
                          "]")
  de_plot <- ggplot() + 
    geom_ribbon(data = sa_de_data,
                aes(
                  x = kappa,
                  ymin = ci_low,
                  ymax = ci_high,
                  fill = spec
                ),
                alpha = 0.3 
    ) +
    geom_line(
      data = sa_de_data,
      aes(
        x = kappa,
        y = de_rd,
        color = spec
      ),
      linewidth = 1.8
    ) +
    # Add the Naive point estimates
    geom_point(
      data = sa_de_naive_data,
      aes(
        x = kappa,
        y = de_rd
        # shape = aug,
        # group = aug
      ),
      color = "black",
      size = 4,
      # position = position_dodge(12)
    ) +
    # Add the Naive error bars
    geom_errorbar(
      data = sa_de_naive_data,
      aes(
        x = kappa,
        ymin = ci_low,
        ymax = ci_high
        # group = aug
      ),
      color = "black",
      width = 0.05,
      linewidth = 0.7,
      # position = position_dodge(12)
    ) +
    # Horizontal line at y=0 (no effect)
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey48"
    ) +
    # Use the custom color palette for lines and fills
    scale_color_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    scale_fill_manual(
      name = "SA Model:",
      values = palette,
      labels = plot_labels
    ) +
    # Customize the shape legend
    # scale_shape_manual(
    #   name = "Naive Model:",
    #   values = c("TRUE" = 17, "FALSE" = 16), # Triangle and Circle
    #   labels = c("TRUE" = "Augmented", "FALSE" = "Unaugmented")
    # ) +
    # Add labels (using latex2exp)
    labs(
      title = TeX(paste0("Direct Effect; ", "$m^e$ = ", m_e)),
      subtitle = subtitle_text,
      y = "Estimated DE",
      x = TeX("$\\kappa$")
    ) +
    
    scale_x_continuous(breaks = seq(min(sa_de_data$kappa), max(sa_de_data$kappa), 0.25)) +
    
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 13),
      axis.text = element_text(size = 13)
    )
  
  return(de_plot)
}



pba_ie_plot <- function(pba_ie_data,
                        null_results,
                        prior_levels,
                        prior_colors,
                        model_labels,
                        model_shapes,
                        aug){
  # prepare null data
  null_results[,`:=`(
    model = spec,
    augmented = aug,
    ie_rd_mean = ie_rd,
    ie_rd_q_low = ci_low,
    ie_rd_q_high = ci_high,
    uncertainty_type = "total",
    prior = "Naive"
    
  )]
  null_results <- null_results[,.(
    uncertainty_type, ie_rd_mean, ie_rd_q_low, ie_rd_q_high, augmented, prior, model)
  ]
  
  # combine with pba results
  plot_data <- rbindlist(list(pba_ie_data, null_results))
  plot_data$prior <- factor(plot_data$prior, levels = prior_levels)
  
  # Create the figure
  pba_ie_plot_aug <- ggplot(
    # pba_ie_test,
    # pba_ie_results[augmented == TRUE & uncertainty_type == "total",],
    data = plot_data[model != "Naive",],
    aes(x = prior,
        y = ie_rd_mean,
        ymin = ie_rd_q_low,
        ymax = ie_rd_q_high,
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # Naive data
    geom_errorbar(
      data = plot_data[model == "Naive",],
      color = "black",
      width = 0.15,     
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_data[model == "Naive",],
      color = "black",
      size = 4,
      shape = 15
    ) +
    geom_errorbar(
      aes(color = prior, group = model),
      position = position_dodge(.3),
      width = 0.3, 
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = prior, shape = model, group = model),
      position = position_dodge(.3),
      size = 4
    ) +
    scale_color_manual(values = prior_colors) +
    scale_shape_manual(
      name = "Model", # Set legend title
      values = model_shapes,
      labels = model_labels
    ) +
    scale_y_continuous(breaks = seq(-0.5,0.1,0.1),
                       limits = c(-0.59,0.1)) +
    guides(
      color = "none" # Hide the color (prior) legend
    ) +
    labs(
      x = "",
      y = "Estimated IE",
      title = "Indirect Effect"
    ) +
  theme_bw(base_size = 14)  +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 14)) 

  
  return(pba_ie_plot_aug)
  
}

pba_de_plot <- function(pba_de_data,
                        null_results,
                        prior_levels,
                        prior_colors,
                        model_labels,
                        model_shapes,
                        aug){
  # prepare null data
  null_results[,`:=`(
    model = spec,
    augmented = aug,
    de_rd_mean = de_rd,
    de_rd_q_low = ci_low,
    de_rd_q_high = ci_high,
    uncertainty_type = "total",
    prior = "Naive"
    
  )]
  null_results <- null_results[,.(
    uncertainty_type, de_rd_mean, de_rd_q_low, de_rd_q_high, augmented, prior, model)
  ]
  
  # combine with pba results
  plot_data <- rbindlist(list(pba_de_data, null_results))
  plot_data$prior <- factor(plot_data$prior, levels = prior_levels)
  
  # Create the figure
  pba_de_plot_aug <- ggplot(
    # pba_de_test,
    # pba_de_results[augmented == TRUE & uncertainty_type == "total",],
    data = plot_data[model != "Naive",],
    aes(x = prior,
        y = de_rd_mean,
        ymin = de_rd_q_low,
        ymax = de_rd_q_high,
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    # Naive data
    geom_errorbar(
      data = plot_data[model == "Naive",],
      color = "black",
      width = 0.15,     
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_data[model == "Naive",],
      color = "black",
      size = 4,
      shape = 15
    ) +
    geom_errorbar(
      aes(color = prior, group = model),
      position = position_dodge(.3),
      width = 0.3, 
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = prior, shape = model, group = model),
      position = position_dodge(.3),
      size = 4
    ) +
    scale_color_manual(values = prior_colors) +
    scale_shape_manual(
      name = "Model", # Set legend title
      values = model_shapes,
      labels = model_labels
    ) +
    scale_y_continuous(breaks = seq(-0.5,0.1,0.1),
                       limits = c(-0.59,0.1)) +

    guides(
      color = "none" # Hide the color (prior) legend
    ) +
    labs(
      x = "",
      y = "Estimated DE",
      title = "Direct Effect"
    ) +
    theme_bw(base_size = 14)  +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 15),
          axis.text.y = element_text(size = 14)) 
  
  return(pba_de_plot_aug)
  
}






