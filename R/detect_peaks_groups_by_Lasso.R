

gcc_phase_data_for_lasso_source_reconstruction_plot <-
  function(gcc,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lambda,
           lag_window_in_meters,
           max_points,
           keep_if_z_is_greater_than = 5,
           keep_the_best_n = 5,
           velocity_of_sound = 334,
           freq_filter = F)
  {


    if (!"gcc_phase_std" %in% names(gcc))
      gcc <- gcc_phase_data(
        rec = gcc,
        t_start = t_start,
        t_end = t_end,
        min_freq = min_freq,
        max_freq = max_freq,
        remove_if_z_is_greater_than = keep_if_z_is_greater_than,
        freq_filter = freq_filter
      )
    nsamples = floor(length(gcc$gcc_phase[[1]][[1]])/2)*2

    lag_max=nsamples/2
    if (!is.null(lag_window_in_meters))
         lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
    gcc$lag_max = lag_max

    stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 <= nsamples)

    lags_1 = ((-lag_max):(lag_max - 1)) / gcc$fs[1]

    factor = ceiling(lag_max / max_points)

    lag_max_f = ceiling(lag_max / factor)

    lags_1f = floor(lags_1 * gcc$fs[1] / factor) * factor / gcc$fs[1]


    i_lags_1 = c((-lag_max):0 + nsamples, 1:(lag_max  - 1))


    n_1_lags = length(unique(lags_1f))
    n_receptors = length(gcc$labels)



    i_lags = 1:n_1_lags

    lags = ((-lag_max_f):(lag_max_f - 1)) / gcc$fs[1] * factor

    lags_j = lags
    lags_k = lags


    i_j_lags = Reduce(function(prev, x)
      c(prev, rep(x, length(i_lags))), x = i_lags, init = c())
    i_k_lags = rep(i_lags, length(i_lags))

    lags_ij = Reduce(function(prev, x)
      c(prev, rep(x, length(lags))), x = lags, init = c())
    lags_ik = rep(lags, length(lags))

    j_k_lags = ((lags_ik - lags_ij) * gcc$fs[1] / factor) + n_1_lags / 2 +
      1

    j_k_lags[(j_k_lags < 1) | (j_k_lags > n_1_lags)] <- NA


    gcc$plot_peak_group_lasso_data = Reduce(function(p, i)
      rbind(p,
            Reduce(
              function(pi, j)
                rbind(pi,
                      Reduce(
                        function(pij, k)
                        {
                          d <- data.frame(
                            lag_value_j = gcc$gcc_phase_std[[i]][[j - i]]$z[i_lags_1],
                            lag_value_k = gcc$gcc_phase_std[[i]][[k - i]]$z[i_lags_1],
                            lag_value_jk = gcc$gcc_phase_std[[j]][[k - j]]$z[i_lags_1],
                            lags = lags_1f
                          ) %>% dplyr::group_by(lags) %>% dplyr::summarise(
                            lag_value_j = max(lag_value_j),
                            lag_value_k = max(lag_value_k),
                            lag_value_jk = max(lag_value_jk)
                          )
                          threshold_j = sort(d$lag_value_j, decreasing = T)[keep_the_best_n]
                          threshold_k = sort(d$lag_value_k, decreasing = T)[keep_the_best_n]
                          threshold_jk = sort(d$lag_value_jk, decreasing = T)[keep_the_best_n]


                          d$lag_value_j[(d$lag_value_j < keep_if_z_is_greater_than) &
                                          (d$lag_value_j < threshold_j)] <-
                            0
                          d$lag_value_k[(d$lag_value_k < keep_if_z_is_greater_than) &
                                          (d$lag_value_k < threshold_k)] <-
                            0
                          d$lag_value_jk[(d$lag_value_jk < keep_if_z_is_greater_than) &
                                           (d$lag_value_jk < threshold_jk)] <-
                            0
                          Y <-
                            c(d$lag_value_j, d$lag_value_k, d$lag_value_jk)
                          build_ABCij_to_ABi_transform_matrix(length(d$lag_value_j)) ->
                            A
                          cv_model <-
                            gglasso::gglasso(
                              x = A$A,
                              y =  Y,
                              group = A$xy_group,
                              lambda = lambda
                            )
                          x <- as.vector(cv_model$beta[, 1])

                          x[x == 0] <- NA
                          n2 = length(i_j_lags)
                          ABij <- x[seq(from = 1,
                                        by = 3,
                                        to = n2 * 3)]
                          ACij <- x[seq(from = 2,
                                        by = 3,
                                        to = n2 * 3)]
                          BCij <- x[seq(from = 3,
                                        by = 3,
                                        to = n2 * 3)]




                          rbind(
                            pij,
                            data.frame(
                              lag_i_j = lags_ij,
                              lag_i_k = lags_ik,
                              lag_values = (d$lag_value_j[i_j_lags] +
                                              d$lag_value_k[i_k_lags] +
                                              d$lag_value_jk[j_k_lags]),
                              source_values_AB = ABij,
                              source_values_AC = ACij,
                              source_values_BC = BCij,
                              source_nvalues_AB = ABij / max(ABij, na.rm = T),
                              source_nvalues_AC = ACij / max(ACij, na.rm = T),
                              source_nvalues_BC = BCij / max(BCij, na.rm = T),
                              ijk = paste0(gcc$labels[i], "_",
                                           gcc$labels[j], "_",
                                           gcc$labels[k]),
                              i = gcc$labels[i],
                              j = gcc$labels[j],
                              k = gcc$labels[k],
                              i_j = paste0(i, "_", j, "_", gcc$labels[i], "_",
                                           gcc$labels[j]),
                              i_k = paste0(i, "_", k, "_", gcc$labels[i], "_",
                                           gcc$labels[k])

                            )
                          )
                        },
                        (j + 1):n_receptors,
                        data.frame()
                      )),
              (i + 1):(n_receptors - 1),
              data.frame()
            )),
      1:(n_receptors - 2),
      data.frame())

    return(gcc)
  }


plot_peaks_groups_by_Lasso <-
  function(gcc,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           max_points = NULL,
           lambda,
           keep_if_z_is_greater_than = 4,
           keep_the_best_n = 5,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    if (!"plot_peak_group_lasso_data" %in% names(gcc))
      gcc = gcc_phase_data_for_lasso_source_reconstruction_plot(
        gcc = gcc,
        t_start = t_start,
        t_end = t_end,
        min_freq = min_freq,
        max_freq = max_freq,
        lag_window_in_meters = lag_window_in_meters,
        lambda=lambda,
        max_points = max_points,
        keep_if_z_is_greater_than = keep_if_z_is_greater_than,
        keep_the_best_n = keep_the_best_n,
        velocity_of_sound = velocity_of_sound,
        freq_filter = freq_filter
      )



    print(
      ggplot2::ggplot(gcc$plot_peak_group_lasso_data) +
        ggplot2::geom_raster(ggplot2::aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        ggplot2::geom_point(
          data = dplyr::filter(
            gcc$plot_peak_group_lasso_data,
            !is.na(source_nvalues_AB)
          ),
          ggplot2::aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_AB,
            color = source_nvalues_AB
          ),
          alpha = 0.2
        ) +
        ggplot2::facet_grid(i_k ~ i_j) + ggplot2::scale_fill_distiller(palette = "Spectral") +
        ggplot2::scale_color_distiller(palette = "Spectral")
    )
    print(
      ggplot2::ggplot(gcc$plot_peak_group_lasso_data) +
        ggplot2::geom_raster(ggplot2::aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        ggplot2::geom_point(
          data = dplyr::filter(
            gcc$plot_peak_group_lasso_data,
            !is.na(source_nvalues_AC)
          ),
          ggplot2::aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_AC,
            color = source_nvalues_AC
          ),
          alpha = 0.2
        ) +
        ggplot2::facet_grid(i_k ~ i_j) + ggplot2::scale_fill_distiller(palette = "Spectral") +
        ggplot2::scale_color_distiller(palette = "Spectral")

    )
    print(
      ggplot2::ggplot(gcc$plot_peak_group_lasso_data) +
        ggplot2::geom_raster(ggplot2::aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        ggplot2::geom_point(
          data = dplyr::filter(
            gcc$plot_peak_group_lasso_data,
            !is.na(source_nvalues_BC)
          ),
          ggplot2::aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_BC,
            color = source_nvalues_BC
          ),
          alpha = 0.2
        ) +
        ggplot2::facet_grid(i_k ~ i_j) + ggplot2::scale_fill_distiller(palette = "Spectral") +
        ggplot2::scale_color_distiller(palette = "Spectral")

    )

  }


