gcc_phase_data_for_tri_plot <-
  function(gcc,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
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

    nsamples = length(gcc$gcc_phase[[1]][[1]])
    lag_max=nsamples
    if (!is.null(lag_window_in_meters))
      lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
    gcc$lag_max = lag_max

    stopifnot("the sampled frame is wider than the lag window" = lag_max * 4 < nsamples)

    lags_2 = ((-lag_max * 2 - 1):(lag_max * 2)) / gcc$fs[1]

    factor = ceiling(2 * lag_max / max_points)

    lag_max_f = ceiling(lag_max / factor)

    lags_2f = floor(lags_2 * gcc$fs[1] / factor) * factor / gcc$fs[1]


    i_lags_2 = c((-lag_max * 2):0 + nsamples, 1:(lag_max * 2 + 1))


    n_2_lags = length(unique(lags_2f))
    n_receptors = get_number_of_receptors(gcc)

    i_lags = (n_2_lags / 4):(n_2_lags * 3 / 4 - 1)
    lags = ((-lag_max_f):(lag_max_f - 1)) / gcc$fs[1] * factor

    lags_j = lags
    lags_k = lags


    i_j_lags = Reduce(function(prev, x)
      c(prev, rep(x, length(i_lags))), x = i_lags, init = c())
    i_k_lags = rep(i_lags, length(i_lags))

    lags_ij = Reduce(function(prev, x)
      c(prev, rep(x, length(lags))), x = lags, init = c())
    lags_ik = rep(lags, length(lags))

    j_k_lags = ((lags_ik - lags_ij) * gcc$fs[1] / factor) + n_2_lags / 2 +
      1

    gcc$lag_max = lag_max
    gcc$gcc_data_for_tri_plot = Reduce(function(p, i)
      rbind(p,
            Reduce(
              function(pi, j)
                rbind(pi,
                      Reduce(
                        function(pij, k)
                        {
                          d <- data.frame(
                            lag_value_j = gcc$gcc_phase_std[[i]][[j - i]]$z[i_lags_2],
                            lag_value_k = gcc$gcc_phase_std[[i]][[k - i]]$z[i_lags_2],
                            lag_value_jk = gcc$gcc_phase_std[[j]][[k - j]]$z[i_lags_2],
                            lags = lags_2f
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

                          rbind(
                            pij,
                            data.frame(
                              lag_i_j = lags_ij,
                              lag_i_k = lags_ik,
                              lag_values = (d$lag_value_j[i_j_lags] +
                                              d$lag_value_k[i_k_lags] +
                                              d$lag_value_jk[j_k_lags]),
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


gcc_phase_tri_plot <-
  function(gcc,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           max_points = NULL,
           keep_if_z_is_greater_than = 4,
           keep_the_best_n = 5,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    if (!"gcc_data_for_tri_plot" %in% names(gcc))
      gcc = gcc_phase_data_for_tri_plot(
        gcc = gcc,
        t_start = t_start,
        t_end = t_end,
        min_freq = min_freq,
        max_freq = max_freq,
        lag_window_in_meters = lag_window_in_meters,
        max_points = max_points,
        keep_if_z_is_greater_than = keep_if_z_is_greater_than,
        keep_the_best_n = keep_the_best_n,
        velocity_of_sound = velocity_of_sound,
        freq_filter = freq_filter
      )


    ggplot2::ggplot(gcc$gcc_data_for_tri_plot) +
      ggplot2::geom_raster(ggplot2::aes(x = lag_i_j, y = lag_i_k, fill = lag_values)) +
      ggplot2::facet_grid(i_k ~ i_j) + ggplot2::scale_fill_distiller(palette = "Spectral")
  }


Tri_Delay_to_gcc_phat_matrix <- function(number_of_lags)
{
  # bi-delay-tripower X array
  #
  n_receptors = 3

  Xa = array(
    dim = c(number_of_lags, number_of_lags, n_receptors),
    dimnames = list("i_lag", "k_lag", receptor_pair = c("AB", "AC", "BC"))
  )

  X = matrix(nrow = 1,
             ncol = number_of_lags * number_of_lags * n_receptors)

  Ya = array(
    dim = c(number_of_lags, n_receptors),
    dimnames = list("lag", receptor_pair = c("AB", "AC", "BC"))
  )

  Y = matrix(nrow = number_of_lags * n_receptors, ncol = 1)

  Aa = array(dim = c(
    number_of_lags,
    n_receptors,
    number_of_lags,
    number_of_lags,
    n_receptors
  ))

  A = matrix(nrow =)

}

