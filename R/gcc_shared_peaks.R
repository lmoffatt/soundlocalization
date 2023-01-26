

obtain_lag_max <-
  function(rec,
           frame,
           lag_window_in_meters,
           velocity_of_sound=334)
  {
    lag_max = floor(frame$nsamples[1] / 4) * 2
    if (is.null(lag_max) || length(lag_max) == 0)
      lag_max = floor(min(rec$raw_nsamples) / 5) * 2

    if (!is.null(lag_window_in_meters))
      lag_max = ceiling(lag_window_in_meters / velocity_of_sound * rec$fs[1])

    return(lag_max)

  }






group_peaks_by_common_lag <- function(shared_peaks)
{
  find_duplicates_per_column <- function(lag_values)
  {
    d = data.frame(index = seq_along(lag_values), values = lag_values)
    d = dplyr::arrange(d, values)
    d$next_index = lag_values
    if (length(lag_values) > 1)
      d$next_index = dplyr::if_else(d$values == c(d$values[2:nrow(d)], d$values[1]),
                                    c(d$index[2:nrow(d)], d$index[1]),
                                    d$index)
    d = dplyr::arrange(d, index)
    return (d$next_index)
  }[]

  find_duplicates_indexes <- function(shared_peaks)
  {
    Reduce(
      f = function(d, n) {
        ni = paste0('next_index_', shared_peaks$i[n], '_', shared_peaks$j[n])

        vi = paste0('values_', shared_peaks$i[n], '_', shared_peaks$j[n])
        d[ni] = find_duplicates_per_column(shared_peaks$d[[vi]])
        return (d)
      },
      x = 1:length(shared_peaks$i),
      data.frame(index = 1:nrow(shared_peaks$d))
    )
  }

  du = find_duplicates_indexes(shared_peaks)
  ci = lapply(1:nrow(du), function(i)
    unique(as.numeric(du[i,])))
  du$connected_indexes = ci

  get_all_duplicated_indexes <-
    function(other_set, current_set, new_set)
    {
      Reduce(
        f = function(p, i) {
          new_indexes = setdiff(ci[[i]], c(other_set, p))
          if (!length(new_indexes) == 0)
            return (get_all_duplicated_indexes(other_set, c(p, new_indexes), new_indexes))
          else
            return (p)
        },
        new_set,
        current_set
      )
    }


  groups = Reduce(
    f = function(p, i) {
      if (!i %in% p$all)
      {
        new_group = get_all_duplicated_indexes(p$all, c(), ci[[i]])
        p$groups = c(p$groups, list(new_group))
        p$all = c(p$all, new_group)
      }
      return (p)
    },
    1:nrow(du),
    list(all = c(), groups = list())
  )
  stopifnot('somehow the grouping failed' = all(sort(groups$all) == 1:nrow(du)))

  return (groups$groups)

}










obtain_shared_peaks_per_frame <-

  function(rec,
           frame,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           keep_if_z_is_greater_than = 5,
           keep_the_best_n = 50,
           velocity_of_sound = 334,
           signal_filter = F,
           freq_filter = F)

  {
    lag_max = obtain_lag_max(rec,frame,lag_window_in_meters = lag_window_in_meters,
                             velocity_of_sound = velocity_of_sound)

    if (has_element_and_matches_this(
      frame,
      "shared_peaks",
      list_to_match = list(
        min_freq = min_freq,
        max_freq = max_freq,
        signal_filter = signal_filter,
        freq_filter = freq_filter,
        keep_if_z_is_greater_than = keep_if_z_is_greater_than,
        keep_the_best_n = keep_the_best_n,
        lag_max = lag_max
      )
    ))
    return(frame)



    frame <- gphase_filter_peaks_data_per_frame(
      rec,
      frame,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      lag_window_in_meters = lag_window_in_meters,
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      keep_the_best_n = keep_the_best_n,
      velocity_of_sound = velocity_of_sound,
      signal_filter = signal_filter,
      freq_filter = freq_filter
    )

    expand <- function(d_matrix, d_vector)
    {
      d = data.frame(lags = d_vector$d$lags, values = d_vector$d$lag_value)
      colnames(d) <-
        c(
          paste0("lag_", d_vector$i, "_", d_vector$j),
          paste0("values_", d_vector$i, "_", d_vector$j)
        )
      d <- dplyr::full_join(d_matrix$d, d, by = character())
      list(
        i = c(d_matrix$i, d_vector$i),
        j = c(d_matrix$j, d_vector$j),
        d = d
      )
    }

    contract <- function(d_matrix, d_vector)
    {
      d = data.frame(lags = d_vector$d$lags, values = d_vector$d$lag_value)
      colnames(d) <- c(
        paste0("lag_", d_vector$i, "_", d_vector$j),
        paste0("values_", d_vector$i, "_", d_vector$j)
      )

      dd = d_matrix$d
      lagij = paste0("lag_", d_vector$i, "_", d_vector$j)
      dd[[lagij]] =
        dd[[paste0("lag_", 1, "_", d_vector$j)]] - dd[[paste0("lag_", 1, "_", d_vector$i)]]
      d <- dplyr::inner_join(dd, d, by = lagij)
      list(
        i = c(d_matrix$i, d_vector$i),
        j = c(d_matrix$j, d_vector$j),
        d = d
      )
    }



    n_receptors = get_number_of_receptors(rec)

    frame$shared_peaks =
      Reduce (
        function (d_matrix, j)
          Reduce(
            function (d_matrixi, i)
              if (frame$gcc_peaks_data[[i]][[j - i]]$j %in% d_matrixi$j)
                contract(d_matrixi, frame$gcc_peaks_data[[i]][[j - i]])
            else
              expand(d_matrixi, frame$gcc_peaks_data[[i]][[j - i]]),
            1:(j - 1),
            d_matrix
          ),
        2:n_receptors,
        list(d = data.frame())
      )

    frame$shared_peaks$d$sum_of_values = frame$shared_peaks$d %>%
      dplyr::select(starts_with("values")) %>% rowSums()
    frame$shared_peaks$d$sum_of_log_values =
      frame$shared_peaks$d %>% dplyr::select(starts_with("values")) %>% Reduce(function(x, y)
        x + log(y), ., init = 0)
    frame$shared_peaks$d = dplyr::arrange(frame$shared_peaks$d, desc(sum_of_log_values))


    n_receptors = get_number_of_receptors(rec)


    d = frame$shared_peaks$d
    frame$shared_peaks$delays = d[paste0("lag_1_", 2:n_receptors)]
    logv = d %>% dplyr::select(!dplyr::contains("sum") &
                                 dplyr::contains("values")) %>% log()

    logR = vapply(1:n_receptors,
                  function(i) {
                    (rowSums(dplyr::select(logv, dplyr::contains(paste0(
                      i
                    )))) -
                      rowSums(dplyr::select(logv, !dplyr::contains(paste0(
                        i
                      )))) / (n_receptors - 2)) / (n_receptors - 1)
                  }, numeric(length = nrow(d)))

    frame = Reduce(function (fr, i)
    {
      if (nrow(d) > 0)
        fr$shared_peaks$delays[[paste0('logR_', i)]] = logR[((i - 1) * nrow(d)) +
                                                              (1:nrow(d))]
      return(fr)
    }, seq_len(n_receptors), frame)
    #  frame$shared_peaks$groups = numeric(0)
    #  if (nrow(frame$shared_peaks$d) > 0)
    #    frame$shared_peaks$groups = group_peaks_by_common_lag(frame$shared_peaks)
    return (frame)
  }


obtain_shared_peaks <- function(rec,
                                t_start = NULL,
                                t_end = NULL,
                                min_freq = NULL,
                                max_freq = NULL,
                                lag_window_in_meters = NULL,
                                keep_if_z_is_greater_than = 5,
                                keep_the_best_n = 50,
                                velocity_of_sound = 334,
                                signal_filter = F,
                                freq_filter = F)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames,
                        function(frame)
                          return(
                            obtain_shared_peaks_per_frame(
                              rec = rec,
                              frame = frame,
                              t_start = t_start,
                              t_end = t_end ,
                              min_freq = min_freq,
                              max_freq = max_freq,
                              lag_window_in_meters = lag_window_in_meters,
                              keep_if_z_is_greater_than = keep_if_z_is_greater_than,
                              keep_the_best_n = keep_the_best_n,
                              velocity_of_sound = velocity_of_sound,
                              signal_filter = signal_filter ,
                              freq_filter = freq_filter
                            )
                          ))
    return(rec)
  }
  else
    return(
      obtain_shared_peaks_per_frame(
        rec = rec,
        frame = rec,
        t_start = t_start,
        t_end = t_end ,
        min_freq = min_freq,
        max_freq = max_freq,
        lag_window_in_meters = lag_window_in_meters,
        keep_if_z_is_greater_than = keep_if_z_is_greater_than,
        keep_the_best_n = keep_the_best_n,
        velocity_of_sound = velocity_of_sound,
        signal_filter = signal_filter ,
        freq_filter = freq_filter
      )
    )
}



shared_peaks_data_for_plot <- function(gcc,
                                       t_start = NULL,
                                       t_end = NULL,
                                       min_freq = NULL,
                                       max_freq = NULL,
                                       lag_window_in_meters,
                                       show_greatest_n_peaks = 10,
                                       keep_if_z_is_greater_than = 5,
                                       keep_the_best_n = 5,
                                       velocity_of_sound = 334,
                                       signal_filter = T,
                                       freq_filter = F)

{
  gcc <- obtain_shared_peaks(
    rec = gcc,
    t_start = t_start,
    t_end = t_end,
    min_freq = min_freq,
    max_freq = max_freq,
    lag_window_in_meters = lag_window_in_meters,
    keep_if_z_is_greater_than = keep_if_z_is_greater_than,
    keep_the_best_n = keep_the_best_n,
    velocity_of_sound = velocity_of_sound,
    signal_filter = signal_filter,
    freq_filter = freq_filter
  )
  peak_threshold = gcc$shared_peaks$d$sum_of_log_values[min(show_greatest_n_peaks, nrow(gcc$shared_peaks$d))]

  shared_peaks = gcc$shared_peaks$d %>%
    dplyr::filter(sum_of_log_values >= peak_threshold)

  max_lag = gcc$lag_max / gcc$fs[1]
  min_lag = -max_lag

  n_receptors = get_number_of_receptors(gcc)

  only_in_range <- function(min, max, x)
  {
    ifelse((x >= min) & (x <= max), x, NA)
  }




  diag_lag <-
    function(min_lag_x,
             max_lag_x,
             min_lag_y,
             max_lag_y,
             lag_jk)
    {
      Reduce(function(p, lag)
      {
        x_min = min_lag_x
        y_min = min_lag_y
        x_max = max_lag_x
        y_max = max_lag_y

        # the diagonal is outside the inner square
        if ((min_lag_x + lag > max_lag_y) |
            (max_lag_x + lag < min_lag_y))
        {
          x_min = NA
          y_min = NA
          x_max = NA
          y_max = NA
        }
        else
        {
          if (min_lag_x + lag > min_lag_y)
            y_min = min_lag_x + lag
          else
            x_min = min_lag_y - lag
          if (max_lag_x + lag < max_lag_y)
            y_max = max_lag_x + lag
          else
            x_max = max_lag_y - lag
        }
        return (list(
          x = c(p$x, x_min, x_max, NA),
          y = c(p$y, y_min, y_max, NA)
        ))
      },
      lag_jk,
      list(x = c(), y = c()))
    }



  get_xy_lags <-
    function(min_lag_x,
             max_lag_x,
             min_lag_y,
             max_lag_y,
             lag_j,
             lag_k,
             lag_jk)
    {
      diagxy = diag_lag(min_lag_x, max_lag_x, min_lag_y, max_lag_y,  lag_jk)
      return(list(
        x = c(rep(
          only_in_range(min_lag_x, max_lag_x, lag_j), each = 3
        ),
        rep(
          c(min_lag_x, max_lag_x, NA), length(lag_k)
        ),
        diagxy$x),
        y =
          c(rep(
            c(min_lag_y, max_lag_y, NA), length(lag_j)
          ),
          rep(
            only_in_range(min_lag_y, max_lag_y, lag_k), each = 3
          ),
          diagxy$y)
      ))
    }

  gcc$shared_peaks_for_plot = Reduce(
    function(p, i)
      Reduce(
        function(pi, j)
          Reduce(function(pij, k)
          {
            lag_j = gcc$gcc_peaks_data[[i]][[j - i]]$d$lags
            lag_k = gcc$gcc_peaks_data[[i]][[k - i]]$d$lags
            lag_jk = gcc$gcc_peaks_data[[j]][[k - j]]$d$lags

            lag_value_j = gcc$gcc_peaks_data[[i]][[j - i]]$d$lag_value
            lag_value_k = gcc$gcc_peaks_data[[i]][[k - i]]$d$lag_value
            lag_value_jk = gcc$gcc_peaks_data[[j]][[k - j]]$d$lag_value

            lag_shared_i_j = shared_peaks[[paste0("lag_", i, "_", j)]]
            lag_shared_i_k = shared_peaks[[paste0("lag_", i, "_", k)]]

            lag_value_shared_i_j = shared_peaks[[paste0("values_", i, "_", j)]]
            lag_value_shared_i_k = shared_peaks[[paste0("values_", i, "_", k)]]
            lag_value_shared_j_k = shared_peaks[[paste0("values_", j, "_", k)]]

            min_lag_x = min(lag_j)
            max_lag_x = max(lag_j)

            min_lag_y = min(lag_k)
            max_lag_y = max(lag_k)

            min_shared_lag_x = min(lag_shared_i_j)
            max_shared_lag_x = max(lag_shared_i_j)

            min_shared_lag_y = min(lag_shared_i_k)
            max_shared_lag_y = max(lag_shared_i_k)





            label_ijk = paste0(gcc$labels[i], "_",
                               gcc$labels[j], "_",
                               gcc$labels[k])
            label_i = gcc$labels[i]
            label_j = gcc$labels[j]
            label_k = gcc$labels[k]
            label_i_j = paste0(i, "_", j, "_", gcc$labels[i], "_",
                               gcc$labels[j])
            label_i_k = paste0(i, "_", k, "_", gcc$labels[i], "_",
                               gcc$labels[k])
            label_j_k = paste0(j, "_", k, "_", gcc$labels[j], "_",
                               gcc$labels[k])


            axis = c(rep("i_j", length(lag_j) * 3),
                     rep("i_k", length(lag_k) * 3),
                     rep("j_k", length(lag_jk) * 3))

            axis_labels = c(rep(label_i_j, length(lag_j) * 3),
                            rep(label_i_k, length(lag_k) *
                                  3),
                            rep(label_j_k, length(lag_jk) *
                                  3))

            # x_lags = c(rep(lag_j, each = 3),
            #
            #            rep(c(min_lag, max_lag, NA), length(lag_k)),
            #            diag_lag_x(min_lag, max_lag, lag_jk))
            #
            # y_lags = c(rep(c(min_lag, max_lag, NA), length(lag_j)),
            #            rep(lag_k, each = 3),
            #            diag_lag_y(min_lag, max_lag, lag_jk))
            value_lags = c(rep(lag_value_j, each = 3),
                           rep(lag_value_k, each = 3),
                           rep(lag_value_jk, each = 3))

            xy_lags_all = get_xy_lags(min_lag, max_lag, min_lag, max_lag, lag_j , lag_k, lag_jk)
            xy_lags_set = get_xy_lags(min_lag_x,
                                      max_lag_x,
                                      min_lag_y,
                                      max_lag_y,
                                      lag_j ,
                                      lag_k,
                                      lag_jk)
            xy_lags_shared = get_xy_lags(
              min_shared_lag_x,
              max_shared_lag_x,
              min_shared_lag_y,
              max_shared_lag_y,
              lag_j ,
              lag_k,
              lag_jk
            )



            d_lines = data.frame(
              axis = axis,

              axis_labels = axis_labels,

              x_lags_all = xy_lags_all$x ,
              y_lags_all = xy_lags_all$y,

              x_lags_set = xy_lags_set$x,
              y_lags_set = xy_lags_set$y,

              x_lags_shared = xy_lags_shared$x,
              y_lags_shared = xy_lags_shared$y,

              value_lags = value_lags,
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

            d_shared = data.frame(
              lag_shared_i_j = lag_shared_i_j,
              lag_shared_i_k = lag_shared_i_k,

              lag_value_shared_i_j = lag_value_shared_i_j,
              lag_value_shared_i_k = lag_value_shared_i_k,
              lag_value_shared_j_k = lag_value_shared_j_k,

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

            return(list(
              d_lines = rbind(pij$d_lines, d_lines),
              d_shared = rbind(pij$d_shared, d_shared)
            ))
          },
          (j + 1):n_receptors,
          pi),
        (i + 1):(n_receptors - 1),
        p
      ),
    1:(n_receptors - 2),
    list(d_lines = data.frame(), d_shared = data.frame())
  )

  return (gcc)
}


plot_shared_peaks <- function(gcc,
                              t_start = NULL,
                              t_end = NULL,
                              min_freq = NULL,
                              max_freq = NULL,
                              lag_window_in_meters,
                              show_greatest_n_peaks = 10,
                              keep_if_z_is_greater_than = 5,
                              keep_the_best_n = 200,
                              velocity_of_sound = 334,
                              freq_filter = F)
{
  if (!"shared_peaks_for_plot" %in% names(gcc))
    gcc <- shared_peaks_data_for_plot(
      gcc,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      lag_window_in_meters = lag_window_in_meters,
      show_greatest_n_peaks = show_greatest_n_peaks,
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      keep_the_best_n = keep_the_best_n,
      velocity_of_sound = velocity_of_sound,
      freq_filter = freq_filter
    )
  print(
    ggplot2::ggplot(gcc$shared_peaks_for_plot$d_lines) +
      ggplot2::geom_path(ggplot2::aes(x_lags_all, y_lags_all), alpha = 0.1) +
      ggplot2::facet_grid(i_k ~ i_j, scales = "free") + ggplot2::scale_color_distiller(palette = "Spectral") +
      ggplot2::geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        ggplot2::aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )
  print(
    ggplot2::ggplot(gcc$shared_peaks_for_plot$d_lines) +
      ggplot2::geom_path(ggplot2::aes(x_lags_set, y_lags_set), alpha = 0.1) +
      ggplot2::facet_grid(i_k ~ i_j, scales = "free") + ggplot2::scale_color_distiller(palette = "Spectral") +
      ggplot2::geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        ggplot2::aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )
  print(
    ggplot2::ggplot(gcc$shared_peaks_for_plot$d_lines) +
      ggplot2::geom_path(ggplot2::aes(x_lags_shared, y_lags_shared), alpha = 0.1) +
      ggplot2::facet_grid(i_k ~ i_j, scales = "free") + ggplot2::scale_color_distiller(palette = "Spectral") +
      ggplot2::geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        ggplot2::aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )

}
