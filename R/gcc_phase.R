calculate_list_of_gcc_phase <-
  function(ffted_signal,
           remove_if_z_is_greater_than = 5,
           freq_filter = FALSE) {
    stopifnot(length(unique(ffted_signal$fs)) == 1)
    if (!"fft" %in% names(ffted_signal))
      ffted_signal <- calculate_ffts(ffted_signal)

    n = length(ffted_signal$fft)

    nsamples = length(ffted_signal$fft[[1]])
    fs = ffted_signal$fs[1]
    freq = c((0:(nsamples / 2 - 1)), ((-nsamples / 2):0) - 1) * fs / nsamples

    filter = 1
    if (freq_filter)
      filter = ifelse((abs(freq) > ffted_signal$min_freq) &
                        (abs(freq) < ffted_signal$max_freq),
                      1,
                      0)
    #
    ffted_signal$gcc_phase =
      lapply(1:(n - 1),
             function(i)
               lapply((i + 1):n,
                      function(j)
                        Re(
                          fft(
                            Conj(ffted_signal$fft[[i]]) *
                              ffted_signal$fft[[j]] /
                              Mod(Conj(ffted_signal$fft[[i]]) *
                                    ffted_signal$fft[[j]]) *
                              filter
                            ,
                            inverse = T
                          )
                        ) / length(ffted_signal$fft[[i]])))

    ffted_signal$gcc_phase_std =
      lapply(ffted_signal$gcc_phase,
             function(xi)
               lapply(xi,
                      function(xij)
                        standarization_blind_to_ouliers(x = xij,
                                                        remove_if_z_is_greater_than = remove_if_z_is_greater_than)))
    return(ffted_signal)
  }

gcc_phase_data <-
  function(rec,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           remove_if_z_is_greater_than = 5,
           freq_filter = F)

  {
    if (!"duration"%in% names(rec))
      rec=truncate_to_same_length(rec)

    if (!(is.null(t_start) & is.null(t_end)))
    {
      if (is.null(t_start))
        t_start = 0
      if (is.null(t_end))
        t_end = rec$duration
      i_start = floor(t_start * rec$fs[1]) + 1
      nsamples = floor((t_end - t_start) * rec$fs[1] / 2) * 2
      rec <- rec %>%
        apply_frame(framed_by = "gcc_phase_data",
                    i_start = i_start, nsamples = nsamples)
    }
    if (!(is.null(min_freq) & is.null(max_freq)))
    {
      if (is.null(min_freq))
        min_freq = 0
      if (is.null(max_freq))
        max_freq = rec$fs[1] / 2
      rec <- rec %>%
        apply_filter(min_freq = min_freq, max_freq = max_freq)
    }
    rec %>%
      calculate_list_of_gcc_phase(freq_filter = freq_filter,
                                  remove_if_z_is_greater_than = remove_if_z_is_greater_than) ->
      gcc
    return(gcc)
  }






gcc_phase_data_for_plot <-
  function(gcc,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters,
           max_points,
           remove_if_z_is_greater_than = 5,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    data_for_plot_gcc_phase <-
      function(gcc,
               lag_max=NULL,
               lag_window_in_meters =NULL,
               max_points,
               decimating_operation = "max",
               velocity_of_sound = 334)
      {
        nsamples = length(gcc$gcc_phase[[1]][[1]])
        if (!is.null(lag_window_in_meters))
          lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
        else if(is.null(lag_max))
          lag_max=nsamples
        gcc$lag_max = lag_max
        stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 < nsamples)
        lags = ((-lag_max):(lag_max - 1)) / gcc$fs[1]
        i_lags = c((-lag_max):0 + nsamples, 1:(lag_max - 1))

        nsamples = length(gcc$gcc_phase[[1]][[1]])
        stopifnot("the sampled frame is wider than the lag window" = lag_max * 4 < nsamples)

        factor = ceiling((2 * lag_max + 1) / max_points)

        lag_max_f = lag_max / factor

        n_receptors = length(gcc$labels)

        lags_f = round(lags * gcc$fs[1] / factor) * factor / gcc$fs[1]

        lags = unique(lags_f)


        gcc$gcc_phase_data_for_plot = Reduce(function(p, i)
          rbind(p,
                Reduce(
                  function(pi, j)
                  {
                    d <- data.frame(lag_value = gcc$gcc_phase_std[[i]][[j]]$z[i_lags],
                                    lags = lags_f) %>% dplyr::group_by(lags) %>%
                      dplyr::summarise(lag_value = match.fun(decimating_operation)(lag_value))

                    rbind(
                      pi,
                      data.frame(
                        lags = lags,
                        std_lag_values = d$lag_value ,

                        i = gcc$labels[i],
                        j = gcc$labels[i + j]
                      )
                    )
                  },
                  1:length(gcc$gcc_phase_std[[i]]),
                  data.frame()
                )),
          1:length(gcc$gcc_phase_std),
          data.frame())

        return(gcc)

      }
    if (!"gcc_phase_std" %in% names(gcc))
      gcc <- gcc_phase_data(
        rec = gcc,
        t_start = t_start,
        t_end = t_end,
        min_freq = min_freq,
        max_freq = max_freq,
        freq_filter = freq_filter
      )
    s = min_freq

    gcc <-
      data_for_plot_gcc_phase(
        gcc,
        max_points = max_points,
        lag_window_in_meters = lag_window_in_meters,
        velocity_of_sound = velocity_of_sound
      )
    return(gcc)
  }




gcc_phase_plot <-
  function(rec,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           max_points = 1000,
           remove_if_z_is_greater_than = 5,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    if (!"gcc_phase_data_for_plot" %in% names(rec))
    {
      rec = gcc_phase_data_for_plot(
        gcc = rec,
        t_start = t_start,
        t_end = t_end,
        min_freq = min_freq,
        max_freq = max_freq,
        lag_window_in_meters = lag_window_in_meters,
        max_points = max_points,
        remove_if_z_is_greater_than = remove_if_z_is_greater_than,
        velocity_of_sound = velocity_of_sound,
        freq_filter = freq_filter
      )
    }
    ggplot2::ggplot(rec$gcc_phase_data_for_plot) + ggplot2::geom_line(ggplot2::aes(
      x = lags,
      y = std_lag_values,
      color = paste0(i, j)
    )) + ggplot2::facet_grid(i ~ j, scales = "free", space = "free")
  }



























gphase_filter_peaks_data <- function(gcc,
                                     t_start = NULL,
                                     t_end = NULL,
                                     min_freq = NULL,
                                     max_freq = NULL,
                                     lag_window_in_meters,
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
  lag_max=floor(nsamples/4)*2
  if (!is.null(lag_window_in_meters))
    lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])

  gcc$lag_max = lag_max

  nsamples = length(gcc$gcc_phase[[1]][[1]])
  stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 < nsamples)

  gcc$lags = ((-lag_max):(lag_max - 1)) / gcc$fs[1]

  n_receptors = length(gcc$labels)
  i_lags = c((-lag_max):0 + nsamples, 1:(lag_max  - 1))
  gcc$gcc_peaks_data = lapply(1:(n_receptors - 1),
                              function (i)
                                lapply((i + 1):n_receptors, function(j)
                                {
                                  d = data.frame(lags = gcc$lags,
                                                 lag_value = gcc$gcc_phase_std[[i]][[j - i]]$z[i_lags])
                                  threshold = min(sort(d$lag_value, decreasing = T)[keep_the_best_n],
                                                  keep_if_z_is_greater_than)
                                  d <-
                                    d %>% dplyr::filter(lag_value > threshold)

                                  return(list(i = i,
                                              j = j,
                                              d = d))
                                }))

  return(gcc)
}


