



has_element_and_matches_this <-
  function(frame, element, list_to_match)
  {
    if (!all(c(element, names(list_to_match)) %in% names(frame)))
      return (FALSE)
    for (name in names(list_to_match))
      if (!frame[[name]] == list_to_match[[name]])
        return (FALSE)
    return (TRUE)

  }





calculate_list_of_gcc_phase_per_frame <-
  function(rec,
           frame,
           remove_if_z_is_greater_than = 5,
           freq_filter = FALSE)
  {
    stopifnot(length(unique(rec$fs)) == 1)

    n = length(frame$fft)

    nsamples = length(frame$fft[[1]])
    fs = rec$fs[1]
    freq = c((0:(nsamples / 2 - 1)), ((-nsamples / 2 + 1):0) - 1) * fs / nsamples

    filter = 1
    if (freq_filter)
      filter = ifelse((abs(freq) > max(0,frame$min_freq)) &
                        (abs(freq) < min(fs/2,frame$max_freq)),
                      1,
                      0)

    frame$freq_filter = freq_filter
    frame$gcc_phase =
      lapply(1:(n - 1),
             function(i)
               lapply((i + 1):n,
                      function(j)
                        Re(
                          fft(
                            Conj(frame$fft[[i]]) *
                              frame$fft[[j]] /
                              Mod(Conj(frame$fft[[i]]) *
                                    frame$fft[[j]]) *
                              filter
                            ,
                            inverse = T
                          )
                        ) / length(frame$fft[[i]])))

    frame$gcc_phase_std =
      lapply(frame$gcc_phase,
             function(xi)
               lapply(xi,
                      function(xij)
                        standarization_blind_to_ouliers(x = xij,
                                                        remove_if_z_is_greater_than = remove_if_z_is_greater_than)))
    return(frame)
  }



calculate_list_of_gcc_phase <- function(rec,
                                        remove_if_z_is_greater_than = 5,
                                        freq_filter = FALSE)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames,
                        function(frame)
                          calculate_list_of_gcc_phase_per_frame(rec,
                                                                frame,
                                                                remove_if_z_is_greater_than,
                                                                freq_filter))
    return (rec)
  }
  else
    return (
      calculate_list_of_gcc_phase_per_frame(rec, rec,
                                            remove_if_z_is_greater_than,
                                            freq_filter)
    )
}







gcc_phase_data_per_frame <-
  function(rec,
           frame,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           remove_if_z_is_greater_than = 5,
           signal_filter = T,
           freq_filter = F)

  {
    if (has_element_and_matches_this(
      frame,
      element = "gcc_phase_std",
      list_to_match = list(
        min_freq = min_freq,
        max_freq = max_freq,
        signal_filter= signal_filter,
        freq_filter = freq_filter
      )
    ))
      return (frame)

    if (!"duration" %in% names(frame))
      frame = truncate_to_same_length(frame)

    if (!(is.null(t_start) & is.null(t_end)))
    {
      if (is.null(t_start))
        t_start = 0
      if (is.null(t_end))
        t_end = frame$duration
      i_start = floor(t_start * rec$fs[1]) + 1
      nsamples = floor((t_end - t_start) * rec$fs[1] / 2) * 2
      frame <- frame %>%
        apply_frame(framed_by = "gcc_phase_data",
                    i_start = i_start,
                    nsamples = nsamples)
    }
    frame$signal_filter=signal_filter
    if (signal_filter && !(is.null(min_freq) & is.null(max_freq)))
    {
      if (is.null(min_freq))
        min_freq = 0
      if (is.null(max_freq))
        max_freq = rec$fs[1] / 2
      frame$min_freq = min_freq
      frame$max_freq = max_freq
      frame <-
        apply_filter_per_frame(rec, frame, min_freq = min_freq, max_freq = max_freq)
      frame <- calculate_ffts_per_frame(rec, frame)

    }
    if (! signal_filter)
    {
      frame$min_freq = min_freq
      frame$max_freq = max_freq

    }
    if (!"fft" %in% names(frame))
      frame <- calculate_ffts_per_frame(rec, frame)

    frame <- calculate_list_of_gcc_phase_per_frame(rec,
                                                   frame,
                                                   freq_filter = freq_filter,
                                                   remove_if_z_is_greater_than = remove_if_z_is_greater_than)
    return(frame)
  }





gcc_phase_data <- function(rec,
                           t_start = NULL,
                           t_end = NULL,
                           min_freq = NULL,
                           max_freq = NULL,
                           remove_if_z_is_greater_than = 5,
                           signal_filter=F,
                           freq_filter = F)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames,
                        function(frame)
                          return(
                            gcc_phase_data_per_frame(
                              rec = rec,
                              frame = frame,
                              t_start = t_start,
                              t_end = t_end ,
                              min_freq = min_freq,
                              max_freq = max_freq,
                              remove_if_z_is_greater_than = remove_if_z_is_greater_than,
                              signal_filter = signal_filter,
                              freq_filter = freq_filter
                            )
                          ))
    return (rec)
  }
  else
    return (
      gcc_phase_data_per_frame(
        rec,
        rec,
        t_start,
        t_end ,
        min_freq,
        max_freq,
        remove_if_z_is_greater_than,
        signal_filter,
        freq_filter
      )
    )
}



gcc_phase_data_for_plot_per_frame <-
  function(gcc,
           frame,
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
    data_for_plot_gcc_phase_per_fame <-
      function(gcc,
               frame,
               lag_max = NULL,
               lag_window_in_meters = NULL,
               max_points,
               decimating_operation = "max",
               velocity_of_sound = 334)
      {
        nsamples = length(frame$nsamples)
        if (!is.null(lag_window_in_meters))
          lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
        else if (is.null(lag_max))
          lag_max = nsamples / 2
        frame$lag_max = lag_max
        stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 <= nsamples)
        lags = ((-lag_max):(lag_max - 1)) / gcc$fs[1]
        i_lags = c((-lag_max):0 + nsamples, 1:(lag_max - 1))

        nsamples = length(frame$gcc_phase[[1]][[1]])
        stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 <=
                    nsamples)

        factor = ceiling((2 * lag_max + 1) / max_points)

        lag_max_f = lag_max / factor

        n_receptors = length(gcc$labels)

        lags_f = round(lags * gcc$fs[1] / factor) * factor / gcc$fs[1]

        lags = unique(lags_f)


        frame$gcc_phase_data_for_plot = Reduce(function(p, i)
          rbind(p,
                Reduce(
                  function(pi, j)
                  {
                    d <- data.frame(lag_value = frame$gcc_phase_std[[i]][[j]]$z[i_lags],
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
                  1:length(frame$gcc_phase_std[[i]]),
                  data.frame()
                )),
          1:length(frame$gcc_phase_std),
          data.frame())

        return(frame)

      }
    if (all(c("gcc_phase_data_for_plot", "min_freq", "max_freq") %in% names(frame)) &&
        (frame$min_freq == min_freq) &&
        (frame$max_freq == max_freq))
      return (frame)

    frame <- gcc_phase_data_per_frame(
      rec = gcc,
      frame = frame,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      freq_filter = freq_filter
    )
    s = min_freq

    frame <-
      data_for_plot_gcc_phase_per_fame(
        gcc,
        frame,
        max_points = max_points,
        lag_window_in_meters = lag_window_in_meters,
        velocity_of_sound = velocity_of_sound
      )
    return(frame)
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
               lag_max = NULL,
               lag_window_in_meters = NULL,
               max_points,
               decimating_operation = "max",
               velocity_of_sound = 334)
      {
        nsamples = length(gcc$gcc_phase[[1]][[1]])
        if (!is.null(lag_window_in_meters))
          lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
        else if (is.null(lag_max))
          lag_max = nsamples
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

gcc_phase_plot_per_frame <-
  function(rec,
           frame,
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
    if (!"gcc_phase_data_for_plot" %in% names(frame))
    {
      frame = gcc_phase_data_for_plot_per_frame(
        gcc = rec,
        frame = frame,
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
    ggplot2::ggplot(frame$gcc_phase_data_for_plot) + ggplot2::geom_line(ggplot2::aes(
      x = lags,
      y = std_lag_values,
      color = paste0(i, j)
    )) + ggplot2::facet_grid(i ~ j, scales = "free", space = "free")
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



























gphase_filter_peaks_data_per_frame <- function(rec,
                                               frame,
                                               t_start = NULL,
                                               t_end = NULL,
                                               min_freq = NULL,
                                               max_freq = NULL,
                                               lag_window_in_meters,
                                               keep_if_z_is_greater_than = 5,
                                               keep_the_best_n = 5,
                                               velocity_of_sound = 334,
                                               signal_filter=T,
                                               freq_filter = F)
{
  lag_max = obtain_lag_max(rec,frame,lag_window_in_meters, velocity_of_sound)

  if ((has_element_and_matches_this(
    frame,
    "gcc_peaks_data",
    list_to_match = list(
      min_freq = min_freq,
      max_freq = max_freq,
      signal_filter=signal_filter,
      freq_filter = freq_filter,
      keep_the_best_n = keep_the_best_n,
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      lag_max = lag_max
    )
  )))
  return (frame)

  frame <- gcc_phase_data_per_frame(
    rec = rec,
    frame = frame,
    t_start = t_start,
    t_end = t_end,
    min_freq = min_freq,
    max_freq = max_freq,
    remove_if_z_is_greater_than = keep_if_z_is_greater_than,
    signal_filter = signal_filter,
    freq_filter = freq_filter
  )

  frame$lag_max = lag_max
  nsamples=frame$nsamples[1]
  stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 <= nsamples)

  frame$lags = ((-lag_max):(lag_max - 1)) / rec$fs[1]

  n_receptors = length(rec$labels)
  i_lags = c((-lag_max + 1):0 + nsamples, 1:(lag_max  ))
  frame$gcc_peaks_data =
    lapply(1:(n_receptors - 1),
           function (i)
             lapply((i + 1):n_receptors, function(j)
             {
               d = data.frame(lags = frame$lags,
                              lag_value = frame$gcc_phase_std[[i]][[j - i]]$z[i_lags])
               threshold = min(sort(d$lag_value, decreasing = T)[keep_the_best_n],
                               keep_if_z_is_greater_than)
               d <-
                 d %>% dplyr::filter(lag_value > threshold)

               return(list(i = i,
                           j = j,
                           d = d))
             }))

  frame$keep_the_best_n = keep_the_best_n
  frame$keep_if_z_is_greater_than = keep_if_z_is_greater_than
  return(frame)
}


gphase_filter_peaks_data <- function(rec,
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
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames,
                        function(frame)
                          return(
                            gphase_filter_peaks_data_per_frame(
                              rec,
                              frame,
                              t_start,
                              t_end ,
                              min_freq,
                              max_freq,
                              lag_window_in_meters,
                              keep_if_z_is_greater_than,
                              keep_the_best_n,
                              velocity_of_sound,
                              freq_filter
                            )
                          ))
    return (rec)
  }
  else
    return (
      gphase_filter_peaks_data_per_frame(
        rec,
        rec,
        t_start,
        t_end ,
        min_freq,
        max_freq,
        lag_window_in_meters,
        keep_if_z_is_greater_than,
        keep_the_best_n,
        velocity_of_sound,
        freq_filter
      )
    )
}
