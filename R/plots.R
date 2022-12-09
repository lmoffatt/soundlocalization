frequency_plot <- function(signal,
                           f_start = 0,
                           f_end = Inf,
                           maxpoints = 1e4,
                           alpha = 0.7,
                           scale_x_log = F,
                           dB = F)
{
  n = length(signal$fft)
  fs = signal$fs
  f_end = min(fs / 2, f_end)
  nsamples = length(signal$fft[[1]])
  df = fs / nsamples
  i_start = f_start / df + 1
  i_end = f_end / df + 1

  d = data.frame(
    fft = Reduce(function(v, i)
      c(
        v, downsample_signal(signal$fft[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c()),
    freq = Reduce(function(v, i)
      c(
        v,
        downsample_frequency(df, signal$signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c()),
    file = Reduce(function(v, i)
      c(
        v,
        downsample_label(signal$filename[i],
                         signal$signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())


  )
  d$dB = 10 * log10(abs(d$fft) / max(abs(d$fft)))
  if (ggplot2::scale_x_log)
  {
    ggplot2::ggplot(d) +
      ggplot2::geom_line(ggplot2::aes(freq, dB, color = file), alpha = alpha) +
      ggplot2::scale_x_log10() +
      ggplot2::facet_wrap(~ file, scales = "free")
  }
  else if (dB)
    ggplot2::ggplot(d) +
    ggplot2::geom_line(ggplot2::aes(freq, dB, color = file), alpha = alpha) +
    ggplot2::facet_wrap(~ file, scales = "free")
  else
    ggplot2::ggplot(d) +
    ggplot2::geom_line(ggplot2::aes(freq, abs(fft), color = file), alpha = alpha) +
    ggplot2::facet_wrap(~ file, scales = "free")

}



#' Title
#'
#' @param recordings recordings object
#' @param t_start  start of the plot measured from the beginning of the recordings
#' @param t_end   end of the plot measured from the beginning of the recordings
#' @param maxpoints  number of points to be shown
#' @param alpha  transparency (1= no transparency)
#'
#' @return
#' @export
#'
#' @examples
sound_plot <-
  function(recordings,
           t_start = 0,
           t_end = Inf,
           maxpoints = 1e4,
           alpha = 0.7)
  {

    signal=recordings$raw_signal
    if ("signal_for_fft"%in% names(recordings))
      signal=recordings$signal_for_fft
    else  if ("framed_raw_signal"%in% names(recordings))
      signal=recordings$framed_raw_signal

    n=length(signal)
    nsamples=recordings$raw_nsamples
    fs = recordings$fs
    i_start = t_start * fs
    i_end = t_end * fs
    if (is.infinite(t_end))
    {
      t_end=n/fs
      i_end=nsamples
    }

    s = Reduce(function(v, i)
      c(
        v,
        downsample_signal(signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())
    t = Reduce(function(v, i)
      c(
        v,
        downsample_time(fs, signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())
    receptor = Reduce(function(v, i)
      c(
        v,
        downsample_label(recordings$filename[i],
                         signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())

    d = data.frame(signal = s,
                   t = t,
                   receptor = receptor)

    ggplot2::ggplot(d) +
      ggplot2::geom_line(ggplot2::aes(t, signal), alpha = alpha) +
      ggplot2::facet_wrap(~ receptor)
  }







cross_plot <-
  function(i,
           j,
           rec,
           t_start,
           t_end,
           min_freq,
           max_freq,
           lag_window_in_meters,
           remove_if_z_is_greater_than = 5,
           lag_error = 10,
           lag_filter_factor = 0.9,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    fs = rec$fs[1]
    rec = select_receptors(rec, c(i, j))
    gcc = gcc_phase_data(
      rec,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      remove_if_z_is_greater_than = remove_if_z_is_greater_than,
      lag_window_in_meters = lag_window_in_meters,
      velocity_of_sound = velocity_of_sound,
      freq_filter = freq_filter
    )
    d = gcc$data_for_plot
    best_lags = d$lags[d$r_lag_values > lag_filter_factor]
    best_lag_value = d$r_lag_values[d$r_lag_values > lag_filter_factor]
    lapply(best_lags, function(lag)
      filter_by_lag(
        fft1 = gcc$fft[[1]],
        fft2 = gcc$fft[[2]],
        lag = lag,
        lag_error = lag_error / fs,
        fs = fs
      ))

  }



sound_plot_frames <-
  function(recordings,
           frame,
           maxpoints = 1e4,
           alpha = 0.7,
           for_fft=FALSE)
  {

    signal=frame$framed_raw_signal
    if (for_fft)
      signal=frame$signal_for_fft

    n=length(signal)
    nsamples=frame$nsamples
    fs = recordings$fs


    s = Reduce(function(v, i)
      c(
        v,
        downsample_signal(signal[[i]], maxpoints)
      )
      , 1:n, c())
    t = Reduce(function(v, i)
      c(
        v,
        downsample_time(fs, signal[[i]], maxpoints)
      )
      , 1:n, c())
    receptor = Reduce(function(v, i)
      c(
        v,
        downsample_label(recordings$filename[i],
                         signal[[i]], maxpoints)
      )
      , 1:n, c())

    d = data.frame(signal = s,
                   t = t,
                   receptor = receptor)

    ggplot2::ggplot(d) +
      ggplot2::geom_line(ggplot2::aes(t, signal), alpha = alpha) +
      ggplot2::facet_wrap(~ receptor)



  }



frequency_plot_frames <- function(signal,
                                  frame,
                           f_start = 0,
                           f_end = Inf,
                           maxpoints = 1e4,
                           alpha = 0.7,
                           scale_x_log = F,
                           dB = F)
{
  n = length(frame$fft)
  fs = signal$fs
  f_end = min(fs / 2, f_end)
  nsamples = length(frame$fft[[1]])
  df = fs / nsamples
  i_start = f_start / df + 1
  i_end = f_end / df + 1

  d = data.frame(
    fft = Reduce(function(v, i)
      c(
        v, downsample_signal(frame$fft[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c()),
    freq = Reduce(function(v, i)
      c(
        v,
        downsample_frequency(df, frame$framed_raw_signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c()),
    file = Reduce(function(v, i)
      c(
        v,
        downsample_label(signal$label[i],
                         frame$framed_raw_signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())


  )
  d$dB = 10 * log10(abs(d$fft) / max(abs(d$fft)))
  if (scale_x_log)
  {
    ggplot2::ggplot(d) +
      ggplot2::geom_line(ggplot2::aes(freq, dB, color = file), alpha = alpha) +
      ggplot2::scale_x_log10() +
      ggplot2::facet_wrap(~ file, scales = "free")
  }
  else if (dB)
    ggplot2::ggplot(d) +
    ggplot2::geom_line(ggplot2::aes(freq, dB, color = file), alpha = alpha) +
    ggplot2::facet_wrap(~ file, scales = "free")
  else
    ggplot2::ggplot(d) +
    ggplot2::geom_line(ggplot2::aes(freq, abs(fft), color = file), alpha = alpha) +
    ggplot2::facet_wrap(~ file, scales = "free")

}


