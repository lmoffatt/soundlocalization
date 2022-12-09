apply_hanning_window_per_frame <- function(rec, frame)
{
  #  hn = e1071::hanning.window(nsamples)
  stopifnot(length(unique(frame$nsamples)) == 1)

  hn = seewave::ftwindow(frame$nsamples[1], wn = "hanning")
  frame$filtered_by = "apply_hanning_window"
  signal_for_fft =
    lapply(1:length(rec$raw_signal), function (i_receptor)
    {
      is = frame$i_start[i_receptor]
      ie = is + frame$nsamples[min(length(frame$nsamples), i_receptor)] - 1
      s = rec$raw_signal[[i_receptor]][is:ie] * hn
      return (s)
    })
  frame$signal_for_fft = signal_for_fft
  return (frame)
}



#' applies a hanning window
#'
#' @param equalized_signal  the signal should be already in matrix form
#'
#' @return it modifies the signal in place
#' @export
#'
#' @examples
apply_hanning_window <- function(rec)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames, function(frame)
      apply_hanning_window_per_frame(rec, frame))
    return (rec)
  }
  else
    return (apply_hanning_window_per_frame(rec, rec))

}




downsample_frequency <-
  function(df, signal, maxpoints, i_start, i_end)
  {
    ns = length(signal)
    i_s = max(1, i_start)
    i_e = min(i_end, ns)
    nsamples = i_e - i_s + 1
    factor = ceiling(nsamples / maxpoints)
    i = seq(from = i_s, to = i_e, by = factor)
    return ((i - 1) * df)
  }


downsample_time <-
  function(fs,
           signal,
           maxpoints,
           i_start = NULL,
           i_end = NULL)
  {
    ns = length(signal)
    i_s = max(1, i_start)
    i_e = min(i_end, ns)
    nsamples = i_e - i_s + 1
    factor = ceiling(nsamples / maxpoints)
    i = seq(from = i_s, to = i_e, by = factor)
    return (i / fs)
  }

downsample_label <-
  function(label,
           signal,
           maxpoints,
           i_start = NULL,
           i_end = NULL)
  {
    ns = length(signal)
    i_s = max(1, i_start)
    i_e = min(i_end, ns)
    nsamples = i_e - i_s + 1
    factor = ceiling(nsamples / maxpoints)
    i = seq(from = i_s, to = i_e, by = factor)
    return (rep(label, length(i)))
  }

downsample_signal <-
  function(signal,
           maxpoints,
           i_start = NULL,
           i_end = NULL)
  {
    ns = length(signal)
    i_s = max(1, i_start)
    i_e = min(i_end, ns)
    nsamples = i_e - i_s + 1
    factor = ceiling(nsamples / maxpoints)
    i = seq(from = i_s, to = i_e, by = factor)
    return (signal[i])
  }


standarize <- function(x, x_mean = NULL, x_sd = NULL)
{
  if (is.null(x_mean))
    x_mean = mean(x)
  if (is.null(x_sd))
    x_sd = sd(x)
  return((x - x_mean) / x_sd)
}


standarization_blind_to_ouliers <-
  function(x, remove_if_z_is_greater_than)
  {
    m_x = mean(x)
    s_x = sd(x)
    z_x = standarize(x, m_x, s_x)
    m_x = mean(x[abs(z_x) < remove_if_z_is_greater_than])
    s_corr = sd(x[abs(z_x) < remove_if_z_is_greater_than])
    if (s_corr < s_x)
    {
      s_x = s_corr
      z_x = standarize(x, m_x, s_x)
      m_x = mean(x[abs(z_x) < remove_if_z_is_greater_than])
      s_corr = sd(x[abs(z_x) < remove_if_z_is_greater_than])
    }
    return(list(mean = m_x
                , sd = s_corr, z = z_x))
  }


apply_filter_per_frame <-
  function(rec,
           frame,
           min_freq = NULL,
           max_freq = NULL)
  {
    f = rec$fs[1]
    if (is.null(min_freq))
      min_freq = frame$min_freq
    else
      frame$min_freq = min_freq
    if (is.null(max_freq))
      max_freq = frame$max_freq
    else
      frame$max_freq = max_freq

    if (is.null(min_freq) && is.null(max_freq))
      return (apply_hanning_window_per_frame(rec, frame))
    else
    {
      frame$filtered_by = "apply_filter"
      if (is.null(min_freq))
        min_freq = 0
      if (is.null(max_freq))
        max_freq = rec$fs / 2
      frame$signal_for_fft =
        lapply(1:length(rec$raw_signal),
               function (i)
               {
                 is = frame$i_start[i]
                 ie = is + frame$nsamples[i] - 1
                 s = rec$raw_signal[[i]][is:ie]
                 if (is.null(min_freq))
                   return(seewave::ffilter(
                     s,
                     f = rec$fs[1],
                     to = max_freq,
                     wl = min(1024, length(s))
                   )[, 1])
                 else  if (is.null(max_freq))
                   return (seewave::ffilter(
                     s,
                     f = rec$fs,
                     from = min_freq,
                     wl = min(1024, length(s))
                   )[, 1])
                 else
                   return (seewave::ffilter(
                     s,
                     f = rec$fs,
                     from = min_freq,
                     to = max_freq,
                     wl = min(1024, length(s))
                   )[, 1])

               })
      return(frame)
    }
  }


#' Title calculates all the cross correlation pairs
#'
#' @param ffted_signal the fft should be already calculated
#'
#' @return inserts the cross correlation function
#' @export
#'
#' @examples


apply_filter <- function(rec,
                         min_freq = NULL,
                         max_freq = NULL)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames, function(frame)
      apply_filter_per_frame(rec, frame, min_freq = min_freq, max_freq = max_freq))
    return (rec)
  }
  else
    return (apply_filter_per_frame(rec, rec, min_freq = min_freq, max_freq = max_freq))


}




filter_by_lag <- function(fft1, fft2, lag, lag_error, fs)
{
  stopifnot(length(fft1) == length(fft2))
  nsamples = length(fft1)

  w = c(0:(nsamples / 2 - 1), (-nsamples / 2):-1) / nsamples * fs * 2 *
    pi
  phi = exp(-lag * w * 1i)
  angle_err = abs(Arg(exp((-lag - lag_error) * w * 1i) / phi))

  phi_delta = abs(Arg(fft2 / (fft1 * phi)))
  phi_filter = phi_delta < angle_err
  fft1_fil = fft1 * phi_filter
  fft2_fil = fft2 * phi_filter

  return(
    list(
      f = w / 2 / pi,
      lag = lag,
      phi_delta = phi_delta,
      angle_err = angle_err,
      phi_filter = phi_filter,
      fft1 = fft1_fil,
      fft2 = fft2_fil,
      signal1 = (fft(fft1_fil, inverse = T)),
      signal2 = fft(fft2_fil, inverse = T)
    )
  )
}
