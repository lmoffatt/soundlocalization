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


select_receptors <- function(rec, receptorList)
{
  d <- lapply(names(rec), function(n)
    rec[[n]][receptorList])
  names(d) <- names(rec)
  return (d)
}




downsample_time <- function(fs, signal, maxpoints, i_start, i_end)
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
  function(label, signal, maxpoints, i_start, i_end)
  {
    ns = length(signal)
    i_s = max(1, i_start)
    i_e = min(i_end, ns)
    nsamples = i_e - i_s + 1
    factor = ceiling(nsamples / maxpoints)
    i = seq(from = i_s, to = i_e, by = factor)
    return (rep(label, length(i)))
  }

downsample_signal <- function(signal, maxpoints, i_start, i_end)
{
  ns = length(signal)
  i_s = max(1, i_start)
  i_e = min(i_end, ns)
  nsamples = i_e - i_s + 1
  factor = ceiling(nsamples / maxpoints)
  i = seq(from = i_s, to = i_e, by = factor)
  return (signal[i])
}

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
  if (scale_x_log)
  {
    ggplot2::ggplot(d) +
      ggplot2::geom_line(aes(freq, dB, color = file), alpha = alpha) +
      scale_x_log10() +
      facet_wrap(~ file, scales = "free")
  }
  else if (dB)
    ggplot2::ggplot(d) +
    ggplot2::geom_line(aes(freq, dB, color = file), alpha = alpha) +
    facet_wrap(~ file, scales = "free")
  else
    ggplot2::ggplot(d) +
    ggplot2::geom_line(aes(freq, abs(fft), color = file), alpha = alpha) +
    facet_wrap(~ file, scales = "free")

}



sound_plot <-
  function(signal,
           t_start = 0,
           t_end = Inf,
           maxpoints = 1e4,
           alpha = 0.7,
           type = "line")
  {
    n = length(signal$signal)
    fs = signal$fs[1]
    i_start = t_start * fs
    i_end = t_end * fs

    s = Reduce(function(v, i)
      c(
        v,
        downsample_signal(signal$signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())
    t = Reduce(function(v, i)
      c(
        v,
        downsample_time(fs, signal$signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())
    receptor = Reduce(function(v, i)
      c(
        v,
        downsample_label(signal$filename[i],
                         signal$signal[[i]], maxpoints, i_start, i_end)
      )
      , 1:n, c())

    d = data.frame(signal = s,
                   t = t,
                   receptor = receptor)

    ggplot2::ggplot(d) +
      ggplot2::geom_line(aes(t, signal), alpha = alpha) +
      ggplot2::facet_wrap(~ receptor)



  }
to_Wave <- function(recordings, i, start, end)
{
  fs = recordings$fs[1]
  en = min(length(recordings$signal[[i]]), end * fs)
  tuneR::Wave(
    left = recordings$signal[[i]][(start * fs + 1):en],
    samp.rate = recordings$fs[i],
    bit = recordings$bit[i],
    pcm = recordings$pcm[i]
  )
}


spectro_plot <-
  function(recordings,
           i,
           start,
           end,
           min_freq,
           max_freq,
           listen = F,
           listen_freq = 1,
           complex = F)
  {
    w = to_Wave(recordings = recordings, i, start, end)
    f = recordings$fs[1]
    w = seewave::ffilter(w, f, from = min_freq, to = max_freq)

    fs = recordings$fs[i]
    if (listen)
    {
      tuneR::setWavPlayer('/usr/bin/aplay')
      seewave::listen(wave = w, f = listen_freq * fs)
    }
    seewave::spectro(
      w,
      f = fs,
      ovlp = 50,
      complex = complex,
      collevels = seq(-50, 0, 2),
      flim = c(0, max_freq * 1e-3 * 1.2)
    )
  }
spectro_data <-
  function(rec,
           start,
           end,
           min_freq,
           max_freq,
           window_length_meters = NULL,
           ovlp,
           min_dB = -50,
           velocity_of_sound = 334)
  {
    fs = rec$fs[1]
    ns = floor((end - start) * fs / 2) * 2
    i_s = floor(start * fs) + 1
    rec = apply_frame(rec, list(i_start = i_s,
                                i_end = i_s + ns + 1))
    wl_points = 512
    if (!is.null(window_length_meters))
      wl_points = 2 ^ ceiling(log(window_length_meters / velocity_of_sound *
                                    rec$fs[1]) / log(2))
    nsamples = rec$nsamples[1]
    fs = rec$fs[1]
    n_frames = floor(nsamples / wl_points * 2) - 1
    indexes = lapply(1:n_frames, function (i)
      list(
        i_start = 1 + (i - 1) * wl_points / 2,
        i_end = (i + 1) * wl_points / 2
      ))

    rec %>%
      apply_filter(min_freq = min_freq, max_freq = max_freq) -> rec

    l_fft = lapply(indexes, function(ind)
      rec %>% apply_frame(ind) %>% apply_hanning_window() %>% calculate_ffts())
    d = Reduce(function(p, i_receptor)
    {
      d_receptor = Reduce(function(pr, i_frame)
      {
        Amplitude = Mod(l_fft[[i_frame]]$fft
                        [[i_receptor]][1:(wl_points / 2)])

        rbind(
          pr,
          data.frame(
            t = indexes[[i_frame]]$i_start / fs,
            f = 0:(wl_points / 2 - 1) * fs / wl_points,
            Amplitude = Amplitude,
            receptor = rec$labels[i_receptor]
          )
        )
      },
      1:length(indexes),
      data.frame())
      d_receptor$dB = 20 * log10(d_receptor$Amplitude / max(d_receptor$Amplitude))

      rbind(p, d_receptor)
    },
    1:length(rec$labels),
    data.frame())
    return(d)

  }


spectro_plot_luc <-
  function(rec,
           start,
           end,
           min_freq,
           max_freq,
           window_length_meters = NULL,
           ovlp,
           min_dB = -50,
           velocity_of_sound = 334)
  {
    d = spectro_data(
      rec,
      start,
      end,
      min_freq,
      max_freq,
      window_length_meters,
      ovlp,
      min_dB,
      velocity_of_sound
    )
    d %>% dplyr::filter(f < max_freq) %>%
      #  dplyr::filter(dB > min_dB) %>%
      ggplot() +
      geom_raster(aes(x = t, y = f, fill = dB), interpolate = T) +
      facet_wrap(~ receptor) +
      scale_fill_distiller(palette = "Spectral", limits =
                             c(min_dB, 0))


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



#' Title calculates all the cross correlation pairs
#'
#' @param ffted_signal the fft should be already calculated
#'
#' @return inserts the cross correlation function
#' @export
#'
#' @examples
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


apply_filter <- function(rec, min_freq, max_freq)
{
  f = rec$fs[1]
  rec$min_freq = min_freq
  rec$max_freq = max_freq
  rec$signal_for_fft =
    lapply(rec$signal,
           function(x)
             seewave::ffilter(x, f, from = min_freq, to = max_freq)[, 1])
  return(rec)
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
    if (!(is.null(t_start) & is.null(t_end)))
    {
      if (is.null(t_start))
        t_start = 0
      if (is.null(t_end))
        t_end = rec$duration
      i_start = floor(t_start * rec$fs[1]) + 1
      nsamples = floor((t_end - t_start) * rec$fs[1] / 2) * 2
      i_end = i_start + nsamples - 1
      rec <- rec %>%
        apply_frame(list(i_start = i_start, i_end = i_end))
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



data_for_plot_gcc_phase <-
  function(gcc,
           lag_window_in_meters,
           max_points,
           decimating_operation = "max",
           velocity_of_sound = 334)
  {
    nsamples = length(gcc$gcc_phase[[1]][[1]])
    lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
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
                                lags = lags_f) %>% group_by(lags) %>%
                  summarise(lag_value = match.fun(decimating_operation)(lag_value))

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
        remove_if_z_is_greater_than = remove_if_z_is_greater_than,
        velocity_of_sound = velocity_of_sound,
        freq_filter = freq_filter
      )
    }
    ggplot(rec$gcc_phase_data_for_plot) + geom_line(aes(
      x = lags,
      y = std_lag_values,
      color = paste0(i, j)
    )) + facet_grid(i ~ j, scales = "free", space = "free")
  }




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

    lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
    gcc$lag_max = lag_max

    nsamples = length(gcc$gcc_phase[[1]][[1]])
    stopifnot("the sampled frame is wider than the lag window" = lag_max * 4 < nsamples)

    lags_2 = ((-lag_max * 2 - 1):(lag_max * 2)) / gcc$fs[1]

    factor = ceiling(2 * lag_max / max_points)

    lag_max_f = ceiling(lag_max / factor)

    lags_2f = floor(lags_2 * gcc$fs[1] / factor) * factor / gcc$fs[1]


    i_lags_2 = c((-lag_max * 2):0 + nsamples, 1:(lag_max * 2 + 1))


    n_2_lags = length(unique(lags_2f))
    n_receptors = length(gcc$labels)

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
                          ) %>% group_by(lags) %>% summarise(
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


    ggplot(gcc$gcc_data_for_tri_plot) +
      geom_raster(aes(x = lag_i_j, y = lag_i_k, fill = lag_values)) +
      facet_grid(i_k ~ i_j) + scale_fill_distiller(palette = "Spectral")
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

    lag_max = ceiling(lag_window_in_meters / velocity_of_sound * gcc$fs[1])
    gcc$lag_max = lag_max

    nsamples = length(gcc$gcc_phase[[1]][[1]])
    stopifnot("the sampled frame is wider than the lag window" = lag_max * 2 < nsamples)

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


    gcc$gcc_data_for_source_plot = Reduce(function(p, i)
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
                          ) %>% group_by(lags) %>% summarise(
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
                            gglasso(
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



gcc_phase_source_plot <-
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
    if (!"gcc_data_for_source_plot" %in% names(gcc))
      gcc = gcc_phase_data_for_source_reconstruction_plot(
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



    print(
      ggplot(r_source$gcc_data_for_source_plot) +
        geom_raster(aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        geom_point(
          data = filter(
            r_source$gcc_data_for_source_plot,
            !is.na(source_nvalues_AB)
          ),
          aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_AB,
            color = source_nvalues_AB
          ),
          alpha = 0.2
        ) +
        facet_grid(i_k ~ i_j) + scale_fill_distiller(palette = "Spectral") +
        scale_color_distiller(palette = "Spectral")
    )
    print(
      ggplot(r_source$gcc_data_for_source_plot) +
        geom_raster(aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        geom_point(
          data = filter(
            r_source$gcc_data_for_source_plot,
            !is.na(source_nvalues_AC)
          ),
          aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_AC,
            color = source_nvalues_AC
          ),
          alpha = 0.2
        ) +
        facet_grid(i_k ~ i_j) + scale_fill_distiller(palette = "Spectral") +
        scale_color_distiller(palette = "Spectral")

    )
    print(
      ggplot(r_source$gcc_data_for_source_plot) +
        geom_raster(aes(
          x = lag_i_j, y = lag_i_k, fill = lag_values
        )) +
        geom_point(
          data = filter(
            r_source$gcc_data_for_source_plot,
            !is.na(source_nvalues_BC)
          ),
          aes(
            x = lag_i_j,
            y = lag_i_k,
            size = source_nvalues_BC,
            color = source_nvalues_BC
          ),
          alpha = 0.2
        ) +
        facet_grid(i_k ~ i_j) + scale_fill_distiller(palette = "Spectral") +
        scale_color_distiller(palette = "Spectral")

    )

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
                                    d %>% filter(lag_value > threshold)

                                  return(list(i = i,
                                              j = j,
                                              d = d))
                                }))

  return(gcc)
}


obtain_shared_peaks <- function(gcc,
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
  if (!"gcc_peaks_data" %in% names(gcc))
    gcc <- gphase_filter_peaks_data(
      gcc,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      lag_window_in_meters = lag_window_in_meters,
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      keep_the_best_n = keep_the_best_n,
      velocity_of_sound = velocity_of_sound,
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
    d <- full_join(d_matrix$d, d, by = character())
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
    dd[[paste0("lag_", d_vector$i, "_", d_vector$j)]] =
      dd[[paste0("lag_", 1, "_", d_vector$j)]] - dd[[paste0("lag_", 1, "_", d_vector$i)]]
    d <- inner_join(dd, d)
    list(
      i = c(d_matrix$i, d_vector$i),
      j = c(d_matrix$j, d_vector$j),
      d = d
    )
  }
  n_receptors = length(gcc$labels)

  gcc$shared_peaks =
    Reduce (
      function (d_matrix, j)
        Reduce(
          function (d_matrixi, i)
            if (gcc$gcc_peaks_data[[i]][[j - i]]$j %in% d_matrixi$j)
              contract(d_matrixi, gcc$gcc_peaks_data[[i]][[j - i]])
          else
            expand(d_matrixi, gcc$gcc_peaks_data[[i]][[j - i]]),
          1:(j - 1),
          d_matrix
        ),
      2:n_receptors,
      list(d = data.frame())
    )

  return (gcc)
}


shared_peaks_data_for_plot <- function(gcc,
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
  if (!"shared_peaks" %in% names(gcc))
    gcc <- obtain_shared_peaks(
      gcc,
      t_start = t_start,
      t_end = t_end,
      min_freq = min_freq,
      max_freq = max_freq,
      lag_window_in_meters = lag_window_in_meters,
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      keep_the_best_n = keep_the_best_n,
      velocity_of_sound = velocity_of_sound,
      freq_filter = freq_filter
    )


  max_lag = gcc$lag_max / gcc$fs[1]
  min_lag = -max_lag

  n_receptors = length(gcc$labels)

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
      return (list(
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

            lag_shared_i_j = gcc$shared_peaks$d[[paste0("lag_", i, "_", j)]]
            lag_shared_i_k = gcc$shared_peaks$d[[paste0("lag_", i, "_", k)]]

            lag_value_shared_i_j = gcc$shared_peaks$d[[paste0("values_", i, "_", j)]]
            lag_value_shared_i_k = gcc$shared_peaks$d[[paste0("values_", i, "_", k)]]
            lag_value_shared_j_k = gcc$shared_peaks$d[[paste0("values_", j, "_", k)]]

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
                              keep_if_z_is_greater_than = 5,
                              keep_the_best_n = 5,
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
      keep_if_z_is_greater_than = keep_if_z_is_greater_than,
      keep_the_best_n = keep_the_best_n,
      velocity_of_sound = velocity_of_sound,
      freq_filter = freq_filter
    )
  print(
    ggplot(gcc$shared_peaks_for_plot$d_lines) +
      geom_path(aes(x_lags_all, y_lags_all), alpha = 0.1) +
      facet_grid(i_k ~ i_j, scales = "free") + scale_color_distiller(palette = "Spectral") +
      geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )
  print(
    ggplot(gcc$shared_peaks_for_plot$d_lines) +
      geom_path(aes(x_lags_set, y_lags_set), alpha = 0.1) +
      facet_grid(i_k ~ i_j, scales = "free") + scale_color_distiller(palette = "Spectral") +
      geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )
  print(
    ggplot(gcc$shared_peaks_for_plot$d_lines) +
      geom_path(aes(x_lags_shared, y_lags_shared), alpha = 0.1) +
      facet_grid(i_k ~ i_j, scales = "free") + scale_color_distiller(palette = "Spectral") +
      geom_point(
        data = gcc$shared_peaks_for_plot$d_shared,
        aes(x = lag_shared_i_j, y = lag_shared_i_k)
      )
  )

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






listen <- function(recordings, i, start, end)
{
  fs = recordings$fs[i]

  seewave::listen(to_Wave(recordings = recordings, i, start, end), f = fs)

}
