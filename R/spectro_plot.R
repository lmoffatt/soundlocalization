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

    fs = recordings$fs[1]
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
    rec = apply_frame(rec,framed_by = "spectro_data,",i_start = i_s,nsamples = ns)
    wl_points = 512
    if (!is.null(window_length_meters))
      wl_points = 2 ^ ceiling(log(window_length_meters / velocity_of_sound *
                                    rec$fs[1]) / log(2))
    nsamples = rec$nsamples[1]
    fs = rec$fs[1]
    n_frames = floor(nsamples / wl_points * 2) - 1
    indexes = lapply(1:n_frames, function (i)
      list(
        framed_by="spectro_data",
        i_start = 1 + (i - 1) * wl_points / 2,
        nsamples =  wl_points
      ))

    rec %>%
      apply_filter(min_freq = min_freq, max_freq = max_freq) -> rec

    l_fft = lapply(indexes, function(ind)
      rec %>% apply_frame(ind$framed_by,ind$i_start ,ind$nsamples) %>% apply_hanning_window() %>% calculate_ffts())
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
      ggplot2::ggplot() +
      ggplot2::geom_raster(ggplot2::aes(x = t, y = f, fill = dB), interpolate = T) +
      ggplot2::facet_wrap(~ receptor) +
      ggplot2::scale_fill_distiller(palette = "Spectral", limits =
                                      c(min_dB, 0))


  }

