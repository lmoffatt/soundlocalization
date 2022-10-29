distance <- function(x, y)
{
  return ((y$x - x$x) ^ 2 + (y$y - x$y) ^ 2) ^ 0.5
}


set_source_list <- function(names,
                            x_pos_min,
                            x_pos_max,
                            y_pos_min,
                            y_pos_max,
                            min_freq,
                            max_freq)
{
  return (
    data.frame(
      names = names,
      x_min = x_pos_min,
      x_max = x_pos_max,

      y_min = y_pos_min,
      y_max = y_pos_max,
      min_freq = min_freq,
      max_freq = max_freq
    )
  )
}



set_frame_list <- function(t_start, t_end, sourceList)
{
  d <- data.frame(t_start = t_start, t_end = t_end)
  d$sources <- list(sourceList)
  d

}

get_receptorList <- function(labels,
                             geo_locations,
                             origin,
                             x_axis,
                             pos_error,
                             Gain_error_factor,
                             lag_error,
                             all_frames)
{
  n_frames = length(all_frames$signal)
  p = geo_to_local(origin = origin,
                   x_axis = x_axis,
                   points = geo_locations)
  data.frame(
    names = labels,
    x_min = p[, 1] - pos_error,
    x_max = p[, 1] + pos_error,
    y_min = p[, 2] - pos_error,
    y_max = p[, 2] + pos_error,
    lag_min = -lag_error,
    lag_max = +lag_error
  ) -> d
  d$Gain_min = rep(list(rep(1 / Gain_error_factor, n_frames)), nrow(p))

  d$Gain_max = rep(list(rep(Gain_error_factor, n_frames)), nrow(p))

  d

}


frame_to_signal_fft <- function(rec, all_frames)
{
  fs = rec$fs[1]
  all_frames$aligned_signal <-

    lapply(1:length(all_frames$t_start)
           , function (i)
           {
             i_start = all_frames$t_start[i] * fs + 1
             nsamples = min((all_frames$t_end[i] - all_frames$t_start[i]) * fs,
                            length(rec$signal[[i]]) - i_start + 1)
             nsamples = 2 * floor(nsamples / 2)

             i_end = i_start + nsamples - 1
             hn = e1071::hanning.window(nsamples)
             lapply(rec$framed_raw_signal, function(s) s[i_start:i_end] * hn)
           })
  all_frames$w <- lapply(1:length(all_frames$t_start), function (i)
  {
    i_start = all_frames$t_start[i] * fs + 1
    nsamples = min((all_frames$t_end[i] - all_frames$t_start[i]) * fs,
                   nrow(rec$signal[[i]]) - i_start + 1)
    nsamples = 2 * floor(nsamples / 2)
    c(0:(nsamples / 2 - 1), (-nsamples / 2):-1) * fs / nsamples * 2 * pi
  })
  all_frames$fft <- lapply(1:length(all_frames$t_start), function (i)
       lapply(all_frames$signal[[i]],fft))
  all_frames
}

select_juicy_freq <- function(all_frames, fs, min_freq, fraction)
{
  sel = all_frames[c("t_start", "t_end", "sources")]
  s = lapply(1:length(all_frames$signal), function (i)
  {
    fft = rowSums(Mod(all_frames[i, ]$fft[[1]]))
    nsamples = length(fft)
    i_min = min_freq / fs * nsamples

    sort.int(fft[i_min:(nsamples / 2)],
             decreasing = T,
             index.return = T)$ix + i_min - 1
  })
  sel$fft = lapply(1:length(all_frames$signal),
                   function (i)
                     all_frames$fft[[i]][s[[i]][1:(length(s[[i]]) *
                                                     fraction)], ])
  sel$w = lapply(1:length(all_frames$signal),
                 function (i)
                   all_frames$w[[i]][s[[i]][1:(length(s[[i]]) * fraction)]])
  sel

}



optimize_localizations <-
  function(rec,
           all_frames,
           all_receptors,
           all_sources,
           lambda,
           vsound = 343,
           beta_optimal = NULL,
           algorithm = "sbplx",
           maxeval = 10000)
  {
    n_receptors = nrow(all_receptors)
    n_frames = length(all_frames$signal)
    n_sources = nrow(all_sources)


    n_par_receptor = 3

    vector_to_receptor <- function(beta, i, offset, n_frames)
    {
      lag = (i - 1) * (n_par_receptor + n_frames) + offset
      data.frame(x = beta[1 + lag],
                 y = beta[2 + lag],
                 lag = beta[3 + n_frames + lag]) -> d
      d$Gain = list(vapply(1:n_frames, function(i_frame)
        10 ^ beta[2 + i_frame + lag], 0))
      d

    }

    receptor_to_vector_min <- function(receptor, n_frames)
    {
      c(
        receptor$x_min,
        receptor$y_min,
        vapply(1:n_frames, function(i_frame)
          log10(receptor$Gain_min[[1]][i_frame]), 0),
        receptor$lag_min
      )
    }
    receptor_to_vector_max <- function(receptor, n_frames)
    {
      c(
        receptor$x_max,
        receptor$y_max,
        vapply(1:n_frames, function(i_frame)
          log10(receptor$Gain_max[[1]][i_frame]), 0),
        receptor$lag_max
      )
    }
    vector_to_source <- function(all_sources, beta, i, offset)
    {
      lag = (i - 1) * 2 + offset
      name = all_sources$names[i]
      data.frame(names = name,
                 x = beta[1 + lag],
                 y = beta[2 + lag])

    }



    source_to_vector_min <- function(source)
    {
      c(source$x_min, source$y_min)
    }
    source_to_vector_max <- function(source)
    {
      c(source$x_max, source$y_max)
    }



    receptors_sources_to_beta_min <-
      function(all_receptors, all_sources, n_frames)
      {
        return(c(
          Reduce(
            function(v, i)
              c(v, receptor_to_vector_min(all_receptors[i, ], n_frames)),
            1:nrow(all_receptors),
            c()
          ),
          Reduce(
            function(v, i)
              c(v, source_to_vector_min(all_sources[i, ])),
            1:nrow(all_sources),
            c()
          )
        ))
      }


    receptors_sources_to_beta_max <-
      function(all_receptors, all_sources, n_frames)
      {
        return(c(
          Reduce(
            function(v, i)
              c(v, receptor_to_vector_max(all_receptors[i, ], n_frames)),
            1:nrow(all_receptors),
            c()
          ),
          Reduce(
            function(v, i)
              c(v, source_to_vector_max(all_sources[i, ])),
            1:nrow(all_sources),
            c()
          )
        ))
      }




    beta_to_receptors <- function(beta, n_receptors, n_frames)
    {
      Reduce(function(d, i)
        rbind(
          d,
          vector_to_receptor(
            beta = beta,
            i = i,
            offset = 0,
            n_frames = n_frames
          )
        ),
        1:n_receptors,
        data.frame())
    }

    beta_to_sources <-
      function(all_sources,
               beta,
               n_receptors,
               n_sources,
               n_frames)
      {
        Reduce(function(d, i)
          rbind(
            d,
            vector_to_source(
              all_sources = all_sources,
              beta = beta,
              i = i,
              offset = n_receptors * (n_par_receptor + n_frames)
            )
          ),
          1:n_sources,
          data.frame())
      }


    sqr_sum_per_frequency <-
      function(w,
               Yfft,
               dist,
               receptor_lags,
               receptor_Gain,
               lambda)
      {
        A = t(receptor_Gain / dist *
                exp(-1i * w *
                      (dist / vsound - receptor_lags / 1000)))
        XX = Conj(t(A)) %*% A
        XX = XX + lambda * diag(XX)
        Xfft = solve(XX, Conj(t(A)) %*% Yfft)
        Yfit = A %*% Xfft
        err = Yfft - Yfit
        sqr_err = Conj(t(err)) %*% err
        return(Re(sqr_err))
      }

    X_per_frequency <-
      function(w,
               Yfft,
               dist,
               receptor_lags,
               receptor_Gain,
               lambda)
      {
        A = (receptor_Gain / dist *
               exp(-1i * w * (dist / vsound - receptor_lags / 1000)))
        XX =  A %*% Conj(t(A))
        XX2 = XX + lambda * diag(XX)
        Xfft = Yfft %*% Conj(t(A)) %*% solve(XX2)
        # Xfft = solve(XX, (t(A)) %*% Yfft)
        return (t(Xfft))
      }


    X_per_frame <-
      function(ws,
               fft,
               receptors,
               sources,
               lambda,
               i_frame)
      {
        dist <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            max(0.01, distance(
              receptors[i, ], sources[j, ]
            )), 0),
          numeric(nrow(sources)))
        receptor_lags <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            receptors[i, ]$lag, 0),
          numeric(nrow(sources)))
        receptor_Gain <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            receptors[i, ]$Gain[[1]][i_frame], 0),
          numeric(nrow(sources)))
        vapply(1:length(ws),
               function(i)
                 X_per_frequency(
                   w = ws[i],
                   Yfft = fft[i, ],
                   dist = dist,
                   receptor_lags = receptor_lags,
                   receptor_Gain = receptor_Gain,
                   lambda = lambda
                 ),
               complex(nrow(sources))) -> res
        #    print(list(res=res))
        return(t(res))
      }

    sqr_sum_per_frame <-
      function(ws,
               fft,
               receptors,
               sources,
               lambda,
               i_frame)
      {
        dist <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            max(0.01, distance(
              receptors[i, ], sources[j, ]
            )), 0),
          numeric(nrow(sources)))
        receptor_lags <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            receptors[i, ]$lag, 0),
          numeric(nrow(sources)))
        receptor_Gain <- vapply(1:nrow(receptors), function(i)
          vapply(1:nrow(sources), function (j)
            receptors[i, ]$Gain[[1]][i_frame], 0),
          numeric(nrow(sources)))
        Reduce(
          function(s, i)
            s + sqr_sum_per_frequency(
              w = ws[i],
              Yfft = fft[i, ],
              dist = dist,
              receptor_lags = receptor_lags,
              receptor_Gain = receptor_Gain,
              lambda = lambda
            ),
          #    1:length(selected_ws),
          1:(length(ws)),
          0
        ) -> res
        #    print(list(res=res))
        return(res)
      }

    n_par_zero = n_par_receptor + n_frames + 1

    nn_eval = 1

    localization_sqr_sum <- function(beta)
    {
      # there are some degrees of freedom we have to remove

      # we arbitrarialy set the x,y coordinates of the first
      # receptor to zero, its Gain to 1 (the log10 is zero)
      # and its lag to zero
      # there is also a degree of freedom in the orientation of the
      # x,y coordinates that we have to remove
      # therefore we also set to zero the x coordinate
      # of the second receptor

      # in total there are 5 parameters that are zero
      beta = c(rep(0.0, n_par_zero), beta)

      receptors = beta_to_receptors(beta = beta,
                                    n_receptors = n_receptors,
                                    n_frames = n_frames)
      sources = beta_to_sources(
        all_sources = all_sources,
        beta = beta,
        n_sources = n_sources,
        n_receptors = n_receptors,
        n_frames = n_frames
      )


      sqr <- vapply(1:n_frames, function(i_frame)
        sqr_sum_per_frame(
          ws = all_frames[i_frame,]$w[[1]],
          fft = all_frames[i_frame,]$fft[[1]],
          receptors = receptors,
          sources =
            sources[is.element(sources$names, all_frames[i_frame,]$sources[[1]]),],
          lambda = lambda,
          i_frame = i_frame
        ),
        numeric(1))

      nn_eval <<- nn_eval + 1
      if (nn_eval %% 50 == 0)
        print(list(
          n = nn_eval,
          receptors = receptors,
          sources = sources,
          ss = sum(sqr)
        ))
      return (sum(sqr))
    }

    localization_X <- function(beta)
    {
      beta = c(rep(0.0, n_par_zero), beta)

      receptors = beta_to_receptors(beta = beta,
                                    n_receptors = n_receptors,
                                    n_frames = n_frames)
      sources = beta_to_sources(
        all_sources = all_sources,
        beta = beta,
        n_sources = n_sources,
        n_receptors = n_receptors,
        n_frames = n_frames
      )



      all_frames$fftSources <- lapply(1:n_frames, function(i)
        X_per_frame(
          ws = all_frames$w[[i]],
          fft = all_frames$fft[[i]],
          receptors = receptors,
          sources =
            sources[is.element(sources$names, all_frames[i,]$sources[[1]]),],
          lambda = lambda,
          i_frame = i
        ))
      all_frames$Sources <- lapply(1:n_frames, function(i)
        mvfft(all_frames[i,]$fftSources[[1]], inverse = TRUE))

      return (all_frames)
    }


    if (is.null(beta_optimal))
    {
      nn_eval = 1
      beta_min = receptors_sources_to_beta_min(all_receptors, all_sources, n_frames =
                                                 n_frames)
      beta_max = receptors_sources_to_beta_max(all_receptors, all_sources, n_frames =
                                                 n_frames)

      beta_min = beta_min[(n_par_zero + 1):length(beta_min)]
      beta_max = beta_max[(n_par_zero + 1):length(beta_max)]
      beta_initial = beta_min * 0.5 + beta_max * 0.5

      if (algorithm == "neldermead")
      {
        nloptr::neldermead(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "lbfgs")
      {
        nloptr::lbfgs(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "bobyqa")
      {
        nloptr::bobyqa(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "mlsl")
      {
        nloptr::mlsl(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "mlsl")
      {
        nloptr::mlsl(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "stogo")
      {
        nloptr::stogo(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          maxeval = maxeval
        )
      }
      else if (algorithm == "crs2lm")
      {
        nloptr::crs2lm(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "ccsaq")
      {
        nloptr::ccsaq(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "cobyla")
      {
        nloptr::cobyla(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "crs2lm")
      {
        nloptr::crs2lm(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "direct")
      {
        nloptr::direct(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "directL")
      {
        nloptr::directL(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "isres")
      {
        nloptr::isres(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          maxeval = maxeval
        )
      }
      else if (algorithm == "mma")
      {
        nloptr::mma(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "sbplx")
      {
        nloptr::sbplx(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else if (algorithm == "direct")
      {
        nloptr::direct(
          beta_initial,
          localization_sqr_sum,
          lower = beta_min,
          upper = beta_max,
          nl.info = T,
          control = list(maxeval = maxeval)
        )
      }
      else
        print("no algorithm recognized")
    }
    else
    {
      beta = c(rep(0.0, n_par_zero), beta_optimal)

      receptors = beta_to_receptors(beta = beta,
                                    n_receptors = n_receptors,
                                    n_frames = n_frames)
      sources = beta_to_sources(
        all_sources = all_sources,
        beta = beta,
        n_sources = n_sources,
        n_receptors = n_receptors,
        n_frames = n_frames
      )
      opt_frames <- localization_X(beta_optimal)
      list(
        beta_optimal = beta_optimal,
        receptors = receptors,
        sources = sources,
        opt_frames = opt_frames
      )

    }

  }
