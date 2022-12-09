





optimize_time_lags_to_single_source_position <-
  function(lags,
           receptors,
           sound_velocity)
  {
    n_receptors = ncol(receptors)


    mrec = rowMeans(receptors)

    expected_lags <- function(beta) {
      delays = colSums((beta - receptors) ^ 2) ^ 0.5 / sound_velocity


      return((delays[2:n_receptors] - delays[1] - lags) * 1000)
    }

    localization_sqr_sum <- function(beta)
    {
      sum(expected_lags(beta) ^ 2)
    }

    beta_initial = mrec
    beta_min = c(-10000, 0)
    beta_max = c(10000, 10000)
    maxeval = 1000

    localization_sqr_sum(beta_initial)

    # opd = nloptr::direct(
    #   localization_sqr_sum,
    #   lower = beta_min,
    #   upper = beta_max,
    #   nl.info = T,
    #   control = list(maxeval = maxeval)
    # )
    #
    # beta_opt2 = pracma::fminsearch(localization_sqr_sum, opd$par,maxiter = 10)
    beta_opt3 = pracma::lsqnonlin(fun = expected_lags, x0 = c(10, 10))
    return (list(
      x = beta_opt3$x[1],
      y = beta_opt3$x[2],
      ss = beta_opt3$ssq
    ))

  }


optimize_all_time_lags_to_their_source_positions <-
  function(lags,
           receptors,
           sound_velocity)

  {
    opt_lags = vapply(
      X = seq_len(nrow(lags)),
      FUN = function(i) {
        op = optimize_time_lags_to_single_source_position(
          lags = as.numeric(lags[i,]),
          receptors = receptors,
          sound_velocity = sound_velocity
        )
        return (c(op$x, op$y, op$ss))
      },
      numeric(3)
    )
    d = data.frame(t(opt_lags))
    colnames(d) <- c("x", "y", "ss")
    d = dplyr::arrange(d, ss)
    return (d)

  }

get_receptors <- function(rec)
{
  receptors = vapply(
    1:length(rec$labels),
    FUN = function(i)
      c(rec$x[i], rec$y[i]),
    numeric(2)
  )
}

get_receptors_max <- function(rec)
{
  receptors = vapply(
    1:length(rec$labels),
    FUN = function(i)
      c(rec$x[i] + rec$pos_err[i], rec$y[i] + rec$pos_err[i]),
    numeric(2)
  )
}

get_receptors_min <- function(rec)
{
  receptors = vapply(
    1:length(rec$labels),
    FUN = function(i)
      c(rec$x[i] - rec$pos_err[i], rec$y[i] - rec$pos_err[i]),
    numeric(2)
  )
}


get_offset_max <- function(rec, velocity_of_sound = 334)
{
  max_lag = ((max(rec$x) - min(rec$x)) ^ 2 + (max(rec$y) - min(rec$y)) ^ 2) ^
    0.5 / velocity_of_sound
}


obtain_source_positions_from_lags_per_frame <-
  function (rec,
            frame,
            max_number_of_sources_per_frame,
            fraction_of_sum_of_values,
            receptors = NULL,
            offset = NULL,
            t_start = NULL,
            t_end = NULL,
            min_freq = NULL,
            max_freq = NULL,
            lag_window_in_meters = NULL,
            keep_if_z_is_greater_than = 5,
            keep_the_best_n = 50,
            velocity_of_sound = 334,
            signal_filter = T,
            freq_filter = F)
  {
    lag_max = obtain_lag_max(
      rec,
      frame = frame,
      lag_window_in_meters = lag_window_in_meters,
      velocity_of_sound = velocity_of_sound
    )

    if (has_element_and_matches_this(
      frame,
      "position",
      list_to_match = list(
        min_freq = min_freq,
        max_freq = max_freq,
        signal_filter = signal_filter,
        freq_filter = freq_filter,
        keep_if_z_is_greater_than = keep_if_z_is_greater_than,
        keep_the_best_n = keep_the_best_n,
        lag_max = lag_max,
        max_number_of_sources_per_frame = max_number_of_sources_per_frame,
        fraction_of_sum_of_values = fraction_of_sum_of_values,
        receptors = receptors,
        offset = offset
      )
    ))
    return(frame)


    if (is.null(receptors))
      receptors = get_receptors(rec)



    frame <- obtain_shared_peaks_per_frame(
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
      signal_filter = signal_filter,
      freq_filter = freq_filter
    )
    delays = dplyr::select(frame$shared_peaks$delays, dplyr::starts_with('lag_'))
    min_sum_of_values = max(frame$shared_peaks$d$sum_of_values)*fraction_of_sum_of_values
    n_values = sum(frame$shared_peaks$d$sum_of_values > min_sum_of_values)

    n_delays = min(nrow(delays), max_number_of_sources_per_frame, n_values)
    delays = delays[seq_len(n_delays),]

    if (!is.null(offset) & (n_delays > 0))
      delays = t(t(delays) - offset / 1000)


    frame$shared_peaks$positions <-
      optimize_all_time_lags_to_their_source_positions(delays,
                                                       receptors = receptors,
                                                       sound_velocity = velocity_of_sound)
    return(frame)
  }


obtain_source_positions_from_lags <- function(rec,
                                              max_number_of_sources_per_frame,
                                              fraction_of_sum_of_values,
                                              receptors = NULL,
                                              offset = NULL,
                                              t_start = NULL,
                                              t_end = NULL,
                                              min_freq = NULL,
                                              max_freq = NULL,
                                              lag_window_in_meters = NULL,
                                              keep_if_z_is_greater_than = 5,
                                              keep_the_best_n = 50,
                                              velocity_of_sound = 334,
                                              signal_filter = F,
                                              freq_filter = T)
{
  if ("frames" %in% names(rec))
  {
    rec$frames = lapply(rec$frames,
                        function(frame)
                          return(
                            obtain_source_positions_from_lags_per_frame(
                              rec = rec,
                              frame = frame,
                              receptors = receptors,
                              offset = offset,
                              max_number_of_sources_per_frame = max_number_of_sources_per_frame,
                              fraction_of_sum_of_values = fraction_of_sum_of_values,
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
                          ))
    return (rec)
  }
  else
    return (
      obtain_source_positions_from_lags_per_frame(
        rec,
        rec,
        max_number_of_sources_per_frame,
        fraction_of_sum_of_values = fraction_of_sum_of_values,
        receptors = receptors,
        offset = offset,
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






optimize_receptors_using_gcc_lags_old <-
  function(rec,
           max_number_of_sources_per_frame,
           fraction_of_sum_of_values,
           only_offset = T,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           keep_if_z_is_greater_than = 5,
           keep_the_best_n = 50,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    n_receptors = length(rec$labels)

    receptors_init = get_receptors(rec)

    receptors_max = get_receptors_max(rec)
    receptors_min = get_receptors_min(rec)

    offset_max = rep(get_offset_max(rec), n_receptors - 1) * 1000

    beta_to_offset <- function(beta)
    {
      return (beta[1:(n_receptors - 1)])
    }

    beta_to_receptors <- function(beta)
    {
      return (matrix(c(0, 0, 0, beta[n_receptors + 0:(n_receptors * 2 - 4)]), nrow = 2))
    }

    source_position_from_beta <- function(offset, receptors)
    {
      pos_rec = obtain_source_positions_from_lags(
        rec = rec,
        max_number_of_sources_per_frame = max_number_of_sources_per_frame,
        fraction_of_sum_of_values = fraction_of_sum_of_values,
        offset = offset,
        receptors = receptors,

        velocity_of_sound = velocity_of_sound
      )
      return (pos_rec)
    }

    sqr_sum_per_offset <- function(beta)
    {
      offset = beta_to_offset(beta)
      if (!only_offset)
        receptors = beta_to_receptors(beta)
      else
        receptors = get_receptors(rec)

      pos_rec = source_position_from_beta(offset, receptors)

      ss = Reduce(function(p, frame)
        p + sum(frame$shared_peaks$positions$ss),
        pos_rec$frames,
        0)
      print(list(
        offset = offset,
        receptors = receptors,
        ss = ss
      ))
      return(ss)
    }


    offset_init = numeric(n_receptors - 1)
    receptors_init = as.numeric(receptors_init)[4:(2 * n_receptors)]



    beta_init = c(offset_init, receptors_init)


    receptor_min = as.numeric(receptors_min)[4:(2 * n_receptors)]
    receptor_max = as.numeric(receptors_max)[4:(2 * n_receptors)]

    beta_min = c(-offset_max, receptor_min)

    beta_max = c(offset_max, receptor_max)

    if (only_offset)
    {
      beta_init = offset_init
      beta_min = -offset_max
      beta_max = offset_max
    }

    #  opt_res = pracma::nelder_mead(sqr_sum_per_offset, beta_init)

    #
    #opt_res=nloptr::neldermead(



    opt_res = nloptr::sbplx(
      x0 = beta_init,
      fn = sqr_sum_per_offset,
      lower = beta_min,
      upper = beta_max
      #  nl.info = T,
      #  control = list(maxeval = 1000)
    )
    offsets_opt = opt_res$xmin

    offset = beta_to_offset(offsets_opt)
    receptors = beta_to_receptors(offsets_opt)
    pos_opt = source_position_from_beta(offset, receptors)
    return (list(opt_res = opt_res, pos_opt = pos_opt))
  }



get_lags_dataframe <- function(rec)
{
  return (Reduce(function(p, i) {
    d = rec$frames[[i]]$shared_peaks$d
    if (nrow(d) == 0)
      return (p)
    else
    {
      d$t_start = sp_rec_1$frames[[i]]$i_start[3] / sp_rec_1$fs[1]
      d$i = i
      return(rbind(p, d))
    }
  }, 1:length(sp_rec_1$frames), data.frame()) -> sp_d_1)

}





#' gcc-phase peaks that are likely to be originated by different sound sources
#' (as opposed to noise correlation artifacts)
#'
#'
#'
#'
#' @param sp_rec
#'
#' @return
#' @export
#'
#' @examples
filter_good_gcc_peaks_for_calibration <- function(sp_rec)
{

}








optimize_receptors_using_gcc_lags <-
  function(rec,
           max_number_of_sources_per_frame,
           fraction_of_sum_of_values,
           only_offset = T,
           t_start = NULL,
           t_end = NULL,
           min_freq = NULL,
           max_freq = NULL,
           lag_window_in_meters = NULL,
           keep_if_z_is_greater_than = 5,
           keep_the_best_n = 50,
           velocity_of_sound = 334,
           freq_filter = F)
  {
    n_receptors = length(rec$labels)

    receptors_init = get_receptors(rec)

    receptors_max = get_receptors_max(rec)
    receptors_min = get_receptors_min(rec)

    offset_max = rep(get_offset_max(rec), n_receptors - 1) * 1000

    beta_to_offset <- function(beta)
    {
      return (beta[1:(n_receptors - 1)])
    }

    beta_to_receptors <- function(beta)
    {
      return (matrix(c(0, 0, 0, beta[n_receptors + 0:(n_receptors * 2 - 4)]), nrow = 2))
    }

    source_position_from_beta <- function(offset, receptors)
    {
      pos_rec = obtain_source_positions_from_lags(
        rec = rec,
        max_number_of_sources_per_frame = max_number_of_sources_per_frame,
        fraction_of_sum_of_values = fraction_of_sum_of_values,
        offset = offset,
        receptors = receptors,

        velocity_of_sound = velocity_of_sound
      )
      return (pos_rec)
    }

    sqr_sum_per_offset <- function(beta)
    {
      offset = beta_to_offset(beta)
      if (!only_offset)
        receptors = beta_to_receptors(beta)
      else
        receptors = get_receptors(rec)

      pos_rec = source_position_from_beta(offset, receptors)

      ss = Reduce(function(p, frame)
        p + sum(frame$shared_peaks$positions$ss),
        pos_rec$frames,
        0)
      print(list(
        offset = offset,
        receptors = receptors,
        ss = ss
      ))
      return(ss)
    }


    offset_init = numeric(n_receptors - 1)
    receptors_init = as.numeric(receptors_init)[4:(2 * n_receptors)]



    beta_init = c(offset_init, receptors_init)


    receptor_min = as.numeric(receptors_min)[4:(2 * n_receptors)]
    receptor_max = as.numeric(receptors_max)[4:(2 * n_receptors)]

    beta_min = c(-offset_max, receptor_min)

    beta_max = c(offset_max, receptor_max)

    if (only_offset)
    {
      beta_init = offset_init
      beta_min = -offset_max
      beta_max = offset_max
    }

    #  opt_res = pracma::nelder_mead(sqr_sum_per_offset, beta_init)

    #
    #opt_res=nloptr::neldermead(



    opt_res = nloptr::sbplx(
      x0 = beta_init,
      fn = sqr_sum_per_offset,
      lower = beta_min,
      upper = beta_max
      #  nl.info = T,
      #  control = list(maxeval = 1000)
    )
    offsets_opt = opt_res$xmin

    offset = beta_to_offset(offsets_opt)
    receptors = beta_to_receptors(offsets_opt)
    pos_opt = source_position_from_beta(offset, receptors)
    return (list(opt_res = opt_res, pos_opt = pos_opt))
  }









optimize_localizations_by_gcc_peaks <-
  function(rec,
           n_max_sources = 20,
           vsound = 343,
           beta_optimal = NULL,
           maxeval = 10000)
  {
    n_receptors = length(rec$labels)
    n_frames = length(f_rec_1$frames)

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
