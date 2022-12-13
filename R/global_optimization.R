

freq_limits_to_their_indexes <- function(f_rec, frame, freq_limits)
{
  nsamples = frame$nsamples[1]
  ii = 1:(nsamples / 2)
  f =  (ii - 1) / nsamples * f_rec$fs[1]
  lapply(seq_len(nrow(freq_limits)), function(irow)
  {
    index = ii[(f >= freq_limits[irow, 1]) & (f < freq_limits[irow, 2])]
    return(list(
      fmin = freq_limits[irow, 1],
      fmax = freq_limits[irow, 2],
      index = index,
      freq = f[index]
    ))
  })
}



init_receptors <- function(rec)
{
  n_receptors = length(rec$labels)
  x = rec$x
  y = rec$y
  logG = rep(0, n_receptors)
  offset = rep(0, n_receptors)
  rec$receptors_init = list(
    x = x,
    y = y,
    logG = logG,
    offset = offset
  )
  return(rec)
}

init_frames_gain <- function(rec)
{
  n_receptors = length(rec$labels)
  n_frames = length(rec$frames)


  logG_diff = matrix(0,
                     nrow = n_frames,
                     ncol = n_receptors)

  logG = apply(logG_diff, 2, cumsum)


  rec$frames_gain_init = list(logG_diff = logG_diff,
                              logG = logG)

  return(rec)
}




init_sources_position <- function(rec,
                                  n_sources = NULL,
                                  x_min = NULL,
                                  x_max = NULL,
                                  y_min = NULL,
                                  y_max = NULL,
                                  freq_limits = rbind(c(0, 500),
                                                      c(500, 4000),
                                                      c(4000, 10000)))
  {
  if (is.null(x_min))
    x_min = min(rec$x)
  if (is.null(x_max))
    x_max = max(rec$x)
  if (is.null(y_min))
    y_min = min(rec$y)
  if (is.null(y_max))
    y_max = max(rec$y)

  n_receptors = length(rec$labels)
  if (is.null(n_sources))
    n_sources = n_receptors - 1


  n_frames = length(rec$frames)
  n_freq = nrow(freq_limits)

  total_rows = n_frames * n_freq


  i_frame_pos = rep(seq_len(n_frames), each = n_freq)
  i_freq_pos = rep(seq_len(n_freq), times = n_frames)

  x = matrix(
    x_min + (x_max - x_min) * runif(n_sources * total_rows),
    nrow = total_rows,
    ncol = n_sources
  )
  y = matrix(
    y_min + (y_max - y_min) * runif(n_sources * total_rows),
    nrow = total_rows,
    ncol = n_sources
  )



  rec$sources_pos_init = list(
    i_frame_pos = i_frame_pos,
    i_freq_pos = i_freq_pos,
    x = x,
    y = y
  )


  return(rec)
}


init_sources_signal <- function(rec,
                                n_sources = NULL,
                                freq_limits = rbind(c(0, 500),
                                                    c(500, 4000),
                                                    c(4000, 10000)))
{
  n_receptors = length(rec$labels)
  if (is.null(n_sources))
    n_sources = n_receptors - 1


  n_frames = length(rec$frames)
  n_freq = nrow(freq_limits)

  nfreq = matrix(numeric(n_frames * n_freq), nrow = n_freq, ncol = n_frames)

  frame_freq = lapply(
    rec$frames,
    FUN = function(frame)
      freq_limits_to_their_indexes(
        f_rec = rec,
        frame = frame,
        freq_limits = freq_limits
      )
  )

  for (i_frame in seq_len(n_frames))
    for (i_freq in seq_len(n_freq))
      nfreq[i_freq, i_frame] = length(frame_freq[[i_frame]][[i_freq]]$index)

  n_freq_total = sum(nfreq)
  cum_n_freq = matrix(cumsum(nfreq), nrow = n_freq, ncol = n_frames)

  i_frame = numeric(n_freq_total)
  i_freq = numeric(n_freq_total)
  i_source = numeric(n_freq_total)
  freq =  numeric(n_freq_total)
  X = matrix(complex(n_freq_total * n_sources),
             nrow = n_freq_total,
             ncol = n_sources)

  for (i in seq_len(n_freq))
    for (j in seq_len(n_frames))
    {
      ii = c()
      if (i == 1)
      {
        if (j == 1)
          ii = 1:cum_n_freq[i, j]
        else
          ii = (cum_n_freq[n_freq, j - 1] + 1):cum_n_freq[i, j]
      }
      else
        ii = (cum_n_freq[i - 1, j] + 1):cum_n_freq[i, j]
      i_frame[ii] = j
      i_freq[ii] = i
      i_source[ii] = i + (j - 1) * n_freq
      freq[ii] = frame_freq[[j]][[i]]$freq
    }

  rec$sources_signal_init = list(
    i_frame = i_frame,
    i_freq = i_freq,
    i_source = i_source,
    freq = freq,
    X = X
  )


  return(rec)
}



init_receptor_signal <- function(rec,
                                 freq_limits = rbind(c(0, 500),
                                                     c(500, 4000),
                                                     c(4000, 10000)))
{
  n_receptors = length(rec$labels)
  rec = calculate_ffts(rec)


  n_frames = length(rec$frames)
  n_freq = nrow(freq_limits)

  nfreq = matrix(numeric(n_frames * n_freq), nrow = n_freq, ncol = n_frames)

  frame_freq = lapply(
    rec$frames,
    FUN = function(frame)
      freq_limits_to_their_indexes(
        f_rec = rec,
        frame = frame,
        freq_limits = freq_limits
      )
  )

  for (i_frame in seq_len(n_frames))
    for (i_freq in seq_len(n_freq))
      nfreq[i_freq, i_frame] = length(frame_freq[[i_frame]][[i_freq]]$index)

  n_freq_total = sum(nfreq)
  cum_n_freq = matrix(cumsum(nfreq), nrow = n_freq, ncol = n_frames)

  i_frame = numeric(n_freq_total)
  i_freq = numeric(n_freq_total)
  i_source = numeric(n_freq_total)
  freq =  numeric(n_freq_total)
  Y = matrix(complex(n_freq_total * n_receptors),
             nrow = n_freq_total,
             ncol = n_receptors)

  for (i in seq_len(n_freq))
    for (j in seq_len(n_frames))
    {
      ii = c()
      if (i == 1)
      {
        if (j == 1)
          ii = 1:cum_n_freq[i, j]
        else
          ii = (cum_n_freq[n_freq, j - 1] + 1):cum_n_freq[i, j]
      }
      else
        ii = (cum_n_freq[i - 1, j] + 1):cum_n_freq[i, j]
      i_frame[ii] = j
      i_freq[ii] = i
      i_source[ii] = i + (j - 1) * n_freq
      freq[ii] = frame_freq[[j]][[i]]$freq
      index = frame_freq[[j]][[i]]$index
      for (i_rec in seq_len(n_receptors))
        Y[ii, i_rec] = rec$frames[[1]]$fft[[i_rec]][index]
    }

  rec$receptors_signal_init = list(
    i_frame = i_frame,
    i_freq = i_freq,
    i_source = i_source,
    freq = freq,
    Y = Y
  )


  return(rec)
}






init_parameters <- function(rec,
                            n_sources = NULL,
                            x_min = NULL,
                            x_max = NULL,
                            y_min = NULL,
                            y_max = NULL,
                            freq_limits = rbind(c(0, 500),
                                                c(500, 4000),
                                                c(4000, 10000)))
{
  rec = init_receptors(rec)
  rec = init_frames_gain(rec)
  rec = init_sources_position(
    rec,
    n_sources = n_sources,
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max,
    freq_limits = freq_limits
  )

  rec = init_sources_signal(rec, n_sources = n_sources, freq_limits = freq_limits)

  rec = init_receptor_signal(rec, freq_limits = freq_limits)

  rec$receptors_init -> rec$receptors_step
  rec$frames_gain_init -> rec$frames_gain_step
  rec$sources_pos_init -> rec$sources_pos_step
  rec$sources_signal_init -> rec$sources_signal_step

  rec = get_source_receptor_distance(rec)

  return(rec)

}






parameters_to_beta <- function(f_rec)
{
  n_receptors = length(f_rec$labels)
  n_frames = length(f_rec$frames)
  count_receptors_par = (n_receptors - 1) * 4 - 1
  count_frames_gain_par = (nrow(f_rec$frames_gain_init$logG_diff) - 1) * n_receptors
  count_sources_pos_par = length(f_rec$sources_pos_init$x) + length(f_rec$sources_pos_init$y)
  count_sources_signal_par = length(f_rec$sources_signal_init$X)

  counts = c(
    receptors = count_receptors_par ,
    frames_gain = count_frames_gain_par,
    sources_pos = count_sources_pos_par,
    sources_signal = count_sources_signal_par
  )

  total = sum(counts)
  cumtotal = cumsum(counts)

  beta = complex(total)

  beta[1:cumtotal['receptors']] =
    c(
      f_rec$receptors_step$x[3:n_receptors],
      f_rec$receptors_step$y[2:n_receptors],
      f_rec$receptors_step$logG[2:n_receptors],
      f_rec$receptors_step$offset[2:n_receptors]
    )

  beta[(cumtotal[1] + 1):cumtotal[2]] = f_rec$frames_gain_step$logG_diff[2:n_frames, ]

  beta[(cumtotal[2] + 1):cumtotal[3]] = c(f_rec$sources_pos_step$x,
                                          f_rec$sources_pos_step$y)

  beta[(cumtotal[3] + 1):cumtotal[4]] = f_rec$sources_signal_step$X

  return(beta)

}

get_source_receptor_distance <- function(f_rec)
{
  n_receptors = length(f_rec$labels)
  f_rec$distances = lapply(seq_len(n_receptors),
                           function(i)
                           {
                             r_x = f_rec$receptors_step$x[i]
                             r_y = f_rec$receptors_step$y[i]
                             s_x = f_rec$sources_pos_step$x
                             s_y = f_rec$sources_pos_step$y
                             d = sqrt((r_x - s_x) ^ 2 + (r_y - s_y) ^ 2)
                             return(d)

                           })

  return(f_rec)
}



beta_to_parameters <- function(beta, f_rec)
{
  n_receptors = length(f_rec$labels)
  n_frames = length(f_rec$frames)
  count_receptors_x = (n_receptors - 2)
  count_receptors_y = (n_receptors - 1)
  count_receptors_logG = (n_receptors - 1)
  count_receptors_offset = (n_receptors - 1)


  count_frames_gain_par = (nrow(f_rec$frames_gain_init$logG_diff) - 1) * n_receptors
  count_sources_pos_x = length(f_rec$sources_pos_init$x)
  count_sources_pos_y = length(f_rec$sources_pos_init$y)
  count_sources_signal_par = length(f_rec$sources_signal_init$X)

  counts = c(
    count_receptors_x,
    count_receptors_y,
    count_receptors_logG,
    count_receptors_offset,
    count_frames_gain_par,
    count_sources_pos_x,
    count_sources_pos_y,
    count_sources_signal_par
  )

  total = sum(counts)
  cumtotal = cumsum(counts)


  Re(beta[1:cumtotal[1]]) ->
    f_rec$receptors_step$x[3:n_receptors]
  Re(beta[(cumtotal[1] + 1):cumtotal[2]]) ->
    f_rec$receptors_step$y[2:n_receptors]
  Re(beta[(cumtotal[2] + 1):cumtotal[3]]) ->
    f_rec$receptors_step$logG[2:n_receptors]
  Re(beta[(cumtotal[3] + 1):cumtotal[4]]) ->
    f_rec$receptors_step$offset[2:n_receptors]


  Re(beta[(cumtotal[4] + 1):cumtotal[5]]) ->
    f_rec$frames_gain_step$logG_diff[2:n_frames, ]

  # actualize logG for frames

  f_rec$frames_gain_step$logG = apply(f_rec$frames_gain_step$logG_diff, 2, cumsum)



  Re(beta[(cumtotal[5] + 1):cumtotal[6]]) -> f_rec$sources_pos_step$x[seq_along(f_rec$sources_pos_step$x)]

  Re(beta[(cumtotal[6] + 1):cumtotal[7]]) -> f_rec$sources_pos_step$y[seq_along(f_rec$sources_pos_step$y)]

  beta[(cumtotal[7] + 1):cumtotal[8]] -> f_rec$sources_signal_step$X[seq_along(f_rec$sources_signal_step$X)]


  # actualize distances

  f_rec = get_source_receptor_distance(f_rec)
  return(f_rec)
}

parameters_to_xdata <- function(f_rec)
{
  return(
    cbind(
      seq_along(f_rec$sources_signal_step$i_frame),
      f_rec$sources_signal_step$i_frame,
      f_rec$sources_signal_step$i_freq,
      f_rec$sources_signal_step$i_source,
      f_rec$sources_signal_step$freq

    )
  )

}

parameters_to_ydata <- function(f_rec)
{
  return(f_rec$receptors_signal_init$Y)
}



get_Amplitudes <- function(f_rec, xdata)
{
  n_receptors = length(f_rec$labels)
  Amplitudes =
    lapply(seq_len(n_receptors),
           function(i)
           {
             dist = f_rec$distances[[i]]
             d = dist[xdata[, 4],]
             w = 2 * pi * xdata[, 5]
             logG = f_rec$receptors_step$logG[i] +
               f_rec$frames_gain_step$logG[xdata[, 2], i]
             offset = f_rec$receptors_step$offset[i]
             A = exp(logG - 1i * w *
                       (d / f_rec$velocity_of_sound + offset)) / d
             return(A)
           })

  return(Amplitudes)
}

get_source_signal <- function(f_rec, xdata)
{
  return(f_rec$sources_signal_step$X[xdata[, 1],])
}

get_receptor_signal <- function(f_rec, xdata)
{
  return(f_rec$receptors_signal_init$Y[xdata[, 1],])
}


get_predicted_signal <- function(f_rec, xdata)
{
  A = get_Amplitudes(f_rec, xdata)
  X = get_source_signal(f_rec, xdata)
  n = nrow(X)

  return(vapply(A, function(a)
    rowSums(a * X), complex(n)))

}




sqr_sum_signal <- function(f_rec, beta, xdata)
{
  f_rec = beta_to_parameters(beta, f_rec)
  Yfit = get_predicted_signal(f_rec, xdata)
  Y = get_receptor_signal(f_rec, xdata)
  error = Yfit - Y
  return(sum(Mod(error)))

}





global_optimization <- function(f_rec)
{


}
