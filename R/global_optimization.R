









recording_to_receptors <- function(rec, i_offset = 0)
{
  n_receptors = length(rec$labels)
  logG = rep(0, n_receptors)
  offset = rep(0, n_receptors)
  x = rec$x
  i_x = i_offset + seq_len(length(x) - 2)
  y = rec$y
  i_y = max(i_x) + seq_len(length(y) - 1)
  i_logG = max(i_y) + seq_len(length(logG) - 1)
  i_offset = max(i_logG) + seq_len(length(offset) - 1)
  return(list(
    values = cbind(x, y, logG, offset),
    indexes = list(
      i_x = i_x,
      i_y = i_y,
      i_logG = i_logG,
      i_offset = i_offset
    ),
    offset = max(i_offset)
  ))
}


distribute_sources <-
  function(rec,
           n_sources = NULL,
           x_min = NULL,
           x_max = NULL,
           y_min = NULL,
           y_max = NULL)
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
    x = x_min + (x_max - x_min) * runif(n_sources)
    y = y_min + (y_max - y_min) * runif(n_sources)

    return(cbind(x, y))



  }


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




recording_to_single_frame_sources_position <- function(f_rec,
                                                       frame,
                                                       n_sources = NULL,
                                                       x_min = NULL,
                                                       x_max = NULL,
                                                       y_min = NULL,
                                                       y_max = NULL,
                                                       freq_limits = rbind(c(0, 500),
                                                                           c(500, 4000),
                                                                           c(4000, 10000)))
{
  n_receptors = length(f_rec$labels)
  if (is.null(n_sources))
    n_sources = n_receptors - 1

  fs = freq_limits_to_their_indexes(f_rec, frame, freq_limits)
  sources = Reduce(function(prev, f)
  {
    pos = distribute_sources(f_rec, n_sources, x_min, x_max, y_min, y_max)
    return(rbind(prev, pos))
  }, fs, data.frame())

  return(sources)
}


recording_to_single_frame_sources_signal <- function(f_rec,
                                                     frame,
                                                     offset,
                                                     n_sources = NULL,
                                                     x_min = NULL,
                                                     x_max = NULL,
                                                     y_min = NULL,
                                                     y_max = NULL,
                                                     freq_limits = rbind(c(0, 500),
                                                                         c(500, 4000),
                                                                         c(4000, 10000)))
{
  n_receptors = length(f_rec$labels)
  if (is.null(n_sources))
    n_sources = n_receptors - 1

  fs = freq_limits_to_their_indexes(f_rec, frame, freq_limits)
  sources = Reduce(function(prev, f)
  {
    pos = distribute_sources(f_rec, n_sources, x_min, x_max, y_min, y_max)

    X = matrix(complex(length(f$index) * n_sources),
               nrow = length(f$index),
               ncol = n_sources)
    return(rbind(prev, X))
  }, fs, matrix(complex(0),
                nrow = 0,
                ncol = n_sources))

  return(sources)
}



recording_to_single_frame_sources <- function(f_rec,
                                              frame,
                                              offset,
                                              n_sources = NULL,
                                              x_min = NULL,
                                              x_max = NULL,
                                              y_min = NULL,
                                              y_max = NULL,
                                              freq_limits = rbind(c(0, 500),
                                                                  c(500, 4000),
                                                                  c(4000, 10000)))
{
  n_receptors = length(f_rec$labels)
  if (is.null(n_sources))
    n_sources = n_receptors - 1

  fs = freq_limits_to_their_indexes(f_rec, frame, freq_limits)
  sources = Reduce(function(prev, f)
  {
    pos = distribute_sources(f_rec, n_sources, x_min, x_max, y_min, y_max)
    i_pos = prev$offset + seq_along(pos)

    X = matrix(complex(length(f$index) * nrow(pos)),
               nrow = length(f$index),
               ncol = nrow(pos))
    i_X = max(i_pos) + seq_along(X)
    offset = max(i_X)
    return(list(offset = offset,
                values = c(prev$values,
                           list(
                             list(
                               pos = pos,
                               i_pos = i_pos,
                               X = X,
                               i_X = i_X
                             )
                           ))))
  }, fs, list(offset = offset + 1, values = list()))

  return(list(
    logG = 0,
    i_logG = offset + 1,
    sources = sources$values,
    offset = sources$offset
  ))
}


recording_to_receptors_frame_group <- function(rec,
                                               n_sources = NULL,
                                               i_offset = 0,
                                               x_min = NULL,
                                               x_max = NULL,
                                               y_min = NULL,
                                               y_max = NULL,
                                               freq_limits = rbind(c(0, 500),
                                                                   c(500, 4000),
                                                                   c(4000, 10000)))
{
  receptors = recording_to_receptors(rec, i_offset = i_offset)



  source_positions = Reduce(function(prev, frame)
  {
    pos = recording_to_single_frame_sources_position(rec,
                                                     frame,
                                                     n_sources,
                                                     x_min,
                                                     x_max = x_max,
                                                     y_min ,
                                                     y_max ,
                                                     freq_limits)
    return(rbind(prev, pos))
  }
  ,
  rec$frames,
  data.frame())


  source_signals = Reduce(function(prev, frame)
  {
    pos = recording_to_single_frame_sources_signal(rec,
                                                   frame,
                                                   n_sources,
                                                   x_min,
                                                   x_max = x_max,
                                                   y_min ,
                                                   y_max ,
                                                   freq_limits)
    return(rbind(prev, pos))
  }
  ,
  rec$frames,
  data.frame())

  return(
    list(
      receptors = receptors,
      source_positions = source_positions,
      source_signals = source_signals
    )
  )

}

recording_to_receptors_frame <- function(rec,
                                         n_sources = NULL,
                                         i_offset = 0,
                                         x_min = NULL,
                                         x_max = NULL,
                                         y_min = NULL,
                                         y_max = NULL,
                                         freq_limits = rbind(c(0, 500),
                                                             c(500, 4000),
                                                             c(4000, 10000)))
{
  receptors = recording_to_receptors(rec, i_offset = i_offset)



  frames = Reduce(function(prev, frame)
  {
    sources = recording_to_single_frame_sources(
      rec,
      frame,
      offset = prev$offset,
      n_sources,
      x_min,
      x_max = x_max,
      y_min ,
      y_max ,
      freq_limits
    )
    offset = sources$offset
    return (list(offset = offset, values = c(prev$values, list(sources))))
  }
  ,
  rec$frames,
  list(offset = receptors$offset, values = list()))


}



list_to_beta <- function(parameters)
{
  n = parameters$frames$offset
  out = complex(n)
  recp = parameters$receptors$values
  i_recp = parameters$receptors$indexes
  out[i_recp$i_x] = recp[3:nrow(recp), "x"]
  out[i_recp$i_y] = recp[2:nrow(recp), "y"]
  out[i_recp$i_logG] = recp[2:nrow(recp), "logG"]
  out[i_recp$i_offset] = recp[2:nrow(recp), "offset"]

  for (i in seq_along(parameters$frames$values))
  {
    fr = parameters$frames$values[[i]]
    out[fr$i_logG] <- fr$logG
    for (j in seq_along(fr$sources))
    {
      so = fr$sources[[j]]
      out[so$i_pos] = so$pos
      out[so$i_X] = so$X
    }
  }

  return(out)
}

beta_to_list <- function(beta, parameters)
{
  n = parameters$frames$offset
  i_recp = parameters$receptors$indexes

  parameters$receptors$values[3:nrow(recp), "x"] = beta[i_recp$i_x]
  parameters$receptors$values[2:nrow(recp), "y"] = beta[i_recp$i_y]
  parameters$receptors$values[2:nrow(recp), "logG"] = beta[i_recp$i_logG]
  parameters$receptors$values[2:nrow(recp), "offset"] = beta[i_recp$i_offset]

  for (i in seq_along(parameters$frames$values))
  {
    fr = parameters$frames$values[[i]]
    beta[fr$i_logG] -> parameters$frames$values[[i]]$logG
    for (j in seq_along(fr$sources))
    {
      so = parameters$frames$values[[i]]$sources[[j]]
      beta[so$i_pos] -> parameters$frames$values[[i]]$sources[[j]]$pos
      beta[so$i_X] -> parameters$frames$values[[i]]$sources[[j]]$X
    }
  }

  return(parameters)
}


get_receptors <- function(beta, parameters)
{
  n = parameters$frames$offset
  recep = parameters$receptors$values
  i_recp = parameters$receptors$indexes

  recep[3:nrow(recep), "x"] = beta[i_recp$i_x]
  recep[2:nrow(recep), "y"] = beta[i_recp$i_y]
  recep[2:nrow(recep), "logG"] = beta[i_recp$i_logG]
  recep[2:nrow(recep), "offset"] = beta[i_recp$i_offset]
  return(recep)
}

get_sources_pos <- function(beta, parameters)
{
  for (i in seq_along(parameters$frames$values))
  {
    fr = parameters$frames$values[[i]]
    beta[fr$i_logG] -> parameters$frames$values[[i]]$logG
    for (j in seq_along(fr$sources))
    {
      so = parameters$frames$values[[i]]$sources[[j]]
      beta[so$i_pos] -> parameters$frames$values[[i]]$sources[[j]]$pos
      beta[so$i_X] -> parameters$frames$values[[i]]$sources[[j]]$X
    }
  }

  return(parameters)
}



global_optimization <- function(f_rec)
{
  get_source_Amplitude <- function(beta, xdata)
  {
    i_rec = xdata$i_rec

    i_frame = xdata$i_frame

    w = xdata$w

    i_sources = xdata$i_sources

    sources = get_sources(beta)

    receptors = get_receptors(beta)

    frames = get_frames(beta)

    fG = get_frame_Gain(beta)

    distances = get_source_receptor_distance(sources, receptors)

    dist_recip = 1 / distances

    lags = distances / velocity_of_sound

    A = dist_recip[i_rec,] * exp(-i * lags[i_rec,] * w)
    return(A)
  }


  predict_single_fft <- function(beta, xdata)
  {
    A = get_source_Amplitude(beta, xdata)
    X = get_source_fft(beta, xdata)
    return(A %*% X)

  }




  predict_fft <- function(beta, xdata)
  {
    beta_to_list(beta, parameters = parameters)

  }


  parameters_to_beta <- function(parameters) {

  }
  an_receptors = length(f_rec$labels)

  f_rec <- calculate_ffts(f_rec)


  parameters <- recording_to_receptors_frame(f_rec, n_sources = 2)


}
