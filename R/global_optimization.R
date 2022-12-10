












recording_to_receptors <- function(rec, i_offset = 0)
{
  n_receptors = length(rec$labels)
  logG = rep(0, n_receptors)
  offset = rep(0, n_receptors)
  x = rec$x
  i_x_beta = i_offset + seq_len(length(x) - 2)
  i_x_par = 3:n_receptors
  y = rec$y
  i_y_beta = max(i_x_beta) + seq_len(length(y) - 1)
  i_y_par = max(i_x_par)+ 2:n_receptors
  i_logG_beta = max(i_y_beta) + seq_len(length(logG) - 1)
  i_logG_par = max(i_y_par) + seq_len(length(logG) - 1) +1
  i_offset_beta = max(i_logG_beta) + seq_len(length(offset) - 1)
  i_offset_par = max(i_logG_par) + seq_len(length(offset) - 1) +1


  return(list(
    values = cbind(x, y, logG, offset),
    indexes_par_to_beta = c(i_x_beta, i_y_beta,i_logG_beta, i_offset_beta),
    indexes_beta_to_par = c(i_x_par, i_y_par,i_logG_par, i_offset_par) ,
    offset = max(i_offset_beta)
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

  logG = numeric(n_receptors)
  i_logG_beta = offset + seq_len(n_receptors - 1)
  i_logG_par = seq_len(n_receptors - 1) + 1

  fs = freq_limits_to_their_indexes(f_rec, frame, freq_limits)
  sources = Reduce(function(prev, f)
  {
    pos = distribute_sources(f_rec, n_sources, x_min, x_max, y_min, y_max)
    i_pos_beta = prev$offset + seq_along(pos)
    i_pos_par = seq_along(pos)

    X = matrix(complex(length(f$index) * nrow(pos)),
               nrow = length(f$index),
               ncol = nrow(pos))
    i_X_beta = max(i_pos_beta) + seq_along(X)
    i_X_par = seq_along(X)
    offset = max(i_X_beta)
    return(list(offset = offset,
                values = c(prev$values,
                           list(
                             list(
                               pos = pos,
                               i_pos_beta = i_pos_beta,
                               i_pos_par = i_pos_par,
                               X = X,
                               i_X_beta = i_X_beta,
                               i_X_par = i_X_par

                             )
                           ))))
  }, fs, list(offset = max(i_logG_beta), values = list()))


  return(list(
    logG = logG,
    i_logG_beta = i_logG_beta,
    i_logG_par =i_logG_par,
    sources = sources$values,
    offset = sources$offset
  ))
}


recording_to_parameters <- function(rec,
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

  return(list(
    receptors = receptors,
    frames = Reduce(
      function(prev, frame)
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
      list(offset = receptors$offset, values = list())
    )
  ))
}

parameters_to_beta <- function(parameters)
{
  n = parameters$frames$offset
  out = complex(n)
  recp = parameters$receptors$values
  i_beta = parameters$receptors$indexes_par_to_beta
  i_par =parameters$receptors$indexes_beta_to_par
  out[i_beta] = recp[i_par]

  for (i in seq_along(parameters$frames$values))
  {
    fr = parameters$frames$values[[i]]
    out[fr$i_logG_beta] <- fr$logG[fr$i_logG_par]
    for (j in seq_along(fr$sources))
    {
      so = fr$sources[[j]]
      out[so$i_pos_beta] = so$pos
      out[so$i_X_beta] = so$X
    }
  }

  return(out)
}

beta_to_parameters <- function(beta, parameters)
{
  n = parameters$frames$offset
  i_beta = parameters$receptors$indexes_par_to_beta
  i_par = parameters$receptors$indexes_beta_to_par

  recp = parameters$receptors$values

  parameters$receptors$values[i_par] = Re(beta[i_beta])
  for (i in seq_along(parameters$frames$values))
  {
    fr = parameters$frames$values[[i]]
    Re(beta[fr$i_logG_beta]) -> parameters$frames$values[[i]]$logG[fr$i_logG_par]
    for (j in seq_along(fr$sources))
    {
      so = parameters$frames$values[[i]]$sources[[j]]
      Re(beta[so$i_pos_beta]) -> parameters$frames$values[[i]]$sources[[j]]$pos[so$i_pos_par]
      beta[so$i_X_beta] -> parameters$frames$values[[i]]$sources[[j]]$X[so$i_X_par]
    }
  }

  return(parameters)
}


get_sources<-function(parameters)
{
   ipos
}




global_optimization <- function(f_rec)
{
  get_source_Amplitude <- function(parameters, xdata)
  {
    i_rec = xdata$i_rec

    i_frame = xdata$i_frame



    w = xdata$w

    sources = get_sources(parameters)

    receptors = get_receptors(beta)

    frames = get_frames(beta)

    fG = get_frame_Gain(beta)

    distances = get_source_receptor_distance()

    dist_recip = 1 / distances

    lags = distances / velocity_of_sound

    A = dist_recip[i_rec, ] * exp(-i * lags[i_rec, ] * w)
    return(A)
  }


  predict_single_fft <- function(beta, xdata)
  {
    parameters=beta_to_parameters(beta,parameters)
    A = get_source_Amplitude(parameters, xdata)
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
