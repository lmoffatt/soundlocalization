

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
  i_y_par = max(i_x_par) + 2:n_receptors
  i_logG_beta = max(i_y_beta) + seq_len(length(logG) - 1)
  i_logG_par = max(i_y_par) + seq_len(length(logG) - 1) + 1
  i_offset_beta = max(i_logG_beta) + seq_len(length(offset) - 1)
  i_offset_par = max(i_logG_par) + seq_len(length(offset) - 1) + 1


  return(
    list(
      values = cbind(x, y, logG, offset),
      indexes_par_to_beta = c(i_x_beta, i_y_beta, i_logG_beta, i_offset_beta),
      indexes_beta_to_par = c(i_x_par, i_y_par, i_logG_par, i_offset_par) ,
      offset = max(i_offset_beta)
    )
  )
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
                                              offset_beta,
                                              offset_xdata,
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
  i_logG_beta = offset_beta + seq_len(n_receptors - 1)
  i_logG_par = seq_len(n_receptors - 1) + 1

  fs = freq_limits_to_their_indexes(f_rec, frame, freq_limits)
  sources = Reduce(function(prev, f)
  {
    pos = distribute_sources(f_rec, n_sources, x_min, x_max, y_min, y_max)
    i_pos_beta = prev$offset_beta + seq_along(pos)
    offset_xdata = prev$offset_xdata + length(f$index)
    i_pos_par = seq_along(pos)

    i_pos_xdata =



      X = matrix(complex(length(f$index) * nrow(pos)),
                 nrow = length(f$index),
                 ncol = nrow(pos))
    i_X_beta = max(i_pos_beta) + seq_along(X)
    i_X_par = seq_along(X)





    offset_beta = max(i_X_beta)

    return(list(
      offset_beta = offset_beta,
      offset_xdata = offset_xdata,
      values = c(prev$values,
                 list(
                   list(
                     pos = pos,
                     i_pos_beta = i_pos_beta,
                     i_pos_par = i_pos_par,
                     X = X,
                     freq = f$freq,
                     i_X_beta = i_X_beta,
                     i_X_par = i_X_par,
                     offset_xdata = offset_xdata
                   )
                 ))
    ))
  },
  fs,
  list(
    offset_beta = max(i_logG_beta),
    offset_xdata = offset_xdata,
    values = list()
  ))


  return(
    list(
      logG = logG,
      i_logG_beta = i_logG_beta,
      i_logG_par = i_logG_par,
      sources = sources$values,
      offset_beta = sources$offset_beta,
      offset_xdata = sources$offset_xdata

    )
  )
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
          offset_beta =  prev$offset_beta,
          offset_xdata = prev$offset_xdata,
          n_sources,
          x_min,
          x_max = x_max,
          y_min ,
          y_max ,
          freq_limits
        )
        offset_beta = sources$offset_beta
        offset_xdata = sources$offset_xdata
        return(list(
          offset_beta = offset_beta,
          offset_xdata = offset_xdata,
          values = c(prev$values, list(sources))
        ))
      }
      ,
      rec$frames,
      list(
        offset_beta = receptors$offset,
        offset_xdata = 0,
        values = list()
      )
    )
  ))
}


init_receptors <- function(rec)
{
  n_receptors = length(rec$labels)
  x = rec$x
  y = rec$y
  logG = rep(0, n_receptors)
  offset = rep(0, n_receptors)
  rec$recptors_init = data.frame(
    x = x,
    y = y,
    logG = logG,
    offset = offset
  )
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
  X = matrix(complex(n_freq_total*n_sources), nrow = n_freq_total,ncol = n_sources)
  Y = matrix(complex(n_freq_total*n_receptors), nrow = n_freq_total,ncol = n_receptors)

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
         Y[ii,i_rec] = rec$frames[[1]]$fft[[i_rec]][index]
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
  Y = matrix(complex(n_freq_total*n_receptors), nrow = n_freq_total,ncol = n_receptors)

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
        Y[ii,i_rec] = rec$frames[[1]]$fft[[i_rec]][index]
    }

  rec$recptors_signal_init = list(
    i_frame = i_frame,
    i_freq = i_freq,
    i_source = i_source,
    freq = freq,
    Y = Y
  )


  return(rec)
}






recording_to_seaprate_parameters <- function(rec,
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
  # idea--> to extract the beta (parameter to be optimized),
  #  the ydata (the observed data)
  # and the xdata (the coordinates for each ydata)
  # and how to recover from the beta the parameter for each individual ydata





}





parameters_to_beta <- function(parameters)
{
  n = parameters$frames$offset
  out = complex(n)
  recp = parameters$receptors$values
  i_beta = parameters$receptors$indexes_par_to_beta
  i_par = parameters$receptors$indexes_beta_to_par
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

parameters_to_xdata <- function(parameters)
{
  n_tot = parameters$frames$offset_xdata

  i_frame = integer(n_tot)

  freq = numeric(n_tot)
  i_freq_range = integer(n_tot)

  i2 = 0
  for (iframe in seq_along(parameters$frames$values))
  {
    for (ifreq in seq_along(parameters$frames$values[[iframe]]$sources))
    {
      i1 = i2 + 1
      i2 = parameters$frames$values[[iframe]]$sources[[ifreq]]$offset_xdata
      i_frame[i1:i2] = iframe
      i_freq_range[i1:i2] = ifreq
      freq[i1:i2] = parameters$frames$values[[iframe]]$sources[[ifreq]]$freq
    }
  }



  matrix(c(
    i_frame = i_frame,
    i_freq_range = i_freq_range,
    freq = freq
  ),
  ncol = 3)
}

get_source_receptor_distance <- function(parameters, xdata)
{
  n_receptors = nrow(parameters$receptors$values)
  n_xdata = nrow(xdata)
  nsources = length(parameters$frames$values[[1]]$sources)
  distances = list()
  for (ir in seq_len(nrow(parameters$receptors$values)))
  {
    dist = matrix(0, nrow = n_xdata, ncol = nsources)
    r_x = parameters$receptors$values[ir, 'x']
    r_y = parameters$receptors$values[ir, 'y']
    i_frame = xdata[, 1]
    i_freq = xdata[, 2]
    for (j in seq_len(nsources))
    {
      s_x =  parameters$frames$values[[i_frame]]$sources[[i_freq]]$pos[j, 1]
      s_y =  parameters$frames$values[[i_frame]]$sources[[i_freq]]$pos[j, 2]
      d = sqrt((r_x - s_x) ^ 2 + (r_y - s_y) ^ 2)
      dist[, j] = d
    }
    distances[[ir]] <- dist
  }

  return(distances)
}



global_optimization <- function(f_rec)
{
  get_source_Amplitude <- function(parameters, xdata)
  {
    receptors = get_receptors(beta)

    frames = get_frames(beta)

    fG = get_frame_Gain(beta)

    distances = get_source_receptor_distance(parameters, xdata)

    dist_recip = 1 / distances

    lags = distances / velocity_of_sound

    A = dist_recip[i_rec, ] * exp(-i * lags[i_rec, ] * w)
    return(A)
  }


  predict_single_fft <- function(beta, xdata, paramters)
  {
    parameters = beta_to_parameters(beta, parameters)
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
