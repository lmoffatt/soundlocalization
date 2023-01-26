










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
  n_receptors = get_number_of_receptors(rec)
  x = rec$x
  y = rec$y
  logG = rep(0, n_receptors)
  offset = rep(0, n_receptors)
  return(list(
    x = x,
    y = y,
    logG = logG,
    offset = offset
  ))
}

init_frames_gain <- function(rec)
{
  n_receptors = get_number_of_receptors(rec)
  n_frames = get_number_of_frames(rec)


  logG_diff = matrix(0,
                     nrow = n_frames,
                     ncol = n_receptors)

  logG = apply(logG_diff, 2, cumsum)


  return(list(logG_diff = logG_diff,
              logG = logG))

  return(rec)
}




init_sources <- function(rec,
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

  n_receptors = get_number_of_receptors(rec)
  if (is.null(n_sources))
    n_sources = n_receptors - 1


  n_frames = get_number_of_frames(rec)
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



  return(list(
    i_frame_pos = i_frame_pos,
    i_freq_pos = i_freq_pos,
    x = x,
    y = y
  ))

}


will_be_calculate_sources_signal <- function(rec,
                                             n_sources = NULL,
                                             freq_limits = rbind(c(0, 500),
                                                                 c(500, 4000),
                                                                 c(4000, 10000)))
{
  n_receptors = get_number_of_receptors(rec)
  if (is.null(n_sources))
    n_sources = n_receptors - 1


  n_frames = get_number_of_frames(rec)
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
  X = matrix(
    0 * rnorm(n_freq_total * n_sources) + 0 * rnorm(n_freq_total * n_sources) *
      1i,
    nrow = n_freq_total,
    ncol = n_sources
  )

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



obtain_receptor_signal <- function(rec,
                                   freq_limits = rbind(c(0, 500),
                                                       c(500, 4000),
                                                       c(4000, 10000)))
{
  n_receptors = get_number_of_receptors(rec)
  n_frames = get_number_of_frames(rec)
  n_freq = nrow(freq_limits)
  n_fre_frames = n_frames * n_freq


  rec = calculate_ffts(rec)

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
             nrow = n_receptors,
             ncol = n_freq_total)

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
        Y[i_rec, ii ] = rec$frames[[1]]$fft[[i_rec]][index]
    }

  return(list(
    i_frame = i_frame,
    i_freq = i_freq,
    i_source = i_source,
    freq = freq,
    Y = Y
  ))


  return(rec)
}






init_model <- function(rec,
                       n_sources = NULL,
                       x_min = NULL,
                       x_max = NULL,
                       y_min = NULL,
                       y_max = NULL,
                       freq_limits = rbind(c(0, 500),
                                           c(500, 4000),
                                           c(4000, 10000)),
                       include_frame_gain)
{
  if (include_frame_gain)
    return(
      list(
        receptors = init_receptors(rec),
        frames_gain = init_frames_gain(rec),
        sources = init_sources(
          rec,
          n_sources = n_sources,
          x_min = x_min,
          x_max = x_max,
          y_min = y_min,
          y_max = y_max,
          freq_limits = freq_limits
        ),
        receptor_signal = obtain_receptor_signal(rec = rec, freq_limits = freq_limits)
      )
    )
  else
    return(
      list(
        receptors = init_receptors(rec),
        sources = init_sources(
          rec,
          n_sources = n_sources,
          x_min = x_min,
          x_max = x_max,
          y_min = y_min,
          y_max = y_max,
          freq_limits = freq_limits
        ),
        receptor_signal = obtain_receptor_signal(rec = rec, freq_limits = freq_limits)

      )
    )

}


get_number_of_source_frames <- function(model)
{
  return(nrow(model$sources$x))
}

get_number_of_simultaneous_sources <- function(model)
{
  return(ncol(model$sources$x))
}




get_number_of_sources <- function(model)
{
  return(length(model$sources$x))
}



model_to_beta <- function(f_rec,
                          model)
{
  n_receptors = get_number_of_receptors(f_rec)
  n_frames = get_number_of_frames(f_rec)
  count_receptors_par = (n_receptors - 1) * 4 - 1


  count_frames_gain_par = 0
  include_frame_gain = 'frames_gain' %in% names(model)
  if (include_frame_gain)
    count_frames_gain_par = (nrow(model$frames_gain$logG_diff) - 1) * n_receptors
  count_sources_pos_par = 2 * get_number_of_sources(model)
  counts = c(receptors = count_receptors_par ,
             frames_gain = count_frames_gain_par,
             sources_pos = count_sources_pos_par)

  total = sum(counts)
  cumtotal = cumsum(counts)

  beta = numeric(total)

  beta[1:cumtotal['receptors']] =
    c(
      model$receptors$x[3:n_receptors],
      model$receptors$y[2:n_receptors],
      model$receptors$logG[2:n_receptors],
      model$receptors$offset[2:n_receptors]
    )

  if (include_frame_gain)
  {
    beta[(cumtotal[1] + 1):cumtotal[2]] = model$frames_gain$logG_diff[2:n_frames,]
    model$frames_gain$logG = apply(model$frames_gain$logG_diff, 2, cumsum)


  }
  beta[(cumtotal[2] + 1):cumtotal[3]] = c(model$sources$x[seq_along(model$sources$x)],
                                          model$sources$y[seq_along(model$sources$y)])


  return(beta)

}

get_source_receptor_distance <- function(f_rec, model)
{
  n_receptors = get_number_of_receptors(f_rec)
  n_sim_sources = get_number_of_simultaneous_sources(model)
  n_source_frames = get_number_of_source_frames(model)

  return(vapply(seq_len(n_source_frames),
                function(i_frame)
                  vapply(seq_len(n_sim_sources),
                         function(i_source)
                           vapply(seq_len(n_receptors),
                                  function(i_receptor)
                                  {
                                    r_x = model$receptors$x[i_receptor]
                                    r_y = model$receptors$y[i_receptor]
                                    s_x = model$sources$x[i_frame, i_source]
                                    s_y = model$sources$y[i_frame, i_source]
                                    d = sqrt((r_x - s_x) ^ 2 + (r_y - s_y) ^ 2)
                                    return(d)
                                  }, numeric(1)),
                         numeric(n_receptors)),
                matrix(0, nrow = n_receptors, ncol = n_sim_sources)))
}



beta_to_model <- function(beta, f_rec,
                          model)
{
  n_receptors = get_number_of_receptors(f_rec)
  n_frames = get_number_of_frames(f_rec)
  count_receptors_par = (n_receptors - 1) * 4 - 1


  count_frames_gain_par = 0
  include_frame_gain = 'frames_gain' %in% names(model)
  if (include_frame_gain)
    count_frames_gain_par = (nrow(model$frames_gain$logG_diff) - 1) * n_receptors
  n_sources = get_number_of_sources(model)
  counts = c(count_receptors_par,
             count_frames_gain_par,
             n_sources, n_sources)

  total = sum(counts)
  cumtotal = cumsum(counts)


  cum_receptors = cumsum(c(n_receptors - 2, n_receptors - 1,
                           n_receptors - 1, n_receptors - 1))


  beta[1:cum_receptors[1]] ->
    model$receptors$x[3:n_receptors]
  beta[(cum_receptors[1] + 1):cum_receptors[2]] ->
    model$receptors$y[2:n_receptors]
  beta[(cum_receptors[2] + 1):cum_receptors[3]] ->
    model$receptors$logG[2:n_receptors]
  beta[(cum_receptors[3] + 1):cum_receptors[4]] ->
    model$receptors$offset[2:n_receptors]


  if (include_frame_gain)
  {
    beta[(cumtotal[1] + 1):cumtotal[2]] -> model$frames_gain$logG_diff[2:n_frames,]
    model$frames_gain$logG = apply(model$frames_gain$logG_diff, 2, cumsum)

  }
  beta[(cumtotal[2] + 1):cumtotal[3]] -> model$sources$x[seq_along(model$sources$x)]
  beta[(cumtotal[3] + 1):cumtotal[4]] -> model$sources$y[seq_along(model$sources$y)]


  return(model)

}




model_to_xdata <- function(model)
{
  return(
    t(cbind(
      seq_along(model$receptor_signal$i_frame),
      model$receptor_signal$i_frame,
      model$receptor_signal$i_freq,
      model$receptor_signal$i_source,
      model$receptor_signal$freq

    ))
  )

}

model_to_ydata <- function(model)
{
  return(model$receptor_signal$Y)
}






get_Amplitudes <- function(f_rec, model, distances, xdata)
{
  n_receptors = get_number_of_receptors(f_rec)
  n_sim_sources = get_number_of_simultaneous_sources(model)
  n_source_frames_freq = ncol(xdata)

  n_source_frames = get_number_of_source_frames(model)


  return(vapply(seq_len(n_source_frames_freq),
                function(i_data)
                  vapply(seq_len(n_sim_sources),
                         function(i_source)
                           vapply(seq_len(n_receptors),
                                  function(i_receptor)
                                  {
                                    w = 2 * pi * xdata[5,i_data]
                                    logG = model$receptors$logG[i_receptor] +
                                      model$frames_gain$logG[xdata[2,i_data], i_receptor]
                                    offset = model$receptors$offset[i_receptor]
                                    d = distances[i_receptor, i_source, xdata[4,i_data]]
                                    A = exp(logG - 1i * w *
                                              (d / f_rec$velocity_of_sound + offset)) / d
                                    return(A)
                                  }, complex(1)),
                         complex(n_receptors)),
                matrix(0 + 0i, nrow = n_receptors, ncol = n_sim_sources)))


}







get_source_signal <- function(f_rec, model, A, Ainv, xdata, ydata)
{

}

get_receptor_signal <- function(model, xdata)
{
  return(model$receptor_signal$Y[,xdata[1, ] ])
}





get_inverse_AA <- function(A, lambda = 0.0)
{
  n = dim(A)[3]
  n_sources = dim(A)[2]
  return( vapply(seq_len(n), function(i)
  {
    AHA = Conj(t(A[,,i])) %*% A[,,i]
    AHAla= AHA + diag(diag(AHA))*lambda
    k=rcond(AHAla)
    if (k <1e-7)
    {
       lambda = 1e-6
       AHAla = AHA + diag(diag(AHA))*lambda

    }
    return(solve(AHAla))
  },
  matrix(complex(n_sources*n_sources),n_sources,n_sources)))
}

get_X <- function(AAinv,A, Y)
{
  n = dim(A)[3]
  n_sources = dim(A)[2]
  return( vapply(seq_len(n), function(i)
  {
     AAinv[,,i]%*%Conj(t(A[,,i])) %*% Y[,i]
  },
  complex(n_sources)))
}

get_predicted_Y<-function(A,X)
{
  n = dim(A)[3]
  n_receptors = dim(A)[1]
  return( vapply(seq_len(n), function(i)
  {
    A[,,i] %*% X[,i]
  },
  complex(n_receptors)))

}


get_predicted_signal <-
  function(f_rec, model, xdata)
  {
    ydata= get_receptor_signal(model,xdata)
    distances = get_source_receptor_distance(f_rec, model)
    A = get_Amplitudes(f_rec, model, distances, xdata)

    AAinv = get_inverse_AA(A)

    X = get_X(AAinv,A,ydata)

    Y = get_predicted_Y(A,X)
    return(Y)


  }





predicted_signal <-
  function(f_rec,
           model,
           beta,
           xdata)
  {
    model = beta_to_model(beta, f_rec, model)
    Yfit = get_predicted_signal(f_rec, model, xdata)
    return(Yfit)

  }


gradient_of_predicted_signal <-
  function(f_rec,
           model,
           beta,
           xdata)
  {
    model = beta_to_model(beta, f_rec, model)
    d_beta = d_f__d_beta(f_rec, model, xdata)
    return(d_beta)

  }




global_optimization <- function(f_rec,
                                model,
                                number_of_chunks = 10,
                                maxiter = 100,
                                alpha = 0.1,
                                eps = 1e-16,
                                eps_J = 1e-8)
{
  predicted  <-
    function(beta, xdata) {
      return(predicted_signal(f_rec, model, beta, xdata))
    }

  Jacobian <-
    function(beta, xdata) {
      return(
        d_Y__d_beta(
          f_rec = f_rec,
          model,
          beta = beta,
          xdata = xdata)
      )
    }


  beta_init = model_to_beta(f_rec, model = model)
  xdata = model_to_xdata(model = model)
  ydata = model_to_ydata(model)

  opt = slm(
    f = predicted,

    x0 = beta_init,
    xdata = xdata,
    ydata = ydata,
    number_of_chunks = number_of_chunks,
    maxiter = maxiter,
    Jacobian = Jacobian,
    eps = eps,
    Y_beta
  )
  opt$sqr_diff = opt$sqrsum_tot - opt$sqrsum_prev
  opt$xdata = xdata
  opt$ydata = ydata
  opt$sqr_ydata = sum(Conj(ydata) * ydata)
  opt$sqr_diff_tot =  opt$sqrsum_tot - opt$sqr_ydata
  f_rec$opt = opt
  #model_opt = beta_to_model(opt$xt, f_rec,model)
  #f_rec$predicted_y = predicted(opt$xt, xdata)
  #f_rec$model_opt=model_opt
  return(f_rec)

}


