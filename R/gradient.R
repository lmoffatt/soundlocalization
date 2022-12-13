











sum_by_index <- function(values, indexes, icol)
{
  stopifnot(nrow(values) == nrow(indexes))
  num_indexes = max(indexes[, icol])
  ncols = ncol(values)
  out = c()
  if (is.complex(values))
    out = matrix(complex(num_indexes * ncols),
                 nrow = num_indexes,
                 ncol = ncols)
  else
    out = matrix(0, nrow = num_indexes, ncol = ncols)
  for (i in seq_len(nrow(values)))
    out[indexes[i, icol],] = out[indexes[i, icol],] + values[i,]
  return(out)
}







d_dist_dsource_pos_x <- function(f_rec)
{
  dist = f_rec$distances
  x_s = f_rec$sources_pos_step$x
  x_r = f_rec$receptors_init$x
  n_receptors = length(f_rec$labels)
  dist_sx = lapply(seq_len(n_receptors), function(i)
    (x_s - x_r[i]) / dist[[i]])
  return(dist_sx)
}

d_dist_dsource_pos_y <- function(f_rec)
{
  dist = f_rec$distances
  y_s = f_rec$sources_pos_step$y
  y_r = f_rec$receptors_init$y
  n_receptors = length(f_rec$labels)
  dist_sy = lapply(seq_len(n_receptors), function(i)
    (y_s - y_r[i]) / dist[[i]])
  return(dist_sy)
}





dA_d_distance <- function(f_rec, xdata)
{
  v = f_rec$velocity_of_sound
  w = 2 * pi * xdata[, 5]
  dist = lapply(f_rec$distances, function(d)
    d[xdata[, 4], ])
  A = get_Amplitudes(f_rec, xdata)
  n_receptor = length(f_rec$labels)
  A_dist = lapply(1:n_receptor, function(i)
    - A[[i]] * (1 / dist[[i]] + 1i / v * w))
  return(A_dist)
}

dA_doffset <- function(f_rec, xdata)
{
  w = 2 * pi * xdata[, 5]
  A = get_Amplitudes(f_rec, xdata)
  n_receptor = length(f_rec$labels)
  A_off = lapply(1:n_receptor, function(i)
    - A[[i]] * 1i * w)
  return(A_off)
}

dA_dlogG <- function(f_rec, xdata)
{
  A_logG = get_Amplitudes(f_rec, xdata)
  return(A_logG)
}

dY_dAmplitude <- function(f_rec, xdata)
{
  X = get_source_signal(f_rec, xdata)
  return(X)
}



df_dY <- function(f_rec, xdata)
{
  Y = get_receptor_signal(f_rec, xdata)
  Yfit = get_predicted_signal(f_rec, xdata)

  return(2 * (Yfit - Y))
}




df_dX <- function(f_rec, xdata)
{
  A = get_Amplitudes(f_rec, xdata)
  X = get_source_signal(f_rec, xdata)

  Y = get_receptor_signal(f_rec, xdata)
  n = nrow(X)
  Yfit = (vapply(A, function(a)
    rowSums(a * X), complex(n)))
  Y_dif = Y - Yfit

  f_X =  (vapply(A, function(a)
    rowSums(a * Y_dif), complex(n)))

  return(f_X)
}


df_d_distance <- function(f_rec, xdata)
{
  f_Y = df_dY(f_rec, xdata)
  Y_A = dY_dAmplitude(f_rec, xdata)
  A_d = dA_d_distance(f_rec, xdata)

  n_receptors = length(f_rec$labels)
  f_d = lapply(seq_along(n_receptors),
               function(i)
               {
                 Y_dr = Y_A * A_d[[i]]
                 f_dr = sum_by_index(2 * Re(f_Y[, i]) * Re(Y_dr) +
                                       2 * Im(f_Y[, i]) * Im(Y_dr),
                                     xdata, 4)
                 return(f_dr)
               })
  return(f_d)
}

df_doffset <- function(f_rec, xdata)
{
  f_Y = df_dY(f_rec, xdata)
  Y_A = dY_dAmplitude(f_rec, xdata)
  A_off = dA_doffset(f_rec, xdata)

  n_receptors = length(f_rec$labels)
  f_off = vapply(seq_along(n_receptors),
                 function(i)
                 {
                   Y_offr = Y_A * A_off[[i]]
                   f_offr = sum(2 * Re(f_Y[, i]) * Re(Y_offr) +
                                  2 * Im(f_Y[, i]) * Im(Y_offr))
                   return(f_offr)
                 }, numeric(1))
  return(f_off)
}

df_d_logG <- function(f_rec, xdata)
{
  f_Y = df_dY(f_rec, xdata)
  Y_A = dY_dAmplitude(f_rec, xdata)
  A_logG = dA_dlogG(f_rec, xdata)

  n_receptors = length(f_rec$labels)
  f_logG = lapply(seq_along(n_receptors),
                  function(i)
                  {
                    Y_logGr = Y_A * A_logG[[i]]
                    f_logGr = sum_by_index(2 * Re(f_Y[, i]) * Re(Y_logGr) +
                                             2 * Im(f_Y[, i]) * Im(Y_logGr),
                                           xdata, 2)
                    return(f_logGr)
                  })
  return(f_logG)
}


df_dsource_pos_x <- function(f_rec, xdata)
{
  f_d_dist = df_d_distance(f_rec, xdata)
  dist_spx = d_dist_dsource_pos_x(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_sposx =
    Reduce(
      f = function(p, i)
        p + f_d_dist[[i]] * dist_spx[[i]],
      x = seq_len(n_receptors),
      init = 0
    )
  return(f_sposx)
}

df_dsource_pos_y <- function(f_rec, xdata)
{
  f_d_dist = df_d_distance(f_rec, xdata)
  dist_spy = d_dist_dsource_pos_y(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_sposy =
    Reduce(
      f = function(p, i)
        p + f_d_dist[[i]] * dist_spy[[i]],
      x = seq_len(n_receptors),
      init = 0
    )
  return(f_sposy)
}

df_dreceptor_pos_x <- function(f_rec, xdata)
{
  f_d_dist = df_d_distance(f_rec, xdata)
  dist_spx = d_dist_dsource_pos_x(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_rposx =
    vapply(seq_len(n_receptors),
      function(i) -sum( f_d_dist[[i]] * dist_spx[[i]]),
      numeric(0)
    )
  return(f_rposx)
}

df_dreceptor_pos_y <- function(f_rec, xdata)
{
  f_d_dist = df_d_distance(f_rec, xdata)
  dist_spy = d_dist_dsource_pos_y(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_rposy =
    vapply(seq_len(n_receptors),
           function(i) -sum( f_d_dist[[i]] * dist_spy[[i]]),
           numeric(0)
    )
  return(f_rposy)
}



d_f__d_beta <- function(f_rec, beta, xdata)
{



}
