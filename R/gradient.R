












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
    out[indexes[i, icol], ] = out[indexes[i, icol], ] + values[i, ]
  return(out)
}







d_dist__d_source_pos_x <- function(f_rec)
{
  dist = f_rec$distances
  x_s = f_rec$sources_pos_step$x
  x_r = f_rec$receptors_init$x
  n_receptors = length(f_rec$labels)
  dist_sx = lapply(seq_len(n_receptors), function(i)
    (x_s - x_r[i]) / dist[[i]])
  return(dist_sx)
}

d_dist__d_source_pos_y <- function(f_rec)
{
  dist = f_rec$distances
  y_s = f_rec$sources_pos_step$y
  y_r = f_rec$receptors_init$y
  n_receptors = length(f_rec$labels)
  dist_sy = lapply(seq_len(n_receptors), function(i)
    (y_s - y_r[i]) / dist[[i]])
  return(dist_sy)
}





d_A__d_distance <- function(f_rec, xdata)
{
  v = f_rec$velocity_of_sound
  w = 2 * pi * xdata[, 5]
  dist = lapply(f_rec$distances, function(d)
    d[xdata[, 4],])
  A = get_Amplitudes(f_rec, xdata)
  n_receptor = length(f_rec$labels)
  A_dist = lapply(1:n_receptor, function(i)
    - A[[i]] * (1 / dist[[i]] + 1i / v * w))
  return(A_dist)
}

d_A__d_offset <- function(f_rec, xdata)
{
  w = 2 * pi * xdata[, 5]
  A = get_Amplitudes(f_rec, xdata)
  n_receptor = length(f_rec$labels)
  A_off = lapply(1:n_receptor, function(i)
    - A[[i]] * 1i * w)
  return(A_off)
}

d_A__d_logG <- function(f_rec, xdata)
{
  A_logG = get_Amplitudes(f_rec, xdata)
  return(A_logG)
}

d_Y__d_Amplitude <- function(f_rec, xdata)
{
  X = get_source_signal(f_rec, xdata)
  return(X)
}



d_f__d_Y <- function(f_rec, xdata)
{
  Y = get_receptor_signal(f_rec, xdata)
  Yfit = get_predicted_signal(f_rec, xdata)

  return(2 * (Yfit - Y))
}




d_f__d_X <- function(f_rec, xdata)
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


d_f__d_distance <- function(f_rec, xdata)
{
  f_Y = d_f__d_Y(f_rec, xdata)
  Y_A = d_Y__d_Amplitude(f_rec, xdata)
  A_d = d_A__d_distance(f_rec, xdata)

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

d_f__d_offset <- function(f_rec, xdata)
{
  f_Y = d_f__d_Y(f_rec, xdata)
  Y_A = d_Y__d_Amplitude(f_rec, xdata)
  A_off = d_A__d_offset(f_rec, xdata)

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

d_f__d_logG <- function(f_rec, xdata)
{
  f_Y = d_f__d_Y(f_rec, xdata)
  Y_A = d_Y__d_Amplitude(f_rec, xdata)
  A_logG = d_A__d_logG(f_rec, xdata)

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


d_f__d_source_pos_x <- function(f_rec, xdata)
{
  f_dist = d_f__d_distance(f_rec, xdata)
  dist_spx = d_dist__d_source_pos_x(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_sposx =
    Reduce(
      f = function(p, i)
        p + f_dist[[i]] * dist_spx[[i]],
      x = seq_len(n_receptors),
      init = 0
    )
  return(f_sposx)
}

d_f__d_source_pos_y <- function(f_rec, xdata)
{
  f_dist = d_f__d_distance(f_rec, xdata)
  dist_spy = d_dist__d_source_pos_y(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_sposy =
    Reduce(
      f = function(p, i)
        p + f_dist[[i]] * dist_spy[[i]],
      x = seq_len(n_receptors),
      init = 0
    )
  return(f_sposy)
}

d_f__d_receptor_pos_x <- function(f_rec, xdata)
{
  f_dist = d_f__d_distance(f_rec, xdata)
  dist_spx = d_dist__d_source_pos_x(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_rposx =
    vapply(seq_len(n_receptors),
           function(i)
             - sum(f_dist[[i]] * dist_spx[[i]]),
           numeric(0))
  return(f_rposx)
}

d_f__d_receptor_pos_y <- function(f_rec, xdata)
{
  f_dist = d_f__d_distance(f_rec, xdata)
  dist_spy = d_dist__d_source_pos_y(f_rec = f_rec)
  n_receptors = length(f_rec$labels)


  f_rposy =
    vapply(seq_len(n_receptors),
           function(i)
             - sum(f_dist[[i]] * dist_spy[[i]]),
           numeric(0))
  return(f_rposy)
}

d_f__d_logG0 <- function(f_rec, xdata)
{
  f_logG = d_f__d_logG(f_rec, xdata)
  f_logG0 = vapply(f_logG, function(f_logGi)
    sum(f_logGi), numeric(1))
  return(f_logG0)

}

d_f__d_logGdiff <- function(f_rec, xdata)
{
  f_logG = d_f__d_logG(f_rec, xdata)
  f_logG0 = d_f__d_logG0(f_rec, xdata)
  n_receptors = length(f_rec$labels)
  n_frames = nrow(f_logG[[1]])
  f_logGdiff = vapply(seq_len(n_receptors), function(i)
    colSums(f_logG0[i] - cumsum(f_logG[[i]])), numeric(n_frames))
  return(f_logGdiff)
}

d_f__d_beta <- function(f_rec, beta, xdata)
{

}
