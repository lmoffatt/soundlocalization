








sum_by_index <- function(values, indexes, icol, num_indexes)
{
  stopifnot(nrow(values) == nrow(indexes))
  ncols = ncol(values)
  out = c()
  if (is.complex(values))
    out = matrix(0 + 0i, nrow = num_indexes,
                 ncol = ncols)
  else
    out = matrix(0, nrow = num_indexes, ncol = ncols)
  for (i in seq_len(nrow(values)))
    out[indexes[i, icol],] = out[indexes[i, icol],] + values[i,]
  return(out)
}

sum_array_by_index <-
  function(array,
           indexes,
           icol,
           dimension,
           num_indexes)
  {
    stopifnot(nrow(array) == nrow(indexes))
    mydim = dim(array)
    mydim[dimension] = num_indexes

    out = c()
    if (is.complex(array))
      out = array(0 + 0i, mydim)
    else
      out = array(0, mydim)
    for (i in seq_len(nrow(array)))
      out[indexes[i, icol],] = out[indexes[i, icol],] + array[i,]
    return(out)
  }


d_dist__d_source_pos_x <- function(dist, x_s, x_r)
{
  n_receptors = dim(dist)[1]
  n_sim_sources = dim(dist)[2]
  n_frame_freq = dim(dist)[3]
  stopifnot(nrow(x_s) == n_sim_sources)
  stopifnot(ncol(x_s) == n_frame_freq)
  stopifnot(length(x_r) == n_receptors)


  dist_sx_n =
    vapply(seq_len(n_frame_freq), function(n)
      vapply(seq_len(n_sim_sources), function(j)
        return ((x_s[, j] - x_r) / dist[, j, n]),
        numeric(n_receptors)),
      matrix(0, nrow = n_receptors, ncol = n_sim_sources))

  dist_sx = aperm(dist_sx_n, c(3, 1, 2))
  return(dist_sx)

}





d_A__d_distance <- function(xdata, A, dist, velocity_of_sound)
{
  n_receptor = dim(dist)[1]
  n_sources_sim = dim(dist)[2]
  n_frame_freq = dim(dist)[3]
  n = dim(A)[3]
  stopifnot(n == ncol(xdata))
  w = rep(2 * pi * xdata[5, ], each = n_receptor * n_sources_sim)
  A_dist_n = -A * (1 / dist[, , xdata[3, ]] + 1i / velocity_of_sound * w)
  # A_dist = aperm(A_dist_n, c(3,1,2))
  return(A_dist_n)
}

d_A__d_offset <- function(xdata, A)
{
  n_receptor = dim(A)[1]
  n_sources_sim = dim(A)[2]
  n = dim(A)[3]
  w = rep(2 * pi * xdata[, 5], each = n_receptor * n_sources_sim)
  stopifnot(n == nrow(xdata))
  A_off_n = -A * 1i * w
  #  A_off = aperm(A_off_n, c(3,1,2))
  return(A_off_n)
}

d_A__d_logG <- function(A)
{
  A_logG_n = A
  # A_logG = aperm(A_logG_n, c(3,1,2))
  return(A_logG_n)
}




MT2v <- function(Matrix, Tensor, vector)
{
  n1 = dim(Matrix)[1]
  n2 = dim(Matrix)[2]
  stopifnot(n2 == dim(Tensor)[1])
  n3 = dim(Tensor)[2]
  stopifnot(n3 == dim(vector)[1])
  out = array(dim = c(n1, n2, n3))
  for (i in seq_len(n1))
    for (j in seq_len(n2))
      for (k in seq_len(n3))
        out[i, j, k] = Matrix[i, j] * Tensor[j, k] * vector[k]
  return (out)

}

MT2 <- function(Matrix, Tensor)
{
  n1 = dim(Matrix)[1]
  n2 = dim(Matrix)[2]
  stopifnot(n2 == dim(Tensor)[1])
  n3 = dim(Tensor)[2]
  out = array(dim = c(n1, n2, n3))
  for (i in seq_len(n1))
    for (j in seq_len(n2))
      for (k in seq_len(n3))
        out[i, j, k] = Matrix[i, j] * Tensor[j, k]
  return(out)

}
T2M <- function(Tensor, Matrix)
{
  n1 = dim(Tensor)[1]
  n2 = dim(Tensor)[2]
  stopifnot(n2 == dim(Matrix)[1])
  n3 = dim(Matrix)[2]
  out = array(0 + 0i, dim = c(n1, n2, n3))
  for (i in seq_len(n1))
    for (j in seq_len(n2))
      for (k in seq_len(n3))
        out[i, j, k] = Tensor[i, j] * Matrix[j, k]
  return(out)

}


T3v <- function(Tensor, vector)
{
  n1 = dim(Tensor)[1]
  n2 = dim(Tensor)[2]
  n3 = dim(Tensor)[3]
  stopifnot(n3 == dim(vector)[1])
  out = array(0 + 0i, dim = c(n1, n2))
  for (i in seq_len(n1))
    for (j in seq_len(n2))
      out[i, j] = Tensor[i, j,] %*% vector
  return(out)

}


DT2v <- function(diagonal_tensor2, vector)
{
  n1 = dim(diagonal_tensor2)[1]
  n2 = dim(diagonal_tensor2)[2]
  stopifnot(n2 == dim(vector)[1])
  out = array(dim = c(n1, n1, n2))
  for (i in seq_len(n1))
    out[i, i,] = diagonal_tensor2[i,] %*% vector
  return(out)


}




MT3v <- function(Matrix, Tensor3, vector)
{
  TV= T3v(Tensor3, vector)
  return(MT2(Matrix, TV))
}


I_plus_lambda_diag <- function(Tensor, lambda)
{
  n1 = dim(Tensor)[1]
  n2 = dim(Tensor)[2]
  n3 = dim(Tensor)[3]
  out = Tensor
  for (i in seq_len(n1))
    for (j in seq_len(n2))
      out[i, j, i] = out[i, j, i] * (1 + lambda)
  return(out)
}




d_X__d_A <- function(dA, A, X, AAinv, Yfit, Y0)
{
  n_dbeta = dim(dA)[4]
  n_frames = dim(A)[3]

  n_receptors = dim(A)[1]
  n_sources = dim(A)[2]
  n_dA = dim(dA)[4]
  return(vapply(seq_len(n_dbeta), function(j)
    return(
      vapply(seq_len(n_frames), function(i)
      {
        M1 = -AAinv[, , i]
        T1 = T2M(Conj(t(dA[, , i, j])), A[, , i]) + MT2(Conj(t(A[, , i])), dA[, , i, j])
        v1 = AAinv[, , i] %*% Conj(t(A[, , i])) %*% Y0[, i]

        M2 = AAinv[, , i]
        T2 = Conj(t(dA[, , i, j]))
        v2 = Y0[, i]

        return(MT3v(M1, T1, v1) + MT2v(M2, T2, v2))
      }, array(
        0 + 0i, dim = c(n_sources, n_receptors, n_sources)
      ))
      , array(0 + 0i,
              dim = c(
                n_sources, n_receptors, n_sources, n_frames
              ))
    )))
}



d_Y__d_A <- function(dA, A, X, AAinv, Yfit, Y0)
{
  n_samples = dim(A)[3]

  n_receptors = dim(A)[1]
  n_sources = dim(A)[2]

  Y_An = vapply(seq_len(n_samples), function(i)
  {
    M = A[, , i] %*% AAinv[, , i]
    T1 = T2M(Conj(t(dA[, , i])), A[, , i]) + MT2(Conj(t(A[, , i])), dA[, , i])
    v1 = AAinv[, , i] %*% Conj(t(A[, , i])) %*% Y0[, i]

    T2 = Conj(t(dA[, , i]))
    v2 = Y0[, i]

    DT3 = dA[, , i]
    v3 = X[, i]
    First =MT3v(-M, T1, v1)
    Second = MT2v(M, T2, v2)
    Third = DT2v(DT3, v3)
    return(First + Second + Third)
  }, array(0 + 0i, dim = c(n_receptors, n_receptors, n_sources)))

  Y_A = aperm(Y_An, perm = c(1, 4, 2, 3))
  return (Y_A)
}

d_f__d_Y <- function(Y0, Yfit)
{
  return(2 * (Yfit - Y0))
}




d_f__d_A <- function(dA, A, X, AAinv, Yfit, Y0)
{
  n_frames = dim(A)[3]

  n_receptors = dim(A)[1]
  n_sources = dim(A)[2]
  return(vapply(seq_len(n_frames), function(i)
  {
    M = A %*% AAinv[, , i]
    T1 = T2M(Conj(t(dA[, , i, j])), A[, , i]) + MT2(Conj(t(A[, , i])), dA[, , i])
    v1 = AAinv[, , i] %*% Conj(t(A[, , i])) %*% Y0[, i]

    T2 = Conj(t(dA[, , i]))
    v2 = Y0[, i]

    DT3 = dA[, , i]
    v3 = X[, i]

    dY = MT3v(-M, T1, v1) + MT2v(M, T2, v2) + DT2v(DT3, v3)

    df = colSums(Re(Conj(Yfit[, i] - Y0[, i]) * dY))
    return(df)
  }, array(0 , dim = c(
    n_receptors, n_sources
  ))))

}



d_Y__d_distance <-
  function(xdata,
           A,
           dist,
           velocity_of_sound,
           X,
           AAinv,
           Yfit,
           Y0)
  {
    A_d = d_A__d_distance(xdata, A, dist, velocity_of_sound)

    Y_d = d_Y__d_A(A_d, A, X, AAinv, Yfit, Y0)

    return(Y_d)
  }

d_Y__d_offset <- function(xdata, A, X, AAinv, Yfit, Y0, n_receptors)
{
  A_off = d_A__d_offset(xdata, A)

  Y_offs = d_Y__d_A(A_off, A, X, AAinv, Yfit, Y0)

  Y_off = apply(Y_offs, MARGIN = c(1, 2, 3), sum)
  return(matrix(Y_off, ncol = n_receptors))

}

d_Y__d_logG <- function(xdata, A, X, AAinv, Yfit, Y0, n_receptors)
{
  A_logG = d_A__d_logG(A)

  Y_logGns = d_Y__d_A(A_logG, A, X, AAinv, Yfit, Y0)

  Y_logGn = apply(Y_logGns, MARGIN = c(1, 2, 3), sum)


  return(matrix(Y_logGn, ncol = n_receptors))
}

d_Y__d_logGdiff <- function(xdata, Y_logG, n_receptors, n_frames)
{
  n_samples = nrow(Y_logG)

  nframe = matrix(rep(2:n_frames, each = n_samples * n_receptors),
                  c(n_samples, n_receptors * (n_frames - 1)))

  iframe = rep(xdata[2,], each = n_receptors)

  return(ifelse(iframe > (nframe - 1) , Y_logG, 0 + 0i))

}


d_f__d_distance <-
  function(xdata,
           A,
           dist,
           velocity_of_sound,
           X,
           AAinv,
           Yfit,
           Y0)
  {
    A_d = d_A__d_distance(xdata, A, dist, velocity_of_sound)

    f_d = d_f__d_A(A_d, A, X, AAinv, Yfit, Y0)

    return(f_d)
  }

d_f__d_offset <- function(xdata, A, X, AAinv, Yfit, Y0)
{
  A_off = d_A__d_offset(xdata, A)

  f_off = d_f__d_A(A_off, A, X, AAinv, Yfit, Y0)

  return(rowSums(f_off))
}

d_Y__d_source_pos_x <-
  function(Y_dist,
           dist_spx,
           xdata,
           n_receptors,
           n_simul_sources,
           n_frame_freq_sources)
  {
    Y_sposx_r = Y_dist * rep(dist_spx[xdata[, 3], , ], each = n_receptors)
    Y_sposxn = apply(Y_sposx_r, MARGIN = c(1, 2, 4), sum)

    Y_sposx = matrix(Y_sposxn, ncol = n_simul_sources)
    n_samples = nrow(Y_sposx)

    nsource = matrix(rep(1:n_frame_freq_sources, each = n_samples * n_simul_sources),
                    nrow=n_samples, ncol=n_frame_freq_sources*n_simul_sources)

    isource = rep(xdata[2,], each = n_receptors)

    return(ifelse(nsource == nsource  , Y_sposx, 0 + 0i))
  }




d_f__d_source_pos_x <-
  function(f_dist,
           xdata,
           dist_spx,
           n_frames)
  {
    f_sposx = sum_by_index(values = colSums(f_dist * dist_spx[, , xdata[, 3]]),
                           indexes = xdata,
                           2,
                           n_frames)

    return(f_sposx)
  }


d_Y__d_receptor_pos_x <-
  function(Y_dist,
           dist_spx,
           xdata,
           n_receptors)
  {
    Y_rposx_s = -Y_dist * rep(dist_spx[xdata[, 3], , ], each = n_receptors)
    Y_rposx = apply(Y_rposx_s, MARGIN = c(1, 2, 3), sum)

    return(matrix(Y_rposx, ncol = n_receptors))
  }




d_f__d_receptor_pos_x <-
  function(f_dist,
           dist_spx,
           xdata)
  {
    f_rposx = rowSums(-f_dist * dist_spx[, , xdata[, 3]])
    return(f_rposx)
  }





d_f__d_logG <- function(xdata, A, X, AAinv, Yfit, Y0)
{
  A_logG = d_A__d_logG(A)

  f_logG = d_f__d_A(A_logG, A, X, AAinv, Yfit, Y0)

  return(f_logG)
}


d_f__d_logG0 <- function(f_logG)
{
  f_logG0 = rowSums(f_logG)
  return(f_logG0)

}


d_f__d_logGdiff <-
  function(xdata,
           f_logGi,
           f_logG0,
           n_frames,
           n_receptors)
  {
    f_logG = sum_by_index(
      values = apply(f_logGi, MARGIN = c(1, 3), sum),
      indexes = xdata,
      icol = 3,
      num_indexes = n_frames
    )

    f_logGdiff = vapply(seq_len(n_receptors), function(i)
      f_logG0[i] - cumsum(f_logG[i, ]), numeric(n_frames))
    return(f_logGdiff)
  }

d_f__d_beta <- function(f_rec, model, beta, xdata)
{
  ydata = get_receptor_signal(model, xdata)
  distances = get_source_receptor_distance(f_rec, model)
  A = get_Amplitudes(f_rec, model, distances, xdata)

  AAinv = get_inverse_AA(A)

  X = get_X(AAinv, A, ydata)

  Yfit = get_predicted_Y(A, X)


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
             n_sources,
             n_sources)

  total = sum(counts)
  cumtotal = cumsum(counts)


  cum_receptors = cumsum(c(n_receptors - 2, n_receptors - 1,
                           n_receptors - 1, n_receptors - 1))


  f_beta = numeric(total)

  x_r = model$receptors$x
  x_s = model$sources$x
  y_r = model$receptors$y
  y_s = model$sources$y
  velocity_of_sound = f_rec$velocity_of_sound

  f_dist = d_f__d_distance(xdata, A, dist, velocity_of_sound, X, AAinv, Yfit, ydata)
  dist_spx = d_dist__d_source_pos_x(dist, x_s, x_r)
  dist_spy = d_dist__d_source_pos_x(dist, y_s, y_r)


  f_rposx = d_f__d_receptor_pos_x(f_dist,
                                  dist_spx,
                                  xdata)
  f_rposy = d_f__d_receptor_pos_x(f_dist,
                                  dist_spy,
                                  xdata)

  f_logG = d_f__d_logG(xdata, A, X, AAinv, Yfit, ydata)
  f_logG0 = d_f__d_logG0(f_logG)
  f_offset = d_f__d_offset(xdata, A, X, AAinv, Yfit, ydata)

  f_beta[1:cumtotal[1]] =
    c(f_rposx[3:n_receptors],
      f_rposy[2:n_receptors],
      f_logG0[2:n_receptors],
      f_offset[2:n_receptors])






  if (model$include_receptor_frame_Gain)
  {
    f_logGdiff = d_f__d_logGdiff(xdata,
                                 f_logG,
                                 f_logG0,
                                 n_frames,
                                 n_receptors)
    f_beta[(cumtotal[1] + 1):cumtotal[2]] = f_logGdiff[2:n_frames,]
  }
  f_spx = d_f__d_source_pos_x(f_dist,
                              xdata,
                              dist_spx,
                              n_frames)
  f_spy = d_f__d_source_pos_x(f_dist,
                              xdata,
                              dist_spy,
                              n_frames)

  f_beta[(cumtotal[2] + 1):cumtotal[3]] = c(f_spx[seq_along(f_spx)], f_spy[seq_along(f_spy)])


  return(f_beta)

}

d_Y__d_beta    <- function(f_rec, model, beta, xdata)
{
  ydata = get_receptor_signal(model, xdata)
  distances = get_source_receptor_distance(f_rec, model)
  A = get_Amplitudes(f_rec, model, distances, xdata)

  AAinv = get_inverse_AA(A)
  X = get_X(AAinv, A, ydata)
  Yfit = get_predicted_Y(A, X)
  n_y = length(Yfit)


  n_receptors = get_number_of_receptors(f_rec)
  n_frames = get_number_of_frames(f_rec)
  n_total_sources = get_number_of_sources(model)
  count_receptors = n_receptors * 5
  count_frames_gain_par = 0
  include_frame_gain = 'frames_gain' %in% names(model)
  if (include_frame_gain)
    count_frames_gain_par = (nrow(model$frames_gain$logG_diff) - 1) * n_receptors
  n_sources = get_number_of_sources(model)
  counts = c(count_receptors,
             count_frames_gain_par,
             n_sources * 3)

  total = sum(counts)
  cumtotal = cumsum(counts)

  x_r = model$receptors$x
  x_s = model$sources$x
  y_r = model$receptors$y
  y_s = model$sources$y
  z_r = model$receptors$z
  z_s = model$sources$z

  velocity_of_sound = f_rec$velocity_of_sound
  dist = get_source_receptor_distance(f_rec,model)

  Y_dist = d_Y__d_distance(xdata, A, dist, velocity_of_sound, X, AAinv, Yfit, ydata)
  dist_spx = d_dist__d_source_pos_x(dist, x_s, x_r)
  dist_spy = d_dist__d_source_pos_x(dist, y_s, y_r)
  dist_spz = d_dist__d_source_pos_x(dist, z_s, z_r)


  Y_rposx = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spx,
                                  xdata,
                                  n_receptors)
  Y_rposy = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spy,
                                  xdata,
                                  n_receptors)
  Y_rposz = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spz,
                                  xdata,
                                  n_receptors)


  Y_logG = d_Y__d_logG(xdata, A, X, AAinv, Yfit, ydata, n_receptors)
  Y_logGdiff = d_Y__d_logGdiff(xdata = xdata, Y_logG , n_receptors , n_frames)
  Y_offset = d_Y__d_offset(xdata, A, X, AAinv, Yfit, ydata, n_receptors)






  Y_spx = d_Y__d_source_pos_x(Y_dist,
                              dist_spx,
                              xdata,
                              n_receptors,
                              n_simul_sources,
                              n_frame_freq_sources)
  Y_spy = d_Y__d_source_pos_x(Y_dist,
                              dist_spy,
                              xdata,
                              n_receptors,
                              n_total_sources)
  Y_spz = d_Y__d_source_pos_x(Y_dist,
                              dist_spz,
                              xdata,
                              n_receptors,
                              n_total_sources)

  i_rec = 1:cumtotal[1]

  i_logGdiff = (cumtotal[1] + 1):cumtotal[2]

  i_sources = (cumtotal[2] + 1):cumtotal[3]


  f_beta2[i_rec, i_rec] = t(Conj(Y_rec)) %*% Y_rec





  f_beta2[i_r_px, i_r_px] = t(Conj(Y_rposx)) %*% Y_rposx
  f_beta2[i_r_px, i_r_px] = t(Conj(Y_rposx)) %*% Y_rposx
  f_beta2[i_r_px, i_r_px] = t(Conj(Y_rposx)) %*% Y_rposx



  f_beta[(cumtotal[1] + 1):cumtotal[2]] = f_logGdiff[2:n_frames,]

  f_beta[(cumtotal[2] + 1):cumtotal[3]] = c(f_spx[seq_along(f_spx)], f_spy[seq_along(f_spy)])

  Y_rec = cbind(Y_rposx, Y_rposy, Y_rposz, Y_logG, Y_offset)

  return(f_beta)

}


d2_f__d_beta2    <- function(f_rec, model, beta, xdata)
{
  ydata = get_receptor_signal(model, xdata)
  distances = get_source_receptor_distance(f_rec, model)
  A = get_Amplitudes(f_rec, model, distances, xdata)

  AAinv = get_inverse_AA(A)

  X = get_X(AAinv, A, ydata)

  Yfit = get_predicted_Y(A, X)

  n_y = length(Yfit)


  n_receptors = get_number_of_receptors(f_rec)
  n_frames = get_number_of_frames(f_rec)

  n_simul_sources = get_number_of_simultaneous_sources(model)
  n_sources_frames = get_number_of_source_frames(model)


  count_receptors = n_receptors * 5
  count_frames_gain_par = 0
  include_frame_gain = 'frames_gain' %in% names(model)
  if (include_frame_gain)
    count_frames_gain_par = (nrow(model$frames_gain$logG_diff) - 1) * n_receptors
  n_sources = get_number_of_sources(model)
  counts = c(count_receptors,
             count_frames_gain_par,
             n_sources * 3)

  total = sum(counts)
  cumtotal = cumsum(counts)




  f_beta2 = matrix(0, nrow = total, ncol = total)

  x_r = model$receptors$x
  x_s = model$sources$x
  y_r = model$receptors$y
  y_s = model$sources$y
  z_r = model$receptors$z
  z_s = model$sources$z

  velocity_of_sound = f_rec$velocity_of_sound

  Y_dist = d_Y__d_distance(xdata, A, dist, velocity_of_sound, X, AAinv, Yfit, ydata)
  dist_spx = d_dist__d_source_pos_x(dist, x_s, x_r)
  dist_spy = d_dist__d_source_pos_x(dist, y_s, y_r)
  dist_spz = d_dist__d_source_pos_x(dist, z_s, z_r)


  Y_rposx = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spx,
                                  xdata,
                                  n_receptors)
  Y_rposy = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spy,
                                  xdata,
                                  n_receptors)
  Y_rposz = d_Y__d_receptor_pos_x(Y_dist,
                                  dist_spz,
                                  xdata,
                                  n_receptors)


  Y_logG = d_Y__d_logG(xdata, A, X, AAinv, Yfit, ydata, n_receptors)
  Y_logGdiff = d_Y__d_logGdiff(xdata = xdata, Y_logG , n_receptors , n_frames)
  Y_offset = d_Y__d_offset(xdata, A, X, AAinv, Yfit, ydata, n_receptors)



  Y_rec = cbind(Y_rposx, Y_rposy, Y_rposz, Y_logG, Y_offset)



  Y_spx = d_Y__d_source_pos_x(Y_dist,
                              dist_spx,
                              xdata,
                              n_receptors,
                              n_simul_sources,
                              n_sources_frames)
  Y_spx = d_Y__d_source_pos_x(Y_dist,
                              dist_spx,
                              xdata,
                              n_receptors,
                              n_simul_sources,
                              n_sources_frames)
  Y_spx = d_Y__d_source_pos_x(Y_dist,
                              dist_spx,
                              xdata,
                              n_receptors,
                              n_simul_sources,
                              n_sources_frames)


  Y_beta = cbind(Y_rposx,Y_rposy,Y_rposz, Y_logG, Y_offset,Y_logGdiff,Y_rposx,Y_rposy,Y_rposz)

  return(Y_beta)

}
