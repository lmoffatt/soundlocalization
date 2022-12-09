
shift<-function(v,n)
{
  if (n<0)
    v=v[length(v):1-n]
  else
    v=v[1:length(v)+n]
  v[is.na(v)]<-0
  return (v)
}



build_ABCij_to_ABi_transform_matrix<-function(nlags)
{
  nlags2=nlags*nlags
  # lags=-nlags:(nlags-1)
  # ABij--> ABi
  # ACij--> ACj
  # BCij --> BCj-i

  # vector ABij_vector=matrix(nrow=nlags*nlags, ncol = 1)
  #ABi_vector=matrix(nrow = nlags, ncol=1)
  #ABij_to_ABi_matrix= matrix(nrow=nlags,ncol=nlags*nlags)
  lags=-(nlags/2):(nlags/2-1)
  X_lags=rep(lags,nlags)
  Y_lags=Reduce(function(p,l) c(p,rep(l,nlags)),lags,c())
  X_Y_lags=1:nlags2
  ACij_to_ACj_matrix=matrix(
    rep(diag(nlags),nlags),
    nrow=nlags,ncol=nlags2)
  ABij_to_ABi_matrix= t(ACij_to_ACj_matrix)

  BCij_to_BCj_i_matrix=
    matrix(
    vapply(-(nlags/2):(nlags/2-1),
           function (i) shift(diag(nlags),-i*nlags),diag(nlags)
    ),
    ncol=nlags,nrow = nlags2)

  A=matrix(0,nrow = 3*nlags, ncol=3*nlags2)
  A[1:nlags,seq(from=1, by=3, to=nlags2*3)]=ABij_to_ABi_matrix
  A[1:nlags+nlags,seq(from=2, by=3, to=nlags2*3)]=ACij_to_ACj_matrix
  A[1:nlags+nlags*2,seq(from=3, by=3, to=nlags2*3)]=BCij_to_BCj_i_matrix

  return(list(A=A, x_lags=rep(X_lags,3), y_lags=rep(Y_lags,3),xy_group=rep(X_Y_lags,each=3)))
}











locatiogram <- function(rec,
                        frames,
                        receptors,
                        axis_length = 13,
                        n_points = 100,
                        sound_velocity = 343)
{
  delay_source_receptor<-function(single_source,all_receptors)
  {
    n_receptors = length(all_receptors)
    delays <-
      vapply(1:n_receptors, function(i)
        distance(source, all_receptors[[i]]),
        numeric(1)) / sound_velocity
    return (delays)

  }

  gcc_amplitude_source_receptor_pairs<-function(single_source, all_receptors)
  {
    delays=delay_source_receptor(single_source, all_receptors)
  }

  source_position_to_gcc_pairs <-
    function(single_source,
             all_receptors,
             single_frame)
    {
      n_receptors = nrow(receptors)
      t_dist <-
        vapply(1:n_receptors, function(i)
          distance(source, all_receptors[[i]]),
          numeric(1)) / sound_velocity
      Reduce(function(s, i)
        s + Reduce(function(s, j)
        {
          nsamples = length(frame$ccf[[1]][[i]][[j - i]])
          lag <- ((t_dist[i] - t_dist[j]) * fs) %% nsamples
          s + frame$ccf[[1]][[i]][[j - i]][lag+1]
        }, (i+1):n_receptors, 0),
        1:(n_receptors - 1),
        0)

    }



  generate_coordinates <- function(axis_length, n)
  {
    merge(data.frame(ix = (-2*n+1):(2*n-1),
                     x = c(-n / (1:n), (n - 1):(-n + 1) / n, n / (n:1)) * axis_length),
          data.frame(iy = 0:(2*n-1),
                     y = c((0:(n - 1))/n, n / (n:1)) * axis_length),
          by = NULL)
  }

  sources = generate_coordinates(axis_length = axis_length, n = n_points)
  sources$ccf_density = vapply(1:nrow(sources), function(i)
    source_position_to_ccf_sum(
      source = sources[i,],
      receptors = receptors,
      frame = single_frames,
      fs = fs,
      sound_velocity = sound_velocity
    ), numeric(1))
  sources
}
