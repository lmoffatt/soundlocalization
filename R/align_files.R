


#' applies a hanning window
#'
#' @param equalized_signal  the signal should be already in matrix form
#'
#' @return it modifies the signal in place
#' @export
#'
#' @examples
apply_hanning_window <- function(equalized_signal)
{
  is = equalized_signal$i_start
  ie = equalized_signal$i_end
  nsamples = ie[1] - is[1] + 1

  hn = e1071::hanning.window(nsamples)
  equalized_signal$signal = equalized_signal$signal * hn
  equalized_signal
}


#' Title obtains the fourier transform of all the signals at once
#'
#' @param windowed_signal the hanning window should be already applied
#'
#' @return inserts a field with the fft
#' @export
#'
#' @examples
calculate_ffts <- function(windowed_signal)
{
  ncols = ncol(windowed_signal$signal)
  nrows = nrow(windowed_signal$signal)
  stopifnot("fft on a wide matrix" = ncols < nrows)
  windowed_signal$fft = mvfft(windowed_signal$signal)
  return (windowed_signal)
}





#' Title calculates all the cross correlation pairs
#'
#' @param ffted_signal the fft should be already calculated
#'
#' @return inserts the cross correlation function
#' @export
#'
#' @examples
calculate_list_of_cross_correlations <- function(ffted_signal)
{
  ncols = ncol(ffted_signal$fft)
  nrows = nrow(ffted_signal$fft)
  ffts = ffted_signal$fft

  ccf = lapply(1:(ncols - 1),
               function(i)
                 lapply(i:(ncols - 1),
                        function(j)
                          Re(fft(
                            Conj(ffts[, i]) * ffts[, j + 1], inverse = T
                          ))))


  ffted_signal$ccf <- ccf
  ffted_signal
}


#' Calculates the optimal lag combination for a set of ccfs pairs
#'
#' @param ccf_list should include double indexed list of cross correlation functions
#'
#' @return it inserts the global and local optimal lags
#' @export
#' @examples
optimize_delays_from_ccf <- function(maxs_ccf)
{
  nsamples = nrow(maxs_ccf$signal)
  n = ncol(maxs_ccf$signal)
  ccf = maxs_ccf$ccf
  ccf_for_delay <- function(x) {
    lags = c(0, x)

    sum = 0
    for (i in 1:(n - 1))
    {
      for (j in (i + 1):n)
      {
        lag = lags[j] - lags[i]
        l0 = floor(lag)
        l1 = l0 + 1
        r = lag - l0
        s = 1 - r
        i0 = lag_to_index(l0,nsamples)
        i1 = lag_to_index(l1,nsamples)
        sum = sum + ccf[[i]][[j - i]][i0] * s + ccf[[i]][[j - i]][i1] *
          r
      }
    }
    return (-sum)
  }
  resp = list()
  initial_delays = colMeans(maxs_ccf$max_lag_by_row)
  resp[[1]] <- nloptr::neldermead	(
    initial_delays,
    ccf_for_delay,
    lower = numeric(n - 1) - nsamples / 4,
    upper = numeric(n - 1) + nsamples / 4
  )
  for (i in 1:n)
  {
    initial_delays = maxs_ccf$max_lag_by_row[i,]
    resp[[i + 1]] <- nloptr::neldermead	(
      initial_delays,
      ccf_for_delay,
      lower = numeric(n - 1) - nsamples / 4,
      upper = numeric(n - 1) + nsamples / 4
    )
  }
  optimals = list(
    par = t(vapply(resp, function(r)
      r$par, numeric(n - 1))),
    value = vapply(resp, function(r)
      - r$value, 0),
    iter = vapply(resp, function(r)
      r$iter, 0),
    convergence = vapply(resp, function(r)
      r$convergence, 0),
    message = vapply(resp, function(r)
      r$message, 'text')
  )
  n_opt = which.max(optimals$value)
  opt = resp[[n_opt]]
  list(local_maximal = optimals, maximal = opt)
}


#' calculates the indexes so all the waves are of the same size
#'
#' @param list_of_waves a list of class Wave objects
#'
#' @return a list with the start end indexes vectors
#' @export
#'
#' @examples
#'files<-list(
#' system.file("extdata", "Grabacion_2_andrea.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_joaco.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_cande.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_3_luciano.wav", package = "soundlocalization")
#' )
#' waveList=loadwavfiles(files)
#' indexes=indexes_for_same_length(waveList)
#'
#'
indexes_for_same_length <- function(signals)
{
  nsamples = signals$meta_data[["nsamples"]]
  min_n = min(nsamples)
  diff_n = nsamples - min_n
  i_start = floor(diff_n / 2) + 1
  i_end = i_start + min_n - 1
  list(i_start = i_start, i_end = i_end)
}


#' applies an interval to a set of signals
#'
#' @param signals a signals list object
#' @param indexes an interval index object
#'
#' @return it appends signal interval to the set
#' @export
#'
#' @examples
#'files<-list(
#' system.file("extdata", "Grabacion_2_andrea.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_joaco.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_cande.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_3_luciano.wav", package = "soundlocalization")
#' )
#' waveList=loadwavfiles(files)
#' indexes=indexes_for_same_length(waveList)
#' aligned=apply_indexes(waveList,indexes)
#'
apply_indexes <- function(signals, indexes)
{
  is = indexes$i_start
  ie = indexes$i_end
  s = signals$signal
  nsamples = ie[1] - is[1] + 1
  signal = vapply(1:length(s), function(i)
    s[[i]][is[i]:ie[i]], numeric(nsamples))

  list(
    i_start = is,
    i_end = ie,
    fs = signals$meta_data$fs[1],
    signal = signal
  )
}


#' convert an index of a ccf to the corresponding lag
#'
#' @param i index of a ccf
#' @param nsmaples  number of samples of the ccf
#'
#' @return the lag that correspond to the index i
#' @export
#'
#' @examples
index_to_lag <- function(i, nsamples)
{
  ifelse(i < nsamples / 2, i - 1, i - nsamples - 1)
}

lag_to_index <- function(lag, nsamples)
{
  ifelse(lag < 0, nsamples + lag + 1, lag + 1)
}


#' finds the lag that correspond to the maximal cross correlation
#'
#' @param ccf a list that contains the pairs of cross correlations
#'
#' @return inserts that optimal lags of each pair of microphones
#' @export
#'
#' @examples
max_lag_ccf <- function(ccf)
{
  nsamples <- nrow(ccf$signal)
  max_lag <- sapply(ccf$ccf,
                    function (ccfi)
                      sapply(ccfi,
                             function (ccfij)
                               index_to_lag(which.max(ccfij), nsamples))
                    )


  ccf$lag <- max_lag
  ccf
}

#' calculates optimal lags for each microphone
#'
#' @param max_ccf list containing the paired lags
#'
#' @return it appends the optimal lags for each microphone pair
#' @export
#'
#' @examples
calc_lag_ccf <- function(max_ccf)
{
  max_lag = max_ccf$lag
  n = ncol(max_ccf$signal)
  d = matrix(0, nrow = n, ncol = n - 1)
  for (j in 2:n)
    d[1, j - 1] = max_lag[[1]][[j - 1]]
  for (i in 2:n)
    for (j in 2:n)
      if (i < j)
        d[i, j - 1] = max_lag[[1]][[i - 1]] + max_lag[[i]][[j - i]]
  else if (i == j)
    d[i, j - 1] = max_lag[[1]][[i - 1]]
  else
    d[i, j - 1] = max_lag[[1]][[i - 1]] - max_lag[[j]][[i - j]]
  max_ccf$max_lag_by_row = d
  max_ccf
}







align_files <- function(list_of_waves)
{
  li <- obtain_i_start_i_end(list_of_waves)
  ccf = calculate_list_of_cross_correlations(li$list_of_waves, li$i_start_n, i_end_n)


  neldermead(
    initial_delays,
    ccf_for_delay,
    lower = numeric(3) - wl / 4 / fs,
    upper = numeric(3) + wl / 4 / fs
  )

}
