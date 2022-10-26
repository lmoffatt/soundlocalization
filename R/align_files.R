
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
  stopifnot(nsamples>2)

#  hn = e1071::hanning.window(nsamples)
  hn = seewave::ftwindow(nsamples,wn="hanning")
  equalized_signal$signal_for_fft = lapply(equalized_signal$signal, function(x)
    x * hn)
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
  if (!"signal_for_fft" %in% colnames(windowed_signal))
    windowed_signal = apply_hanning_window(windowed_signal)

  windowed_signal$fft = lapply(windowed_signal$signal_for_fft, fft)
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
  if (! "fft" %in% names(ffted_signal))
    ffted_signal<-calculate_ffts(ffted_signal)
  n = length(ffted_signal$fft)

  ffted_signal$ccf = lapply(1:(n - 1),
                            function(i)
                              lapply(i:(n - 1),
                                     function(j)
                                       Re(
                                         fft(Conj(ffted_signal$fft[[i]]) *
                                               ffted_signal$fft[[j + 1]], inverse = T)
                                       )))

  ffted_signal
}


#' Calculates the optimal lag combination for a set of ccfs pairs
#'
#' @param ccf_list should include double indexed list of cross correlation functions
#'
#' @return it inserts the global and local optimal lags
#' @export
#' @examples
optimize_delays_from_ccf <- function(data, maxs_ccf, max_lag_by_row)
{
  nsamples = length(data$ccf[[1]][[1]])
  n = nrow(max_lag_by_row)
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
        i0 = lag_to_index(l0, nsamples)
        i1 = lag_to_index(l1, nsamples)
        sum = sum + data$ccf[[i]][[j - i]][i0] * s + data$ccf[[i]][[j - i]][i1] *
          r
      }
    }
    return (-sum)
  }
  resp = list()
  initial_delays = colMeans(max_lag_by_row)
  resp[[1]] <- nloptr::neldermead	(
    initial_delays,
    ccf_for_delay,
    lower = numeric(n - 1) - nsamples / 4,
    upper = numeric(n - 1) + nsamples / 4
  )
  for (i in 1:n)
  {
    initial_delays = max_lag_by_row[i, ]
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
  maxs_ccf$local_maximal = optimals
  maxs_ccf$maximal = opt
  data$local_maximal=optimals
  data$maximal=opt
  return(data)
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
  stopifnot(length(unique(signals$fs))==1)
  nsamples = signals[["nsamples"]]
  min_n = floor(min(nsamples)/2)*2
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
#' aligned=apply_frame(waveList,indexes)
#'
apply_frame <- function(signals, indexes)
{
  is = indexes$i_start
  ie = indexes$i_end

  s = signals$signal
  nsamples = ie[1] - is[1] + 1
  stopifnot(nsamples%%2==0)
  signal = lapply(1:length(s), function(i)
    s[[i]][is[min(i,length(is))]:ie[min(i,length(ie))]])

  out  = signals
  out$i_start = is
  out$i_end = ie
  out$nsamples=1+ie-is
  out$signal = signal
  out$duration = out$nsamples/out$fs
  out$interval = interval(out$time,
                          out$time+microseconds(out$duration*1e6) )

  out
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
max_lag_ccf <- function(signal)
{
  nsamples <- length(signal$ccf[[1]][[1]])
  max_lag <- sapply(signal$ccf,
                    function (ccfi)
                      sapply(ccfi,
                             function (ccfij)
                               ifelse(
                                 is.null(ccfij), 0,
                                 index_to_lag(which.max(ccfij), nsamples)
                               )))
}

#' calculates optimal lags for each microphone
#'
#' @param max_ccf list containing the paired lags
#'
#' @return it appends the optimal lags for each microphone pair
#' @export
#'
#' @examples
calc_lag_ccf <- function(max_lag)
{
  n = length(max_lag) + 1
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
  d
}





indexes_for_maximal_alignment <- function ( o_frame)
{
  delays = c(0, round(o_frame$maximal$par))
  i_start = o_frame$i_start + delays
  i_start_new = i_start - min(i_start) + 1
  nsamples = o_frame$nsamples
  n_samples_new = min(nsamples - i_start_new + 1)
  i_end_new = i_start_new + n_samples_new-1
  list(i_start = i_start_new, i_end = i_end_new)
}




align_frame <- function(signal_frame)
{
  signal_frame %>%
    apply_hanning_window() %>%
    calculate_ffts() -> signal_frame

  ccf_signal = calculate_list_of_cross_correlations(signal_frame)
  max_ccf = max_lag_ccf(ccf_signal)
  calc_lag_ccf(max_ccf) -> maxs_ccf
  optimize_delays_from_ccf(ccf_signal, max_ccf, maxs_ccf) -> opt_ccf

}

align_recordings <- function(signals, ncylces = 3, show_iterations=FALSE)
{
  i_same_length = indexes_for_same_length(signals)

  apply_frame(signals, i_same_length) %>%
    align_frame() -> o_frame

  i_opt_frame = indexes_for_maximal_alignment( o_frame = o_frame)
  apply_frame(signals, i_opt_frame) %>%
    align_frame() -> o_frame_2

  i_opt_frame_2 = indexes_for_maximal_alignment( o_frame = o_frame_2)
  apply_frame(signals, i_opt_frame_2) -> opt
  if (show_iterations)
  list(opt = opt,
       frame_1 = o_frame,
       frame_2 = o_frame_2)
  else
    opt
}
