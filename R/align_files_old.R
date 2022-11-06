align_recordings_old <- function(signals, ncylces = 3, show_iterations=FALSE)
{
  align_frame <- function(signal_frame)
  {
    max_lag_ccf <- function(signal)
    {
      index_to_lag <- function(i, nsamples)
      {
        ifelse(i < nsamples / 2, i - 1, i - nsamples - 1)
      }
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
    signal_frame %>%
      apply_hanning_window() %>%
      calculate_ffts() -> signal_frame

    ccf_signal = calculate_list_of_cross_correlations(signal_frame)
    max_ccf = max_lag_ccf(ccf_signal)
    calc_lag_ccf(max_ccf) -> maxs_ccf
    optimize_delays_from_ccf(ccf_signal, max_ccf, maxs_ccf) -> opt_ccf

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
