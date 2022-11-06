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
  n = length(ffted_signal$fft[[1]])

  ffted_signal$ccf = lapply(1:(n - 1),
                            function(i)
                              lapply(i:(n - 1),
                                     function(j)
                                       Re(
                                         fft(Conj(ffted_signal$fft[[1]][[i]]) *
                                               ffted_signal$fft[[1]][[j + 1]], inverse = T)
                                       )))

  ffted_signal
}

