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



