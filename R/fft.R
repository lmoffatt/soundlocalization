


calculate_ffts_per_frame <- function(rec,frame)
{
  if (!"signal_for_fft" %in% names(frame))
    frame = apply_hanning_window_per_frame(rec,frame)
  frame$fft = lapply(frame$signal_for_fft,fft)
  return (frame)
}





#' Title obtains the fourier transform of all the signals at once
#'
#' @param windowed_signal the hanning window should be already applied
#'
#' @return inserts a field with the fft
#' @export
#'
#' @examples
calculate_ffts <- function(rec)
{

  if ("frames" %in% names(rec))
  {
     rec$frames=lapply(rec$frames,
                       function(frame)
                         calculate_ffts_per_frame(rec,frame))
     return (rec)
  }
  else
    return (calculate_ffts_per_frame(rec,rec))
}



