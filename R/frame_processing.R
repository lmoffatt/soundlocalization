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
  nsamples = signals[["raw_nsamples"]]
  min_n = floor(min(nsamples)/2)*2
  diff_n = nsamples - min_n
  i_start = floor(diff_n / 2) + 1
  i_end = i_start + min_n - 1
  list(i_start = i_start, nsamples = min_n)
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
apply_frame <- function(rec, framed_by,i_start,nsamples)
{
  is = i_start
  if (length(i_start)==1)
    is = rec$i_start +i_start

  ie = is +nsamples -1




  s = rec$raw_signal
  stopifnot(nsamples%%2==0)

  signal = lapply(1:length(s), function(i)
    s[[i]][is[i]:ie[i]])



  nsamples=1+ie-is
  duration=nsamples/rec$fs
  interval=interval(rec$time,
                    rec$time+microseconds(duration*1e6) )

  out=extract_raw_recording(rec)
  out$framed_by=framed_by
  out$i_start=is
  out$nsamples=nsamples
  out$framed_raw_signal=signal
  out$duration=duration
  out$iterval=interval

  return(out)
}
truncate_to_same_length<-function(rec)
{
  nsamples = rec[["raw_nsamples"]]
  min_n = floor(min(nsamples)/2)*2
  diff_n = nsamples - min_n
  i_start = floor(diff_n / 2) + 1
  return (apply_frame(rec,framed_by = "truncate_to_same_length",i_start = i_start,nsamples = min_n))

}

