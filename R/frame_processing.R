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
indexes_for_same_length <- function(rec)
{
  nsamples = rec[["raw_nsamples"]]
  min_n = floor(min(nsamples)/2)*2
  diff_n = nsamples - min_n
  rec$i_start = floor(diff_n / 2) + 1
  rec$nsamples = rec$i_start + min_n - 1
  list(i_start = i_start, nsamples = min_n)
}

apply_frame_on_frames <- function(rec, out,framed_by,i_start,nsamples)
{
  is = i_start
  if (length(i_start)==1)
    is = rec$i_start +i_start-1

  ie = is +nsamples -1
  s = rec$raw_signal
  stopifnot(nsamples%%2==0)

  signal = lapply(seq_along(s), function(i)
    s[[i]][is[i]:ie[i]])



  nsamples=1+ie-is
  duration=nsamples/rec$fs
  interval=interval(rec$time,
                    rec$time+microseconds(duration*1e6) )

  out$framed_by=framed_by
  out$i_start=is
  out$t_start=is/rec$fs[1]
  out$t_end=(is+nsamples)/rec$fs[1]

  out$nsamples=nsamples
  out$framed_raw_signal=signal
  out$duration=duration
  out$iterval=interval
  return(out)

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
apply_frame <- function(rec, framed_by,i_start,nsamples)
{
  out=extract_raw_recording(rec)
return(apply_frame_on_frames(rec,out,framed_by,i_start,nsamples))
}

truncate_to_same_length<-function(rec)
{
  nsamples = rec[["raw_nsamples"]]
  min_n = floor(min(nsamples)/2)*2
  diff_n = nsamples - min_n
  i_start = floor(diff_n / 2) + 1
  return (apply_frame(rec,framed_by = "truncate_to_same_length",i_start = i_start,nsamples = min_n))

}



apply_this_frames <- function(rec, framed_by,i_start,nsamples)
{
  rec$frames=lapply(1:length(i_start),
                    function(i)
                      apply_frame_on_frames(rec,list(),framed_by,i_start[i], nsamples[i]))
  return (rec)
}




frame_recordings<-function(rec, framed_by="frame_recordings",
                           i_start_diff=NULL, nsamples=NULL,
                           for_spectrogram=FALSE,
                           window_size_in_seconds=NULL,
                           minimum_frame_size_in_meters=NULL,
                           sound_velocity=334)
{
  if (!"framed_by" %in% names(rec))
    rec <- align_recordings(rec)


  stopifnot("nsamples are not equal"=length(unique(rec$nsamples))==1)
  if (is.null(nsamples))
     nsamples=rec$nsamples[1]
  max_distance=(diff(range(rec$x))^2+diff(range(rec$y))^2)^0.5

  frame_size=512
  if (is.null(minimum_frame_size_in_meters)&& !for_spectrogram)
    minimum_frame_size_in_meters=max_distance

  if (!is.null(minimum_frame_size_in_meters))
     frame_size=minimum_frame_size_in_meters/sound_velocity*8*rec$fs
  if (!is.null(window_size_in_seconds))
    frame_size=window_size_in_seconds*rec$fs

  frame_size=ceiling_to_power_of_2(frame_size)

  if (is.null(i_start_diff))
     i_start_diff=1
  i_start=seq(i_start_diff,i_start_diff+nsamples-frame_size,frame_size/2)
  samples=rep(frame_size,length(i_start))


 return (apply_this_frames(rec = rec,framed_by = 'frame_recordings',i_start = i_start,nsamples = samples))
}





