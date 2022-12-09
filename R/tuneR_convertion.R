to_Wave <- function(recordings, i, start, end)
{
  fs = as.numeric(recordings$fs[1])
  signal=NA
  if ("framed_raw_signal" %in% names(recordings))
    signal=recordings$framed_raw_signal
  else if ("raw_signal" %in% names(recordings))
    signal=recordings$raw_signal

  en = min(length(signal[[i]]), end * fs)
  mysignal=signal[[i]][(start * fs + 1):en]
  tuneR::Wave(
    left = mysignal,
    samp.rate = fs,
    bit = recordings$bit[i],
    pcm = recordings$pcm[i]
  )
}

frame_to_Wave <- function(recordings, frame,i,use_raw=FALSE)
{
  fs = as.numeric(recordings$fs[1])
  signal=NA
  if (("signal_for_fft" %in% names(frame))&!use_raw)
    signal=frame$signal_for_fft
  else if ("framed_raw_signal" %in% names(frame))
    signal=frame$framed_raw_signal

  mysignal=signal[[i]]
  tuneR::Wave(
    left = mysignal,
    samp.rate = fs,
    bit = recordings$bit[i],
    pcm = recordings$pcm[i]
  )
}
