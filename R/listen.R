
listen_frame <- function(recordings,frame, i,use_raw_signal=F,factor=1)
{
  tuneR::setWavPlayer('/usr/bin/aplay')

  fs = recordings$fs

  seewave::listen(frame_to_Wave(recordings = recordings, frame=frame,i = i,use_raw = use_raw_signal),
                  f = fs*factor)

}


listen <- function(recordings, i, start, end)
{
  fs = recordings$fs[i]

  seewave::listen(to_Wave(recordings = recordings, i, start, end), f = fs)

}
