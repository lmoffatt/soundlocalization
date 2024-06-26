

read_wav_file <- function(filename)
{
  w = tuneR::readWave(filename = filename)
  return(list(fs = w@samp.rate, signal = w@left, bit = w@bit))
}



write_wav_file <- function(filename,fs,signal, bit =16)
{
  w = tuneR::Wave(left = signal,samp.rate = fs, bit = bit)
  tuneR::writeWave(w,filename)
}

obtain_signals <- function(session)
{
  s = lapply(session$files,read_wav_file)
  fs = vapply(s,function(w) w$fs,numeric(1))
  fs = unique(fs)
  stopifnot("all frequencies equal" = length(fs) == 1)
  signal = lapply(s,function(w) w$signal)
  return (list(fs = fs,signal = signal))
}


