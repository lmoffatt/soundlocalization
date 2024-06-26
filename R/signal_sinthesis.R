







simulate_record<-function(emitters_position=list(x=c(0,10),y=c(1,10)),
                          emitters_signals,
                          labels=c("Andrea","Joaco","Candela","Luciano"),
                          geo_locations=rbind(geo_coordinates(lat=-34.61410605441556,
                                                              lon= -58.38404404725583),
                                              geo_coordinates(-34.614103295104634, -58.3839950969408),
                                              geo_coordinates(-34.6140999839314, -58.383942123312224),
                                              geo_coordinates(-34.61409722462029, -58.38389920796756)),

                          time_of_recordings=ymd_hm("2022-09-03 10:52",tz="America/Buenos_Aires")
                          )

{

  if (length(time_of_recording)==1)
    time_of_recording=rep(time_of_recording,length(files))

  lapply(files, tuneR::readWave) -> l
  nsamples = vapply(l, function(x)
    length(x@left), 1)
  fs = vapply(l, function(x)
    x@samp.rate, 1)
  duration =nsamples/fs

  d=list(
    labels=labels,
    lat=vapply(1:nrow(geo_locations),
               function(i) geo_locations[i,"lat"],1),
    lon=vapply(1:nrow(geo_locations),
               function(i) geo_locations[i,"lon"],1),

    time=time_of_recording,
    interval = interval(time_of_recording,time_of_recording+microseconds(duration*1e6) ),

    dirname = vapply(files, function(x)
      dirname(x), 'c'),
    filename = vapply(files, function(x)
      basename(x), 'c'),
    nsamples=nsamples,
    fs=fs,
    stereo = vapply(l, function(x)
      x@stereo, TRUE),
    bit = vapply(l, function(x)
      x@bit, 1),
    pcm = vapply(l, function(x)
      x@pcm, TRUE)
  )
  d$signal= lapply(l, function(x)
    x@left)
  return (d)

}
