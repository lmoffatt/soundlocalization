


align_recordings<-function(recordings,lag_window_in_meters=NULL)
{
  round_to_power_of_2<-function(x,d=6)
  {
    n=round(floor(log2(x))/2)
    floor(x/2^n)*2^n
  }
  if (!"shared_peaks" %in% names(recordings))
    recordings <- obtain_shared_peaks(
      recordings,
      lag_window_in_meters = lag_window_in_meters
    )

   n_receptors=length(recordings$labels)
   peak_start=recordings$shared_peaks$d[1,paste0("lag_1_",c(2:n_receptors))]*recordings$fs

   peak_start=c(1,as.numeric(peak_start))

   peak_start=peak_start+recordings$i_start
   peak_start=peak_start-min(peak_start)+1

   nsamples=round_to_power_of_2(min(recordings$raw_nsamples-peak_start+1))



   recordings<-apply_frame(recordings,"shared_peaks",i_start = peak_start, nsamples = nsamples)
   return (recordings)

}
