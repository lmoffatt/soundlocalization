

chop_in_time<-function(recordings,aligned_start,intervals)
{
 fs=recordings$fs

 i_intervals=lappy(1:length(intervals), function (i)
   list(i_start=aligned_start[i]+intervals[[i]]$start/fs,
        i_end=aligned_start[i]+intervals[[i]]$end/fs
   ))

 lapply(1:length(intervals),
                  function (i)
                    apply_indexes(recordings,i_intervals[[i]])
   )

}
