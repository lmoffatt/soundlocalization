




#' Title
#'
#' @param files a list with the files where the quasi-simultaneous recordings are stored
#' @param labels labels for the recordings
#' @param geo_locations location of the recordings
#' @param time_of_recording time of the recordings
#'
#' @return a list object containing the recordings ready for processing
#' @export
#'
#' @examples
recordings <-
  function(files,
           labels,
           geo_locations,
           time_of_recording)
  {
    # lets check that we have the geo_location for each receptor
    stopifnot(
      "the number files and geo_locations is different" =
        length(files) == nrow(geo_locations) &
        length(files) == length(labels)
    )



    if (length(time_of_recording) == 1)
      time_of_recording = rep(time_of_recording, length(files))

    lapply(files, tuneR::readWave) -> l
    nsamples = vapply(l, function(x)
      length(x@left), 1)
    fs = vapply(l, function(x)
      x@samp.rate, 1)
    duration = nsamples / fs

    stopifnot("recordings differ in their sampling rate" = length(unique(fs)) ==
                1)

    d = list(
      labels = labels,
      lat = vapply(1:nrow(geo_locations),
                   function(i)
                     geo_locations[i, "lat"], 1),
      lon = vapply(1:nrow(geo_locations),
                   function(i)
                     geo_locations[i, "lon"], 1),
      time = time_of_recording,
      dirname = vapply(files, function(x)
        dirname(x), 'c'),
      filename = vapply(files, function(x)
        basename(x), 'c'),
      fs = unique(fs),
      stereo = vapply(l, function(x)
        x@stereo, TRUE),
      bit = vapply(l, function(x)
        x@bit, 1),
      pcm = vapply(l, function(x)
        x@pcm, TRUE),
      raw_nsamples = nsamples,
      raw_interval = interval(
        time_of_recording,
        time_of_recording + microseconds(duration * 1e6)
      )
    )
    d$raw_signal = lapply(l, function(x)
      x@left)
    return (d)

  }


extract_raw_recording <- function(somewhat_processed_recording)
{
  field_list = c(
    "labels",
    "lat",
    "lon",
    "time",
    "dirname",
    "filename",
    "stereo",
    "bit",
    "pcm",
    "fs",
    "raw_interval",
    "raw_nsamples",
    "raw_signal"
  )

  return (somewhat_processed_recording[field_list])
}
