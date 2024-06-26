

set_receptor <- function (label,
                          file,
                          geo_location,
                          time_of_recording)
{
  return(
    list(
      label = label,
      file = file,
      geo_location = geo_location,
      time_of_recording = time_of_recording
    )
  )
}

#' Title
#'
#' @param label
#' @param x_pos
#' @param y_pos
#' @param error
#' @param t_recordings
#' @param t_starts
#' @param t_ends
#' @param min_freq
#' @param max_freq
#'
#' @return
#' @export
#'
#' @examples
set_sources <- function(label,
                        x_pos,
                        y_pos,
                        error,
                        t_recordings,
                        t_starts,
                        t_ends,
                        min_freq,
                        max_freq)
{
  return (
    list(
      label = label,
      x = x_pos,
      y = y_pos,
      error = error,
      t_starts = t_starts,
      t_ends = t_ends,
      min_freq = min_freq,
      max_freq = max_freq
    )
  )
}









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
           time_of_recording,
           velocity_of_sound =334,
           position_x =NULL,
           position_y =NULL,
           position_error =NULL,
           origin = NULL,
           x_axis = NULL)
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


  p = geo_to_local(points = geo_locations,
                   origin = origin,
                     x_axis = x_axis)


  if (is.null(position_x))
    position_x =  p[,1]
  if (is.null(position_y))
    position_y =  p[,2]
  if (is.null(position_error))
    position_error =   vapply(1:nrow(geo_locations),
                              function(i)
                                geo_locations[i, "error"], 1)




    d = list(
      labels = labels,
      lat = vapply(1:nrow(geo_locations),
                   function(i)
                     geo_locations[i, "lat"], 1),
      lon = vapply(1:nrow(geo_locations),
                   function(i)
                     geo_locations[i, "lon"], 1),
      x = position_x,
      y = position_y,
      pos_err = position_error,
      velocity_of_sound = velocity_of_sound,
      origin = origin,
      x_axis = x_axis,
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
    "x",
    "y",
    "velocity_of_sound",
    "pos_err",
    "origin",
    "x_axis",
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



simulate_recording <- function(receptor_positions_x,
                               receptor_positions_y,
                               receptor_gain_change_time,
                               receptors_start_logGain,
                               receptor_new_logGain,
                               receptor_offsets,
                              source_files,
                              source_labels,
                              source_start_time,
                              source_end_time,
                              source_positions_x,
                              source_positions_y,
                              velocity_of_sound =334,
                               origin = NULL,
                              x_axis = NULL)

{
  sim=list()
  sim$receptors_sim = list(
    x = receptor_positions_x,
    y = receptor_positions_y,
    logG = receptors_start_logGain,
    offset = receptor_offsets
  )





}



get_number_of_receptors <- function(rec)
{
  return(length(rec$x))
}

