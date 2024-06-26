



#' Title
#'
#' @param label
#' @param receptors
#' @param latitudes
#' @param longitudes
#' @param elevations
#' @param time_of_recording
#' @param files
#'
#' @return
#' @export
#'
#' @examples
build_session <-
  function(label,
           receptors,
           latitudes,
           longitudes,
           elevations,
           time_of_recording,
           files)
  {
    return(
      list(
        label = label,
        receptors = receptors,
        lat = latitudes,
        lon = longitudes,
        elevations = elevations,
        time_of_recording = time_of_recording,
        files = files
      )
    )
  }






