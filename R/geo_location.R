
#' geolocalization to local distance (in meters)
#'
#' @param lat latitude
#' @param lon longitude
#' @param origin latitude and longitude of the (0,0) point
#'
#' @return  a list (x,y,origin) where x is the distance to East or West of the origin
#'                                    y is the distance to the North or South of the origin
#'                                    origin is a vector with the latitude and longitude of the origin
#' @export
#'
#' @examples
geo_to_local <- function(lat,lon,origin=NULL)
{
  if (is.null(origin))
    origin=c(lat=mean(lat),lon=mean(lon))


  max_lat = 1e-3

  max_lon = 1e-3

  dlat = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], origin["lat"] + max_lat))/max_lat

  dlon = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"] + max_lon, origin["lat"]))/max_lon


  x = (lon - origin["lon"])*dlon
  y = (lat - origin["lat"])*dlat

  return(list(origin = origin, x = x, y = y))

  }


#' transform local coordinates to the equivalent geolocation, given the geolocation of the origin
#'
#' @param origin a vector containing the latitude and longitude of the origin
#' @param x  distance from the origin in the W-E axis
#' @param y  distance from the origin in the S-N axis
#'
#' @return  the latitude and longitudes of the coordinates also the origin used
#' @export
#'
#' @examples
local_to_geo <- function(origin,x,y)
{


  max_lat = 1e-3

  max_lon = 1e-3

  dlat = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], origin["lat"] + max_lat))/max_lat

  dlon = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"] + max_lon, origin["lat"]))/max_lon


  lon =   x/dlon + origin["lon"]
  lat = y/dlat +  origin["lat"]

  return(list(lat = lat, lon = lon, origin = origin))

}




obtain_receptor_position <- function(session)
{
  session$position = geo_to_local(lat = session$lat,lon = session$lon)
  session$position$z = session$elevations
  return(session)
}

