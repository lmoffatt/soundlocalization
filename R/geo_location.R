






geo_distance <- function(x, y)
{
  geosphere::distGeo(p1 = c(x["lon"], x["lat"]), p2 = c(y["lon"], y["lat"]))
}


geo_to_local2 <- function(origin, x_axis, points)
{
  dp = t(vapply(1:nrow(points),
              function (i)
                c(lat=geosphere::distGeo(
                  c(origin["lon"], origin["lat"]),
                  c(origin["lon"], points[i,"lat"])
                ),
                lon=geosphere::distGeo(
                  c(origin["lon"], origin["lat"]),
                  c(points[i,"lon"], origin["lat"])
                )),
              numeric(2)))

  print(dp)

  dx = c(geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], x_axis["lat"])),
         geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(x_axis["lon"], origin["lat"])))

  dx = dx / norm(dx, type = "2")
  dy = c(-dx[2], dx[1])
  dy = dy / norm(dy, type = "2")

  print(dx%*%dy)

  cbind(x = dp %*% dx,
       y = dp %*% dy)

}


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

  return(list(lat = lat, lon = lon))

}




#' Title
#'
#' @param lat
#' @param lon
#'
#' @return
#' @export
#'
#' @examples
geo_coordinates <- function(lat, lon,error=4)
{
  c(lon = lon, lat = lat, error=error)
}
