

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
geo_to_local <- function(origin, x_axis, points)
{
  dy = c(geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], x_axis["lat"])),
         geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(x_axis["lon"], origin["lat"])))

  dy = dy / norm(dy, type = "2")
  dx = c(-dy[2], dy[1])
  dlat = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], origin["lat"] + 1e-3))/1e-3

  dlon = geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"] + 1e-3, origin["lat"]))/1e-3



  dp = t(vapply(1:nrow(points),
                function (i)
                  c((points[i,"lat"]- origin["lat"])*dlat,
                    (points[i,"lon"]-origin["lon"])*dlon)
                  ,
                numeric(2)))
    cbind(x = dp %*% dx,
        y = dp %*% dy)

}


local_to_geo <- function(origin, x_axis, x, y)
{
  dx = c(geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(origin["lon"], x_axis["lat"])),
         geosphere::distGeo(c(origin["lon"], origin["lat"]),
                            c(x_axis["lon"], origin["lat"])))

  dx = dx / norm(dx, type = "2")
  dy = c(-dx[2], dx[1])
  dy = dy / norm(dy, type = "2")

  d=rbind(dx,dy)

  dx=d[,1]
  dy=d[,2]



  p = matrix(c(x, y), nrow = length(x), ncol = 2)



  dlat = 1e-3 / geosphere::distGeo(c(origin["lon"], origin["lat"]),
                                   c(origin["lon"], origin["lat"] + 1e-3))

  dlon = 1e-3 / geosphere::distGeo(c(origin["lon"], origin["lat"]),
                                   c(origin["lon"] + 1e-3, origin["lat"]))


  cbind(lat = p %*% dx * dlat + origin["lat"],
                  lon = p %*% dy * dlon + origin["lon"])->c
  colnames(c)<-c("lat","lon")
  return(c)

}




geo_coordinates <- function(lat, lon)
{
  c(lon = lon, lat = lat)
}
