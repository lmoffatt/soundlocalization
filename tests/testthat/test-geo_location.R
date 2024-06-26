test_that("geolocation invertibility", {
  lat = c(-35.100313,-35.099890,-35.1003703,-35.099870,-35.100110)
  lon = c(-58.304676,-58.304550,-58.3041342,-58.304039,-58.304306)

  local = geo_to_local(lat = lat,lon = lon)
  geo = local_to_geo(origin = local$origin, x = local$x, y = local$y)

  expect_equal(geo,list(lat = lat,lon = lon,origin = local$origin))
})
