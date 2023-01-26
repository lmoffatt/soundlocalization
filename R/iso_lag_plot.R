iso_lag_data <-
  function(rec,
           zoom_out = 5,
           n_points = 100,
           sound_velocity = 334)
  {
    dx = max(rec$x) - min(rec$x)
    dy = max(rec$y) - min(rec$y)
    dist = max(dx, dy) * zoom_out
    x_mean = max(rec$x) * 0.5 + min(rec$x) * 0.5
    y_mean = max(rec$y) * 0.5 + min(rec$y) * 0.5

    p = ((-n_points / 2):(n_points / 2)) / n_points * dist

    x = p + x_mean
    y = p + y_mean
    dp = full_join(data.frame(x = x),
                   data.frame(y = y),by=character())

    n_receptors = get_number_of_receptors(rec)
    d_1 = ((dp$x - rec$x[1]) ^ 2 + (dp$y - rec$y[1]) ^ 2) ^ 0.5

    d = Reduce(function (p, i)
    {
      dp = full_join(data.frame(x = x),
                     data.frame(y = y),by=character())

      dp$var = paste0('lag_1_', i)

      dp$lag= (((dp$x - rec$x[i]) ^ 2 +
                                      (dp$y - rec$y[i]) ^ 2) ^ 0.5 - d_1) / sound_velocity

      return (rbind(p, dp))
    },
    2:n_receptors,
    data.frame())
    return (d)
  }
