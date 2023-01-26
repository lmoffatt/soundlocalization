


random_framed_sources <- function(rec,
                          n_sources = NULL,
                          x_min = NULL,
                          x_max = NULL,
                          y_min = NULL,
                          y_max = NULL,
                          freq_limits = rbind(c(0, 500),
                                              c(500, 4000),
                                              c(4000, 10000))
)
  {

      if (is.null(x_min))
        x_min = min(rec$x)
      if (is.null(x_max))
        x_max = max(rec$x)
      if (is.null(y_min))
        y_min = min(rec$y)
      if (is.null(y_max))
        y_max = max(rec$y)

      n_receptors = get_number_of_receptors(rec)
      if (is.null(n_sources))
        n_sources = n_receptors - 1


      n_frames = get_number_of_frames(rec)
      n_freq = nrow(freq_limits)

      total_rows = n_frames * n_freq


      i_frame_pos = rep(seq_len(n_frames), each = n_freq)
      i_freq_pos = rep(seq_len(n_freq), times = n_frames)

      x = matrix(
        x_min + (x_max - x_min) * runif(n_sources * total_rows),
        nrow = total_rows,
        ncol = n_sources
      )
      y = matrix(
        y_min + (y_max - y_min) * runif(n_sources * total_rows),
        nrow = total_rows,
        ncol = n_sources
      )



      return(list(
        i_frame_pos = i_frame_pos,
        i_freq_pos = i_freq_pos,
        x = x,
        y = y
      ))


      return(rec)
}
