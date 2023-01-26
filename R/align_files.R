

floor_to_power_of_2 <- function(x)
{
  n = floor(log2(x))
  floor(x / (2 ^ n)) * 2 ^ n
}

round_to_power_of_2 <- function(x)
{
  n = floor(floor(log2(x)))
  round(x / (2 ^ n)) * 2 ^ n
}

ceiling_to_power_of_2 <- function(x)
{
  n = floor(floor(log2(x)))
  ceiling(x / (2 ^ n)) * 2 ^ n
}


floor_to_semi_power_of_2 <- function(x)
{
  n = floor(floor(log2(x)) / 2) + 1
  floor(x / (2 ^ n)) * 2 ^ n
}

round_to_semi_power_of_2 <- function(x)
{
  n = floor(floor(log2(x)) / 2) + 1
  round(x / (2 ^ n)) * 2 ^ n
}

ceiling_to_semi_power_of_2 <- function(x)
{
  n = floor(floor(log2(x)) / 2) + 1
  ceiling(x / (2 ^ n)) * 2 ^ n
}



align_recordings <- function(recordings, lag_window_in_meters = NULL)
{
  if (!"shared_peaks" %in% names(recordings))
    recordings <- obtain_shared_peaks(recordings,
                                      lag_window_in_meters = lag_window_in_meters)

  n_receptors = get_number_of_receptors(recordings)
  peak_start = recordings$shared_peaks$d[1, paste0("lag_1_", c(2:n_receptors))] *
    recordings$fs

  peak_start = c(0, as.numeric(peak_start))

  peak_start = peak_start + recordings$i_start
  peak_start = peak_start - min(peak_start) + 1

  nsamples = floor_to_semi_power_of_2(min(recordings$raw_nsamples - peak_start +
                                            1))



  recordings <-
    apply_frame(recordings,
                "shared_peaks",
                i_start = peak_start,
                nsamples = nsamples)
  return (recordings)

}
