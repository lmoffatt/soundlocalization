





build_source_type_statistics <-
  function(maximum_altitude,
           number_per_ha,
           fundamental_freq,
           harmonic_amplitudes,
           duration,
           decay_factor)
  {
    return(
      list(
        maximum_altitude = maximum_altitude,
        number_per_ha = number_per_ha,
        fundamental_freq = fundamental_freq,
        harmonic_amplitudes = harmonic_amplitudes,
        duration = duration,
        decay_factor = decay_factor
      )
    )
  }




simulate_harmonic_signal <-
  function(fs,
           nsamples = NULL,
           trace_duration = NULL,
           fundamental_freq,
           harmonic_amplitudes,
           t_start,
           t_end,
           decay_factor = NULL)
  {
    if (is.null(nsamples))
      nsamples = ceiling(trace_duration * fs)
    t = 0:(nsamples - 1) / fs

    A = harmonic_amplitudes
    w = 1:length(harmonic_amplitudes) * fundamental_freq * 2 * pi

    if (!is.null(decay_factor))
      y = vapply(
        X = t,
        FUN = function (x) {
          if ((x > t_start) & (x < t_end))
            sum(A * exp(-(x - t_start) / (t_end - t_start) / decay_factor) *
                  sin(w * (x - t_start)))
          else
            0
        },
        0
      )
    else
      y = vapply(
        X = t,
        FUN = function(x) {
          if ((x > t_start) & (x < t_end))
            sum(A * sin(w * (x - t_start)))
          else
            0
        },
        0
      )

    return(y)
  }


sample_parameter <- function(n, parameter)
{
  return(rnorm(n, mean = parameter$mean, sd = parameter$sd))
}







sample_sources_signals <- function(sources,
                                   fs,
                                   nsamples = NULL,
                                   n_sources,
                                   sources_statistics)
{
  trace_duration = nsamples / fs

  fundamental_freq = sample_parameter(n_sources, sources_statistics$fundamental_freq)
  harmonic_amplitudes = sample_parameter(n_sources, sources_statistics$harmonic_amplitudes)
  t_start = runif(n_sources, 0, trace_duration)
  t_end = t_start + sample_parameter(n_sources, sources_statistics$duration)
  decay_factor = sample_parameter(n_sources, sources_statistics$decay_factor)

  sources$signals = c(lapply(1:n_sources, function(i)
    return(
      simulate_harmonic_signal(
        fs,
        nsamples = nsamples,
        trace_duration = trace_duration,
        fundamental_freq = fundamental_freq[i],
        harmonic_amplitudes = harmonic_amplitudes[i],
        t_start = t_start[i],
        t_end = t_end[i],
        decay_factor = decay_factor[i]
      )
    )))
  return(sources)
}


sample_source_positions <-
  function(sources,
           n_sources,
           max_distance,
           maximum_altitude)
  {
    sources$x = c(sources$x, runif(n_sources, -max_distance, max_distance))
    sources$y = c(sources$y, runif(n_sources, -max_distance, max_distance))
    sources$z  = c(sources$z,
                   runif(n_sources,
                         0, maximum_altitude))

    return (sources)

  }


sample_all_sources_types <- function(fs,
                                     nsamples ,
                                     min_dist,
                                     noise_level,
                                     all_sources_types)

{
  return(Reduce(function(sources, sources_statistics)
  {
    d = 10 ^ (
      20 * sources_statistics$sources_parameters$harmonic_amplitudes[1]$mean - noise_level
    ) + min_dist
    lambda = sources_statistics$number_per_ha * (d / 100) ^ 2

    n_sources = rpois(1, lambda)

    sources = sample_sources_signals(
      sources,
      fs,
      nsamples = nsamples,
      n_sources = n_sources,
      sources_statistics = sources_statistics
    )

    sources = sample_source_positions(
      sources = sources,
      max_distance = d,
      maximum_altitude = sources_statistics$max_altitude,
      n_sources = n_sources
    )
    return(sources)
  }, all_sources_types, list()))

}


obtain_distances <- function(sources, receptors)
{
  ns = length(sources)
  nr = length(receptors)
  return(vapply(
    1:ns,
    FUN = function(is) {
      return(vapply(
        1:nr,
        FUN = function(ir)
          return(sqrt(
            (sources[is]$x - receptors[ir]$x) ^ 2 +
              (sources[is]$y - receptors[ir]$y) ^ 2 +
              (sources[is]$z - receptors[ir]$z) ^ 2
          )),
        FUN.VALUE = numeric(1)
      ))
    },
    FUN.VALUE = numeric(nr)
  ))
}


obtain_A <- function(fs,
                     nsamples,
                     sources,
                     receptors,
                     velocity_of_sound)
{
  w = c(0:(nsamples / 2 - 1), (-nsamples / 2):-1) * fs / nsamples * 2 * pi

  ns = length(sources)
  nr = length(receptors)
  return(
    vapply(
      1:ns,
      FUN = function(is) {
        return(vapply(
          1:nr,
          FUN = function(ir)
          {
            d = sqrt(
              (sources[is]$x - receptors[ir]$x) ^ 2 +
                (sources[is]$y - receptors[ir]$y) ^ 2 +
                (sources[is]$z - receptors[ir]$z) ^ 2
            )
            t0 = d / velocity_of_sound
            A = 1.0 / d * exp(1i * w * t0)
            return(A)
          },
          FUN.VALUE = numeric(nsamples)
        ))
      },
      FUN.VALUE = matrix(ncol = nsamples, nrow = nr)
    ))
}






obtain_fft <- function(nsamples,porter)
{
  hn = seewave::ftwindow(nsamples, wn = "hanning")
  porter$fft = lapply(porter$signal,function(signal)
    fft(signal*hn))
  return(porter)
}


hear_the_sources <-
  function(fs,
           nsamples,
           sources,
           receptors,
           velocity_of_sound)
  {

    if (!"fft" %in% names(sources))
      sources=obtain_fft(nsamples,sources)
    A = obtain_A(fs = fs,
                 nsamples = nsamples,
                 sources = sources,
                 receptors = receptors,
                 velocity_of_sound = velocity_of_sound)






  }







#' synthesize a simulated session
#' for a given set of emitters of known position and signals
#'
#'
#' @return a session
#' @export
#'
#' @examples
simulate_session <- function(label,
                             fs,
                             nsamples,
                             time_of_recording,
                             noise_level,
                             receptors,
                             all_sources_types)
{
  min_dist = max(max(receptors_x) - min(receptors_x),
                 max(receptors_y) - min(receptors_y))
  sources = sample_all_sources_types(
    fs = fs,
    nsamples = nsamples,
    min_dist = min_dist,
    noise_level = noise_level,
    all_sources_types = all_sources_types
  )


}
