#' Calculates the optimal lag combination for a set of ccfs pairs
#'
#' @param ccf_list should include double indexed list of cross correlation functions
#'
#' @return it inserts the global and local optimal lags
#' @export
#' @examples
optimize_delays_from_ccf <- function(data, maxs_ccf, max_lag_by_row)
{
  nsamples = length(data$ccf[[1]][[1]])
  n = nrow(max_lag_by_row)
  ccf_for_delay <- function(x) {
    lags = c(0, x)

    sum = 0
    for (i in 1:(n - 1))
    {
      for (j in (i + 1):n)
      {
        lag = lags[j] - lags[i]
        l0 = floor(lag)
        l1 = l0 + 1
        r = lag - l0
        s = 1 - r
        i0 = lag_to_index(l0, nsamples)
        i1 = lag_to_index(l1, nsamples)
        sum = sum + data$ccf[[i]][[j - i]][i0] * s + data$ccf[[i]][[j - i]][i1] *
          r
      }
    }
    return (-sum)
  }
  resp = list()
  initial_delays = colMeans(max_lag_by_row)
  resp[[1]] <- nloptr::neldermead	(
    initial_delays,
    ccf_for_delay,
    lower = numeric(n - 1) - nsamples / 4,
    upper = numeric(n - 1) + nsamples / 4
  )
  for (i in 1:n)
  {
    initial_delays = max_lag_by_row[i, ]
    resp[[i + 1]] <- nloptr::neldermead	(
      initial_delays,
      ccf_for_delay,
      lower = numeric(n - 1) - nsamples / 4,
      upper = numeric(n - 1) + nsamples / 4
    )
  }
  optimals = list(
    par = t(vapply(resp, function(r)
      r$par, numeric(n - 1))),
    value = vapply(resp, function(r)
      - r$value, 0),
    iter = vapply(resp, function(r)
      r$iter, 0),
    convergence = vapply(resp, function(r)
      r$convergence, 0),
    message = vapply(resp, function(r)
      r$message, 'text')
  )
  n_opt = which.max(optimals$value)
  opt = resp[[n_opt]]
  maxs_ccf$local_maximal = optimals
  maxs_ccf$maximal = opt
  data$local_maximal=optimals
  data$maximal=opt
  return(data)
}

