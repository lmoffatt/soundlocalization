












Jacobian_by_diff <- function(f, x, xdata, eps)
{
  k = length(x)
  n=ncol(xdata)
  Y = f(x, xdata)

  dYdx = vapply(seq_len(k),function(i)
  {
    e = numeric(length(x))
    e[i] = max(1, abs(x[i])) * eps
    Yi = f(x + e, xdata)
    dYi = (Yi - Y) / e[i]
    return(dYi)
  },
  Y[seq_along(Y)])

  f_name = deparse(substitute(f))
  x_names = names(x)
  if (is.null(x_names))
    x_names = 1:length(x)

  return (dYdx)
}


get_stochastic_indexes <-
  function(number_of_samples, number_of_chunks)
  {
    random_indexes = sample.int(number_of_samples)
    return(lapply(1:number_of_chunks, function (i)
    {
      i_start = 1 + (i - 1) / number_of_chunks * number_of_samples
      i_end = i / number_of_chunks * number_of_samples
      return(random_indexes[i_start:i_end])
    }))
  }







#' Stochastic Levenberg Marquardt Algorithm
#' Minimizes the square sum of f(x,xdata), taking samples of xdata for calculating the steps
#'
#' @param f  function to be optimized
#' @param x0 initial guess of the multidimentional parameter to be optimized
#' @param xdata  independent possibly multivariate variable where f is applied
#' @param ydata  rowise monovariate variable given for each xdata
#' @param Jacobian optional, function that produces the Jacobian of f
#' @param number_of_chunks number of samples from the data
#' @param maxiter maximum number of iteration
#' @param return_Fn  whether the function values are returned
#' @param mu  initial mu parameter for Levenberg-Marquardt
#' @param nu  initial nu parameter (it multiplies mu in adaptation)
#' @param eps epsilon value for Jacobian
#' @param max_mu maximum mu value
#'
#' @return a list containing the optimal x and other things
#' @export
#'
#' @examples
#'
#'
#' sum_of_Normals <- function(x, X){
#'   return(x["A1"] * exp(-(X[,"x1"] - x["mean1"])^2/x["sd1"]^2) +
#'          x["A2"] * exp(-(X[,"x2"] - x["mean2"])^2/x["sd2"]^2))
#' }
#'
#' x = c(A1 = 12.3,mean1 = 24, sd1 = 17.3, A2 = 0.23,mean2 = 21.2 ,sd2 = 5)
#' X=cbind(runif(1e5) * 100 - 20, runif(1e5) * 100 - 20)
#' colnames(X)<-c('x1','x2')
#'
#' Yn = sum_of_Normals(x,X)+rnorm(1e5)*0.2
#'
#'
#' opt=slm(f =sum_of_Normals, x0 = x, xdata = X, ydata = Yn,number_of_chunks = 100)
#'
slm <- function(f,
                x0,
                xdata,
                ydata,
                number_of_chunks,
                maxiter = 10000,
                return_Fn = FALSE,
                Jacobian = NULL,
                mu = 1,
                nu = 2,
                eps = 1e-5,
                max_mu = 1e10)
{
  #ydata = matrix(data = as.numeric(ydata), ncol = 1)
  stopifnot("same number of rows" = ncol(xdata) == ncol(ydata))

  no_progress <- function(ca, eps)
  {
    if (length(ca$previously) == 0)
      return (FALSE)
    return (all(abs(ca$previously$x - ca$x) < eps * (abs(ca$x) + eps)))
  }

  calc_rho <- function(pr, ca)
  {
    return(2 * (pr$sqrsum_i - ca$sqrsum_0) / (t(ca$hn) %*% (pr$ln %*% ca$hn - pr$gradient_i)))
  }

  initial_step <- function(xn, mu, nu)
  {
    ii = stc_ind[[1]]
    Fn = f(xn, xdata[,ii])
    yn = ydata[,ii]
    fn = sum(Re(Conj(Fn - yn)*(Fn - yn)))
    Jn = Jacobian(xn, xdata[,ii])
    gn = Re(Conj(t(Jn)) %*% (Fn - yn)[seq_along(Fn)])
    A = Re(Conj(t(Jn)) %*% Jn)
    ln = mu * diag(diag(A))
    ca = list(
      previously = list(),
      x = xn,
      sqrsum_tot = fn,
      sqrsum_i = fn,
      gradient_tot = gn,
      gradient_i = gn,
      A = A,
      ln = ln,
      nsteps = 1,
      done = FALSE,
      mu = mu,
      nu = nu
    )
    if (return_Fn)
    {
      Fn_tot  = ydata - ydata
      Fn_tot[ii] <- Fn
      ca$Fn_tot = Fn_tot
    }
    return (ca)

  }


  if (is.null(Jacobian))
    Jacobian <-
    function(x, xdata, ydata)
      Jacobian_by_diff(
        f = f,
        x = x,
        xdata = xdata,
        eps = eps
      )

  n = ncol(xdata)
  stc_ind = get_stochastic_indexes(number_of_samples = n, number_of_chunks = number_of_chunks)
  Reduce(
    f = function(pr, i) {
      if ((pr$done))
        return (pr)
      ca = list(
        previously = pr$previously,
        nsteps = i,
        mu = pr$mu,
        nu = pr$nu,
        done = FALSE
      )
      ca$ln = ca$mu * diag(diag(pr$A+eps))

      while ((rcond(pr$A + ca$ln) < eps) & (ca$mu < max_mu/2))
      {
        ca$mu = min(ca$mu * ca$nu, max_mu)
        ca$nu = min(2 * ca$nu, 2 ^ 5)
        ca$ln = ca$mu * (diag(diag(pr$A))+diag(eps,nrow=nrow(pr$A)))
        print(rcond(pr$A + ca$ln))
        print(ca$mu)
       }
       print('out')
      if ((rcond(pr$A + ca$ln) < eps ))
      {
        ca$ln = ca$mu * diag(nrow = nrow(pr$A))
      }


      ca$hn = solve(pr$A + ca$ln,-pr$gradient_i)
      ca$x = as.numeric(pr$x + ca$hn)
      names(ca$x) <- names(pr$x)
      ii0 = stc_ind[[(i - 2) %% number_of_chunks + 1]]
      Fn0 = f(ca$x, xdata[,ii0])
      yn0 = ydata[,ii0]
      ca$sqrsum_0 = sum(Re(Conj(Fn0 - yn0) * (Fn0 - yn0)))
      rho = calc_rho(pr, ca)
      while (rho <= 0 && (ca$mu < max_mu))
      {
        ca$mu = min(ca$mu * ca$nu, max_mu)
        ca$nu = min(2 * ca$nu, 2 ^ 5)
        ca$ln = ca$mu * diag(diag(pr$A))
        ca$hn = solve(pr$A + ca$ln,-pr$gradient_i)
        ca$x = as.numeric(pr$x + ca$hn)
        names(ca$x) <- names(pr$x)
        ii0 = stc_ind[[(i - 2) %% number_of_chunks + 1]]
        Fn0 = f(ca$x, xdata[,ii0])
        yn0 = ydata[,ii0]
        ca$sqrsum_0 = sum(Re(Conj(Fn0 - yn0) * (Fn0 - yn0)))
        rho = calc_rho(pr, ca)
      }
      if (rho > 0)
      {
        ca$mu = pr$mu * max(1 / 3, 1 - (2 * rho - 1) ^ 3)
        ca$nu = 2
        ca$stepcondition = list(condition = 'acceptable', rho = rho)
      }
      else
      {
        ca$mu = min(pr$mu * pr$nu, max_mu)
        ca$nu = min(2 * pr$nu, 2 ^ 6)
        stepcondition = list(condition = 'reject', rho = rho)
      }

      ii = stc_ind[[(i - 1) %% number_of_chunks + 1]]
      Fn = f(ca$x, xdata[,ii])
      yn = ydata[,ii]
      ca$sqrsum_i = sum(Re(Conj(Fn - yn)*(Fn - yn)))
      Jn = Jacobian(ca$x, xdata[,ii],yn)
      ca$A = Re(Conj(t(Jn)) %*% Jn)
      ca$gradient_i = Re(Conj(t(Jn)) %*% (Fn - yn)[seq_along(Fn)])
      ca$sqrsum_tot = pr$sqrsum_tot + ca$sqrsum_i
      ca$gradient_tot = pr$gradient_tot + ca$gradient_i
      if (return_Fn)
      {
        ca$Fn_tot = pr$Fn_tot
        ca$Fn_tot[ii] = Fn
      }
      if (i %% number_of_chunks != 0)
        return (ca)
      else if (max(abs(ca$gradient_tot)) < eps)
      {
        ca$stopcondition = list(reason = 'small gradient norm', gradient_norm =
                                  max(abs(ca$gradient_tot)))
        ca$done = TRUE
        return(ca)
      }
      else if (no_progress(ca, eps))
      {
        hn = ca$previously$x - ca$x
        ca$stopcondition = list(reason = 'no progress',
                                delta_x = hn)
        ca$done = TRUE
        return (ca)
      }
      else
      {
        print(ca$sqrsum_tot)
        ca$previously = within(ca, rm(previously))
        ca$sqrsum_tot = 0
        ca$gradient_tot = 0

        return (ca)
      }
    },
    x = 2:maxiter,
    init = initial_step(xn = x0, mu = mu, nu = nu)
  )

}
