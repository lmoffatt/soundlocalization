#' Load several wav files at once
#'
#' @param list of files to be loaded
#'
#' @return a list of wave objects
#' @export
#'
#' @examples
#'
#'
#'
loadwavfiles<-function(fileList)
{
   lapply(fileList,readWave)
}
