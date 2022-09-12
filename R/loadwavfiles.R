#' Load several wav files at once
#'
#' @param list of files to be loaded
#'
#' @return a list of wave objects
#' @export
#'
#' @examples
#' files<-list("~/papers/ranas/recordings/venezuela 0/Grabacion_3_luciano.wav")
#' d<-loadwavfiles(files)
#'
#'
loadwavfiles<-function(fileList)
{
   lapply(fileList,tuneR::readWave)
}
