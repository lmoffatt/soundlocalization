#' Load several wav files at once
#'
#' @param filelist the list of files to be loaded
#'
#' @return list with a dataframe of file names and metadata and a list of signals.
#' @export
#'
#' @examples
#' files<-list(
#' system.file("extdata", "Grabacion_2_andrea.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_joaco.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_1_cande.wav", package = "soundlocalization"),
#' system.file("extdata", "Grabacion_3_luciano.wav", package = "soundlocalization")
#' )
#' waveList<-loadwavfiles(files)
#'
#'
loadwavfiles <- function(fileList)
{
  lapply(fileList, tuneR::readWave) -> l
  list(
    meta_data = data.frame(
      dirname = vapply(fileList, function(x)
        dirname(x), 'c'),
      filename = vapply(fileList, function(x)
        basename(x), 'c'),
      nsamples = vapply(l, function(x)
        length(x@left), 1),
      stereo = vapply(l, function(x)
        x@stereo, TRUE),
      fs = vapply(l, function(x)
        x@samp.rate, 1),
      bit = vapply(l, function(x)
        x@bit, 1),
      pcm = vapply(l, function(x)
        x@pcm, TRUE)
    ),
    signal = lapply(l, function(x)
      x@left)
  )
}






savewavfiles <- function(meta_data,
                         signal,
                         suffix)
{
  addsuffix_to_file<-function (x,suffix)
  {
    paste0(substr(x,0,nchar(x)-4),suffix,".wav")
  }


  for (i in 1:ncol(signal))
    tuneR::writeWave(
      object = to_Wave(meta_data, signal, i),
      filename = addsuffix_to_file(meta_data$filename[i], suffix)
    )

}

