test_that("wave_IO duality", {
  fname = system.file("extdata", "/campo/Grabacion 10_andrea.wav", package = "soundlocalization")
  test_fname = "test.wav"
  w = read_wav_file(fname)
  write_wav_file(filename = test_fname ,fs = w$fs,signal = w$signal)
  w2 = read_wav_file(test_fname)

  expect_equal(w,w2)

})
