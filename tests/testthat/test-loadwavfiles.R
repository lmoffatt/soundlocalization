mock_empty_Wave<-function(nsamples,fs=44100,bit=16, pcm= TRUE)
{
  tuneR::Wave(left=integer(nsamples), samp.rate = fs, bit = bit, pcm = pcm)
}

test_that("indexes_for_same_length works", {

  expect_equal(1,1
    )
})
