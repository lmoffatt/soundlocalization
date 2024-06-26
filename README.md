
<!-- README.md is generated from README.Rmd. Please edit that file -->

# soundlocalization

<!-- badges: start -->
<!-- badges: end -->

The goal of soundlocalization is to isolate and localize the sources of
sound in a set of simultaneous recordings, with emphasis in birds and
frogs.

## Installation

You can install the development version of soundlocalization from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lmoffatt/soundlocalization")
```

## Example

``` r
library(soundlocalization)
library(ggplot2)
library(lubridate)
#> 
#> Adjuntando el paquete: 'lubridate'
#> The following objects are masked from 'package:base':
#> 
#>     date, intersect, setdiff, union
```

The purpose of this package is to provide tools for isolating and
localizing the source of sounds using an array of receptors on a not
particularly regular order.

We define as a recording session, a list of recordings taken
approximately simultaneously on an array of receptors.

For each session we need an identifying arbitrary label for the session,
the approximate time of the start of the recordings and a list of
recordings each with a unique identifying label, the file name of the
audio archive and the approximate location of the receptor (latitude and
longitude and height above the ground).

``` r
files <- list(
  system.file("extdata", "/campo/Grabacion 10_andrea.wav", package = "soundlocalization"),
  system.file("extdata", "/campo/Grabacion 6_candela.wav", package = "soundlocalization"),
  system.file("extdata", "/campo/Grabacion 3_joaquin.wav", package = "soundlocalization"),
  system.file("extdata", "/campo/Voz 006_julieta.wav", package = "soundlocalization"),
  system.file("extdata", "/campo/Grabacion 18_luciano.wav", package = "soundlocalization")
)
labels = c("Andrea", "Candela", "Joaco", "Julieta","Luciano")
lat= c(-35.100313,-35.099890,-35.1003703,-35.099870,-35.100110)
lon= c(-58.304676,-58.304550,-58.3041342,-58.304039,-58.304306)
elevations =c(1.4,1.42,1.48,1.4,1.54)
time_of_recording = lubridate::ymd_hm("20221225 1735")
```

``` r
session <- build_session(label = "detras_de_la_casa",receptors = labels,latitudes = lat,longitudes = lon,elevations = elevations,time_of_recording = time_of_recording, files = files)
```

Lets see the organization of rec

``` r
str(session)
#> List of 7
#>  $ label            : chr "detras_de_la_casa"
#>  $ receptors        : chr [1:5] "Andrea" "Candela" "Joaco" "Julieta" ...
#>  $ lat              : num [1:5] -35.1 -35.1 -35.1 -35.1 -35.1
#>  $ lon              : num [1:5] -58.3 -58.3 -58.3 -58.3 -58.3
#>  $ elevations       : num [1:5] 1.4 1.42 1.48 1.4 1.54
#>  $ time_of_recording: POSIXct[1:1], format: "2022-12-25 17:35:00"
#>  $ files            :List of 5
#>   ..$ : chr "/tmp/RtmpdhiCus/temp_libpath2f4a61d486441/soundlocalization/extdata//campo/Grabacion 10_andrea.wav"
#>   ..$ : chr "/tmp/RtmpdhiCus/temp_libpath2f4a61d486441/soundlocalization/extdata//campo/Grabacion 6_candela.wav"
#>   ..$ : chr "/tmp/RtmpdhiCus/temp_libpath2f4a61d486441/soundlocalization/extdata//campo/Grabacion 3_joaquin.wav"
#>   ..$ : chr "/tmp/RtmpdhiCus/temp_libpath2f4a61d486441/soundlocalization/extdata//campo/Voz 006_julieta.wav"
#>   ..$ : chr "/tmp/RtmpdhiCus/temp_libpath2f4a61d486441/soundlocalization/extdata//campo/Grabacion 18_luciano.wav"
```

Now the idea of soundlocalization is to locate the time, frequency and
point of origin of the sounds of a recording.

There are two elements in soundlocalization:

1.  Sound receptors
2.  Sound sources (or emitters)

We have an initial measure of the position of the receptors and possibly
of some sources.

We have the data of the time of each of the simultaneous recordings done
on each receptor.

With purpose of a clear analysis of the recordings it is convenient to
define frames as intervals of time and frequency.

The main tool for the separation of the sound sources is the generalized
cross correlation -phase, a normalization of the cross correlation
function. The idea is that separate sources would generate separate
peaks in the cross correlation between the receptor pairs. The possition
of the peaks is given by the time difference in the arrival of the sound
to different receptors, therefore by analyzing all the possible pairs,
it is possible to assign the peaks of the pair A-B to the peaks in B-C
and A-C. One of the algorithms of soundlocalization identifies the
shared peaks that arrive at given times at all the receptors. A second
algorithm assign the shared peaks to source positions. However those
positions depend on a proper sincronization of the recordings, something
that is usually not feasible, so we need a third algorithm to
sinchronize the recordings, that is determine the relative offset of
each receptor. This is possible if we know with precision the position
of one or more sound sources or if we have a good number of them.

The whole idea of soundlocalization consist in optimize everything at
the same time: the offset of each receptor, the precise relative
position of each receptor and sound source and the gain of each receptor
at each time (sometimes the recordings have automatic gain adjustments)

For doing that it is necessary to do three things: 1. Indicate the
possible position of the receptors with some error estimates 2. indicate
the possible time, frequency and position coordinates of some sounds
sources. 3. Organize the recording on different frames defined by time
and frequency intervals (might be also regions of space¿?).

# Indicate the position of receptors with their errors

The idea is to use geocoordinates for the position (including altitude
for height) and UTC for time of the recording. The frequecy is simpler
to understand.

Now we define a list of some identified sound sources

``` r
sources = set_sources(
  label = c(
    'Andrea_voice',
    'Joaquin_voice',
    'Candela_voice',
    'Luciano_voice'
  ),
  x_pos = c(0, 4.5, 9.3, 13.3),
  y_pos = c(0.1, 0.1, 0.1, 0.1),
  error = 2,
  t_recordings = c('Andrea',
                   'Joaquin',
                   'Candela',
                   'Luciano'),
  t_starts = c(123.221, 121.876, 122.556,120.460),
  t_ends = c(124.109, 122.433,123.249,122.969),
  min_freq = 500,
  max_freq = 6000)
sources
#> $label
#> [1] "Andrea_voice"  "Joaquin_voice" "Candela_voice" "Luciano_voice"
#> 
#> $x
#> [1]  0.0  4.5  9.3 13.3
#> 
#> $y
#> [1] 0.1 0.1 0.1 0.1
#> 
#> $error
#> [1] 2
#> 
#> $t_starts
#> [1] 123.221 121.876 122.556 120.460
#> 
#> $t_ends
#> [1] 124.109 122.433 123.249 122.969
#> 
#> $min_freq
#> [1] 500
#> 
#> $max_freq
#> [1] 6000
```

Now, with the information on the sound sources, we can organize the
analysis of the recordings in different frames. Each frame has a start
and end time (according to a given recording)
