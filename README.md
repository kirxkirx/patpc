### photon arrival time periodicity code (patpc)
This code will search for a period in photon arrival times (as opposed to "countrate as a function of time" = "lightcurve") data in the context of X-ray and gamma-ray astronomy. This is a C implementation of the Hm test [(de Jager, Raubenheimer & Swanepoel, 1989)](https://ui.adsabs.harvard.edu/abs/1989A&A...221..180D). The chance occurance probability for a periodogram peak is calculated as `Prob(> H) = exp (−0.4H)` following equation (5) of [de Jager & Büsching (2010)](https://ui.adsabs.harvard.edu/abs/2010A&A...517L...9D) and corrected for multiple trials using an estimate of the number of independent frequencies as suggested by [Schwarzenberg-Czerny (2003)](https://adsabs.harvard.edu/abs/2003ASPC..292..383S). While writing the code I was following the explanation of the Hm test by [Kerr (2011)](http://adsabs.harvard.edu/abs/2011ApJ...732...38K). The code makes use of multiple processing cores thanks to [OpenMP](https://en.wikipedia.org/wiki/OpenMP).

### Requirements
In order to compile and run the code you'll need:
* A C compiler (tested with `gcc`)
* `GNU make` (you'll need to install it if you are on Mac or something)
* [CFITSIO library](https://heasarc.gsfc.nasa.gov/fitsio/) (optional, needed to read [OGIP FITS event files](https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/events/ogip_94_003/ogip_94_003.html))
* [GNU Scientific Library](https://en.wikipedia.org/wiki/GNU_Scientific_Library) (optional)
* [gnuplot](https://en.wikipedia.org/wiki/Gnuplot) and [ImageMagick](https://en.wikipedia.org/wiki/ImageMagick) `covert` for plotting (optional)

### Installation
Download and unpack the archive or clone this repository, then run `make`.

### Usage
Specify the OGIP FITS event file or a single-column ASCII file lisitng the photon arrival times in seconds as the command line argument:
````
./patpc event_file.evt
````
or
````
./patpc phton_arrival_times.txt
````
This will print the parameters of the highest peak of the Hm periodogram and the power spectrum to the terminal. The code will also try to create plots of the power spectrum, Hm periodogram, binned lighcurve as a function of time, binned lightcurve as a functon of phase (folded with the period corresponding to the highest Hm periodogram peak).

You may also manually specify the search parameters:
````
./patpc event_file.evt 5000 5 0.1
````
where 5000 is the maximum trial period in seconds, 5 is the minimum trial period in seconds, 0.1 is the maximum phase shift between the first and the last point of the lightcurve (determines the width of the step between the trial frequencies).


Bug reports, pull requests and feature suggestions are warmly welcome!
