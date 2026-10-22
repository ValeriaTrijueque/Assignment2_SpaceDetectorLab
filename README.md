# Assignment2_SpaceDetectorLab
This assignment is divided in two parts: Part1(HST Calibration Routine), Part2 (Hunga Tonga Analysis Sentinel 2)

# Part1: HST Calibration Routine
In this part of the assignment, raw data from the Hubble Space Telescope are to be calibrated. The processing needs to be done with the aid of some calibration files provided in the assignment:

- MASKFILE
- BLEVFILE
- BIASFILE
- DARKFILE
- FLATFILE

The following formula is needed to apply the calibration:

stars + noise = (raw - bias - dark) / flat - sky

After that, the calibration files are to be compared against the science-ready files provided in the assignment using quantitative analysis to see how accurate the files have been calibrated compared to the reference.
In this case 

# Part2: HST Hunga Tonga Analysis Sentinel 2
In this part of the assignment, Google Earth images of the Hunga Tonga volcano in ... where used to make a plot of the land area as a function of time and then Sentinel-2 images of the same volcano to make plots of the land area covered in vegetation as a function of time. The pixel size had to be determined both in the Google Earth and Sentinel-2 images to calculate the area required.

After that the spectral index images were compared with the Sentinel-2 data reduction pipeline.