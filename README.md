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

IMPORTANT
There are three Jupyter files inside the folder called: HSTubai2505m, HSTubai2507m, HSTubai2508m. In each of these files the respective exposures are calibrated and compared with the science-ready images. This was done to prioritise organisation and provide a clean and easy to understand code. 
In addition, a python script named "Funcitions" will also be found in the folder. This script containes all the functions used to calibrate the data and do all the quantitative process and analysis. 

# Part2: HST Hunga Tonga Analysis Sentinel 2
In this part of the assignment, Google Earth images of the Hunga Tonga volcano from years 2014, 2016, 2020, 2021, 2022 and 2023 were provided. The main goal is to plot the land area as a function of time, to analyse the variation due to volcano eruption. However, because of high cloud density the land area in 2020 was not calculated satisfactory.
The Sentinel-2 images of the volcano in 2017, 2018 and 2019 were also provided. In this case the area covered in vegetation had to be calculated and ploted as a function of time.
After obtaining all this information, the results were compared with the vegetation image provided by Sentinel-2 "Sentinel-2_L2A_NDVI"

IMPORTANT
As in Part1, there is a Jupyter notebook with all the processing done in this part of the assignment for both sentinel-2 and google earth. There is also a python script where the all functions used can be found. 
