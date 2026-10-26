Valeria Trijueque Cerda

# Assignment2_SpaceDetectorLab
This assignment is divided in two parts: Part1(HST Calibration Routine), Part2 (Hunga Tonga Analysis Sentinel 2)

# Part1: HST Calibration Routine
The purpose of this part of the assignment is to calibrate raw data from three exposure of the Hubble Space Telescope. This process needs to be done with the aid of some calibration files provided in the assignment:

- MASKFILE
This file provides information about the bad pixels, so they can be masked from the raw data. 0 represents good pixels, 2 defective pixels and 256 questionable pixels. Therefore, a mask was performed to remove the defective pixels from the raw data. The questionable pixels were not removed, because some information was getting lost and a bad quality image produced.

- BLEVFILE
This file is applied to the mask data to ensure the values stored in the pixels are positive and, in this way, avoiding saturation of the pixels and losing data. This was accomplished by making an overscan region of the pixels that are always covered and then apply this to the masked data by calculating separately the mean of the odd and even overscan regions and subtracting it from the odd and even columns of the masked image.

- BIASFILE
The data in the bias file compensates for the readout noise. Consequently, it is subtracted from the image being calibrated.

- DARKFILE
The data in the dark file accounts for the thermal noise. Consequently, it is multiplied by the exposure time and subtracted from the bias corrected image.

- FLATFILE
The data in the flat file accounts for the pixel sensitivity, which is multiplied by the dark corrected image.

- Once all these files were applied to the raw image for calibration, an estimation of the sky had to be computed to also remove it from the flat corrected image. The sky estimate was obtained by calculating the mean and standard deviation of the flat corrected image, and then getting the median of pixels with a value greater than (mean-3*standev) and smaller than (mean+3*standev).

-To finish, photometry calibration was performed with PHOTFLAM, PHOTZPT, and ZP_CORR. The energy flux density was first calculated by multiplying the count rate by the PHOTOFLAM (Inverse Sensitivity), and the astronomical magnitude by performing the following calculation:
magnitude = -2.5*np.log10(fluxdata) + PHOTZPT + ZP_CORR

After that, the computed calibration files are compared against the science-ready files provided in the assignment. Some quantitative analysis was used to see how accurate the files have been calibrated with respect to the reference. 
The comparition was achieved by using histograms, calculating possible sources from both, fitting a 2d gaussian to three astronomical sources and performing some statistical analysis.


IMPORTANT
There are three Jupyter notebooks inside the folder called: HSTubai2505m, HSTubai2507m, HSTubai2508m. In each of these notebooks the respective exposures from Huble Space Telescope are calibrated and compared with the science-ready images. (All exposures files can be found inside the "calibration" folder of the hst-data, and the science-ready files are found inside the "F336W" folder). This was done to prioritise organisation and provide a clean and easy to understand code. 
In addition, a python script named "Functions" can also be found in the same place as the Jupyter notebooks. This script contains all the functions used to calibrate the data and do all the quantitative process and analysis.


# Part2: HST Hunga Tonga Analysis Sentinel 2
In this part of the assignment the Hunga Tonga volcano was analysed using sentinel-2 and Google Earth images:

SENTINEL-2:
The Sentinel-2 provides images of the volcano in 2017, 2018 and 2019 years. From these images the area covered in vegetation was calculated for each year and plotted in a graph as a function of time. Given that no volcano eruption occurred during those 3 years, sentinel-2 images are used to analyse the growth in vegetation during that time. This was achieved by applying a mask of the vegetation to the images and calculating the pixel size. 

The pixel size is first converted into meters from angular units. Once the width and height of the pixels are obtained, the area of a pixel can be calculated. Knowing the pixel size and the pixels containing the vegetation information, the area covered in vegetation can be obtained by multiplying both values together.

After obtaining all this information, the results were compared with the vegetation image provided by Sentinel-2 "Sentinel-2_L2A_NDVI", performing some statistical analysis of the mean, median, standard deviation and minimum and maximum of both.

GOOGLE EARTH:
The Google Earth images show the evolution of Hunga Tonga volcano in 2014, 2016, 2020, 2021, 2022 and 2023 years. The main goal is to plot the volcano's land area as a function of time, to analyse its variation due to volcano eruption. However, because of high cloud density, the land area in 2020 was not calculated satisfactory and therefore it is rejected from the graph.

To obtain the land area, a mask of the land was performed. The mask was achieved by converting the rgb images into hsv images and then applying the mask. As in Sentinel-2, the pixel size was calculated but, in this case, using the scale bar provided in the Google Earth images. Once the pixel size and the pixels containing the land information were acquired, both are multiplied to get the land area at specific year.

IMPORTANT
As in Part1, there is a Jupyter notebook with all the processing done in this part of the assignment for both sentinel-2 (vegetation analysis) and Google Earth (land area analysis). There is also a python script named "Sentinel" where all functions used can be found. The sentinel-2 files for each year are obtained from the folders "S2_2017-07-12", "S2_2018-07-12", "S2_2019-07-12" in the sentinel-data, and the google earth images for each year are obtained from the folder "White_light" in the sentinel-data.



