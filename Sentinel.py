import rasterio
import numpy as np
from matplotlib.colors import rgb_to_hsv
import matplotlib.image as mpimg

##THIS SCRIPT CONTAINS THE FUNCTIONS USED FOR THE ANALYSIS AND PROCESSING OF THE SENTINEL-2 AND GOOGLE EARTH FILES FOR HUNGA TONGA

def obtain_files(name, folder): 
    """
    Obtain the files from a specific folder

    Parameters:
    -----------
    name : string - name of the file
    folder: path to the folder with the file

    Returns:
    -----------
    files : 'list' - list of file paths
    """
    library= {
        "2017" : list((folder).rglob("*.tiff")),
        "2018" : list((folder).rglob("*.tiff")),
        "2019" : list((folder).rglob("*.tiff"))
    }
    files= library[name]
    return files

def readata(file, number):
    """
    Read the data in the file

    Parameters:
    -----------
    file : 'array' - file
    number: number of the file from the list you can to read the data 

    Returns:
    -----------
    readdata : 'array' data
    """
    with rasterio.open(file[number]) as src:
        readdata = src.read(1)
        profile = src.profile
    return readdata

def veggie(Band1, Band2):
    """
    take two bands and make a normalized difference image - vegetation with band 4 and 8

    Parameters:
    -----------
    Band1 : 'array' 
    Band2: 'array'

    Returns:
    -----------
    vegetation : 'array' normalized difference image of vegetation
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        vegetation = (Band1 - Band2) / (Band1 + Band2)
        vegetation = np.nan_to_num(vegetation, nan=0.0, posinf=0.0, neginf=0.0)
    return vegetation

    
def calculate_veggie_area(vegyear, pixelarea):
    """
    Calculate the area covered in vegetation in m^2

    Parameters:
    -----------
    vegyear : 'array'  vegetation in a specific year
    pixelarea: 'float' pixel size in m^2

    Returns:
    -----------
    total_vegetation_area_m2 : 'float' area covered in vegetation in m^2
    """
    mask = vegyear>=0.45
    vegetation_pixel_count = np.sum(mask)
    total_vegetation_area_m2 = vegetation_pixel_count * pixelarea
    return total_vegetation_area_m2

# FUNCTION OBTAINED FROM HST WORKSHOP 
def comparition_stats(data):
    "prints out statistics to compare quantitatively the vegetation image provided from sentinel and the vegetation computed"
    stats = {
        "Mean": np.mean(data),
        "Median": np.median(data),
        "Std Dev": np.std(data),
        "Min": np.min(data),
        "Max": np.max(data)
    }
    return stats

# FUNCTION OBTAINED FROM HST WORKSHOP 
def hsv_image(image):
    "converts rgb(red, green, blue) image array into hsv image array(hue, saturation, value)"
    rgb_image = mpimg.imread(image)
    rgb_norm = rgb_image / 255
    hsv_image = rgb_to_hsv(rgb_norm)
    return hsv_image, rgb_image

def calculate_land_area(image, area_pixel):
    """
    Calculate land area in m^2

    Parameters:
    -----------
    image : 'array'  
    area_pixel: 'float' pixel size in m^2

    Returns:
    -----------
    land_area : 'float' land area in m^2
    """
    mask = image < 0.45
    pixel_count = np.sum(mask)
    land_area = pixel_count * area_pixel
    return land_area