import rasterio
import numpy as np

def obtain_files(name, folder): #Name of the file you want to obtain, Folder from where you want to obtain file
    library= {
        "2017" : list((folder).rglob("*.tiff")),
        "2018" : list((folder).rglob("*.tiff")),
        "2019" : list((folder).rglob("*.tiff"))
    }
    files= library[name]
    return files

def readata(file, number):
    with rasterio.open(file[number]) as src:
        readdata = src.read(1)
        profile = src.profile
    return readdata, profile

def veggie(Band1, Band2):
    with np.errstate(divide="ignore", invalid="ignore"):
        vegetation = (Band1 - Band2) / (Band1 + Band2)
        vegetation = np.nan_to_num(vegetation, nan=0.0, posinf=0.0, neginf=0.0)
    return vegetation

    

def calculate_veggie_area(vegyear, pixelarea):
    mask = vegyear>=0.45
    vegetation_pixel_count = np.sum(mask)
    total_vegetation_area_m219 = vegetation_pixel_count * pixelarea
    return total_vegetation_area_m219

def comparition_stats(data):
    "print out stats with nice formatting"
    stats = {
        "Mean": np.mean(data),
        "Median": np.median(data),
        "Std Dev": np.std(data),
        "Min": np.min(data),
        "Max": np.max(data)
    }
    return stats