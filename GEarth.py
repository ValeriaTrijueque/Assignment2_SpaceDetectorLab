from matplotlib.colors import rgb_to_hsv
import matplotlib.image as mpimg
import numpy as np

def hsv_image(image):
    rgb_image = mpimg.imread(image)
    rgb_norm = rgb_image / 255
    hsv_image = rgb_to_hsv(rgb_norm)
    return hsv_image, rgb_image

def calculate_land_area(image, area_pixel):
    mask = image < 0.45
    pixel_count = np.sum(mask)
    land_area = pixel_count * area_pixel
    return land_area