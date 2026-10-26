from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter, convolve 
from scipy.ndimage import maximum_filter

##THIS SCRIPT CONTAINS THE FUNCTIONS USED FOR THE PROCESSING AND CALIBRATION OF THE HST RAW FILES

def obtain_files(name, folder): 
    """
    Obtain the files from a specific folder

    Parameters:
    -----------
    name : string - name of the file
    folder: path to the folder with the file

    Returns:
    -----------
    files : path to file
    """

    library= {
        "raw" : next(folder.rglob("*_raw.fits")),
        "bias" : next(folder.rglob("*_biasfile.fits")),
        "blev" : next(folder.rglob("*_blevfile.fits")),
        "dark" : next(folder.rglob("*_darkfile.fits")),
        "flat" : next(folder.rglob("*_flatfile.fits")),
        "mask" : next(folder.rglob("*_maskfile.fits")),
    }
    files= library[name]
    return files


def obtain_data(name, file):
    """
    Extract the required data from the files to compute calibration 

    Parameters:
    -----------
    name : string - this is done for the raw and blev files of HSTubai2505m and  
    HSTubai2508m since data is not obtained in hdul[1].data, but in hdul[0].data
    file: path to file

    Returns:
    -----------
    data :  `array` - data from the file
    """
    if name != "raw" and name!= "blev":
        with fits.open(file) as hdul:
            data = hdul[0].data
    else:
        with fits.open(file) as hdul:
            data = hdul[1].data
    return data


def exp_time(file_raw):
    """
    Extract the exposure time from the raw file
    
    Parameters:
    -----------
    file_raw : path to raw file

    Returns:
    -----------
    exp_time : 'float' - exposure time
    """
    with fits.open(file_raw) as hdul:
        exp_time = hdul[0].header.get('EXPTIME',0.0)
    return exp_time


def masking(data_mask, raw_data): 
    """
    Mapping the bad pixels: 0 represents good pixel, 2 bad/defective pixel and 256 questionable pixel. - the questionable pixels are not mapped, because a lot of information is lost, producing bad quality images. It was decided to map just the bad/defective pixels.
    
    Parameters:
    -----------
    data_mask : 'array' mask data
    raw_data : 'array' raw data

    Returns:
    -----------
    masked_raw : 'array' - masked raw data
    """
    condition = (data_mask == 0) | (data_mask == 256) #remove bad pixels, to provide better data, but not questionable pixels, since it is not necessary and it would remove a lot of data from the image, creating vertical lines at the top of the image
    masked_raw = np.where(condition, raw_data, np.nan)
    return masked_raw

    

def apply_blev(data_blev, mask,difference):
    """
    Ensure positive bits to avoid saturation of pixels and loss of data. Use overscan regions of chips, make the mean of the odd and even regions and substract them to the odd and even columns of the masked data.
    
    Parameters:
    -----------
    data_blev : 'array' blev data
    mask : 'array' masked raw data
    difference : string - this is used for the blev data in HSTubai2507m since the shape of the array is shape=(4, 800, 14), so data_blev[0] needs to be applied.

    Returns:
    -----------
    corrected_raw_data : 'array' - data with applied blev file
    """
    if difference == "blev7":
        overscan_region = data_blev[0][9:790, 8:14]
    else:
        overscan_region = data_blev[9:790, 8:14]
    odd = overscan_region[:, 0::2]
    even = overscan_region[:, 1::2]
    mean_odd = np.mean(odd)
    mean_even =  np.mean(even)
    corrected_raw_data = mask.copy()
    corrected_raw_data[:, 0::2] = mask[:, 0::2] - mean_odd
    corrected_raw_data[:, 1::2] = mask[:, 1::2] - mean_even
    return corrected_raw_data


def apply_bias(data_bias, corrected_raw_data):
    """
    Substract the bias, which compensates for readout noise, from the corrected_raw_data
    
    Parameters:
    -----------
    data_bias : 'array' bias data
    corrected_raw_data : 'array' masked raw data with applied blev
    
    Returns:
    -----------
    bias_correct : 'array' - bias corrected image
    """
    bias_correct = corrected_raw_data-data_bias
    return bias_correct


def apply_dark(data_dark, bias_correct, exp_time):
    """
    Multiply the dark data, which accounts for thermal noise, by the eposure time to then substract it from the bias corrected image.
    
    Parameters:
    -----------
    data_dark : 'array' dark data
    bias_correct : 'array' bias corrected image.
    exp_time: 'float' emposure time
    
    Returns:
    -----------
    dark_correct : 'array' - dark or thermal noise corrected image
    """
    dark_correct = bias_correct-(data_dark*exp_time)
    return dark_correct


def apply_flat(data_flat, dark_correct):
    """
    Flat data accounts for pixel sensitivity variation across detector. Multiply the flat data to the dark corrected image.
    
    Parameters:
    -----------
    data_flat : 'array' flat data
    dark_correct : 'array' dark corrected image
    
    Returns:
    -----------
    flat_correct : 'array' - flat corrected image
    """
    flat_correct = dark_correct*data_flat
    return flat_correct


def sky(flat_correct): 
    """
    Estimate the sky contribution to substract it from the flat corrected image. The sky contribution is estimated by calculating the mean and standard deviation, to then make the median of pixels with a value greater than (mean-3*standev) and smaller than (mean+3*standev).
    
    Parameters:
    -----------
    flat_correct : 'array' flat corrected image
    
    Returns:
    -----------
    sky_corr : 'float' - sky corrected image
    """
    sum = 0
    sum2 = 0
    for i in flat_correct:
        n = len(flat_correct)
        sum +=i
    mean = sum/n
    for i in flat_correct:
        sum2 += (i - mean)**2
    standev = ((1/(n-1))*(sum2))**0.5
    condition = ((mean-3*standev)< flat_correct) &  (flat_correct<(mean+3*standev))
    sky = np.median(condition)
    sky_corr = flat_correct - sky
    return sky_corr


def flux(exp_time, PHOTOFLAM, skycorr_data):
    """
    Converts image data from pixel, into energy flux density units.
    
    Parameters:
    -----------
    exp_time : 'float' exposure time
    PHOTOFLAM : 'float' Inverse Sensitivity
    skycorr_data : 'array' sky corrected image
    
    Returns:
    -----------
    image : 'array' - energy flux image
    """
    count = skycorr_data/exp_time #calculate count rate
    image = count * PHOTOFLAM 
    return image

def magnitude(fluxdata,PHOTZPT,ZP_CORR):
    """
    Calculate astronomical magnitude or brightness.

    Parameters:
    -----------
    fluxdata : 'array' flat corrected image
    PHOTZPT : 'float' Zero point
    ZP_CORR : 'float' Delta magnitude correction for contamination
    
    Returns:
    -----------
    sky_corr : 'array' - magnitude
    """
    magnitude = -2.5*np.log10(fluxdata) + PHOTZPT+ZP_CORR
    return magnitude

# FUNCTION OBTAINED FROM HST WORKSHOP 
def histogram(ax, x, bins=100, xrange=None, density=False,
              normed=True, histtype='step', **kwargs):
    """
    compute and plot a histogram from some given data, correctly handling normalization

    Parameters:
    -----------
    ax: matplotlib.axes.Axes instance
        where to plot the data
    x: array
        input values
    bins: int
        number of bins to use
    range: tuple or None
        bin range
    density: bool (default False)
        if True, return a noramlized probability density
    normed: bool, default True
        if True, return a normalized histogram
    histtype: one of 'bar', 'stacked', 'step', 'stepfilled'; default 'step'
        what type of histogram to draw

    Returns:
    ----------
    n, bins, patches:
        see Axes.stairs for more information
    """
    counts, bin_edges = np.histogram(x, bins, range=xrange)
    widths = bin_edges[1:] - bin_edges[:-1]

    if density:
        counts = counts / (widths * x.size)
    elif normed:
        counts = counts / widths
    
    return ax.hist(
        bin_edges[:-1], bin_edges, weights=counts, histtype=histtype, **kwargs
    )

# FUNCTION OBTAINED FROM HST WORKSHOP
def get_local_max(data, nsize=10):
    """ apply a maximum filter to extract coordinates """
    local_max = maximum_filter(data, size=nsize)
    local_max_mask = (local_max == data)
    coordinates = np.argwhere(local_max_mask)
    print(f"The maximum_filter found {len(coordinates)} possible sources")
    return coordinates

# FUNCTION OBTAINED FROM HST WORKSHOP 
def normalize_data(x, pmin=None, pmax=None):
    """
    normalize data between 0 and 1, optionally clipping the original data between user-specified min/max values, or percentiles

    Parameters:
    -----------
    x: np.ndarray
        The data to normalize
    dmin: float or None
        The minimum data value to use for normalization
    dmax: float or None)
        The maximum data value to use for normalization
    percentiles: list, tuple, int, or None
        Instead of using manual dmin/dmax values, use percentiles instead

    Returns:
    -----------
    norm: np.ndarray
        data normalized to [0, 1]
    dmin, dmax: float, float
        min/max values used for normalization
    """
    res = np.asarray(x, dtype=np.float64)
    # quick fix: np.percentile will complain if pmin, pmax are None!
    # if no pmin/pmax values are passed, we'll just assume they want
    # the global min and max values
    if not pmin:
        pmin = 0
    if not pmax:
        pmax = 100
    
    vmin, vmax = np.percentile(res, [pmin, pmax])

    # basic min-max normalization here
    norm = (res - vmin) / (vmax - vmin)
    
    # make sure it's definitely between 0 and 1
    return np.clip(norm, 0, 1)


# FUNCTION OBTAINED FROM HST WORKSHOP 
def transform_by(func, x, pmin=None, pmax=None):
    """
    transform/stretch the data by a named function, normalizing the data if necessary
    
    Parameters:
    -----------
    func: str
        one of 'linear', 'sqrt' or 'log'
    x: np.ndarray
        data to transform
    pmin: int, float, or None
        lower percentile to use for min/max normalization
    pmax: int, float, or None
        upper percentile to use for min/max normalization
    Returns:
    ----------
    out : np.array
        the transformed data
    """

    x_norm = normalize_data(x, pmin=pmin, pmax=pmax)
    if func.lower() == 'sqrt':
        return np.sqrt(x_norm)
    elif func.lower() == 'log':
        scale = 1000
        log_scaled = np.log10(1 + scale + x_norm) / np.log10(1 + scale)
        return log_scaled
    elif func.lower() == 'linear':
        return x_norm
    else:
        raise ValueError("Unsupported function: ", func)

# FUNCTION OBTAINED FROM HST WORKSHOP 
def gauss_2d(coords, amp, x0, y0, sigx, sigy):
    """ return a 2d gaussian evaluation 
    
    Parameters:
    -----------
    coords : `list` coordinates
    amp, x0, y0, sigx, sigy: best fit parameters 

    Returns:
    -----------
    g.ravel() : `np.ndarray` fitted gaussian
    """
    y, x = coords
    term1 = (x - x0)**2 / (2 * sigx**2)
    term2 = (y - y0)**2 / (2 * sigy**2)
    g = amp * np.exp(-(term1 + term2))
    return g.ravel()

# FUNCTION OBTAINED FROM HST WORKSHOP
def square_aperture(data, center, box_size):
    """
    Extract a square cutout from the input data

    Parameters:
    -----------
    data : `np.ndarray`
        the input data
    center : list or tuple
        the center of the box to extract, in y-x order
    box_size : int
        the size of the box in pixels

    Returns:
    -----------
    res : `np.ndarray`
        a cutout of the data
    """

    y, x = map(int, center)

    # ensure the box_size is odd
    if box_size % 2 == 0:
        box_size += 1

    half = box_size // 2

    # use slices
    sly = slice(y - half, y + half + 1)
    slx = slice(x - half, x + half + 1)

    return data[sly, slx]

# FUNCTION OBTAINED FROM HST WORKSHOP
def fit_gauss_2d(data, init_par=None):
    """
    Compute the optimal or best fit parameters from the 2d gaussian fitting

    Parameters:
    -----------
    data : `np.ndarray`
        the input data
    init_par : initial parameters

    Returns:
    -----------
    popt : 'float' best fit or optimal parameters
    pcov: 'float' uncertinties
    """

    res = data.copy()
    res[res < 0] = 0.
    if init_par is None:
        amp = res.max() - res.min()
        x0, y0 = moment_1(res)
        sigx, sigy = moment_2(res)
        init_par = [amp, x0, y0, sigx, sigy]

    y, x = np.indices(res.shape)
    popt, pcov = curve_fit(
        gauss_2d,
        (y, x),
        res.ravel(),
        p0=init_par
    )

    return popt, pcov

# FUNCTION OBTAINED FROM HST WORKSHOP 
def moment_1(data):
    """ get the intensity weighted centroid of an image """
    y, x = np.indices(data.shape)
    x0 = np.sum(x * data) / data.sum()
    y0 = np.sum(y * data) / data.sum()
    return x0, y0

# FUNCTION OBTAINED FROM HST WORKSHOP 
def moment_2(data):
    """ get the second moment to estimate the peak width """
    y, x = np.indices(data.shape)
    x0, y0 = moment_1(data)
    sigx = np.sqrt(np.sum(((x - x0)**2) * data) / data.sum())
    sigy = np.sqrt(np.sum(((y - y0)**2) * data) / data.sum())
    return sigx, sigy

# FUNCTION OBTAINED FROM HST WORKSHOP 
def comparition_stats(data):
    "prints out statistics to compare quantitatively"
    stats = {
        "Mean": np.mean(data),
        "Median": np.median(data),
        "Std Dev": np.std(data),
        "Min": np.min(data),
        "Max": np.max(data)
    }
    return stats