from astropy.io import fits
import numpy as np
from scipy.optimize import curve_fit



#Function to obtain files from the specific folder: inputs: name of the data and folder #from where you want the data
def obtain_files(name, folder): #Name of the file you want to obtain, Folder from where you want to obtain file
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

#Function to obtain data from the desire file
def obtain_data(name, file):
    if name != "raw" and name!= "blev":
        with fits.open(file) as hdul:
            data = hdul[0].data
    else:
        with fits.open(file) as hdul:
            data = hdul[1].data
    return data

#Calculate exposure time
def exp_time(file_raw):
    with fits.open(file_raw) as hdul:
        exp_time = hdul[0].header.get('EXPTIME',0.0)
    return exp_time

#Function to mask raw data
def masking(data_mask, raw_data): 
    condition = (data_mask == 0) | (data_mask == 256) #remove bad pixels, to provide better data, but not questionable pixels, since it is not necessary and removes a lot of data, creating lines 
    masked_raw = np.where(condition, raw_data, np.nan)
    return masked_raw

#Function to apply the blev data
def apply_blev(data_blev, raw_data,difference):
    if difference == "blev7":
        overscan_region = data_blev[0][9:790, 8:14]
    else:
        overscan_region = data_blev[9:790, 8:14]
    odd = overscan_region[:, 0::2]
    even = overscan_region[:, 1::2]
    mean_odd = np.mean(odd)
    mean_even =  np.mean(even)
    corrected_raw_data = raw_data.copy()
    corrected_raw_data[:, 0::2] = raw_data[:, 0::2] - mean_odd
    corrected_raw_data[:, 1::2] = raw_data[:, 1::2] - mean_even
    return corrected_raw_data

#Function to apply bias data
def apply_bias(data_bias, corrected_raw_data):
    bias_correct = corrected_raw_data-data_bias
    bias_correct
    return bias_correct

#Function to apply dark data
def apply_dark(data_dark, bias_correct, exp_time):
    dark_correct = bias_correct-(data_dark*exp_time)
    dark_correct
    return dark_correct

#Function to apply flat data
def apply_flat(data_flat, dark_correct):
    flat_correct = dark_correct*data_flat
    return flat_correct

#Function to remove the sky
def calculate_mean_and_stndev(flat_correct): # It's okay (Check)
    sum = 0
    sum2 = 0
    for i in flat_correct:
        n = len(flat_correct)
        sum +=i
    mean = sum/n
    for i in flat_correct:
        sum2 += (i - mean)**2
    standev = ((1/(n-1))*(sum2))**0.5
    condition = ((mean-3*standev)> flat_correct) &  (flat_correct<(mean+3*standev))
    sky = np.median(condition)
    sky_corr = flat_correct - sky
    return sky_corr


def flux(exp_time, PHOTOFLAM, skycorr_data):
    count = skycorr_data/exp_time
    image = count * PHOTOFLAM
    return image

def magnitude(fluxdata,PHOTZPT):
    magnitude = -2.5*np.log(fluxdata) + PHOTZPT
    return magnitude

def histogram(ax, x, bins=100, xrange=None, density=False,
              normed=True, histtype='step', **kwargs):
    """
    compute and plot a histogram from some given data, correctly
    handling normalization

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

def normalize_data(x, pmin=None, pmax=None):
    """
    normalize data between 0 and 1, optionally clipping the original data
    between user-specified min/max values, or percentiles

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

    
def transform_by(func, x, pmin=None, pmax=None):
    """
    transform/stretch the data by a named function,
    normalizing the data if necessary
    
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

def gauss_2d(coords, amp, x0, y0, sigx, sigy):
    """ return a 2d gaussian evaluation """
    y, x = coords
    term1 = (x - x0)**2 / (2 * sigx**2)
    term2 = (y - y0)**2 / (2 * sigy**2)
    g = amp * np.exp(-(term1 + term2))
    return g.ravel()


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

def fit_gauss_2d(data, init_par=None):
    """ a fitting routine for a 2d gaussian function """

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

def moment_1(data):
    """ get the intensity weighted centroid of an image """
    y, x = np.indices(data.shape)
    x0 = np.sum(x * data) / data.sum()
    y0 = np.sum(y * data) / data.sum()
    return x0, y0

def moment_2(data):
    """ get the second moment to estimate the peak width """
    y, x = np.indices(data.shape)
    x0, y0 = moment_1(data)
    sigx = np.sqrt(np.sum(((x - x0)**2) * data) / data.sum())
    sigy = np.sqrt(np.sum(((y - y0)**2) * data) / data.sum())
    return sigx, sigy