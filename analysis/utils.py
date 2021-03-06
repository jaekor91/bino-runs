from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from scipy.signal import medfilt 
import pandas as pd

# ---- Global variables
SN_thres = 25

def OII_model_simple(x, red_z, A, sig0=0.75, A_ratio = 1, side=None):
	"""
	Simple OII model as double gaussians with the amplitude (in flux)
	equal to A and A * A_ratio
	and sig in wavelength (AA).
	"""
	x0 = 3727.5 # Rest frame peak center [O II] 3726, 3729
	x_center = x0 * (1 + red_z)
	d0 = 3 # The restframe peak distance in AA.
	d = d0 * (1 + red_z)
	sig = (1+red_z) * sig0
	N1 = np.exp(-(x - x_center - d/2.)**2/ (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
	N2 = np.exp(-(x - x_center + d/2.)**2/ (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
	if side is not None:
		if side == 0: 
			return N1
		else:
			return N2
#     print(x_center)
	return A * (N1 + A_ratio * N2)

def find_nearest_idx(arr, x):
	return np.argmin(np.abs(arr-x))


def parse_tanveer_file():
    """
    Return the result of visual inspection completed by Tanveer.
    """

    #--- St82-1hr Side A 2D spectra analyzed by Tanveer.
    tanveer = pd.read_excel("./st82-1-Tanveer.xlsx")

    Nobjs = 64
    z_270, z_600 = np.ones(Nobjs) * -1, np.ones(Nobjs) * -1
    objnums = 64 - np.arange(0, 64, 1)

    counter = 0
    for i, x in tanveer["z"].iteritems():
        if not np.isnan(x):
            z_270[counter] = x
        counter +=1

    counter = 0
    for i, x in tanveer["z.1"].iteritems():
        if not np.isnan(x):
            z_600[counter] = x
        counter +=1


    # ---- Check whether the redshifts agree
    # If unresolved disagreement, then there should be print(.)
    counter = 0
    z_union = np.ones_like(z_600)
    for i in range(Nobjs):
        if (z_600[i] > 0) or (z_270[i] >0):
            counter += 1
            if (z_600[i] > 0) and (z_270[i] >0) and np.abs(z_600[i] - z_270[i]) > 1e-2:
                print(z_270[i], z_600[i])
                
    return z_270, z_600, objnums

def load_train_data_files():
    Nobjs = 64

    #---- Loading in the 2D data for St82-1hr
    data2D_270 = fits.open("/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-270/obj_abs_slits_lin.fits")
    data2D_600 = fits.open("/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-600/obj_abs_slits_lin.fits")
    data2D_270_err = fits.open("/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-270/obj_abs_err_slits_lin.fits")
    data2D_600_err = fits.open("/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-600/obj_abs_err_slits_lin.fits")
   

    #---- Confirm that all spectra are on the same linear grid scale.
    # CRVAL1  =        231.216992188
    # CDELT1  =       0.128999996185

    i = 1
    crval1_270 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
    cdelt1_270 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])

    for i in range(2, Nobjs):
        tmp1 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
        tmp2 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])    
        if np.abs(tmp1-crval1_270)>1e-6:
            print(tmp1)
        if np.abs(tmp2-cdelt1_270)>1e-6:
            print(tmp2)
    wave_grid_270 = crval1_270 + cdelt1_270 * np.arange(data2D_270[1].shape[1])
    wave_grid_270 *= 10

    i = 1
    crval1_600 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
    cdelt1_600 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])

    for i in range(2, Nobjs):
        tmp1 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
        tmp2 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])    
        if np.abs(tmp1-crval1_600)>1e-6:
            print(tmp1)
        if np.abs(tmp2-cdelt1_600)>1e-6:
            print(tmp2)

    wave_grid_600 = crval1_600 + cdelt1_600 * np.arange(data2D_600[1].shape[1])
    wave_grid_600 *= 10

    fname_270 = "/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-270/obj_abs_slits_lin.fits"
    fname_err_270 = "/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-270/obj_abs_err_slits_lin.fits"
    fname_600 = "/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-600/obj_abs_slits_lin.fits"
    fname_err_600 = "/Users/jaehyeon/Documents/Research/binospec/data/v1/st82-1hr-600/obj_abs_err_slits_lin.fits"

    data2D_270, header_270 = bino_data_preprocess(fname_270, fname_err_270)
    data2D_600, header_600 = bino_data_preprocess(fname_600, fname_err_600)

    return data2D_270, data2D_600, header_270, header_600, wave_grid_270, wave_grid_600


def bino_data_preprocess(data_fname, err_fname):
    """
    Import data and corresponding errors and return a list of length (Nobjs+1), where each element
    corresponds to an image of (num_rows, num_cols, 2). +1 is for the first empty HDU.
    num_rows and num_cols can vary.
    
    Note that error is set to infinity wherever NaN appears in the image.

    """
    #---- Loading in the 2D data for St82-1hr
    data2D = fits.open(data_fname)
    data2D_err = fits.open(err_fname)

    #---- Generate the place holder list for numpy array of both data and errors
    Nobjs = len(data2D)-1
    data = [None] #
    header = [None]

    #---- Loop through the imported data. Multiply by a common factor 1e15 to adjust the overall scale.
    for i in range(1, Nobjs+1):
        im = np.copy(data2D[i].data) * 1e15
        err = np.copy(data2D_err[i].data) * 1e15
        num_rows, num_cols = im.shape
        
        # Find pixels in the image or error that has NaN values
        iNaN = np.logical_or(np.isnan(im), np.isnan(err), err==0)
        
        # Set the NaN values equal to zero in the image and error to the infinity.
        im[iNaN] = 0
        err[iNaN] = 1e30
        
        # Properly save
        data_tmp = np.zeros((num_rows, num_cols, 2))
        data_tmp[:, :, 0] = im
        data_tmp[:, :, 1] = err
        
        data.append(data_tmp)
        header.append(data2D[i].header)
        
    return data, header



def post_stamp_from_HDU(HDU, objnum, idx, width=32, row_min=5, row_max=25, m = 20, remove_outlier = False):
    """
    Create a post stamp of size (32, 32) where the relevant spectrum is placed in the middle.
    """
    post_stamp = np.zeros((32, 32))
    row_range = row_max - row_min
    center = 16 # In the post stamp
    row_low = center - row_range // 2
    row_high = center + row_range // 2
    col_low = center - width//2
    col_high = center + width//2
    
    if ((idx-width//2) < 0) or ((idx+width//2) >= HDU[objnum].data.shape[1]):
        return post_stamp
    
    post_stamp[row_low:row_high, col_low:col_high] \
        = np.copy(HDU[objnum].data)[row_min:row_max, idx-width//2:idx+width//2]            
    
    # ---- Perform additional pre-processing
    # Zero-out NaN values
    post_stamp[np.isnan(post_stamp)] = 0

    ibool = np.abs(post_stamp) < 1e-25 # Zero regions

    # if (ibool.sum() / float(width**2)) < 1e-1: # If not even 10 percent has real coverage, then return immediately.
    #     return post_stamp
    # This one is commented out because otherwise it will mess up down stream analysis.
    
    # Eliminate extreme outliers -- within non-zero footprint.
    if remove_outlier:
        post_stamp_strip = post_stamp[~ibool].ravel()
        val_median = np.median(post_stamp_strip) # Calculate median based on non-zero region
        val_std = np.std(post_stamp_strip)
        post_stamp[abs(post_stamp - val_median) > (m * val_std)] = 0
        post_stamp[ibool] = 0 # Set equal to zero originally zero positions.
    
    return post_stamp


def post_stamp_from_imerr_arr(HDU, objnum, idx, width=32, row_min=5, row_max=25, m = 10, remove_outlier = True):
    """
    HDU is a list of spectra formatted.
    Create a post stamp of size (32, 32), one for image and one for error, where the relevant spectrum is placed in the middle.
    """
    post_stamp = np.zeros((32, 32, 2))
    post_stamp[:, :, 1] = 1e30
    row_range = row_max - row_min
    center = 16 # In the post stamp
    row_low = center - row_range // 2
    row_high = center + row_range // 2
    col_low = center - width//2
    col_high = center + width//2
    
    if ((idx-width//2) < 0) or ((idx+width//2) >= HDU[objnum].shape[1]):
        im = post_stamp[:, :, 0]
        err = post_stamp[:, :, 1]
        return im, err
    
    post_stamp[row_low:row_high, col_low:col_high, :] \
        = np.copy(HDU[objnum])[row_min:row_max, idx-width//2:idx+width//2, :]            

    im = post_stamp[:, :, 0]
    err = post_stamp[:, :, 1]
    
    # ---- Perform additional pre-processing    
    # Eliminate extreme outliers -- within non-zero footprint.
    if remove_outlier:
        ibool = np.abs(im) > 1e-20
        if ibool.sum() < 5:
            pass
        else:
            im_strip = im[ibool].ravel()
            val_median = np.median(im_strip) # Calculate median based on non-zero region
            val_std = np.std(im_strip)
            im[abs(im - val_median) > (m * val_std)] = 0
            # im[ibool] = 0 # Set equal to zero originally zero positions.
    
    return im, err


def post_stamp_from_imerr(imerr, idx, width=32, row_min=5, row_max=25, m = 10, remove_outlier = True):
    """
    Create a post stamp of size (32, 32), one for image and one for error, where the relevant spectrum is placed in the middle.
    """
    post_stamp = np.zeros((32, 32, 2))
    post_stamp[:, :, 1] = 1e30
    row_range = row_max - row_min
    center = 16 # In the post stamp
    row_low = center - row_range // 2
    row_high = center + row_range // 2
    col_low = center - width//2
    col_high = center + width//2
    
    if ((idx-width//2) < 0) or ((idx+width//2) >= imerr.shape[1]):
        im = post_stamp[:, :, 0]
        err = post_stamp[:, :, 1]
        return im, err
    
    post_stamp[row_low:row_high, col_low:col_high, :] \
        = np.copy(imerr[row_min:row_max, idx-width//2:idx+width//2, :])

    im = post_stamp[:, :, 0]
    err = post_stamp[:, :, 1]
    
    # ---- Perform additional pre-processing    
    # Eliminate extreme outliers -- within non-zero footprint.
    if remove_outlier:
        ibool = np.abs(im) > 1e-20
        if ibool.sum() < 5:
            pass
        else:
            im_strip = im[ibool].ravel()
            val_median = np.median(im_strip) # Calculate median based on non-zero region
            val_std = np.std(im_strip)
            im[abs(im - val_median) > (m * val_std)] = 0
            # im[ibool] = 0 # Set equal to zero originally zero positions.
    
    return im, err    
     

def idx_peaks(wavegrid, redz):
    """
    Given a wavelength grid and a redshift, return the indices corresponding to
    the following emission line peaks: OII, Ha, Hb, OIII (1, 2)
    """
    names = ["OII", "Ha", "Hb", "OIII1", "OIII2"]
    OII = 3727
    Ha = 6563
    Hb = 4861
    OIII1 = 4959
    OIII2 = 5007
    peak_list = [OII, Ha, Hb, OIII1, OIII2]
    
    # Compute redshifted location
    peak_redshifted_list = []            
    for pk in peak_list:
        peak_redshifted_list.append(pk * (1+redz))
        
    # Compute wavegrid index corresponding to the location. Return -1 if outside the bound.
    index_list = []
    for pk in peak_redshifted_list:
        idx = find_nearest_idx(wavegrid, pk)
        if (idx >=1) and (idx < wavegrid.size-1):
            index_list.append(idx)
        else:
            index_list.append(-1)
    
    return names, index_list




def num_other_matches(detection_list, peak_wavelength, peak_proposed="OII", tol=20):
    OII = 3727
    Ha = 6563    
    Hb = 4861
    OIII1 = 4959
    OIII2 = 5007
    
    if peak_proposed == "OII":
        redz = (peak_wavelength/OII)-1
    elif peak_proposed == "Hb":
        redz = (peak_wavelength/Hb)-1
    elif peak_proposed == "OIII1":
        redz = (peak_wavelength/OIII1)-1
    elif peak_proposed == "OIII2":
        redz = (peak_wavelength/OIII2)-1        
    else: 
        assert False
        
    counter = 0 # Number of hits
    expected_list = [OII * (1+redz), Ha * (1+redz), Hb * (1+redz), OIII1 * (1+redz), OIII2 * (1+redz)]
    # Check if the expected is actually detetcted.
    for x in expected_list:
        # Find the closest wavelength position 
        y = detection_list[find_nearest_idx(detection_list, x)]
        if (np.abs(x-y) < tol):
            counter += 1
#             print(x,y)
        
    return counter

def header2info(header):
    """
    Given a header return the following limited information
    in dictionary format.
    RA, DEC, SLITOBJ (bitcode)
    
    header must be a string.
    """
    RA = float(header.split("RA")[1].split("=")[1].split("/")[0])
    DEC = float(header.split("DEC")[1].split("=")[1].split("/")[0])
    BIT = (int(header.split("SLITOBJ")[1].split("=")[1].split("/")[0][2:-2]))
#     info = {"RA": RA, "DEC": DEC, "BITNUM": BIT}
    
    return RA, DEC, BIT
    


def generalized_gauss_PSF(num_rows, num_cols, x, y, FWHM, rho=0, num_comps=10, scatter = 0):
    """
    Given num_rows x num_cols of an image, generate a generalized PSF
    at location x, y.
    
    - Rho: covariance element of the 2D covariance matrix with sigma as diagonal std.
    - num_comps: Number of components overwhich to divide up the components.
    - scatter: The scatter in row and column direction.
    
    bivariate formula: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    """
    sigma = FWHM / 2.354
    im = np.zeros((num_rows, num_cols))
    xv = np.arange(0.5, num_rows)
    yv = np.arange(0.5, num_cols)
    yv, xv = np.meshgrid(xv, yv) # In my convention xv corresponds to rows and yv to columns
    
    for _ in range(num_comps):
        dx, dy = np.random.randn(2) * scatter
        dy *= 2
        PSF = np.exp(-(np.square(xv-x-dx) + np.square(yv-y-dy) - 2 * rho * (yv-y-dy) * (xv - x -dx))/(2*sigma**2 * (1-rho**2))) \
            /(np.pi * 2 * sigma**2 * np.sqrt(1 - rho**2))

    return PSF / num_comps    

def poisson_realization(D0):
    """
    Given a truth image D0, make a Poisson realization of it.
    """
    D = np.zeros_like(D0)
    for i in range(D0.shape[0]):
        for j in range(D0.shape[1]):
            D[i, j] = np.random.poisson(lam=D0[i, j], size=1)
    return D

def wave_grid_extractor(header, grid_length):
    """
    Given the header output of the function bino_data_preprocess, 
    return the corresponding linear grid.
    """

    i = 1
    crval1 = float(str(header[i]).split("CRVAL1")[1].split("=")[1].split("/")[0])
    cdelt1 = float(str(header[i]).split("CRVAL1")[1].split("=")[2].split("/")[0])

    for i in range(2, len(header)):
        tmp1 = float(str(header[i]).split("CRVAL1")[1].split("=")[1].split("/")[0])
        tmp2 = float(str(header[i]).split("CRVAL1")[1].split("=")[2].split("/")[0])    
        if np.abs(tmp1-crval1)>1e-6:
            print(tmp1)
        if np.abs(tmp2-cdelt1)>1e-6:
            print(tmp2)

    wave_grid = crval1 + cdelt1 * np.arange(grid_length)
    wave_grid *= 10
    
    return wave_grid