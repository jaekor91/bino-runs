from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from scipy.signal import medfilt 
import pandas as pd

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

    return data2D_270, data2D_270_err, data2D_600, data2D_600_err, wave_grid_270, wave_grid_600



def post_stamp_from_HDU(HDU, objnum, idx, width=32, row_min=5, row_max=25, m = 20):
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
    
    # # Eliminate extreme outliers
    # post_stamp[abs(post_stamp - np.median(post_stamp)) > m * np.std(post_stamp)] = 0            
    
    return post_stamp

def idx_peaks(wavegrid, redz):
    """
    Given a wavelength grid and a redshift, return the indices corresponding to
    the following emission line peaks: OII, Ha, Hb, OIII (1, 2)
    """
    names = ["OII", "Ha", "Hb", "OIII1", "OIII3"]
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