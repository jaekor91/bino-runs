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
#     print x_center
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
    # If unresolved disagreement, then there should be prints.
    counter = 0
    z_union = np.ones_like(z_600)
    for i in xrange(Nobjs):
        if (z_600[i] > 0) or (z_270[i] >0):
            counter += 1
            if (z_600[i] > 0) and (z_270[i] >0) and np.abs(z_600[i] - z_270[i]) > 1e-2:
                print z_270[i], z_600[i]
                
    return z_270, z_600, objnums

def load_train_data_files():
    Nobjs = 64
	
    #---- Loading in the 2D data for St82-1hr
    data2D_270 = fits.open("../../data/from_vulcan/st82-1hr-270/obj_counts_slits_lin.fits")
    data2D_600 = fits.open("../../data/from_vulcan/st82-1hr-600/obj_counts_slits_lin.fits")
    data1D_270 = fits.open("../../data/from_vulcan/st82-1hr-270/obj_counts_slits_extr.fits")[1].data
    data1D_600 = fits.open("../../data/from_vulcan/st82-1hr-600/obj_counts_slits_extr.fits")[1].data
    
    #---- Confirm that all spectra are on the same linear grid scale.
    # CRVAL1  =        231.216992188
    # CDELT1  =       0.128999996185

    i = 1
    crval1_270 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
    cdelt1_270 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])

    for i in xrange(2, Nobjs):
        tmp1 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
        tmp2 = float(str(data2D_270[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])    
        if np.abs(tmp1-crval1_270)>1e-6:
            print tmp1
        if np.abs(tmp2-cdelt1_270)>1e-6:
            print tmp2
    wave_grid_270 = crval1_270 + cdelt1_270 * np.arange(data2D_270[1].shape[1])
    wave_grid_270 *= 10

    i = 1
    crval1_600 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
    cdelt1_600 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])

    for i in xrange(2, Nobjs):
        tmp1 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[1].split("/")[0])
        tmp2 = float(str(data2D_600[i].header).split("CRVAL1")[1].split("=")[2].split("/")[0])    
        if np.abs(tmp1-crval1_600)>1e-6:
            print tmp1
        if np.abs(tmp2-cdelt1_600)>1e-6:
            print tmp2

    wave_grid_600 = crval1_600 + cdelt1_600 * np.arange(data2D_600[1].shape[1])
    wave_grid_600 *= 10

    return data2D_270, data2D_600, data1D_270, data1D_600, wave_grid_270, wave_grid_600



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
    
    post_stamp[row_low:row_high, col_low:col_high] \
        = np.copy(HDU[objnum].data)[row_min:row_max, idx-width//2:idx+width//2]            
    
    # ---- Perform additiona pre-processing
    # Zero-out NaN values
    post_stamp[np.isnan(post_stamp)] = 0
    
    # Eliminate extreme outliers
    post_stamp[abs(post_stamp - np.median(post_stamp)) > m * np.std(post_stamp)] = 0            
    
    return post_stamp
