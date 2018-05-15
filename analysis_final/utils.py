from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage.filters import median_filter
import os

def preprocess_bino(fname_data, fname_err, data_dir):
	"""
	Preprocessor goals.
	Data: 
		- If NaN: Data ZERO and Error infinity.
	Error:
		- If NaN: Error infinity and Data ZERO.
	Output: 
		- A numpy array of shape (Ntargets+1, 2, 32, num_cols). 
			- Though the native data have different number of rows, we use a single fixed number here.
			- Channel 0: Data
			- Channel 1: Error
		- List of headers

	The native data unit is ergs/cm^2/s/nm. Preprocessor changes this to
	10^-17 ergs/cm^2/s/Angstrom.

	First spectrum in native data is saved in loc "1". We follow the same convention.
	"""
	infinity = 1e60
	unit_conversion = 10**18

	# ---- Output variables
	data_err = None
	list_headers = [None]

	# ---- Python list of spectral data/err
	# Zeroth element is a blank.
	data = fits.open(data_dir + fname_data)
	err = fits.open(data_dir + fname_err)

	# ---- Place holder for the output array
	Ncols = data[1].data.shape[1]
	data_err = np.zeros((len(data), 2, 32, Ncols))
	data_err[:, 1, :, :] = infinity  # All data/errors are initially set to zero and infinity.

	for i in range(1, len(data)):
		# ---- Import data
		data_tmp = data[i].data * unit_conversion
		err_tmp = err[i].data * unit_conversion

		# ---- Apply preprocessing
		ibool = np.logical_or(np.isnan(err_tmp), np.isnan(data_tmp), err_tmp <=0.)
		data_tmp[ibool] = 0
		err_tmp[ibool] = infinity

		# ---- Trim the data
		idx_min, idx_max = index_edges(data_tmp)
		L_trim = 50
		data_tmp[:, :idx_min+L_trim] = 0
		data_tmp[:, idx_max-L_trim:] = 0
		err_tmp[:, :idx_min+L_trim] = infinity
		err_tmp[:, idx_max-L_trim:] = infinity


		# ---- Save data
		# Nrows = min(32, data_tmp.shape[0])
		# Data is usually crap outside this range
		data_err[i, 0, 4:25] = data_tmp[4:25]
		data_err[i, 1, 4:25] = err_tmp[4:25]

		# ---- Save header
		header_tmp = data[i].header
		list_headers.append(header_tmp)

	# # ---- Plot for bebugging 
	# vmin = -4
	# vmax = 4
	# plt.close()
	# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 5))
	# # Plot data
	# data_plot = ax1.imshow(data_tmp, aspect="auto", cmap="gray", interpolation="none", vmin=vmin, vmax=vmax)
	# plt.colorbar(data_plot, ax = ax1)
	# # Plot err
	# err_plot = ax2.imshow(err_tmp, aspect="auto", cmap="gray", interpolation="none", vmin=0.02, vmax=0.05)
	# plt.colorbar(err_plot, ax = ax2)

	# plt.show()
	# plt.close()        
		
	return data_err, list_headers

def bit_from_header(header):
	name = header["SLITOBJ"]
	if name == "stars":
		name = 2**1
	elif name == "gal":
		name = 2**2
	return int(name)

def extract_single_data(data_err, list_headers, specnum):
	"""
	Extract single spectrum data, err, header from the list, arr provided by Preprocessor.
	- specnum: Target object number in a file. Ranges from 1 through approx. 140.
	"""
	header = list_headers[specnum]
	data = data_err[specnum, 0]
	err = data_err[specnum, 1]
	
	return data, err, header

def ivar_from_err(err):
	return 1./np.square(err)

def naive_profile(data, ivar, idx_min=0, idx_max=-1, L_trim = 500):
	"""
	Assumes D_ij ~ Norm(f_j * k_i, sig_ij) where f_j = f.
	"""
	if L_trim > 0:
		data = data[:, idx_min+L_trim:idx_max-L_trim]
		ivar = ivar[:, idx_min+L_trim:idx_max-L_trim]
	K = np.sum(data * ivar, axis = 1) / np.sum(ivar, axis = 1)
	K /= np.sum(K) # Normalization step.
	return K 

def produce_spec1D(data_err, list_headers, sig_K):
	"""
	Given 2D spectrum and the extraction kernel width sig_K,
	produce 1D spectra (Ntargets+1, 2, Ncols) and their inverse variance.
	"""
	data_ivar_1D = np.zeros((data_err.shape[0], 2, data_err.shape[3]))
	for specnum in range(1, len(list_headers)):
		data, err, header = extract_single_data(data_err, list_headers, specnum)
		ivar = ivar_from_err(err)

		spec1D_ivar = np.sum(np.square(K_T) * ivar, axis=0)
		spec1D = np.sum(K_T * data * ivar, axis=0) / spec1D_ivar
		
		data_ivar_1D[specnum, 0] = spec1D
		data_ivar_1D[specnum, 1] = spec1D_ivar
	return data_ivar_1D
		
def SIDE_from_header(header):
	return header["SIDE"]

def index_edges(data, num_thres=20):
	"""
	Given long postage stamp of data, return the edges.
	"""
	idx_min = 0
	idx_max = data.shape[1]-1
	tally = np.sum(data == 0., axis=0)
	while tally[idx_min] > num_thres:
		idx_min += 1
	while tally[idx_max] > num_thres:
		idx_max -=1
	return idx_min, idx_max

def gauss_fit2profile(K):
	# ---- Grid search for mu and sigma for the best Gaussian representation of the empirical kernel.
	Nrows = 32 
	mu_arr = np.arange(5, 20, 0.1)
	sig_arr = np.arange(1., 3., 0.05)
	chi_arr = np.zeros((mu_arr.size, sig_arr.size))
	x_arr = np.arange(0, Nrows, 1)
	for i, mu in enumerate(mu_arr):
		for j, sig in enumerate(sig_arr):
			A = np.exp(-np.square(x_arr-mu) / (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
			chi_arr[i, j] = np.sum(np.square(K - A))
	# ---- Best fit
	idx = np.unravel_index(np.argmin(chi_arr), chi_arr.shape)
	mu_best = mu_arr[idx[0]] 
	sig_best = sig_arr[idx[1]]
	
	return mu_best, sig_best

def extract_stellar_profiles(data_err, list_headers):
	K_collection = []
	for specnum in range(1, len(list_headers)):
		data, err, header = extract_single_data(data_err, list_headers, specnum)
		ivar = ivar_from_err(err)

		BIT = bit_from_header(header)
		# ---- Perform optimal extraction 1D spectrum from 2D
		if (BIT == 2):
			idx_min, idx_max = index_edges(data)
			K = naive_profile(data, ivar, idx_min, idx_max, L_trim = 500)
			K_collection.append(K) # Collect K into a list.
	return K_collection
		
def remove_outlier(arr, std_thres = 2):
	"""
	Remove outliers in 1D array by sigma clipping.
	"""
	std = np.std(arr)
	mu = np.median(arr)
	
	return arr[(arr - mu) < (std_thres * std)]


def extraction_kernel_sig(K_collection):
	"""
	Based on the extracted stellar profiles, 
	compute a reasonable gaussian extraction kernal
	width (sigma).
	"""
	# Array of estimates gaussian means and widths
	K_gauss_mus = np.zeros(len(K_collection))
	K_gauss_sig = np.zeros(len(K_collection))
	for i in range(len(K_collection)):
		mu_best, sig_best = gauss_fit2profile(K_collection[i])    
		K_gauss_mus[i] = mu_best
		K_gauss_sig[i] = sig_best

	return np.median(K_gauss_sig)

def K_gauss_profile(mu, sig, Nrows = 32):
	"""
	Generate gaussian extraction profile of length Nrows
	given mu and sig.
	"""
	x_arr = np.arange(0, Nrows, 1)
	K_gauss = np.exp(-(x_arr - mu)**2 / (2 * sig**2))
	K_gauss /= np.sum(K_gauss)
	
	return K_gauss

def plot_kernels(K_collection, K_extract, fname):
	"""
	Plot the collection of stellar kernels and the ultimate
	extraction kernel at the center.
	"""

	fig, ax = plt.subplots(1, figsize=(10, 5))
	for i in range(len(K_collection)):
		ax.plot(K_collection[i], c="red", lw=0.5)    
	ax.plot(K_extract, c="blue", lw=1.5)
	ax.set_ylim([-0.03, 0.3])
	ax.axhline(y=0, c="black", ls="--", lw=1.)
	plt.savefig(fname, dpi=200, bbox_inches="tight")
#     plt.show()
	plt.close()

	return


def K_gauss_profile(mu, sig, Nrows = 32):
	"""
	Generate gaussian extraction profile of length Nrows
	given mu and sig.
	"""
	x_arr = np.arange(0, Nrows, 1)
	K_gauss = np.exp(-(x_arr - mu)**2 / (2 * sig**2))
	K_gauss /= np.sum(K_gauss)
	
	return K_gauss


def produce_spec1D(data_err, list_headers, sig_K, fname_prefix=None, verbose=True):
    """
    Given 2D spectrum and the extraction kernel width sig_K,
    produce 1D spectra (Ntargets+1, 2, Ncols) and their inverse variance.
    """
    Ncols = 32    
    data_ivar_1D = np.zeros((data_err.shape[0], 2, data_err.shape[3]))

    for specnum in range(1, len(list_headers)):
        if verbose and ((specnum % 10) == 0):
            print("Processing spec num: %d" % specnum)

        # ---- Extract the individual data
        data, err, header = extract_single_data(data_err, list_headers, specnum)
        ivar = ivar_from_err(err)

        # ---- Compute the center of the extraction
        idx_min, idx_max = index_edges(data)

        # ---- Algorithmically determine the row location of the spectra.
        # Note that I assume the center of the spectrum falls between 10 and 20.
        row_centers = []
        data_tmp = np.copy(data)
        ivar_tmp = np.copy(ivar)
        data_tmp[:10, :] = 0.
        data_tmp[20:, :] = 0.        
        ivar_tmp[:10, :] = 1e-120
        ivar_tmp[20:, :] = 1e-120        

        for idx in range(idx_min, idx_max-Ncols, Ncols//2):
            # Compute naive profile based on clipped 2D (32, 32) post stamps
            K = naive_profile(data_tmp[:, idx:idx+Ncols], ivar_tmp[:, idx:idx+Ncols], L_trim=-1)
            # Median filtering to reduce noise
            K = median_filter(K, size=5)
    #             # Savgol filtering of the naive profile         
    #             K_filtered = savgol_filter(K, window_length=9, polyorder=3)
            # Compute the center
            row_centers.append(np.argmax(K))

        # Compute the extraction profile using the above computed center and using extraction width.
        row_centers = np.asarray(row_centers)
        row_centers = row_centers[(row_centers > 9) & (row_centers < 21)]
        mu = np.round(np.median(row_centers))
        K_T = K_gauss_profile(mu, sig_K).reshape((K.size, 1))

        # ---- 1D extraction performed here
        spec1D_ivar = np.sum(np.square(K_T) * ivar, axis=0)
        spec1D = np.sum(K_T * data * ivar, axis=0) / spec1D_ivar

        # ---- Save the extracted spectrum
        data_ivar_1D[specnum, 0] = spec1D
        data_ivar_1D[specnum, 1] = spec1D_ivar        

        if fname_prefix is not None:
            plt.close()
            # ---- Spec figures
            fname = fname_prefix + "spec%d-2D.png" %specnum
            fig, ax = plt.subplots(1, figsize=(17, 1))
            ax.imshow(data, aspect="auto", cmap="gray", interpolation="none", vmin=-0.5, vmax=0.5)
            ax.axhline(y=mu+0.5, c="red", ls="--", lw=0.4)
            plt.savefig(fname, dpi=200, bbox_inches="tight")
            plt.close()

            # ---- Histogram of centers determined
            fname = fname_prefix + "spec%d-centers.png" %specnum
            fig, ax = plt.subplots(1, figsize=(7, 3))
            ax.hist(row_centers, bins=np.arange(0.5, 32.5, 1), histtype="step", color="black", normed=True)
            ax.plot(K_T, c="red", label="K_stellar")
            ax.axvline(x=mu, c="red", ls="--", lw=0.4)
            plt.savefig(fname, dpi=200, bbox_inches="tight")
            plt.close()

    return data_ivar_1D

