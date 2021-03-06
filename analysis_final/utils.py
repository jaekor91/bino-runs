from astropy.io import fits # , ascii
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.ndimage.filters import median_filter, uniform_filter
import os
import time
from astropy.coordinates import SkyCoord
from astropy import units as u
deg2arcsec=3600

import matplotlib as mpl
mpl.rcParams['xtick.major.size'] = 15
mpl.rcParams['xtick.major.width'] = 1.
mpl.rcParams['ytick.major.size'] = 15
mpl.rcParams['ytick.major.width'] = 1.
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

def flux2mag(flux):
    return 22.5-2.5*np.log10(flux)    

def mag2flux(mag):
    return 10**(0.4*(22.5-mag))

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
		ibool = np.logical_or.reduce((np.isnan(err_tmp), np.isnan(data_tmp), err_tmp <=0., data_tmp > 10**4))
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

def gauss_fit2profile(K, mu_min=2., mu_max=20., sig_min=1., sig_max=3., dsig=0.05):
	# ---- Grid search for mu and sigma for the best Gaussian representation of the empirical kernel.
	Nrows = 32 
	mu_arr = np.arange(mu_min, mu_max, 0.1)
	sig_arr = np.arange(sig_min, sig_max, dsig)
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


def extraction_kernel_sig(K_collection, all_sigs=False, dsig=0.01, sig_min=0.5, sig_max=5):
	"""
	Based on the extracted stellar profiles, 
	compute a reasonable gaussian extraction kernal
	width (sigma).
	"""
	# Array of estimates gaussian means and widths
	K_gauss_mus = np.zeros(len(K_collection))
	K_gauss_sig = np.zeros(len(K_collection))
	for i in range(len(K_collection)):
		mu_best, sig_best = gauss_fit2profile(K_collection[i], dsig=dsig, sig_min=sig_min, sig_max=sig_max)
		K_gauss_mus[i] = mu_best
		K_gauss_sig[i] = sig_best

	if all_sigs:
		return K_gauss_sig
	else:
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

def wavegrid_from_header(header, Ncols):
	"""
	Construct a linear grid based on the header
	and a user specified number of columns.
	"""
	x0 = header["CRVAL1"] * 10
	dx = header["CDELT1"] * 10
	return x0 + np.arange(0, Ncols, 1.) * dx

def extract_post_stamp(data, err):
	"""
	Given data and err of a single spectrum, extract a 32 x 32 post stamp
	at a random location.
	"""
	num_rows = 32
	idx_min, idx_max = index_edges(data)
	idx_start = np.random.randint(idx_min, idx_max-num_rows/2, size=1)[0]

	# Store the selected stamp
	post_stamp = np.zeros((2, 32, 32))
	post_stamp[0] = data[:, idx_start:idx_start+num_rows]
	post_stamp[1] = err[:, idx_start:idx_start+num_rows]
	
	return post_stamp

def idx_peaks(wavegrid, redz, idx_min=0, idx_max=None):
	"""
	Given a wavelength grid and a redshift, return the indices corresponding to
	the following emission line peaks: OII, Ha, Hb, OIII (1, 2)
	"""
	names = ["OII", "Hb", "OIII1", "OIII2", "Ha"]
	OII = 3727
	Ha = 6563
	Hb = 4861
	OIII1 = 4959
	OIII2 = 5007
	peak_list = [OII, Hb, OIII1, OIII2, Ha]
	
	if idx_max is None:
		idx_max = wavegrid.size-1
	
	# Compute redshifted location
	peak_redshifted_list = []            
	for pk in peak_list:
		peak_redshifted_list.append(pk * (1+redz))
		
	# Compute wavegrid index corresponding to the location. Return -1 if outside the bound.
	index_list = []
	for pk in peak_redshifted_list:
		idx = find_nearest_idx(wavegrid, pk)
		if (idx >=idx_min) and (idx < idx_max):
			index_list.append(idx)
		else:
			index_list.append(-1)
	
	return names, index_list

def find_nearest_idx(arr, x):
	return np.argmin(np.abs(arr-x))


def plot_post_stamps(stamps, num_start=0, fname="test.png", N_stamps2plot = 100, figsize_per_stamp = 3, vmin = -3, vmax = +5):
	"""
	Plot post stamps input in (Nstamps, 32, 32) format.
	"""
	N_stamps_per_row = int(np.ceil(np.sqrt(100)))
	figsize = figsize_per_stamp * N_stamps_per_row
	plt.close()
	fig, ax_list = plt.subplots(N_stamps_per_row, N_stamps_per_row, figsize=(figsize, figsize))
	for i in range(N_stamps2plot):
		idx_row = i // N_stamps_per_row
		idx_col = i % N_stamps_per_row
		ax_list[idx_row, idx_col].imshow(stamps[i], aspect="auto", cmap="gray", interpolation="none")# , vmin=vmin, vmax=vmax)
		ax_list[idx_row, idx_col].axis("off")        
		ax_list[idx_row, idx_col].set_title(i+num_start, fontsize=20)
	plt.savefig(fname, dpi=100, bbox_inches="tight")
	plt.close()

def generalized_gauss_PSF(num_rows, num_cols, x, y, sigma_x, sigma_y, rho=0, num_comps=50, scatter=0.5):
	"""
	Given num_rows x num_cols of an image, generate a generalized PSF
	at location x, y.

	- Rho: covariance element of the 2D covariance matrix with sigma as diagonal std.
	- num_comps: Number of components overwhich to divide up the components.
	- scatter: The scatter in row and column direction. 

	bivariate formula: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	"""
	xv = np.arange(0.5, num_rows)
	yv = np.arange(0.5, num_cols)
	yv, xv = np.meshgrid(xv, yv) # In my convention xv corresponds to rows and yv to columns

	im = np.zeros((32, 32))
	for _ in range(num_comps):
		dx, dy = np.random.randn(2) * scatter
		PSF = np.exp(-( (np.square(xv-x-dx) / sigma_x**2) + (np.square(yv-y-dy) / sigma_y**2) \
					- (2 * rho * (yv-y-dy) * (xv - x-dx) /(sigma_x * sigma_y)) )/ (2 * (1-rho**2))) \
			/ (np.pi * 2 * sigma_x * sigma_y * np.sqrt(1 - rho**2))
		im += PSF
		
	return im/num_comps 

def gen_SN_train_example(data, err, double = True):
	"""
	Given a blank image, generate a blob image by injection.
	"""
	# --- Sample row width from est. distribution
	sig_sig_x = 0.382
	mu_sig_x = 2.
	sig_x = max(np.random.randn() * sig_sig_x + mu_sig_x, 1.5)

	# --- Col width distribution
	sig_sig_y = mu_sig_x * 1.25
	mu_sig_y = mu_sig_x * 1.25
	sig_y = max(np.random.randn() * sig_sig_y + mu_sig_y, 1.5)

	# --- Y up and down scatter
	x = 10 * (np.random.random()-0.5) + 15
	y = 15 + np.random.random() # Fixed at center

	# --- Generate random angle
	sig_rho = 0.2
	rho = np.random.randn() * sig_rho
	if rho > 0.7:
		rho = 0.7
	if rho < -0.7:
		rho = -0.7

	# --- Generate a blob image
	if double:
		sep_min = 3
		sep_max = 10
		sep = np.random.random() * (sep_max - sep_min) + sep_min
		PSF1 = generalized_gauss_PSF(32, 32, x, y-sep/2., sig_x, sig_y, rho=rho)
		PSF2 = generalized_gauss_PSF(32, 32, x, y+sep/2., sig_x, sig_y, rho=rho)    
		PSF = PSF1 + PSF2
	else:
		PSF = generalized_gauss_PSF(32, 32, x, y, sig_x, sig_y, rho=rho)

	fudge_factor = 50 * sig_x * sig_y
	im = data + fudge_factor * np.percentile(data[4:25], 80) * PSF

#     plt.close()
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3))
#     ax1.imshow(data, aspect="auto", cmap="gray", interpolation="none") #, vmin=vmin, vmax=vmax)
#     ax2.imshow(data/err, aspect="auto", cmap="gray", interpolation="none") #, vmin=vmin, vmax=vmax)
#     plt.show()
#     plt.close()

#     plt.close()
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3))
#     ax1.imshow(im, aspect="auto", cmap="gray", interpolation="none") #, vmin=vmin, vmax=vmax)
#     ax2.imshow(im/err, aspect="auto", cmap="gray", interpolation="none") #, vmin=vmin, vmax=vmax)
#     plt.show()
#     plt.close()
	SN = im/err
	return SN

def gen_post_stamp_arr(SN, query):
    """
    Given (32, Ncol) image, generate Nquery post stamps
    centered at where query vector is True.
    """
    Nquery = query.sum()
    stamps = np.zeros((Nquery, 32, 32))
    counter = 0
    i = 0
    while i < query.size:
        if query[i]:
            stamps[counter] = SN[:, i-16:i+16]        
            counter +=1
        i+=1
    return stamps.reshape((Nquery, 32, 32, 1))

def reject_heuristics(stamps):
    reject = np.zeros(stamps.shape[0], dtype=bool)
    for i in range(reject.size):
        SN = stamps[i, 4:25, :, :]
        too_many_zeros = (np.abs(SN) < 1e-20).sum() > 32 * 20 * 0.5 
        too_many_high_SN = (np.abs(SN) > 50).sum() > 32 * 20 * 0.5
        if too_many_zeros or too_many_high_SN:
            reject[i] = True
    return reject

def z_candidates(wavegrid1, wavegrid2, detection1, detection2, z_min = 0., z_max = 2.,\
    idx_min1=0, idx_max1=None, idx_min2=0, idx_max2=None):
    """
    Based on detection/wavegrid arrays, determine candidates.
    """
    z_arr = np.arange(z_min, z_max, 1e-4)
    hits_arr = np.zeros(z_arr.size, dtype=int)

    for i, z in enumerate(z_arr):
        # Count the number of hits from first data
        _, idx_list = idx_peaks(wavegrid1, z, idx_min=idx_min1, idx_max=idx_max1)
        for idx in idx_list:
            if (idx > 10) and (detection1[idx-2:idx+3].sum()>0):
                    hits_arr[i] += 1

        # Second
        _, idx_list = idx_peaks(wavegrid2, z, idx_min=idx_min2, idx_max=idx_max2)
        for idx in idx_list:
            if (idx > 10) and (detection2[idx-2:idx+3].sum()>0):
                    hits_arr[i] += 1    
    # Concentrate on those that have detections
    z_arr = z_arr[hits_arr>0]
    hits_arr = hits_arr[hits_arr>0]

    # ---- Perform cluster reduction
    z_final = [] 
    hits_final = []
    cluster_sizes = []
    # Neighbors are defined as any other zs that are separated within 1.5 * dz =1.5e-4
    i = 0 
    while i < z_arr.size:
        cluster_z = [z_arr[i]]
        cluster_hits = [hits_arr[i]]
        cluster_size = 1 
        j = 1
        while ((i+j) < z_arr.size) and ((z_arr[i+j] - z_arr[i+j-1]) < 1.5e-4):
            cluster_z.append(z_arr[i+j])
            cluster_hits.append(hits_arr[i+j])            
            cluster_size += 1
            j+=1
        # Save the mean of the clusters
        z_final.append(np.mean(cluster_z))
        hits_final.append(np.mean(cluster_hits))        
        cluster_sizes.append(cluster_size)
        
        # Move on to the next cluster
        i += j
        
    zs = np.asarray(z_final)
    num_hits = np.asarray(hits_final)
    cluster_sizes = np.asarray(cluster_sizes)    

    # Do not sort by number of hits    
    # idx_sort = np.argsort(hits_final)[::-1]
    # zs = z_final[idx_sort]
    # num_hits = hits_final[idx_sort]
    return zs, num_hits, cluster_sizes


def crossmatch_cat1_to_cat2(ra1, dec1, ra2, dec2, tol=1./(deg2arcsec+1e-12)):
	"""
	Return indices of cat1 (e.g., DR3) and cat2 (e.g., DEE2) cross matched to tolerance. 

	Note: Function used to cross-match DEEP2 and DR3 catalogs in each field 
	and test for any astrometric discrepancies. That is, for every object in 
	DR3, find the nearest object in DEEP2. For each DEEP2 object matched, 
	pick DR3 object that is the closest. The surviving objects after these 
	matching process are the cross-matched set.
	"""
	
	# Match cat1 to cat2 using astropy functions.
	idx_cat1_to_cat2, d2d = match_cat1_to_cat2(ra1, dec1, ra2, dec2)
	
	# Indicies of unique cat2 objects that were matched.
	cat2matched = np.unique(idx_cat1_to_cat2)
	
	# For each cat2 object matched, pick cat1 object that is the closest. 
	# Skip if the closest objects more than tol distance away.
	idx1 = [] # Place holder for indices
	idx2 = []
	tag = np.arange(ra1.size,dtype=int)
	for e in cat2matched:
		ibool = (idx_cat1_to_cat2==e)
		candidates = tag[ibool]
		dist2candidates = d2d[ibool]
		# Index of the minimum distance cat1 object
		if dist2candidates.min()<tol:
			idx1.append(candidates[np.argmin(dist2candidates)])
			idx2.append(e)
	
	# Turning list of indices into numpy arrays.
	idx1 = np.asarray(idx1)
	idx2 = np.asarray(idx2)
	
	# Return the indices of cat1 and cat2 of cross-matched objects.
	return idx1, idx2



def match_cat1_to_cat2(ra1, dec1, ra2, dec2):
	"""
	"c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
	catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
	idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

	idx are indices into catalog that are the closest objects to each of the coordinates in c, d2d are the on-sky distances between them, and d3d are the 3-dimensional distances." -- astropy documentation.  

	Fore more information: http://docs.astropy.org/en/stable/coordinates/matchsep.html#astropy-coordinates-matching 
	"""    
	cat1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
	cat2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
	idx, d2d, d3d = cat1.match_to_catalog_sky(cat2)
	
	return idx, d2d.degree

def closest_idx(arr, val):
	return np.argmin(np.abs(arr-val))  

def bit_true(bit1, bit2):
    """
    Compare if the two ints have any common number
    """
    return np.bitwise_and(bit1, bit2).sum() > 0

def FDR_cut(grz):
    """
    Given a list [g,r,z] magnitudes, apply the cut and return an indexing boolean vector.
    """
    g,r,z=grz; yrz = (r-z); xgr = (g-r)
    ibool = (r<23.4) & (yrz>.3) & (yrz<1.6) & (xgr < (1.15*yrz)-0.15) & (xgr < (1.6-1.2*yrz))
    return ibool


def emission_line_model(wg, w, z, peak_name="OII"):
    """
    Return a model normalized to 1 evaluated over a grid wg, width w, and redshift z.
    """
    line_dict = {"Ha": 6563, "Hb": 4861, "OIII1": 4959, "OIII2": 5007}
    
    w *= (1+z)

    # Construct a model first
    if peak_name == "OII":                
    # ---- Double gaussian
        mu1 = 3726.1 * (1+z) # Restframe wavelength position
        mu2 = 3728.8 * (1+z)
        model = 0.73 * np.exp(-(wg - mu1)**2 / (2 * w**2)) / (np.sqrt(2 * np.pi) * w)
        model += np.exp(-(wg - mu2)**2 / (2 * w**2)) / (np.sqrt(2 * np.pi) * w)
        model /= 1.73
    else:
        # ---- Single gaussian
        mu = line_dict[peak_name] * (1+z)
        model = np.exp(-(wg - mu)**2 / (2 * w**2)) / (np.sqrt(2 * np.pi) * w)
    return model

def fit_emission_line(wg, spec1D, spec1D_ivar, z_initial=0., peak_name="OII", continuum=0.):
    """
    Given a 1D spectrum and the corresponding spec1D_ivar, estimate the best fitting
    Gaussian model by a brute force calculation.
    - Given wavegrid, wg, compute wavegrid resolution and the corresponding redshift increment dz.
    - Redshift search range is equal to z_initial pm 10 * dz with redshift grid defined over dz/10.
    - The fit parambers are A_best, sig_best, and z_best. Given the last two, A_best can be estimated and
    so can its uncertainty.
    """
    # Subtract the continuum level from the spectrum
    spec1D = np.copy(spec1D) - continuum

    # Construct redshift grid to fit
    dz = (wg[1] - wg[0]) / 3727.
    zrange = np.arange(z_initial-10*dz, z_initial+10*dz, dz/10.)

    # Construct width grid to fit
    widths = np.arange(0.5, 2.5, 0.05)

    # ---- Make a fit     
    # Best is judged by chi-sq (the smaller the better)    
    chi_sq_best = np.infty
    A_best = 0.
    sig_A_best = np.infty
    width_best = 5.
    z_best = z_initial
    SN_best = 0.
    for z in zrange:
        for w in widths:
            model = emission_line_model(wg, w, z, peak_name=peak_name) # Normalized to 1
            invar_A = np.sum(model**2 * spec1D_ivar)
            A = np.sum(spec1D * model * spec1D_ivar) / invar_A
            sig_A = 1./np.sqrt(invar_A)
            chi_sq = np.sum((spec1D - A * model)**2 * spec1D_ivar)
            SN = A/sig_A
            if SN_best < SN:
                A_best = A
                sig_A_best = sig_A
                width_best = w
                z_best = z
                SN_best = SN
                chi_sq_best = chi_sq

    if (chi_sq_best - np.sum(spec1D**2 * spec1D_ivar)) > -10: # Some minimum chi-square difference
        A_best = 0.
        sig_A_best = 999.
        width_best = 999.
        z_best = -999.
    return A_best, sig_A_best, width_best, z_best
