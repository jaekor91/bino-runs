from utils import *

# ----- Image generator functions
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
        im += np.exp(-(np.square(xv-x-dx) + np.square(yv-y-dy) - 2 * rho * (yv-y-dy) * (xv - x -dx))/(2*sigma**2 * (1-rho**2))) \
            /(np.pi * 2 * sigma**2 * np.sqrt(1 - rho**2))

    return im / float(num_comps)

def blob_im_generator(nrows=32, ncols=32, double=False, fdensity=0):
    """
    Generate an emission line peak blob.
    """
    FWHM_min = 3
    FWHM_max = 6
    FWHM = (FWHM_max - FWHM_min) * np.random.random() + FWHM_min # FWHM selection    
    fmin = 12 * FWHM**2 * fdensity * (7 + 3 * (np.random.random()))
    fmax = fmin * 1.01
    rho_min = -0.6
    rho_max = 0.6
    num_comps_min = 10
    num_comps_max = 25 
    # For double peaks only
    sep_peaks_min = 3
    sep_peaks_max = 10
    scatter_max = .75
    scatter_min = 0
    
    im = np.zeros((nrows, ncols))
    x = nrows//2 + (np.random.random() - 0.5) * 8
    y = ncols//2 + (np.random.random() - 0.5) * 1
    scatter = np.random.random() * scatter_max        
    f = (fmax - fmin) * np.random.random() + fmin # Random flux selection
    f *= (1 + scatter)**2
    rho = (rho_max - rho_min) * np.random.random() + rho_min # Covarince selection
    num_comps = np.random.randint(num_comps_min, num_comps_max, 1)[0]
    if not double:
        peak = f * generalized_gauss_PSF(nrows, ncols, x,  y, FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        im += peak # Add peak
    else:
        sep_peaks = np.random.random() * (sep_peaks_max - sep_peaks_min) + sep_peaks_min # Separation in peaks in pixels
        peak1 = f * generalized_gauss_PSF(nrows, ncols, x,  y-sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        peak2 = f * generalized_gauss_PSF(nrows, ncols, x,  y+sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak    
        im += (peak1 + peak2) # Add peak    
    
    return im # Do not add any noise at this point.

# ----- Training data (equal number of blanks and blobs)
blanks = np.load("./train_data/blanks_filtered1_blacklist.npz")
iblack_list_blanks = blanks["black_list"]
blanks = blanks["data"]

N_blanks = blanks.shape[0]
Nsample = 128 * 40
SN_train = np.zeros((Nsample, 32, 32), dtype=float)
target_train = np.zeros(Nsample, dtype=bool)

i = 0
while i < Nsample:
    if (i % 1000) == 0:
        print("Generating %d example" % i)

    idx = np.random.randint(0, N_blanks, size=1)[0]
    
    # Data and err to be passed
    data = np.copy(blanks[idx, 0]) * np.sqrt(0.95) # Take the blank  
    err = np.copy(blanks[idx, 1])

    ibool = np.logical_or((data[4:25] == 0), (err[4:25] > 1e15)) # Identify all pixels that are either zero or have crazy errors.
    frac_zero = ibool.sum() / float(32*20)
    SN = data[4:25] / err[4:25] # Identify blanks whose SN is out of wack
    frac_SN_high = (SN > 10).sum() / float(32*20)

    # Location of zero pixels
    izero = np.abs(data) < 1e-15

    if (frac_zero > 0.5) or (frac_SN_high > 0.5) or (iblack_list_blanks[idx]): 
        # If more than 50% of pixels are zero or high SN
        # Then consider it a blank 
        data += np.random.randn(32, 32) * err * np.sqrt(0.05)
        data[izero] == 0.
        SN_train[i] = data/err
        i += 1
    else:    
        im_strip = data[4:25, :].flatten()
        im_strip_low = np.percentile(im_strip, 80) # Noise level
        im_strip_high = np.percentile(im_strip, 20)
        if im_strip_low == im_strip_high:
            pass
        else:
            ibool2 = np.logical_and((im_strip > im_strip_low), (im_strip < im_strip_high))
            if ibool2.sum() > 0:
                B = min(np.std(im_strip), np.std(im_strip[ibool]))
            else:
                B = np.std(im_strip)

            r = np.random.random()

            if (r >= 0.333) and (r < 0.6666):
                im_blob = blob_im_generator(double=False, fdensity = B / 15.)
                data += im_blob
                target_train[i] = True                
            elif (r >= 0.6666):
                im_blob = blob_im_generator(double=True, fdensity = B / 15.)
                data += im_blob
                target_train[i] = True                
                
            # Add Gaussian noise to the image. Note that this effectively changes the error but let's not be too concerned about it.
            im_noise = np.random.randn(32, 32) * err * np.sqrt(0.05)
            data += im_noise
            data[izero] = 0 # Restore zero positions
            SN_train[i] = data/err

            # Sample nuber increment   
            i += 1

np.savez("./train_data/train_data.npz", SN=SN_train, target=target_train)
# Finished generated train data.

# ---- Plot five panels of blank/blob        
N_panels = 5 
blobs_train = SN_train[target_train]
for i in range(N_panels):
    plot_post_stamps(blobs_train[i*100:(i+1)*100], fname="./train_data/fig_sim_data/blobs-%d.png" % i)
    
blanks_train = SN_train[~target_train]
for i in range(N_panels):
    plot_post_stamps(blanks_train[i*100:(i+1)*100], fname="./train_data/fig_sim_data/blanks-%d.png" % i)

