from utils import *

im_arr_blobs = np.load("./blob_finder_training_data/real-blobs.npz")["image"]
im_arr_blanks = np.load("./blob_finder_training_data/blanks-filtered.npz")["image"]

# Basic parameters
# - Total flux could range 2,500 to 25,000
# - Diagonal 2D covariance 16 = r^2 of varying from 0.1 to 0.9 randomly


def blob_im_generator(nrows=32, ncols=32, double=False):
    """
    Generate an emission line peak blob.
    """
    fmax = 25000
    fmin = 20000
    FWHM_min = 5
    FWHM_max = 7
    rho_min = -0.8
    rho_max = 0.8
    num_comps_max = 10
    # For double peaks only
    sep_peaks_min = 5
    sep_peaks_max = 8
    scatter_max = 0.5
    scatter_min = 0
    
    im = np.zeros((nrows, ncols))
    x = nrows//2
    y = ncols//2
    if not double:
        f = (fmax - fmin) * np.random.random() + fmin # Random flux selection
        rho = (rho_max - rho_min) * np.random.random() + rho_min # Covarince selection
        FWHM = (FWHM_max - FWHM_min) * np.random.random() + FWHM_min # FWHM selection
        num_comps = np.random.randint(1, num_comps_max, 1)[0]
        scatter = np.random.random() * scatter_max
        peak = f * generalized_gauss_PSF(nrows, ncols, x,  y, FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        im += peak # Add peak
    else:
        sep_peaks = np.random.random() * (sep_peaks_max - sep_peaks_min) + sep_peaks_min # Separation in peaks in pixels
        f = (fmax - fmin) * np.random.random() + fmin # Random flux selection
        rho = (rho_max - rho_min) * np.random.random() + rho_min # Covarince selection
        FWHM = (FWHM_max - FWHM_min) * np.random.random() + FWHM_min # FWHM selection
        num_comps = np.random.randint(1, num_comps_max, 1)[0]
        scatter = np.random.random() * scatter_max
        peak1 = f * generalized_gauss_PSF(nrows, ncols, x,  y-sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        peak2 = f * generalized_gauss_PSF(nrows, ncols, x,  y+sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak    
        im += (peak1 + peak2) # Add peak    
    
    return poisson_realization(im)

Nsample = 10000
im_sim_training = np.zeros((Nsample, 32, 32))
label_training = np.zeros(Nsample, dtype=int)

idx = 0 
num_blanks = xrange(im_arr_blanks.shape[0]) # Number of blanks
B = 100
while idx < Nsample:
    idx_blank = np.random.choice(num_blanks) # Choose a random blank
    im = im_arr_blanks[idx_blank] # Take the blank
    r = np.random.random()
#     r = 0
    if (r < 0.5): # Add the blank to the training data
        pass
    elif (r > 0.5) and (r < 0.75):
        im_blob = blob_im_generator(double=False)
        im += im_blob
        label_training[idx] = 1
    else:
        im_blob = blob_im_generator(double=True)
        im += im_blob
        label_training[idx] = 2
        
    # Add poisson noise to the image
    # Generate poisson noise image, add the to the blank and subtract
    im_poisson = poisson_realization(np.ones((32, 32)) * B) - B    
    im += im_poisson
    
    # Save the image
    im_sim_training[idx] = im
    idx += 1
        
# # ---- View a sample of images.
# plt.close()
# fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))

# i_start = 0
# i_end = i_start + 81
# for i in range(i_start, i_end):
#     idx_row = (i-i_start) // 9
#     idx_col = (i-i_start) % 9
#     ax_list[idx_row, idx_col].imshow(im_sim_training[i, :, :], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
#     title_str = "%4d" % (label_training[i])
#     ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
#     ax_list[idx_row, idx_col].axis("off")    

# # plt.savefig("blob_sim_training_examples.png", dpi=200, bbox_inches="tight")
# plt.show()
# plt.close()

np.savez("./blob_finder_training_data/blob_finder_train_data.npz", image=im_sim_training, label=label_training)