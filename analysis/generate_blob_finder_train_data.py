from utils import *

print("/---- Train set")
im_arr_blobs = np.copy(np.load("./blob_finder_training_data/real-blobs.npz")["image"])
err_arr_blobs = np.copy(np.load("./blob_finder_training_data/real-blobs.npz")["err"])
im_arr_blanks = np.copy(np.load("./blob_finder_training_data/blanks-filtered.npz")["image"])
err_arr_blanks = np.copy(np.load("./blob_finder_training_data/blanks-filtered.npz")["err"])

# --- Basic parameters
# - Multiply everything by 1e16 factor in the beginning.

# --- Important considerations
# - Input both noise and values.
# - Visual inspection of positives
# - Zero out the original zero points

# --- During testing
# - If more than 80% is blank, then pass.


def blob_im_generator(nrows=32, ncols=32, double=False, fdensity=0):
    """
    Generate an emission line peak blob.
    """
    FWHM_min = 3
    FWHM_max = 6
    FWHM = (FWHM_max - FWHM_min) * np.random.random() + FWHM_min # FWHM selection    
    fmin = 10 * FWHM**2 * fdensity * (5 + 2 * (np.random.random()))
    fmax = fmin * 1.01
    rho_min = -0.6
    rho_max = 0.6
    num_comps_min = 10
    num_comps_max = 25 
    # For double peaks only
    sep_peaks_min = 3
    sep_peaks_max = 10
    scatter_max = 0.5
    scatter_min = 0
    
    im = np.zeros((nrows, ncols))
    x = nrows//2 + (np.random.random() - 0.5) * 3
    y = ncols//2 + (np.random.random() - 0.5) * 3
    
    f = (fmax - fmin) * np.random.random() + fmin # Random flux selection
    rho = (rho_max - rho_min) * np.random.random() + rho_min # Covarince selection
    num_comps = np.random.randint(num_comps_min, num_comps_max, 1)[0]
    scatter = np.random.random() * scatter_max    
    if not double:
        peak = f * generalized_gauss_PSF(nrows, ncols, x,  y, FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        im += peak # Add peak
    else:
        sep_peaks = np.random.random() * (sep_peaks_max - sep_peaks_min) + sep_peaks_min # Separation in peaks in pixels
        peak1 = f * generalized_gauss_PSF(nrows, ncols, x,  y-sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak
        peak2 = f * generalized_gauss_PSF(nrows, ncols, x,  y+sep_peaks/2., FWHM=FWHM, rho=rho, scatter=scatter, num_comps=num_comps) # Double peak    
        im += (peak1 + peak2) # Add peak    
    
    return poisson_realization(im * 100) / 100 # Use arbitrary counts-to-flux conversion.

Nsample = 256 * 10
# im_sim_training = np.zeros((Nsample, 32, 32, 1)) 
im_sim_training = np.zeros((Nsample, 32, 32, 2)) # Image and SN information provided.
label_training = np.zeros(Nsample, dtype=bool)

idx = 0 
num_blanks = range(im_arr_blanks.shape[0]) # Number of blanks
while idx < Nsample:
    if (idx % 5000) == 0:
        print(idx)
    idx_blank = np.random.choice(num_blanks) # Choose a random blank
    im = np.copy(im_arr_blanks[idx_blank]) # Take the blank and err
    err = np.copy(err_arr_blanks[idx_blank])

    ibool = (im == 0)
    if (((ibool).sum() / float(32**2)) > 0.70): # If more than 70% of pixels are zeros
        pass
    else: # Otherwise continue
        im_strip = im[6:25, :].flatten()
        im_strip_low = np.percentile(im_strip, 80) # Noise level
        im_strip_high = np.percentile(im_strip, 20)
        if im_strip_low == im_strip_high:
            pass
        else:
            # print(im_strip_low, im_strip_high)x
            ibool2 = np.logical_and((im_strip > im_strip_low), (im_strip < im_strip_high))
            if ibool2.sum() > 0:
                B = min(np.std(im_strip), np.std(im_strip[ibool]))
            else:
                B = np.std(im_strip)

            r = np.random.random()
        #     r = 0
            if (r < 0.5): # Add the blank to the training data
                pass
            elif (r > 0.5) and (r < 0.55):
                im_blob = blob_im_generator(double=False, fdensity = B)
                im += im_blob
                label_training[idx] = True
            else:
                im_blob = blob_im_generator(double=True, fdensity = B)
                im += im_blob
                label_training[idx] = True
                
            # Add poisson noise to the image
            # Generate poisson noise image, add the to the blank and subtract
            im_poisson = (poisson_realization(np.ones((32, 32)) * B * 100) - B * 100) / 100
            im += im_poisson
            im[ibool] = 0 # Restore zero positions
            
            # Save the SN image
            im_sim_training[idx, :, :, 0] = im
            im_sim_training[idx, :, :, 1] = im / err
            idx += 1


# ---- View only the positives
# -- Image
plt.close()
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))

num_panels = 5
idx = 0
for l in range(num_panels):
    counter = 0
    while counter < 81:
        if label_training[idx]: 
            idx_row = counter // 9
            idx_col = counter % 9
            ax_list[idx_row, idx_col].imshow(im_sim_training[idx, :, :, 0], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
        #     title_str = "%4d" % (label_training[i])
            title_str = label_training[idx]
            ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
            ax_list[idx_row, idx_col].axis("off") 
            counter += 1
        idx +=1

    plt.savefig("blob_sim_training_examples_%d.png" % l, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()
# -- Error
plt.close()
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))

num_panels = 5
idx = 0
for l in range(num_panels):
    counter = 0
    while counter < 81:
        if label_training[idx]: 
            idx_row = counter // 9
            idx_col = counter % 9
            ax_list[idx_row, idx_col].imshow(im_sim_training[idx, :, :, 1], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
        #     title_str = "%4d" % (label_training[i])
            title_str = label_training[idx]
            ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
            ax_list[idx_row, idx_col].axis("off") 
            counter += 1
        idx +=1

    plt.savefig("blob_sim_training_examples_%d_SN.png" % l, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    

# --- Blank images
# -- Image
plt.close()
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))

num_panels = 5
idx = 0
for l in range(num_panels):
    counter = 0
    while counter < 81:
        if not label_training[idx]: 
            idx_row = counter // 9
            idx_col = counter % 9
            ax_list[idx_row, idx_col].imshow(im_sim_training[idx, :, :, 0], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
        #     title_str = "%4d" % (label_training[i])
            title_str = label_training[idx]
            ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
            ax_list[idx_row, idx_col].axis("off") 
            counter += 1
        idx +=1

    plt.savefig("blanks_sim_training_examples_%d.png" % l, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    

# Err
plt.close()
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))

num_panels = 5
idx = 0
for l in range(num_panels):
    counter = 0
    while counter < 81:
        if not label_training[idx]: 
            idx_row = counter // 9
            idx_col = counter % 9
            ax_list[idx_row, idx_col].imshow(im_sim_training[idx, :, :, 1], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
        #     title_str = "%4d" % (label_training[i])
            title_str = label_training[idx]
            ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
            ax_list[idx_row, idx_col].axis("off") 
            counter += 1
        idx +=1

    plt.savefig("blanks_sim_training_examples_%d_SN.png" % l, dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    


np.savez("./blob_finder_training_data/blob_finder_train_data.npz", sample=im_sim_training, label=label_training)

print("Completed")





print("/---- Test set")
im_arr_blobs = np.copy(np.load("./blob_finder_training_data/real-blobs.npz")["image"])
err_arr_blobs = np.copy(np.load("./blob_finder_training_data/real-blobs.npz")["err"])
im_arr_blanks = np.copy(np.load("./blob_finder_training_data/blanks-filtered.npz")["image"])
err_arr_blanks = np.copy(np.load("./blob_finder_training_data/blanks-filtered.npz")["err"])


N_blanks = im_arr_blanks.shape[0]
N_blobs = im_arr_blobs.shape[0]
N_total = N_blanks + N_blobs
label = np.zeros(N_total, dtype=bool)
label[:N_blanks] = False # Blank is False
label[N_blanks:] = True
im_arr = np.zeros((N_total, 32, 32, 2))

# ---- Blanks first
for i in range(N_blanks):
    im = im_arr_blanks[i]
    err = err_arr_blanks[i]
    im_arr[i, :, :, 0] = im
    im_arr[i, :, :, 1] = im / err

# ---- Blobs second
for i in range(N_blobs):
    im = im_arr_blobs[i]
    err = err_arr_blobs[i]
    im_arr[i + N_blanks, :, :, 0] = im
    im_arr[i + N_blanks, :, :, 1] = im / err 

np.savez("./blob_finder_training_data/blob_finder_test_data.npz", sample=im_arr, label=label)
print("Completed")