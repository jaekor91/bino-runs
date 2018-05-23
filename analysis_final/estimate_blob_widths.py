from utils import *
blobs = np.load("./train_data/blobs_filtered2.npy")



# Estimate the width distribution
sig_x = [] # Width in row direction
sig_y = [] # Width in column direction

for i in range(blobs.shape[0]):
    data = blobs[i, 0]
    err = blobs[i, 1]
    ivar = ivar_from_err(err)
    
    # Compute extract kernel 
    # Dispersion in x direction for sig_x
    K_x = np.sum(data * ivar, axis = 1) / np.sum(ivar, axis = 1)
    K_x /= np.sum(K_x) # Normalization step.
    _, sig_x_tmp = gauss_fit2profile(K_x, sig_min=0.5, sig_max=3.5)
    sig_x.append(sig_x_tmp)

# ---- Plot histogram and estimate 
plt.hist(sig_x, bins=np.arange(0.5, 3.5, 0.1))
plt.show()
plt.close()
# My estimates based on eyeballing.
mean = 1.9
FWHM = 0.9
sig = FWHM / 2.355
sig


sig_y = [] # Width in column direction

for i in range(blobs.shape[0]):
    data = blobs[i, 0]
    err = blobs[i, 1]
    ivar = ivar_from_err(err)
    
    # Compute extract kernel 
    # Dispersion in y direction for sig_y
    K_y = np.sum(data * ivar, axis = 0) / np.sum(ivar, axis = 0)
    K_y /= np.sum(K_y) # Normalization step.
    
    _, sig_y_tmp = gauss_fit2profile(K_y, sig_min=1., sig_max=6., dsig=0.1)
    sig_y.append(sig_y_tmp)
    

    # ---- Plot histogram and estimate 
plt.hist(sig_y, bins=np.arange(0.5, 6., 0.1))
plt.show()
plt.close()
# My estimates based on eyeballing.
# Very noisy.
# Mean and sig should be double of the sig_x case. 