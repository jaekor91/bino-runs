# Note: We convert all model spectrum fluxes to f_lambda.
from utils import *
import simple_model
import h5py

lib = h5py.File("./fstars/c3k_v1.3_R5K.Fstars.h5", "r")
libwave = lib["wavelengths"][:]
params = lib["parameters"][:]

# plt.close()
# # Plot the coverage in Teff-logg space at solar Fe/H and \alpha/Fe
# sel = (params['feh'] == -.25) & (params['afe'] == 0)
# fig, ax = plt.subplots()
# ax.plot(10**params[sel]['logt'], params[sel]['logg'], 'o')
# ax.set_xscale('log')
# ax.set_xlabel('Teff')
# ax.set_ylabel('log g')
# ax.invert_xaxis()
# ax.invert_yaxis()
# plt.show()
# plt.close()

# # choose a particular temperature and gravity and plot the spectrum
# logt = np.unique(params['logt'])[15]
# sel = sel & (params['logt'] == logt) & (params['logg'] == 4.5)
# fig, ax = plt.subplots()
# ax.plot(libwave, np.squeeze(lib['spectra'][sel, :]))
# ax.set_xlabel('lambda (angstroms)')
# ax.set_ylabel('f_nu (erg/s/cm^2/Hz/sr)')
# ax.set_title('Teff={}'.format(10**logt))
# plt.show()
# plt.close()

def res_experimental(mask_name):
    if "270" in mask_name:
        return 1500
    else:
        return 3000

# Save directory for flux calibration vector
save_dir = "./flux_calibration/"

def sky_mask(wavegrid):
    """
    Given wavegrid return a boolean mask where True if fall into the sky line region
    """
    sky_lines = [(7570, 7725), (8600, 8700)]
    mask = np.zeros(wavegrid.size, dtype=bool)
    for lines in sky_lines:
        mask |= (wavegrid > lines[0]) & (wavegrid < lines[1])

    return mask


def bad_pixel_mask(sepc1d, err1d, ratio_low=0.25, ratio_high=4.):
    """
    Flag bad pixels by making them false.
    """
    mask = (err1d < 10.) # Unusually high errors
    # Extremely high variation in signal to noise
    SN = spec1d / err1d
    SN_median = median_filter(SN, size=301)
    SN2SN_med = SN/SN_median
    mask &=  np.logical_or((SN2SN_med > ratio_low), (SN2SN_med < ratio_high))
    
    # Extremely high spectral values compare to the median
    spec1d_med = median_filter(spec1d, size=301)
    spec2spec_med = spec1d/spec1d_med
    mask &= np.logical_or((spec2spec_med> ratio_low), (spec2spec_med < ratio_high))

    # Extremely high error values compare to the median
    err1d_med = median_filter(err1d, size=301)
    err2err_med = err1d/err1d_med
    mask &= np.logical_or((err2err_med> ratio_low), (err2err_med < ratio_high))
    
    return mask

# ----- Load all the stars
all_stars =np.load("stars/all-stars.npy").item()
mask_names_sorted = sorted(list(all_stars.keys()))
for mask_name in mask_names_sorted:
    print(mask_name)
    
    plot_dir = save_dir + mask_name + "/"
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    # Get stars corresponding to a mask
    stars = all_stars[mask_name]

    # Resolution
    res = res_experimental(mask_name)

    # Get wave grid
    wavegrid = stars["wavegrid"]

    # Get the keys corresponding to all the stars corresponding to a mask
    star_indices = list(stars.keys())
    del star_indices[star_indices.index("wavegrid")]
    star_indices.sort()

    # Different points on the spectrum is masked.
    cal_list = []
    mask_list = []
    model_index_list = []
    iflux_list = []

    for star_num in star_indices:
        # Get a star
        star = stars[star_num]

        # Save iflux 
        iflux_list.append(star["SDSS_iflux"])

        # Retrieve mask and construct spectrum
        spec1d = star["spec1D"][0]
        err1d = 1./np.sqrt(star["spec1D"][1])
        mask = ~np.logical_or.reduce((np.abs(spec1d) < 1e-6, np.abs(err1d) > 1e60, sky_mask(wavegrid), spec1d <= 0.))
        # Additional bad masks
        mask &= bad_pixel_mask(spec1d, err1d)

        plt.close()
        fig, ax = plt.subplots(1, figsize=(15, 7))
        ax.scatter(wavegrid, spec1d, c="black", edgecolors="none", s=5)
        # ----- Iterate through each model spectrum and find the best fit.
        # Save the best model_index
        model_index_best = None
        spec1d_masked = spec1d[mask]
        ivar_masked = 1. / err1d[mask]**2
        chi_sq_best = np.infty
        model_spec_best = None
        cal_best = None
        for model_index in range(params.size):
            libspec = np.squeeze(lib['spectra'][model_index, :])

            # Compute the model spectrum
            cal, model_spec = simple_model.make_model(wavegrid, spec1d, err1d, libwave, libspec/libwave**2, \
                                                          mask=mask, order=11, resolution=res)
            if model_spec is not None:    
                model_spec_calibrated = cal * model_spec
                chi_sq = np.sum(np.square(spec1d_masked - model_spec_calibrated[mask]) * ivar_masked)
                if chi_sq < chi_sq_best:
                    chi_sq_best = chi_sq
                    model_index_best = model_index
                    cal_best = cal
                    model_spec_best = model_spec_calibrated
        ax.scatter(wavegrid[mask], model_spec_best[mask], edgecolors="none", s=5)
        ax.set_ylim([-5, np.max(model_spec_best[mask]) * 1.5])
        plt.savefig(plot_dir+"star-specnum%d.png" % star_num, dpi=200, bbox_inches="tight")
#         plt.show()
        plt.close()

        model_index_list.append(model_index_best)
        cal_list.append(cal_best)
        mask_list.append(mask)
        
    all_stars[mask_name]["iflux_list"] = iflux_list
    all_stars[mask_name]["cal_list"] = cal_list
    all_stars[mask_name]["mask_list"] = mask_list
    all_stars[mask_name]["model_index_list"] = model_index_list

    
# ---- Compute and save the calibration vector
cal_dict = {} 
for mask_name in mask_names_sorted:
    print(mask_name)
    iflux_list = all_stars[mask_name]["iflux_list"]
    cal_list = all_stars[mask_name]["cal_list"]
    mask_list = all_stars[mask_name]["mask_list"]
    model_index_list = all_stars[mask_name]["model_index_list"]
    
    # ---- stars
    # Trimming the edges of cal estimators
    # Also mask
    L_trim = 500 
    trimmed_mask_list = []
    for mask_tmp in mask_list:
        mask = np.copy(mask_tmp)
        # Find min, max
        idx_min = 0
        idx_max = mask.size-1
        while ~mask[idx_min]:
            idx_min +=1
        while ~mask[idx_max]:
            idx_max -=1
        mask[:idx_min + L_trim] = False
        mask[idx_max - L_trim:] = False    
        mask[idx_min + L_trim:idx_max - L_trim] = True
        trimmed_mask_list.append(mask)

    iflux_ratio_list = []
    plt.close()
    fig, ax = plt.subplots(1, figsize=(15, 7))
    for i, cal in enumerate(cal_list):
        mask = trimmed_mask_list[i]
        model_index = model_index_list[i]
        libspec = np.squeeze(lib['spectra'][model_index, :])
        imag_model = simple_model.project_filter(libwave, libspec/libwave**2, bands=["sdss_i0"]).item()
        iflux_model = mag2flux(imag_model)
        iflux_obs = iflux_list[i]
        iflux_ratio = (iflux_model / iflux_obs)
        iflux_ratio_list.append(iflux_ratio)
        cal_tmp = cal[mask] * iflux_ratio
        if (np.abs(cal_tmp) < 3e18).all(): # Exclude bad lines
            ax.scatter(wavegrid[mask], cal_tmp, s=5, edgecolors="none")
    # ---- Overall calibration vector
    # Take the average and smooth many times
    cal_vec = np.zeros(mask.size)
    cal_list_tmp = []
    for k in range(len(cal_list)):
        tmp = cal_list[k] * iflux_ratio_list[k]
        if (np.abs(tmp[mask[k]]) < 3e18).all(): # Exclude bad lines
            cal_list_tmp.append(tmp)   
    cal_arr = np.vstack(cal_list_tmp)
    mask_arr = np.vstack(trimmed_mask_list)
    tally = mask_arr.sum(axis=0)
    idx_min = 0
    idx_max = tally.size-1
    while tally[idx_min] < 1:
        idx_min += 1
    while tally[idx_max] < 1:
        idx_max -=1
    for idx in range(idx_min, idx_max+1):
        cal_vec[idx] = np.median(cal_arr[:, idx][mask_arr[:, idx]])
    cal_vec[:idx_min] = np.median(cal_vec[idx_min:idx_min+500])
    cal_vec[idx_max:] = np.median(cal_vec[idx_max-500:idx_max])
    # cal_vec = savgol_filter(cal_vec, window_length=101, polyorder=5)
    cal_vec = median_filter(cal_vec, size=501)
    for _ in range(2):
        cal_vec = uniform_filter(cal_vec, size=501)

    ax.scatter(wavegrid, cal_vec, s=20, edgecolors="none", c="black")
    plt.savefig(save_dir + "calvec-" + mask_name+".png", dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()
    
    # --- Save
    cal_dict[mask_name] = cal_vec
np.save(save_dir + "calibration-vectors.npy", cal_dict)