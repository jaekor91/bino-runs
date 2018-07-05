# ----- Union catalog polishing
# -- Goals
# A more precise estimation of redshift location (keep the same if redshift deviates more than by 1e-3)
# Fluxing of individual peaks
# -- Extraction strategy:
# 1) Compute the center of peak by weighted summation (naive kernel determination)
# 2) Extract the peak according to the center with a box extraction kernel
# Display the 1D and 2D images of the peaks and a doublet fit for the former.

from utils import *

N_peaks = 5
peak2int = {"OII": 0, "Hb": 1, "OIII1": 2, "OIII2":3 , "Ha":4}
name = ["OII", "Hb", "OIII1", "OIII2", "Ha"]

# ---- Import calibration vectors
cal_dict = np.load("flux_calibration/calibration-vectors.npy").item()

# ---- Import union catalog
union_catalog = np.load("union-catalog-results.npy").item()
REDZ = union_catalog["REDZ"]
MASK_NUM = union_catalog["MASK_NUM"]
OBJ_NUM = union_catalog["OBJ_NUM"]
# ---- Create placeholders for fluxes
# First axis: Objnum
# Second: Blue (0) or red (1) grating
# Third: Five fluxes. peak2int
# Fourth: Continuum, A, sig_A, width, REDZ_new
FLUXES = np.ones((REDZ.size, 2, 5, 5), dtype=float) * -999

data_dir = "../../data/"
mask_dirs = os.listdir(data_dir)
fname_data = "obj_abs_slits_lin.fits" # File names
fname_err = "obj_abs_err_slits_lin.fits"
L_fit = 18
L_continuum = 64


black_listed = [8, 9, 10, 14] # Skip masks in the black list.
# For each mask 
for mask_num in range(16):
    # Check whether mask is black listed
    if mask_num not in black_listed:
        # Get the data directory names
        mask_name = union_catalog["MASK_NAMES"][mask_num]
        
        # Focus on the ones that matter
        ibool = (REDZ > 0) & (MASK_NUM == mask_num)
        objnums = OBJ_NUM[ibool]
        redzs = REDZ[ibool]
        # mini-place holder
        fluxes = np.ones((redzs.size, 2, 5, 5), dtype=float) * -999
                
        # ---- Load the data 
        fnames = []
        for q in mask_dirs:
            if mask_name in q:
                fnames.append(q)
        print(fnames)      
        data_err_1, list_headers_1 = preprocess_bino(fname_data, fname_err, data_dir+fnames[0]+"/")
        data_err_2, list_headers_2 = preprocess_bino(fname_data, fname_err, data_dir+fnames[1]+"/")

        save_dir = "review_panels_final/" + mask_name
        if not os.path.isdir(save_dir):
            os.mkdir(save_dir)

        ##### ----- Params used repeatedly
        Nobjs = data_err_1.shape[0] # Including the blank image at index 0. 
        data, err, header = extract_single_data(data_err_1, list_headers_1, 1)
        wavegrid1 = wavegrid_from_header(header, data.shape[1])
        cal1 = cal_dict[fnames[0]]
        data, err, header = extract_single_data(data_err_2, list_headers_2, 1)
        wavegrid2 = wavegrid_from_header(header, data.shape[1])
        cal2 = cal_dict[fnames[1]]
        
        counter = 0
        for z, specnum in zip(redzs, objnums):
            print(specnum, z)

            # Import 2D spectra
            data1, err1, header1 = extract_single_data(data_err_1, list_headers_1, specnum)
            data2, err2, header2 = extract_single_data(data_err_2, list_headers_2, specnum)

            # Compute edges where there is intensity in the data
            idx_min1, idx_max1 = index_edges(data1)
            idx_min2, idx_max2 = index_edges(data2)


            plt.close()
            fig, ax_list = plt.subplots(4, 5, figsize=(17, 10))
            # ------ First data
            peaks, indices = idx_peaks(wavegrid1, z, idx_min=idx_min1, idx_max=idx_max1)
            # Determine the center trace
            mu_centers = []
            for j in range(len(peaks)):
                idx = indices[j]
                if idx > 0: # If peak is within data region.
                    # Data post stamp
                    ps = median_filter((median_filter(data1[:, idx-L_fit:idx+L_fit], size=(5,5))), size=(5,5))
                    ps_ivar = 1./np.square(err1[:, idx-L_fit:idx+L_fit])
                    # Wavegrid snippet
                    wg = wavegrid1[idx-L_fit:idx+L_fit]
                    # Cal snippet
                    cl = cal1[idx-L_fit:idx+L_fit]
                    # Find the slit direction center by computing the naive kernel and fitting Gaussian to it.
                    K_naive = naive_profile(ps, ps_ivar, L_trim=0)
        #             mu_best = np.argmax(K_naive)
                    mu_best, _ = gauss_fit2profile(K_naive, mu_min=4., mu_max=22., sig_min=0.5, sig_max=3., dsig=0.1)   
                    if (mu_best < 4.2) or (mu_best > 21.8):
                        mu_best = 15.
        #             mu_centers.append(mu_best)
        #     mu_best = np.median(mu_centers)
        #     for j in range(len(peaks)):
        #         idx = indices[j]
        #         if idx > 0: # If peak is within data region.
                    # Continuum estimation
                    ps_continuum = data1[:, idx-L_continuum:idx+L_continuum]
                    ps_continuum_ivar = 1./np.square(err1[:, idx-L_continuum:idx+L_continuum])
                    # Cal snippet
                    cl = cal1[idx-L_continuum:idx+L_continuum]
                    # ---- Perform 1D extraction using Gaussian kernel
                    K = K_gauss_profile(mu_best, 2.257, Nrows = 32)
                    K_T = K.reshape((K.size, 1))            
                    spec1D_ivar = np.sum(np.square(K_T) * ps_continuum_ivar, axis=0)
                    spec1D = np.sum(K_T * ps_continuum * ps_continuum_ivar, axis=0) / spec1D_ivar
                    # Multiply by the flux calibration
                    spec1D_ivar *= cl**2
                    spec1D /= cl
                    # Estimate continuum as the median
                    continuum = np.median(spec1D)

                    # --- Perform Gaussian fitting to a smaller region
                    # Note that the Gaussians are fit independently
                    ps = data1[:, idx-L_fit:idx+L_fit]
                    ps_ivar = 1./np.square(err1[:, idx-L_fit:idx+L_fit])
                    wg = wavegrid1[idx-L_fit:idx+L_fit] # Wavegrid snippet
                    cl = cal1[idx-L_fit:idx+L_fit] # Cal snippet
                    # Perform 1D extraction and multiply by the flux calibration
                    K = K_gauss_profile(mu_best, 2.257, Nrows = 32)
                    K_T = K.reshape((K.size, 1))            
                    spec1D_ivar = np.sum(np.square(K_T) * ps_ivar, axis=0)
                    spec1D = np.sum(K_T * ps * ps_ivar, axis=0) / spec1D_ivar
                    spec1D_ivar *= cl**2
                    spec1D /= cl

                    # Make a gaussian fit.
                    A_best, sig_A_best, width_best, z_best = fit_emission_line(wg, spec1D, spec1D_ivar, z_initial=z, peak_name=peaks[j], continuum=continuum)
                    model_spec = A_best * emission_line_model(wg, width_best, z_best, peak_name=peaks[j]) + continuum
                    print(("%s %.2f %.2f %.2f %.2f")  % (peaks[j], A_best * 1e17, sig_A_best * 1e17, width_best, z_best))
                    # Save the result of fitting here
                    fluxes[counter, 0, peak2int[peaks[j]]] = np.array([continuum, A_best, sig_A_best, width_best, z_best])
                        
                    # -- Plot the fit
                    # 2D 
                    ax_list[0, peak2int[peaks[j]]].imshow(ps, cmap="gray", interpolation="None")                    
                    ax_list[0, peak2int[peaks[j]]].axvline(x=L_fit, c="red", ls="--", lw=0.5)
                    ax_list[0, peak2int[peaks[j]]].axhline(y=mu_best, c="red", ls="--", lw=0.5)                        
                    # 1D
                    ax_list[1, peak2int[peaks[j]]].plot(spec1D, c="black", lw=1)
                    ax_list[1, peak2int[peaks[j]]].plot(model_spec, c="blue", lw=1)            
                    ax_list[1, peak2int[peaks[j]]].set_title("%.2f e17" % (A_best * 1e17))            
                    ax_list[1, peak2int[peaks[j]]].axvline(x=L_fit, c="red", ls="--", lw=0.5)
                    ax_list[1, peak2int[peaks[j]]].axhline(y=continuum, c="blue", ls="--", lw=0.5)
                else:
                    ax_list[1, j].axis("off")        
                ax_list[0, j].axis("off")            

            # ------ Second data
            peaks, indices = idx_peaks(wavegrid2, z, idx_min=idx_min2, idx_max=idx_max2)
            # Determine the center trace
            mu_centers = []
            for j in range(len(peaks)):
                idx = indices[j]
                if idx > 0: # If peak is within data region.
                    # Data post stamp
                    ps = median_filter((median_filter(data2[:, idx-L_fit:idx+L_fit], size=(5,5))), size=(5,5))
                    ps_ivar = 1./np.square(err2[:, idx-L_fit:idx+L_fit])
                    # Wavegrid snippet
                    wg = wavegrid2[idx-L_fit:idx+L_fit]
                    # Cal snippet
                    cl = cal2[idx-L_fit:idx+L_fit]
                    # Find the slit direction center by computing the naive kernel and fitting Gaussian to it.
                    K_naive = naive_profile(ps, ps_ivar, L_trim=0)
                    mu_best, _ = gauss_fit2profile(K_naive, mu_min=4., mu_max=22., sig_min=0.5, sig_max=3., dsig=0.1)   
                    if (mu_best < 4.2) or (mu_best > 21.8):
                        mu_best = 15.
        #             mu_centers.append(mu_best)
        #     mu_best = np.median(mu_centers)
        #     for j in range(len(peaks)):
        #         idx = indices[j]
        #         if idx > 0: # If peak is within data region.
                    # Continuum estimation
                    ps_continuum = data2[:, idx-L_continuum:idx+L_continuum]
                    ps_continuum_ivar = 1./np.square(err2[:, idx-L_continuum:idx+L_continuum])
                    # Cal snippet
                    cl = cal2[idx-L_continuum:idx+L_continuum]
                    # ---- Perform 1D extraction using Gaussian kernel
                    K = K_gauss_profile(mu_best, 2.257, Nrows = 32)
                    K_T = K.reshape((K.size, 1))            
                    spec1D_ivar = np.sum(np.square(K_T) * ps_continuum_ivar, axis=0)
                    spec1D = np.sum(K_T * ps_continuum * ps_continuum_ivar, axis=0) / spec1D_ivar
                    # Multiply by the flux calibration
                    spec1D_ivar *= cl**2
                    spec1D /= cl
                    # Estimate continuum as the median
                    continuum = np.median(spec1D)

                    # --- Perform Gaussian fitting to a smaller region
                    # Note that the Gaussians are fit independently
                    ps = data2[:, idx-L_fit:idx+L_fit]
                    ps_ivar = 1./np.square(err2[:, idx-L_fit:idx+L_fit])
                    wg = wavegrid2[idx-L_fit:idx+L_fit] # Wavegrid snippet
                    cl = cal2[idx-L_fit:idx+L_fit] # Cal snippet
                    # Perform 1D extraction and multiply by the flux calibration
                    K = K_gauss_profile(mu_best, 2.257, Nrows = 32)
                    K_T = K.reshape((K.size, 1))            
                    spec1D_ivar = np.sum(np.square(K_T) * ps_ivar, axis=0)
                    spec1D = np.sum(K_T * ps * ps_ivar, axis=0) / spec1D_ivar
                    spec1D_ivar *= cl**2
                    spec1D /= cl

                    # Make a gaussian fit.
                    A_best, sig_A_best, width_best, z_best = fit_emission_line(wg, spec1D, spec1D_ivar, z_initial=z, peak_name=peaks[j], continuum=continuum)
                    model_spec = A_best * emission_line_model(wg, width_best, z_best, peak_name=peaks[j]) + continuum
                    print(("%s %.2f %.2f %.2f %.2f")  % (peaks[j], A_best * 1e17, sig_A_best * 1e17, width_best, z_best))
                    # Save the result of fitting here
                    fluxes[counter, 1, peak2int[peaks[j]]] = np.array([continuum, A_best, sig_A_best, width_best, z_best])                    

                    # -- Plot the fit
                    # 2D 
                    ax_list[2, peak2int[peaks[j]]].imshow(ps, cmap="gray", interpolation="None")                    
                    ax_list[2, peak2int[peaks[j]]].axvline(x=L_fit, c="red", ls="--", lw=0.5)
                    ax_list[2, peak2int[peaks[j]]].axhline(y=mu_best, c="red", ls="--", lw=0.5)                        
                    # 1D
                    ax_list[3, peak2int[peaks[j]]].plot(spec1D, c="black", lw=1)
                    ax_list[3, peak2int[peaks[j]]].plot(model_spec, c="blue", lw=1)
                    ax_list[3, peak2int[peaks[j]]].set_title("%.2f e17" % (A_best * 1e17))
                    ax_list[3, peak2int[peaks[j]]].axvline(x=L_fit, c="red", ls="--", lw=0.5)
                    ax_list[3, peak2int[peaks[j]]].axhline(y=continuum, c="blue", ls="--", lw=0.5)
                else:
                    ax_list[3, j].axis("off")        
                ax_list[2, j].axis("off")      

            for m in range(N_peaks):
                ax_list[0, m].set_title(name[m], fontsize=15)
            title_str = "specnum%03d-z%.3f" % (specnum, z)
            plt.suptitle(title_str, fontsize=15)
            plt.savefig(save_dir+"/" + title_str + ".png", dpi=100, bbox_inches="tight")
        #     plt.show()
            plt.close()
            print("\n")
            
            # Save the estimated fluxes.
            FLUXES[ibool] = fluxes
            counter += 1
        print("\n")

union_catalog["FLUXES"] = FLUXES
np.save("union-catalog-results-fluxed.npy", union_catalog)