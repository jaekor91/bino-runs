from utils import *
from resnet_clf import ResnetBuilder
from keras import optimizers
from keras import losses
from keras.models import load_model
from sklearn.utils import shuffle
import os

N_peaks = 5
peak2int = {"OII": 0, "Hb": 1, "OIII1": 2, "OIII2":3 , "Ha":4}
name = ["OII", "Hb", "OIII1", "OIII2", "Ha"]


# ---- Load the model to use.
model = load_model('../analysis/david/images2/classification_model_SN.h5') # From david




# ----- Perform detections
# ------------------ Threshold value to use
thres = 0.9999671 # David 80% recall rate
# Recall rate: 0.8000000
# Precision rate: 0.8877246
# Estimated thres: 0.9999671
# Total number of sample: 39614
# Number of failures: 224
# Val -- Failed
# Recall rate: 0.8500000
# Precision rate: 0.8246073
# Estimated thres: 0.9995236
# Total number of sample: 39614
# Number of failures: 244
# Val -- Failed
# Recall rate: 0.9000000
# Precision rate: 0.6834016
# Estimated thres: 0.7355986
# Total number of sample: 39614
# Number of failures: 382
# Val -- Failed
# Estimated recall rate: 0.8893
# Estimated precision rate: 0.7086
# Total number of sample: 39614
# Number of failures: 353
# Val -- Failed
# ------------------

# ---- Mask directory names
data_dir = "../../data/"
dirs = os.listdir(data_dir)

# --- Review panel directory
if not os.path.isdir("review_panel"):
    os.mkdir("review_panel")

fname_data = "obj_abs_slits_lin.fits"
fname_err = "obj_abs_err_slits_lin.fits"

# ---- For each set
for q in range(1, len(dirs)//2):
    print(q, dirs[2*q:2*(q+1)])

    # ---- Load the data 
    data_err_1, list_headers_1 = preprocess_bino(fname_data, fname_err, data_dir+dirs[2*q]+"/")
    data_err_2, list_headers_2 = preprocess_bino(fname_data, fname_err, data_dir+dirs[2*q+1]+"/")

    save_dir = "review_panel/" + dirs[2*q]
    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    ##### ----- Params used repeatedly
    Nobjs = data_err_1.shape[0] # Including the blank image at index 0. 
    data, err, header = extract_single_data(data_err_1, list_headers_1, 1)
    wavegrid1 = wavegrid_from_header(header, data.shape[1])
    data, err, header = extract_single_data(data_err_2, list_headers_2, 1)
    wavegrid2 = wavegrid_from_header(header, data.shape[1])

    # ---- Analyze each spectrum
    for specnum in range(1, Nobjs):
        start = time.time()
        print("Working on spec %d" % specnum)
        data1, err1, header1 = extract_single_data(data_err_1, list_headers_1, specnum)
        data2, err2, header2 = extract_single_data(data_err_2, list_headers_2, specnum)

        # Access bit number
        bit = bit_from_header(header1)

        if (not bit_true(bit, 2**1)) and (not bit_true(bit, 2**2)): # If the object is a target
            # Generate SN image
            SN1 = data1/err1
            SN2 = data2/err2    

            # Compute edges where there is intensity in the data
            idx_min1, idx_max1 = index_edges(data1)
            idx_min2, idx_max2 = index_edges(data2)

            # Boolean vctors that indicate where (center) to take the postage stamps from.
            query1 = np.zeros(wavegrid1.size, dtype=bool)
            query2 = np.zeros(wavegrid2.size, dtype=bool)
            query1[idx_min1+16:idx_max1-16:2] = True
            query2[idx_min2+16:idx_max2-16:2] = True

            # Array that records results of NN query
            detection1 = np.zeros_like(query1)
            detection2 = np.zeros_like(query2)    

            # ---- Generate postage stamps for vectorized query
            stamps1 = gen_post_stamp_arr(SN1, query1)
            stamps2 = gen_post_stamp_arr(SN2, query2)    

            # ----- Reject certain cases using heuristics
            reject1 = reject_heuristics(stamps1)
            reject2 = reject_heuristics(stamps2)
            # Do not query the points that are pre-rejected.
            query1[query1] &= ~reject1
            query2[query2] &= ~reject2
            # Re-generate post satmps
            stamps1 = gen_post_stamp_arr(SN1, query1)
            stamps2 = gen_post_stamp_arr(SN2, query2)    

            # ---- Use neural network to predict.
            P1 = model.predict(stamps1).ravel()
            P2 = model.predict(stamps2).ravel()            
            detection1[query1] = P1 > thres
            detection2[query2] = P2 > thres
            
            # ----- Find redshift candidate based on the detections
            zs, num_hits = z_candidates(wavegrid1, wavegrid2, detection1, detection2, \
                            idx_min1=idx_min1, idx_max1=idx_max1, idx_min2=idx_min2, idx_max2=idx_max2)

            # ---- For each redshift candidate create a review panel
            if zs.size == 0: 
            	print("No candidate detected")
            else:
	            print("Number of candidates %d" % len(zs))
	            print("Max hits (cluster average): %.3f" % num_hits[0])
	            # Only look at top five
	            for i, z in enumerate(zs):
	                print(i, z, num_hits[i])
	                plt.close()
	                fig, ax_list = plt.subplots(4, N_peaks, figsize=(7, 7))
	                # Location of peaks in data 1 given the redshift.
	                peaks, indices = idx_peaks(wavegrid1, z, idx_min=idx_min1, idx_max=idx_max1)
	                for j in range(len(peaks)):
	                    idx = indices[j]
	                    if idx > 0: # If peak is within data region.
	                        ax_list[0, peak2int[peaks[j]]].imshow(data1[:, idx-16:idx+16], cmap="gray", interpolation="None")                    
	                        ax_list[0, peak2int[peaks[j]]].axvline(x=15.5, c="red", ls="--", lw=1)            
	                        ax_list[1, peak2int[peaks[j]]].imshow(SN1[:, idx-16:idx+16], cmap="gray", interpolation="None")                                
	                        ax_list[1, peak2int[peaks[j]]].axvline(x=15.5, c="red", ls="--", lw=1)            
	                    ax_list[0, j].axis("off")
	                    ax_list[1, j].axis("off")        

	                # Location of peaks in data 2 given the redshift.
	                peaks, indices = idx_peaks(wavegrid2, z, idx_min=idx_min2, idx_max=idx_max2)
	                for j in range(len(peaks)):
	                    idx = indices[j]
	                    if idx > 0: # If peak is within data region.
	                        ax_list[2, peak2int[peaks[j]]].imshow(data2[:, idx-16:idx+16], cmap="gray", interpolation="None")  
	                        ax_list[2, peak2int[peaks[j]]].axvline(x=15.5, c="red", ls="--", lw=1)                        
	                        ax_list[3, peak2int[peaks[j]]].imshow(SN2[:, idx-16:idx+16], cmap="gray", interpolation="None")  
	                        ax_list[3, peak2int[peaks[j]]].axvline(x=15.5, c="red", ls="--", lw=1)                        
	                    ax_list[2, j].axis("off")
	                    ax_list[3, j].axis("off")        

	                for m in range(N_peaks):
	                    ax_list[0, m].set_title(name[m], fontsize=15)
	                title_str = "specnum%03d-cand%d-hits%.3f-z%.5f" % (specnum, i, num_hits[i], z)
	                plt.suptitle(title_str, fontsize=15)
	                plt.savefig(save_dir+"/" + title_str + ".png", dpi=100, bbox_inches="tight")
	                plt.close()
	            end = time.time()
	            print("time taken: %d sec" % (end-start))
	            print("\n")
	        print("\n")