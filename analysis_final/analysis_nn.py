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
thres = 0.9999670 # David # Set based on recall threshold of 0.83
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
for q in range(4, len(dirs)//2):
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
		BIT = bit_from_header(header1)

		if (BIT != 2) and (BIT !=4): # If the object is a true target.
			# Generate SN image
			SN1 = data1/err1
			SN2 = data2/err2    

			# Compute edges where there is intensity in the data
			idx_min1, idx_max1 = index_edges(data1)
			idx_min2, idx_max2 = index_edges(data2)

			# Boolean vctors that indicate where (center) to take the postage stamps from.
			query1 = np.zeros(wavegrid1.size, dtype=bool)
			query1[idx_min1+16:idx_max1-16:10] = True
			query2 = np.zeros(wavegrid2.size, dtype=bool)
			query2[idx_min2+16:idx_max2-16:10] = True

			# Array that records results of NN query
			detection1 = np.zeros_like(query1)
			detection2 = np.zeros_like(query2)    

			# ---- Generate postage stamps for vectorized query
			stamps1 = gen_post_stamp_arr(SN1, query1)
			stamps2 = gen_post_stamp_arr(SN2, query2)    

			# ---- Use neural network to predict.
			detection1[query1] = model.predict(stamps1).ravel() > thres
			detection2[query2] = model.predict(stamps2).ravel() > thres

			# ----- Reject certain cases using heuristics
			reject1 = reject_heuristics(stamps1)
			reject2 = reject_heuristics(stamps2)
			detection1[query1] &= ~reject1
			detection2[query2] &= ~reject2

			# ----- Find redshift candidate based on the detections
			zs, num_hits = z_candidates(wavegrid1, wavegrid2, detection1, detection2, \
							idx_min1=idx_min1, idx_max1=idx_max1, idx_min2=idx_min2, idx_max2=idx_max2)
				
			# # ---- For each redshift candidate create a review panel
			if (len(zs) > 0) and (len(zs) < 50):
				
				# --- Look at only top fives
	#             if num_hits[0] >=3 :
	#                 zs = np.asarray(zs[:5])
	#                 num_hits = np.asarray(num_hits[:5])
	#             else:
				zs = np.asarray(zs)
				num_hits = np.asarray(num_hits)
					
				# Sort by redshift
				idx_sort = np.argsort(zs)
				zs = zs[idx_sort][:50]
				num_hits[idx_sort][:50]
				
				print("Number of candidates %d" % len(zs))
				print("Max hits: %d" % num_hits[0])
				for i, z in enumerate(zs):
					print(i, z, num_hits[i])
					plt.close()
					fig, ax_list = plt.subplots(2, N_peaks, figsize=(7, 3))
					# Location of peaks in data 1 given the redshift.
					peaks, indices = idx_peaks(wavegrid1, z, idx_min=idx_min1, idx_max=idx_max1)
					for j in range(len(peaks)):
						idx = indices[j]
						if idx > 0: # If peak is within data region.
							ax_list[0, peak2int[peaks[j]]].imshow(SN1[:, idx-16:idx+16], cmap="gray", interpolation="None")                    
						ax_list[0, j].axis("off")

					# Location of peaks in data 2 given the redshift.
					peaks, indices = idx_peaks(wavegrid2, z, idx_min=idx_min2, idx_max=idx_max2)
					for j in range(len(peaks)):
						idx = indices[j]
						if idx > 0: # If peak is within data region.
							ax_list[1, peak2int[peaks[j]]].imshow(SN2[:, idx-16:idx+16], cmap="gray", interpolation="None")                    
						ax_list[1, j].axis("off")

					for m in range(N_peaks):
						ax_list[0, m].set_title(name[m], fontsize=15)

					plt.savefig(save_dir+"/specnum%03d-cand%d-hits%d-z%.5f.png" % (specnum, i, num_hits[i], z), dpi=50, bbox_inches="tight")

				#     plt.show()
					plt.close()
		end = time.time()
		print("time taken: %d sec" % (end-start))	    
		print("\n")
	print("\n")