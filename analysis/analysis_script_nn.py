from utils import *
from resnet_clf import ResnetBuilder
from keras import optimizers
from keras import losses
from keras.models import load_model
from sklearn.utils import shuffle
import os

# ---- Load the model to use.
model = load_model('./david/images2/classification_model_SN.h5') # From david
# model = load_model('./classification_model_SN.h5') # From daviz

# ----- Perform detections
# ------------------ Threshold value to use
thres = 0.98 # David
# thres = 0.14 # Jae
# ------------------

# ---- File names
data_dir = "/Users/jaehyeon/Documents/Research/binospec/data/v1/"
dirs = os.listdir(data_dir)[1:]



# --- Review panel directory
if not os.path.isdir("review_panel"):
	os.mkdir("review_panel")

for mask_num in range(1, len(dirs)//2):
	dir_name = dirs[2 * mask_num]
	save_dir = "./review_panel/" + dir_name
	#---- Create the directory afresh
	if os.path.isdir(save_dir):
		os.rmdir(save_dir)
	os.mkdir(save_dir)
	
	print("/----------- %s" % dir_name)
	
	fname1 = data_dir + dirs[2 * mask_num] + "/obj_abs_slits_lin.fits"
	fname_err1 = data_dir + dirs[2 * mask_num] + "/obj_abs_err_slits_lin.fits"    
	fname2 = data_dir + dirs[2 * mask_num+1] + "/obj_abs_slits_lin.fits"
	fname_err2 = data_dir + dirs[2 * mask_num+1] + "/obj_abs_err_slits_lin.fits"
	
	fnames = [fname1, fname2]
	fnames_err = [fname_err1, fname_err2]

	# ---- Data extraction and preprocessing
	data2Ds = []
	headers = []
	wave_grids = []
	for i in range(2):
		data2D, header = bino_data_preprocess(fnames[i], fnames_err[i])
		wave_grid = wave_grid_extractor(header, data2D[1].shape[1])
		data2Ds.append(data2D)
		headers.append(header)
		wave_grids.append(wave_grid)
		

	detection_arr_list = []
	Nsample = len(data2D) # Number of objects
	detection_wavegrid_list = [[]]
	for k in range(2): # For each of the data
		print("/---- Working on dataset %d" % k)
		Ntargets = 0
		# Select the dataset
		header = headers[k]
		data2D = data2Ds[k]
		wave_grid = wave_grids[k]
		
		# Set up the detection bit array
		detection_arr = np.zeros((Nsample, wave_grid.size), dtype=bool)
		for i in range(1, Nsample):
			wave_tmp = []
			# Obtain the bitcode
			bitcode = int(str(header[i]).split("SLITOBJ")[1].split(" = ")[1].split("'")[1])

			# If neither a standard star or galaxy
			if (bitcode != 2) and (bitcode != 4):
				Ntargets += 1
				print("Working on spectrum %d" % i)        
				# Extract the relevant data.
				imerr = np.copy(data2D[i])

				# At each wavegrid position, perform detection.
				for idx in range(0, wave_grid.size, 10):
					# Extract the appropriate post stamps
					im, err = post_stamp_from_imerr(imerr, idx, remove_outlier = True)
					# Construct signal to noise image
					SN = im/err

					if ((np.abs(SN) < 1e-6).sum() > (0.7 * 32**2)):
						pass
					else:
						detected=False
						# Predict on each post stamp image and save. 
						detected = model.predict(SN.reshape((1, 32, 32, 1)))[0, 0] > thres
						detection_arr[i, idx] = detected
						if detected:
							wave_tmp.append(wave_grid[idx])
			else:
				print("Skipped spectrum %d" % i)
			if k == 0:
				detection_wavegrid_list.append(wave_tmp)
			else:
				detection_wavegrid_list[i]= wave_tmp + detection_wavegrid_list[i]
				
		detection_arr_list.append(detection_arr)

	# ---- Based on first order detections above, create panels. 
	peak_names = ["OII", "Hb", "OIII1", "OIII2"]
	peak_locations = [3727, 4861, 4959, 5007]        
	vmin = -5
	vmax = 5
	search_range = 60

	for num_obj in range(1, len(data2Ds[0]), 1):
		# All the detections collected.
		test_list = detection_wavegrid_list[num_obj]
		test_list.sort()
		num_detects = len(test_list)
		
		if (num_detects > 0): # If there is even a single detection
			print("/---- %d" % num_obj)        
			# ---- Repeat this process for OII, Hb, OIII1, OIII2
			for q, lam in enumerate(peak_names):  
				peak_loc = peak_locations[q]
				peak_name = lam
	#             print(peak_name)
				
				# Initially set the redshift to None
				redz = None 

				# ---- Find best redshift candidate
				TIE = True

				# ---- Find best OII candidate
				# Propose OII peak at each detection points.
				num_hits_max = 0
				best_peak_model = 0
				for peak in test_list:
					if (peak > peak_loc):
	#                     print(peak)
						num_hits = num_other_matches(test_list, peak, peak_proposed=peak_name, tol=100)
						if (num_hits > num_hits_max) and (num_hits > 0):
							num_hits_max = num_hits
							best_peak_model = peak
							TIE = False # No need to search further        
				# If multiple hits with a proposed OII
				if not TIE:
					redz = (best_peak_model/peak_loc) - 1 # Best redz estimate.
				else: # If tied
					redz = (test_list[0]/peak_loc) - 1 # Best redz estimate.

				# Given the redshift perform neural network search around the detection point.
				# Select the point with the highest probabiity.
				# ---- Data to scan
				idx_redz0 = find_nearest_idx(wave_grids[0], best_peak_model)
				idx_redz1 = find_nearest_idx(wave_grids[1], best_peak_model)
				detect0 = detection_arr_list[0][num_obj, idx_redz0]
				detect1 = detection_arr_list[1][num_obj, idx_redz1]

				if detect0: # If detected in the first data, then use this to find the best.
					idx_redz = idx_redz0
					k = 0# Data index 
				else: 
					idx_redz = idx_redz1
					k = 1

				# Probability grid over which to test. 
				idx_arr = np.arange(idx_redz - search_range, idx_redz + search_range, 3, dtype=int)
				prob_arr = np.zeros(idx_arr.size, dtype=float)

				# Extract the relevant data.
				imerr = np.copy(data2Ds[k][num_obj])

				for i, idx in enumerate(idx_arr):
					# Extract the appropriate post stamps
					im, err = post_stamp_from_imerr(imerr, idx, remove_outlier = True)
					# Construct signal to noise image
					SN = im/err

					# Predict on each post stamp image and save the probability.
					prob_arr[i] = model.predict(SN.reshape((1, 32, 32, 1)))[0, 0]

				# Obtain the best estimate of redz based on the above run. 
				redz = (wave_grids[k][idx_arr[np.argmax(prob_arr)]]/peak_loc) - 1


				# ----------------- Make the review plot
				# At each detection point, propose that the detection is OII
				name, idx_list0 = idx_peaks(wave_grids[0], redz)
				_, idx_list1 = idx_peaks(wave_grids[1], redz)            

				# Extract the relevant data
				imerr0 = data2Ds[0][num_obj]
				imerr1 = data2Ds[1][num_obj]                    

				plt.close()
				fig, ax_list = plt.subplots(4, 5, figsize = (8, 6))
				# --First data
				for k, idx in enumerate(idx_list0):
					if idx > 0:
						# Extract the appropriate post stamps
						im, err = post_stamp_from_imerr(imerr0, idx, remove_outlier = True)
						if ((im==0).sum()/32**2) < 0.7:
							# Construct signal to noise image
							SN = im/err
							# Plot
							ax_list[0, k].imshow(im, interpolation="none", cmap="gray", vmin=vmin, vmax=vmax)
							ax_list[0, k].axvline(x=16, c="red", lw=1)
							# Plot
							ax_list[1, k].imshow(SN, interpolation="none", cmap="gray", vmin=vmin, vmax=vmax)
							ax_list[1, k].axvline(x=16, c="red", lw=1)

				# --Second data
				for k, idx in enumerate(idx_list1):
					if idx > 0:
						# Extract the appropriate post stamps
						im, err = post_stamp_from_imerr(imerr1, idx, remove_outlier = True)
						if ((im==0).sum()/32**2) < 0.7:
							# Construct signal to noise image
							SN = im/err
							# Plot
							ax_list[2, k].imshow(im, interpolation="none", cmap="gray", vmin=vmin, vmax=vmax)
							ax_list[2, k].axvline(x=16, c="red", lw=0.5)
							# Plot
							ax_list[3, k].imshow(SN, interpolation="none", cmap="gray", vmin=vmin, vmax=vmax)
							ax_list[3, k].axvline(x=16, c="red", lw=0.5)

				# ---- Deocrations
				for m in range(5):
					ax_list[0, m].set_title(name[m], fontsize=15)
				for n in range(4):
					for m in range(5):
						ax_list[n, m].axis("off") # tick_params(bottom=False, top=False, left=False, right=False)                
				plt.suptitle("%3d    Peak/Num hits/Redz:%5s/%2d/%1.4f" % (num_obj, peak_name, \
								num_hits_max, redz),  fontsize=20)
				plt.savefig(save_dir+"/specnum%d-%s-z%.3f.png" % (num_obj, peak_name, redz), dpi=100, bbox_inches="tight")
				# plt.show()
				plt.close()    