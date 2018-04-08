#---- Script used to generate blanks and blobs.
from utils import *

Nobjs = 64
z_270, z_600, objnums = parse_tanveer_file()
data2D_270, data2D_270_err, data2D_600, data2D_600_err, wave_grid_270, wave_grid_600 = load_train_data_files()


#---- For each spectrum, plot 2D region where Tanveer found peaks
center_line = 16 
wave_grids = [wave_grid_270, wave_grid_600]
data_list = [data2D_270, data2D_600]
err_list = [data2D_270_err, data2D_600_err]
z_list = [z_270, z_600]

post_stamp_collection = [] # Collect all images of blobs.
err_post_stamp_collection = [] # Collect all images of blobs.


# for i in xrange(1, 5): 
for i in range(1, Nobjs+1): 
    objnum = objnums[i-1] # Get the object number (used to extract 1D/2D data)
    
    # ---- Declare plot
    for l in range(2): # For each of the two coverages
        wave_grid = wave_grids[l]
        data2D = data_list[l]
        err2D = err_list[l]        
        redz = z_list[l]
        
        if redz[i-1] > 0.: # If the redshift recorded    
            plot =True            
            # ---- Names and indices corresponding to various peaks
            names, indices = idx_peaks(wave_grid, redz[i-1])
            for k in range(5):
                idx = indices[k]
                name = names[k]
                if idx > -1: # If the peak is within the wave_grid range.
                    im = post_stamp_from_HDU(data2D, objnum, idx) # Get the post stamp
                    err = post_stamp_from_HDU(err2D, objnum, idx) # Get the post stamp of error file.
                    if (im==0).all(): # If post stamp is all zero, then turn off the axis.
                        pass
                    else: # Otherwise make a plot.
                        # Collect post stamps for visual inspection 
                        post_stamp_collection.append(im)
                        err_post_stamp_collection.append(err)
                        
# ---- Collect the post stamps into a single numpy array
im_arr = np.zeros((len(post_stamp_collection), 32, 32))
err_arr = np.zeros((len(post_stamp_collection), 32, 32))
SN_arr = np.zeros((len(post_stamp_collection), 32, 32))
for i, x in enumerate(post_stamp_collection):
    im_arr[i] = x
    err_arr[i] = err_post_stamp_collection[i]
    SN_arr[i] = im_arr[i] / err_arr[i] # Error is unit of flux.
    SN_arr[i][np.isnan(SN_arr[i])] = 0

# # --- To explore how outliers should be handled.
# plt.hist(SN_arr[44].ravel())
# plt.show()
# plt.close()

# ---- Plot them for visual inspection
plt.close()

flagged = [5, 45, 48, 52, 53, 58, 68, 71]

# flagged = []
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
for i in range(81):
    idx_row = i // 9
    idx_col = i % 9    
    if i < im_arr.shape[0]:
        ax_list[idx_row, idx_col].imshow(im_arr[i], cmap="gray", interpolation ="none")
        if i in flagged:
            ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
        else:
            ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
    ax_list[idx_row, idx_col].axis("off")
plt.savefig("./st82-1hr-sideA/all-peaks.png", dpi=200, bbox_inches ="tight")
# plt.show()
plt.close()


# ---- Corresponding error image
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
for i in range(81):
    idx_row = i // 9
    idx_col = i % 9    
    if i < err_arr.shape[0]:
        ax_list[idx_row, idx_col].imshow(err_arr[i], cmap="gray", interpolation ="none")
        if i in flagged:
            ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
        else:
            ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
    ax_list[idx_row, idx_col].axis("off")
plt.savefig("./st82-1hr-sideA/all-peaks-err.png", dpi=200, bbox_inches ="tight")
# plt.show()
plt.close()



# ---- Corresponding SN image
SN_thres = 30
fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
for i in range(81):
    idx_row = i // 9
    idx_col = i % 9    
    if i < err_arr.shape[0]:
        im  = im_arr[i]
        im[SN_arr[i] > SN_thres] = 0
        ax_list[idx_row, idx_col].imshow(im, cmap="gray", interpolation ="none")
        if i in flagged:
            ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
        else:
            ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
    ax_list[idx_row, idx_col].axis("off")
plt.savefig("./st82-1hr-sideA/all-peaks-SN-corrected.png", dpi=200, bbox_inches ="tight")
# plt.show()
plt.close()


im_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))
err_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))

idx = 0
for i in range(im_arr.shape[0]):
    if i not in flagged:
        im_arr_filtered[idx] = im_arr[i]
        err_arr_filtered[idx] = err_arr[i]
        idx +=1
        
# ---- Plot them for visual inspection
plt.close()

fig, ax_list = plt.subplots(8, 9, figsize=(9, 10))
for i in range(72):
    idx_row = i // 9
    idx_col = i % 9    
    if i < im_arr_filtered.shape[0]:
        ax_list[idx_row, idx_col].imshow(im_arr_filtered[i], cmap="gray", interpolation ="none")
        ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
    ax_list[idx_row, idx_col].axis("off")
plt.savefig("./st82-1hr-sideA/all-peaks-filtered.png", dpi=200, bbox_inches ="tight")
# plt.show()
plt.close()

np.savez("./blob_finder_training_data/real-blobs.npz", image=im_arr_filtered, err=err_arr_filtered)








# Fix the random seed. Generate 10,000 random samples

np.random.seed(0)

# ---- Find random regions where there could be peaks
center_line = 16 
wave_grids = [wave_grid_270, wave_grid_600]
data_list = [data2D_270, data2D_600]
err_list = [data2D_270_err, data2D_600_err]

post_stamp_collection = [] # Collect all images of potential blanks
err_post_stamp_collection = [] # Collect all images of potential blanks


sample_rate = 100 # Number of samples to draw from each spectrum.

for i in range(1, Nobjs+1): 
    objnum = objnums[i-1] # Get the object number (used to extract 1D/2D data)
    
    # ---- Declare plot
    for l in range(2): # For each of the two coverages
        wave_grid = wave_grids[l]
        data2D =data_list[l]
        err2D = err_list[l]
        indices = np.random.choice(np.arange(16, wave_grid.size-16, 1, dtype=int), sample_rate, replace=False)
        
        # Skip the whole process if 
        num = int(str(data2D[objnum].header).split("SLITOBJ =")[1].split("/ object")[0].split("'")[1])
        if (num == 2) or (num==4):
            pass
        else:
            for idx in indices:
                im = post_stamp_from_HDU(data2D, objnum, idx) # Get the post stamp
                err = post_stamp_from_HDU(err2D, objnum, idx) # Get the post stamp for errors
                if ((im==0).sum() / 32**2) > 0.8: # If post stamp is all zero, then turn off the axis.
                    pass
                else: # Otherwise make a plot.
                    post_stamp_collection.append(im)
                    err_post_stamp_collection.append(err)

# ---- Collect the post stamps into a single numpy array
im_arr = np.zeros((len(post_stamp_collection), 32, 32))
err_arr = np.zeros((len(post_stamp_collection), 32, 32))
for i, x in enumerate(post_stamp_collection):
    im_arr[i] = x
    err_arr[i] = err_post_stamp_collection[i]
# im_arr = im_arr[:4050] # Make hard cut at 4050
# err_arr = err_arr[:4050] # Make hard cut at 4050
np.savez("./blob_finder_training_data/blanks-all.npz", image=im_arr, err=err_arr)






with open("./blank_rejects_indices.txt") as f:
    line = f.readline()
flagged = list(set([int(x) for x in line.split(",")]))
flagged.sort()
# flagged = []


data = np.load("./blob_finder_training_data/blanks-all.npz")
im_arr = data["image"]
err_arr = data["err"]


# ---- Plot them for visual inspection
plt.close()
for l in range(im_arr.shape[0]//81):
    i_start = l * 81
    i_end = i_start + 81
    fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
    for i in range(i_start, i_end):
        idx_row = (i-i_start) // 9
        idx_col = (i-i_start) % 9    
        if i < im_arr.shape[0]:
            im = np.copy(im_arr[i])
            SN = im/err_arr[i]
            SN[np.isnan(SN)] = 0
            im[SN>SN_thres] = 0
            ax_list[idx_row, idx_col].imshow(im, cmap="gray", interpolation ="none")
            if i in flagged:
                ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
            else:
                ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
        ax_list[idx_row, idx_col].axis("off")
    plt.savefig("./st82-1hr-sideA/blanks-%d.png" % l, dpi=200, bbox_inches ="tight")
#     plt.show()
    plt.close()

# ---- Plot them for visual inspection -- err
plt.close()
for l in range(im_arr.shape[0]//81):
    i_start = l * 81
    i_end = i_start + 81
    fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
    for i in range(i_start, i_end):
        idx_row = (i-i_start) // 9
        idx_col = (i-i_start) % 9    
        if i < im_arr.shape[0]:
            ax_list[idx_row, idx_col].imshow(err_arr[i], cmap="gray", interpolation ="none")
            if i in flagged:
                ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
            else:
                ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
        ax_list[idx_row, idx_col].axis("off")
    plt.savefig("./st82-1hr-sideA/blanks-%d-err.png" % l, dpi=200, bbox_inches ="tight")
#     plt.show()
    plt.close()




# ---- Sift out the flagged onese
im_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))
err_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))
counter = 0
for l in range(im_arr.shape[0]):
    if l not in flagged:
        im_arr_filtered[counter] = im_arr[l]
        err_arr_filtered[counter] = err_arr[l]        
        counter += 1
        


# ---- Plot them for visual inspection
plt.close()

for l in range(im_arr_filtered.shape[0]//81):
    i_start = l * 81
    i_end = i_start + 81
    fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
    for i in range(i_start, i_end):
        idx_row = (i-i_start) // 9
        idx_col = (i-i_start) % 9    
        if i < im_arr_filtered.shape[0]:
            im = np.copy(im_arr_filtered[i])
            SN = im/err_arr_filtered[i]
            SN[np.isnan(SN)] = 0
            im[SN>SN_thres] = 0            
            ax_list[idx_row, idx_col].imshow(im, cmap="gray", interpolation ="none")
            ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
        ax_list[idx_row, idx_col].axis("off")
    plt.savefig("./st82-1hr-sideA/blanks-filtered-%d.png" % l, dpi=200, bbox_inches ="tight")
#     plt.show()
    plt.close()

np.savez("./blob_finder_training_data/blanks-filtered.npz", image=im_arr_filtered, err=err_arr_filtered)