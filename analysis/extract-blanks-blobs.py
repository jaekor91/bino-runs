#---- Script used to generate blanks and blobs.
from utils import *

Nobjs = 64
z_270, z_600, objnums = parse_tanveer_file()
data2D_270, data2D_600, data1D_270, data1D_600, wave_grid_270, wave_grid_600 = load_train_data_files()

#---- For each spectrum, plot 2D region where Tanveer found peaks
center_line = 16 
wave_grids = [wave_grid_270, wave_grid_600]
data_list = [data2D_270, data2D_600]
z_list = [z_270, z_600]

post_stamp_collection = [] # Collect all images of blobs.

# for i in xrange(1, 5): 
for i in xrange(1, Nobjs+1): 
    objnum = objnums[i-1] # Get the object number (used to extract 1D/2D data)
    
    # ---- Declare plot
    for l in xrange(2): # For each of the two coverages
        wave_grid = wave_grids[l]
        data2D =data_list[l]
        redz = z_list[l]
        
        if redz[i-1] > 0.: # If the redshift recorded    
            plot =True            
            # ---- Names and indices corresponding to various peaks
            names, indices = idx_peaks(wave_grid, redz[i-1])
            for k in xrange(5):
                idx = indices[k]
                name = names[k]
                if idx > -1: # If the peak is within the wave_grid range.
                    im = post_stamp_from_HDU(data2D, objnum, idx) # Get the post stamp
                    if (im==0).all(): # If post stamp is all zero, then turn off the axis.
                        pass
                    else: # Otherwise make a plot.
                        # Collect post stamps for visual inspection 
                        post_stamp_collection.append(im)
                        
                        
# ---- Collect the post stamps into a single numpy array
im_arr = np.zeros((len(post_stamp_collection), 32, 32))
for i, x in enumerate(post_stamp_collection):
    im_arr[i] = x
                      
        
# ---- Plot them for visual inspection
plt.close()

flagged = [5, 45, 48, 52, 53, 58, 68, 71]

fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
for i in xrange(81):
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


im_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))
idx = 0
for i in xrange(im_arr.shape[0]):
    if i not in flagged:
        im_arr_filtered[idx] = im_arr[i]
        idx +=1
        
# ---- Plot them for visual inspection
plt.close()

fig, ax_list = plt.subplots(8, 9, figsize=(9, 10))
for i in xrange(72):
    idx_row = i // 9
    idx_col = i % 9    
    if i < im_arr_filtered.shape[0]:
        ax_list[idx_row, idx_col].imshow(im_arr_filtered[i], cmap="gray", interpolation ="none")
        ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
    ax_list[idx_row, idx_col].axis("off")
plt.savefig("./st82-1hr-sideA/all-peaks-filtered.png", dpi=200, bbox_inches ="tight")
# plt.show()
plt.close()

np.savez("./blob_finder_training_data/real-blobs.npz", image=im_arr_filtered)
# data = np.load("./blob_finder_training_data/real-blobs.npz")





# !!!!!!!!!!! Do not re-sample blanks as they have already been tagged !!!!!!!!!!!!!!!! #
# # ---- Find random regions where there could be peaks
# center_line = 16 
# wave_grids = [wave_grid_270, wave_grid_600]
# data_list = [data2D_270, data2D_600]

# post_stamp_collection = [] # Collect all images of potential blanks

# sample_rate = 48 # Number of samples to draw from each spectrum.

# for i in xrange(1, Nobjs+1): 
#     objnum = objnums[i-1] # Get the object number (used to extract 1D/2D data)
    
#     # ---- Declare plot
#     for l in xrange(2): # For each of the two coverages
#         wave_grid = wave_grids[l]
#         data2D =data_list[l]
#         indices = np.random.choice(np.arange(16, wave_grid.size-16, 1, dtype=int), sample_rate, replace=False)
        
#         for idx in indices:
#             im = post_stamp_from_HDU(data2D, objnum, idx) # Get the post stamp
#             if (im==0).all(): # If post stamp is all zero, then turn off the axis.
#                 pass
#             else: # Otherwise make a plot.
#                 post_stamp_collection.append(im)

# # ---- Collect the post stamps into a single numpy array
# im_arr = np.zeros((len(post_stamp_collection), 32, 32))
# for i, x in enumerate(post_stamp_collection):
#     im_arr[i] = x
# im_arr = im_arr[:4050] # Make hard cut at 4050
# np.savez("./blob_finder_training_data/blanks-all.npz", image=im_arr)



with open("./blob_finder_training_data/blanks-flag.txt") as f:
    line = f.readline()
flagged = list(set([int(x) for x in line.split(",")]))



data = np.load("./blob_finder_training_data/blanks-all.npz")
im_arr = data["image"]

# ---- Plot them for visual inspection
plt.close()

for l in xrange(50):
    i_start = l * 81
    i_end = i_start + 81
    fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
    for i in xrange(i_start, i_end):
        idx_row = (i-i_start) // 9
        idx_col = (i-i_start) % 9    
        if i < im_arr.shape[0]:
            ax_list[idx_row, idx_col].imshow(im_arr[i], cmap="gray", interpolation ="none")
            if i in flagged:
                ax_list[idx_row, idx_col].set_title(i, fontsize=10, color="red")          
            else:
                ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
        ax_list[idx_row, idx_col].axis("off")
    plt.savefig("./st82-1hr-sideA/blanks-%d.png" % l, dpi=200, bbox_inches ="tight")
#     plt.show()
    plt.close()


# ---- Sift out the flagged onese
im_arr_filtered = np.zeros((im_arr.shape[0] - len(flagged), 32, 32))
counter = 0
for l in xrange(im_arr.shape[0]):
    if l not in flagged:
        im_arr_filtered[counter] = im_arr[l]
        counter += 1
        


# ---- Plot them for visual inspection
plt.close()

for l in xrange(im_arr_filtered.shape[0]//81):
    i_start = l * 81
    i_end = i_start + 81
    fig, ax_list = plt.subplots(9, 9, figsize=(10, 10))
    for i in xrange(i_start, i_end):
        idx_row = (i-i_start) // 9
        idx_col = (i-i_start) % 9    
        if i < im_arr_filtered.shape[0]:
            ax_list[idx_row, idx_col].imshow(im_arr_filtered[i], cmap="gray", interpolation ="none")
            ax_list[idx_row, idx_col].set_title(i, fontsize=5, color="black")                      
        ax_list[idx_row, idx_col].axis("off")
    plt.savefig("./st82-1hr-sideA/blanks-filtered-%d.png" % l, dpi=200, bbox_inches ="tight")
#     plt.show()
    plt.close()

np.savez("./blob_finder_training_data/blanks-filtered.npz", image=im_arr_filtered)