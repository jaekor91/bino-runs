from utils import *

Nobjs = 64
z_270, z_600, objnums = parse_tanveer_file()
data2D_270, data2D_600, data1D_270, data1D_600, wave_grid_270, wave_grid_600 = load_train_data_files()

#---- For each spectrum, plot 2D region where Tanveer found peaks
center_line = 16 
wave_grids = [wave_grid_270, wave_grid_600]
data_list = [data2D_270, data2D_600]
z_list = [z_270, z_600]

# for i in xrange(1, 5): 
for i in xrange(1, Nobjs+1): 
    plt.close()
    plot = False
    objnum = objnums[i-1] # Get the object number (used to extract 1D/2D data)
    
    # ---- Declare plot
    fig, ax_list = plt.subplots(2, 5, figsize=(15, 6))    
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
                        ax_list[l, k].axis("off")
                    else: # Otherwise make a plot.
                        ax_list[l, k].imshow(im, aspect="equal", cmap="gray", interpolation="none")
                        ax_list[l, k].axvline(x=center_line, c="red")
                        ax_list[l, k].set_title("%5s" % name, fontsize=12)
                else: # If peak not found within the wavegrid range.
                    ax_list[l, k].axis("off")                    
                    
        else: # If redshift not recorded
            for k in xrange(5):
                ax_list[l, k].axis("off")

    
    if plot:
        title_str = "Obj Num: %3d \t z_600: %.3f \t z_270: %.3f" % (objnum, z_600[i-1], z_270[i-1])
        plt.suptitle(title_str, fontsize=15, y=1.01)
        plt.savefig("./st82-1hr-sideA/peaks/objnum-%04d-peaks.png" % (objnum), dpi=200, bbox_inches="tight")
        # plt.show()
    plt.close()
