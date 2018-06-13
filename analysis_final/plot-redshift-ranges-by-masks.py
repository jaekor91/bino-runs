from utils import *

mask_dirs = os.listdir("../../data/")


# ---- Histogram of redshift ranges
for k in range(1, len(mask_dirs)//2):
    mask1 = mask_dirs[2 * k]
    mask2 = mask_dirs[2 * k + 1]    
    
    print("/-----", mask1)                
    data_dir1 = "../../data/" + mask1 + "/"
    data_dir2 = "../../data/" + mask2 + "/"    

    # ---- Load 2D data
    data_err1, list_headers1 = preprocess_bino(fname_data, fname_err, data_dir1)
    data_err2, list_headers2 = preprocess_bino(fname_data, fname_err, data_dir2)

    # ---- Extract wavegrid
    wavegrid1 = wavegrid_from_header(list_headers1[1], data_err1[1, 0].shape[1])
    wavegrid2 = wavegrid_from_header(list_headers2[1], data_err2[1, 0].shape[1])    

    # ---- Place holder for redshifts
    z_range = np.zeros((len(list_headers1), 2))

    # ---- Loop over each target and determine redshift range
    for i in range(1, len(list_headers1)):
        data1 = data_err1[i, 0]
        data2 = data_err2[i, 0]        

        # Find min/max index corresponding to the data
        idx_min,  _ = index_edges(data1)
        _,  idx_max = index_edges(data2)        
        z_range[i, 0] = (wavegrid1[idx_min]/3727. - 1)
        z_range[i, 1] = (wavegrid2[idx_max]/3727. - 1)            

    plt.hist(z_range[1:].ravel(), bins=np.arange(-0.01, 2, 0.01))
    plt.savefig("./figures/redshift-range-by-mask/%s.png" % mask1, dpi=200, bbox_inches="tight")
    plt.close()