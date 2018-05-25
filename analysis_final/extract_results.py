from utils import *

import os 
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
    
    # Collect analyzed redshift info    
    dir_current = "review_panel/" + dirs[2*q]
    dir_conf = dir_current + "/multiple/"
    dir_single = dir_current + "/single/"
    dir_error = dir_current + "/multiple/error/"
    
    confidence_list = []
    error_list = []
    single_list = []
    for x in os.listdir(dir_conf):
        if (x.endswith(".png")):
            if (x not in specnum_list):
                confidence_list.append(x)
            
    for x in os.listdir(dir_single):
        if (x.endswith(".png")):
            if (x not in specnum_list):
                single_list.append(x)

    if os.path.exists(dir_error):
        for x in os.listdir(dir_error):
            if (x.endswith(".png")):
                if (x not in specnum_list):
                    error_list.append(x)
    
    # ----- Retrieving data for header info
    mask_name = dirs[2*q]
    print(q, mask_name)
    dir_data = data_dir +  mask_name + "/"

    fname_data = "obj_abs_slits_lin.fits"
    fname_err = "obj_abs_err_slits_lin.fits"

    data_err, list_headers  = preprocess_bino(fname_data, fname_err, dir_data)
    
    Nobjs = len(list_headers)
    # Export RA, DEC, BIT, Redshift, Confidence

    RA = np.ones(Nobjs, dtype=float) * -999
    DEC = np.ones(Nobjs, dtype=float) * -999
    BIT = np.ones(Nobjs, dtype=int) * -999
    REDZ = np.ones(Nobjs, dtype=float) * -999
    CONFIDENCE = np.ones(Nobjs, dtype=int) * -999

    for specnum in range(1, Nobjs):
        data, err, header = extract_single_data(data_err, list_headers, specnum)
        RA[specnum] = header["RA"]
        DEC[specnum] = header["DEC"]
        BIT[specnum] = bit_from_header(header)

        # Iterate through the entire review panel list and check if there is the redshift you want. 
        # "High confidence list"
        for i, x in enumerate(confidence_list):
            num = int(x[7:10])
            if specnum == num:
                REDZ[specnum] = float(x.split("z")[1][:-4])
                CONFIDENCE[specnum] = 2

        # Single 
        for i, x in enumerate(single_list):
            num = int(x[7:10])
            if specnum == num:
                REDZ[specnum] = float(x.split("z")[1][:-4])
                CONFIDENCE[specnum] = 1

        # Error
        for i, x in enumerate(error_list):
            num = int(x[7:10])
            if specnum == num:
                z_wrong = float(x.split("z")[1][:-4])
                REDZ[specnum] = z_wrong
                CONFIDENCE[specnum] = 0
                
    # Save the results
    np.savez("./result/" + mask_name + "-results.npz", RA=RA, DEC=DEC, REDZ=REDZ, CONFIDENCE=CONFIDENCE, BIT=BIT)