from utils import *

import os 
# ---- Mask directory names
data_dir = "../../data/"
dirs = os.listdir(data_dir)

# ---- Bino file names
fname_data = "obj_abs_slits_lin.fits"
fname_err = "obj_abs_err_slits_lin.fits"

# ---- For each set
for q in range(1, len(dirs)//2):
    
    # Collect analyzed redshift info    
    dir_current = "review_panel/" + dirs[2*q]
    
    confidence_list = []
    for x in os.listdir(dir_current):
        if (x.endswith(".png")):
            confidence_list.append(x)

    print(confidence_list)
    print("\n" * 5)
            
    # ----- Retrieving data for header info
    mask_name = dirs[2*q]
    print(q, mask_name)
    dir_data = data_dir +  mask_name + "/"

    data_err, list_headers  = preprocess_bino(fname_data, fname_err, dir_data)
    
    Nobjs = len(list_headers)
    # Export RA, DEC, BIT, Redshift, Confidence

    RA = np.ones(Nobjs, dtype=float) * -999
    DEC = np.ones(Nobjs, dtype=float) * -999
    BIT = np.ones(Nobjs, dtype=int) * -999
    REDZ = np.ones(Nobjs, dtype=float) * -999
    # CONFIDENCE = np.ones(Nobjs, dtype=int) * -999 # Appended later.

    for specnum in range(1, Nobjs):
        data, err, header = extract_single_data(data_err, list_headers, specnum)
        RA[specnum] = header["RA"]
        DEC[specnum] = header["DEC"]
        BIT[specnum] = bit_from_header(header)

        # Iterate through the entire review panel list and check if there is the redshift you want. 
        for i, x in enumerate(confidence_list):
            num = int(x.split("-")[0][-3:])
            z = float(x.split("z")[1].split("png")[0][:-1])
            if specnum == num:
                REDZ[specnum] = float(x.split("z")[1][:-4])
                
    # Save the results
    np.savez("./results/" + mask_name + "-results.npz", RA=RA, DEC=DEC, REDZ=REDZ, BIT=BIT)