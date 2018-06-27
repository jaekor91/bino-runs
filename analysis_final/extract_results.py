# ------ Confidence assignment for non-3 level objects
# '2-8h30m-270',
list0 = [(11, 2), (16, 2), (19, 2), (30, 1), (65, 1), (72, 2), (78, 2), (80, 2), (83, 2), \
 (88, 2), (108, 2), (111, 2), (115, 2), (131, 1), (133, 2)]
#  'Eisenstein-COSMOS-1_325_6500',
list1 = [(32, 2), (45, 1), (46, 1), (47, 1), (72, 1), (103, 2), (158, 2)]
#  'Eisenstein-COSMOS-2_326_6500',
list2 = [(45, 1), (89, 1), (91, 1), (111, 1)]
#  'Eisenstein-DR6-2_334_6500',
list3 = [(42, 2), (46, 1), (102, 1), (109, 1), (112, 1), (140, 1)]
#  'Eisenstein-DR6-3_335_6500',
list4 = [(20, 1), (30, 1), (38, 2), (39, 2), (143, 1), (145, 1), (149, 2), (151, 1)] 
#  'Eisenstein-DR6-5_338_6500',
list5 = []
#  'Eisenstein-DR6-6_339_6500',
list6 = [(120, 2), (124, 2), (126, 1)]
#  'ngc-1-1_600_6300',
list7 = [(14, 1), (33, 2), (34, 2), (122, 2), (136, 1)]
#  'NGC-3-1_317_6500',
list8 = [(34, 1), (68, 1), (119, 2)]
#  'ngc-5-1_600_6300',
list9 = []
#  'NGC-6-1_318_6500',
list10 = []
#  'NGC-7-1-updated_319_6500',
list11 = [(12, 1), (24, 2), (36, 1), (49, 1), (123, 2), (125, 2)]
#  'sgc-0hr-1-1_600_6300',
list12 = [(1, 1), (5, 2), (90, 1), (96, 1), (101, 1)]
#  'sgc-3hr-6-1_600_6300',
list13 = [(10, 2), (11, 1), (12, 1), (16, 2), (20, 1), (32, 1), (52, 1), (77, 2), (95, 1), (98, 1), (102, 1),\
         (104, 1), (106, 2), (108, 1), (129, 2)]
#  'st82-1hr-270',
list14 = [(28, 1), (32, 1), (38, 1), (86, 1), (118, 1)]
#  'st82-3hr-125-270'
list15 = [(24, 1)]

master_list = [None, list0, list1, list2, list3, list4, list5, list6, list7, \
               list8, list9, list10, list11, list12, list13, list14, list15]

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
    list_conf = master_list[q] # Non level 3 confidence
    print("Current directory name", dir_current)
    print("List of non-3 confidence", list_conf)

    
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
    CONFIDENCE = np.ones(Nobjs, dtype=int) * -999 # Appended later.

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
                conf = 3 # Default confidence level
                for tags in list_conf:
                    sn, lvl = tags
                    if specnum == sn:
                        conf = lvl
                        break
                CONFIDENCE[specnum] = conf
                
    # Save the results
    np.savez("./results/" + mask_name + "-results.npz", RA=RA, DEC=DEC, REDZ=REDZ, CONFIDENCE=CONFIDENCE, BIT=BIT)