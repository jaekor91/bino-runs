from utils import *

mask_dirs = os.listdir("../../data/") # Mask for which we have data.
fname_data = "obj_abs_slits_lin.fits"
fname_err = "obj_abs_err_slits_lin.fits"

# ---- Mask set name for which we have redshift
dir_redz = "../analysis/" # Redshift results directory.
mask_analyzed = []
for fname in os.listdir(dir_redz):
    if fname.startswith("results-"):
        mask_analyzed.append(fname)
        
# ---- For each mask set, record the proper name in the data directory
# Ugly piece of code but works.
mask_analyzed_dict = {}
for mask in mask_dirs:
    for fname in mask_analyzed:
        if fname[8:-8] in mask:
            if fname not in mask_analyzed_dict.keys():
                mask_analyzed_dict[fname] = [mask]
            else:
                mask_analyzed_dict[fname].append(mask)
                
# ---- For each mask set that has been analyzed iterate through all masks and collect the blobs. 
counter = 0 
counter_max = 3000
post_stamp_blobs = np.zeros((counter_max, 2, 32, 32))

for fname in mask_analyzed_dict.keys():
    print("Mask set nmae: ", fname)
    data = np.load(dir_redz + fname)
    redz = data["redz"]
    
    for mask in mask_analyzed_dict[fname]:
        print("mask-grating", mask)
        data_dir = "../../data/" + mask + "/"
        
        # ---- Import data
        data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)
        Nspecs = data_err.shape[0]        
        
        # Compute the wavegrid based on the header
        data, err, header = extract_single_data(data_err, list_headers, specnum=1)
        wavegrid = wavegrid_from_header(header, Ncols=data.shape[1])

        specnum = 1
        while (specnum < (Nspecs)) and (counter < counter_max):
            z = redz[specnum]
            if z > 0.: 
#                 print(specnum, z)
                # ---- Show stamps corresponding to various peaks.
                # Extract data
                data, err, header = extract_single_data(data_err, list_headers, specnum)

                # Compute edge index
                idx_min, idx_max = index_edges(data)

                # Compute the location of various peaks based on redshift.
                peak_names, idx_locs = idx_peaks(wavegrid, z, idx_min = idx_min, idx_max = idx_max)

                # Save ones that were used to determine redshift.
                for i, idx in enumerate(idx_locs):
                    if idx > 0: 
                        # ---- Collection happens here
                        post_stamp_blobs[counter, 0] = data[:, idx-16:idx+16]
                        post_stamp_blobs[counter, 1] = err[:, idx-16:idx+16]
#                         plt.close()
#                         fig, ax = plt.subplots(1, figsize=(2, 2))                
#                         ax.imshow(, cmap="gray", interpolation="None")
#                         ax.axis("off")
#                         plt.show()
#                         plt.close()
                        counter += 1
            specnum += 1
        
post_stamp_blobs = post_stamp_blobs[:counter, :, :, :]
# ---- Save
np.save("./train_data/blobs_unfiltered.npy", post_stamp_blobs)

# ---- Plot
for i in range(counter//100):
    stamps = post_stamp_blobs[i*100:(i+1)*100, 0] / post_stamp_blobs[i*100:(i+1)*100, 1]
    plot_post_stamps(stamps, num_start=i*100, fname="./train_data/blobs/blobs-%02d.png" % i)
    
    
# ----- First time filtering

rejects1 = [2, 3, 5, 7, 8, 12, 15, 18, 22, (37, 41), 46, (48, 54), 57, 58, (60, 65), (70, 75), 77, 78, 80, (86, 90), (95, 99), 100, 102, (105, 107), (111, 113), (115, 118), (120, 122), (126, 131), 138, (145, 152), (155, 156), (159, 164), 166, (169, 174), (176, 178), 196, 201, 203, (209, 212), 215, 224, 230, 240, 241, (245, 247), (249, 250), 252, 253, 261, 262, 265, 276, (279, 282), 288, 289, 292, (296, 299), (301, 304), (308, 314), 318, 319, (329, 333), (392, 394), (351, 355), (360, 374), (379, 383), 386, (388, 390), (392, 394), 397, 407, 409, 417, 421, (440, 443), (447, 448), 461, 462, (466, 468), (480, 482), 487, 488, 491, (519, 521), (530, 533), 536, 538, (540, 542), (545, 548), (553, 555), 557, 570, 571, 574, 577, 586, 587, 598, 599, (603, 606), 614, 619, (621, 624), 629, (632, 635), 642, 648, 649, 653, 658, 664, 665, (723, 727), (733, 735), 737, 739, (757, 760), (800, 803), (805, 809), (821, 823), (835, 837), 848, 858, 859, (886, 889), 893, 895, 900, 903, 906, (909, 916), 918, 920, 922, 925, 932, (935, 951), (953, 959), (963, 966), 968, 974, 975, 978, 981, (983, 989), (993, 995), (997, 998), (1000, 1006), 1008, 1010, 1013, (1015, 1024), (1030, 1032), (1034, 1039), (1041, 1043), (1045, 1047), (1052, 1055), (1060, 1062), (1073, 1075), (1077, 1084), (1087, 1089), 1092, 1093, 1098, 1099, 1100, (1106, 1107), 1116, 1119, 1120, 1122, (1128, 1141), 1144, (1148, 1189), (1193, 1224), (1226, 1249), (1260, 1266), 1279, (1289, 1295), (1305, 1307), (1315, 1318), 1324, 1325, (1328, 1331), (1333, 1341), (1347, 1351), 1353, (1360, 1362), (1367, 1370), 1373, 1377, (1380, 1382), (1384, 1391)]
post_stamps_blob = np.load("./train_data/blobs_unfiltered.npy")

iselect = np.ones(post_stamps_blob.shape[0], dtype=bool)

for x in rejects1:
    if type(x) is tuple:
        idx_low, idx_high = x
        assert idx_low < idx_high
        iselect[idx_low:idx_high+1] = False
    else:
        assert type(x) is int
        iselect[x] = False
# Reject uninspected.
iselect[1400:] = False

post_stamps_blob_ft1 = post_stamps_blob[iselect]
np.save("./train_data/blobs_filtered1.npy", post_stamps_blob_ft1)

# ---- Plot
for i in range(post_stamps_blob_ft1.shape[0] // 100):
    stamps = post_stamps_blob_ft1[i*100:(i+1)*100, 0] / post_stamps_blob_ft1[i*100:(i+1)*100, 1]
    plot_post_stamps(stamps, num_start=i*100, fname="./train_data/blobs/blobs-filtered1-%02d.png" % i)
    
    
    
# ----- Second time filtering
rejects2 = [6, 37, 44, 51, 149, (151,156), (161, 162), (186, 194), 223, 225, (230,232), 240, 249, 259, 326, 363, 413, 414, (421,426), 435, 436, 452, 453, (460,464), 512, 514, 518, 519, 530, 539, 553, 564, 580, 640, 641, 645]
iselect = np.ones(post_stamps_blob_ft1.shape[0], dtype=bool)
for x in rejects2:
    if type(x) is tuple:
        idx_low, idx_high = x
        assert idx_low < idx_high
        iselect[idx_low:idx_high+1] = False
    else:
        assert type(x) is int
        iselect[x] = False

post_stamps_blob_ft2 = post_stamps_blob_ft1[iselect]
np.save("./train_data/blobs_filtered2.npy", post_stamps_blob_ft2)

# ---- Plot
for i in range(post_stamps_blob_ft2.shape[0] // 100):
    stamps = post_stamps_blob_ft2[i*100:(i+1)*100, 0] / post_stamps_blob_ft2[i*100:(i+1)*100, 1]
    plot_post_stamps(stamps, num_start=i*100, fname="./train_data/blobs/blobs-filtered2-%02d.png" % i)

# ---- Looking at the last few
counter = 700

for i in range(counter, post_stamps_blob_ft2.shape[0]):
    image = post_stamps_blob_ft2[i, 0]/ post_stamps_blob_ft2[i, 1]
    plt.close()
    fig, ax = plt.subplots(1, figsize=(3, 3))
    ax.imshow(image, aspect="auto", cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
    # plt.show()
    plt.close()