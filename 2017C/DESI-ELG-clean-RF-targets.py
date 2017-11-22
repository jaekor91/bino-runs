from utils import *
from os import listdir
from os.path import isfile, join

paths = ["./DR5-tractor/0hr-files/", "./DR5-tractor/3hr-files/", "./DR5-tractor/NGC-files/"]


def wrap_ra(ra, new_break_point = 180):
    """
    Given ra [0, 360], wrap it such that the resulting objects have ra [-(360-new_break_point), new_break_point]
    """
    ibool = ra>new_break_point
    ra[ibool] = ra[ibool]-360.

    return ra

RF_files = ["List_RF_0hr.fits", "List_RF_3hr.fits", "List_RF_8hr30.fits"]
# Concatenate all RF fits files
RF_targets = []
for f in RF_files:
    RF_targets.append(load_fits_table(f))
RF_targets = np.hstack(RF_targets)
ra_RF, dec_RF = RF_targets["ra"], RF_targets["dec"]

RF_targets_final = []
for path in paths:
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    # For each file
    for f in onlyfiles:
        data = load_fits_table(path+f)

        # mask conditions
        ibool = (data["brick_primary"]>0) & \
        (data["allmask_g"]==0) & (data["allmask_r"]==0) & (data["allmask_z"]==0) &\
        (data["flux_ivar_g"] >0) & (data["flux_ivar_r"] >0) & (data["flux_ivar_z"] >0)

        # Considering only those that pass the condition.
        data = data[ibool]
        
        # For each tractor file, compute min/max ra/dec and select RF selected objects that fall within the box.
        ra, dec = data["ra"], data["dec"]
        ra_min, ra_max, dec_min, dec_max = ra.min(), ra.max(), dec.min(), dec.max()
        
        RF_tmp = RF_targets[(ra_RF > ra_min) & (ra_RF < ra_max) & (dec_RF > dec_min) & (dec_RF < dec_max)]
        RF_tmp_ra = RF_tmp["ra"]
        
        # Cross-match (and there should be ~100% match) the two catalogs and make sure to enforce all_mask, 
        ibool  =np.in1d(RF_tmp_ra, ra)
        
        if RF_tmp_ra.size > 0:
            RF_targets_final.append(RF_tmp[ibool])
        
            print f
            print "# in brick; # selected; # RF matched"
            print ra.size, RF_tmp.size, ibool.sum()
            print "\n"        

            if RF_tmp.size * 0.9 > ibool.sum():
                print "Too many objects dropped."
                assert False

# Concatenate all targets
RF_targets_final = np.hstack(RF_targets_final)
        
# Apply Tycho column and save
RF_targets_final = apply_tycho(RF_targets_final, "./tycho2.fits", galtype="ELG")

np.save("DESI-ELG-RF-target.npy", RF_targets_final)
