from selection_script import *
from utils import *
from os import listdir
from os.path import isfile, join

paths = ["./DR5-tractor/0hr-files/", "./DR5-tractor/3hr-files/", "./DR5-tractor/NGC-files/"]

# Basic strategy
# - Work on one directory at a time.
# - Consolidate resulting targets into a single file.

# counter = 0
targets = []
for path in paths:
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    # For each file
    for f in onlyfiles:
        data = load_fits_table(path+f)
        
        # Make selection corresponding to each option
        iselected1 = apply_selection(path+f, option=1).astype(int)
        iselected2 = apply_selection(path+f, option=2).astype(int)
        iselected3 = apply_selection(path+f, option=3).astype(int)
        
        # Augment the data with three new columns
        data = rec.append_fields(data, ['NDM1'], data=[iselected1], dtypes=iselected1.dtype, usemask=False)
        data = rec.append_fields(data, ['NDM2'], data=[iselected2], dtypes=iselected1.dtype, usemask=False)
        data = rec.append_fields(data, ['NDM3'], data=[iselected3], dtypes=iselected1.dtype, usemask=False)
        
        # Apply Tycho column
        data = apply_tycho(data, "./tycho2.fits", galtype="ELG")
    
        # Finding the union of three selections
        iselected = np.logical_or.reduce((iselected1, iselected2, iselected3))
        
        print f, iselected.sum()
        
        # Take the subset and store      
        targets.append(data[iselected])
        
#         counter += 1
#         if counter == 5:
#             assert False

targets = np.hstack(targets)
save_fits(targets, "DESI-ELG-NDM-targets.fits")