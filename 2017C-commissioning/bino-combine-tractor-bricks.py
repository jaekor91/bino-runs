# Used to combine all tractor bricks relevant for run 1. 

# Loading modules
import numpy as np
from os import listdir
from os.path import isfile, join
from astropy.io import ascii, fits
from astropy.wcs import WCS
import numpy.lib.recfunctions as rec
from xd_elg_utils import *
import sys

large_random_constant = -999119283571
deg2arcsec=3600

dr_v = "5" # Data release version
data_directory = "/Users/jaehyeon/Documents/Research/ELG_target_selection/data-repository/Tractor_binospec_test/run1-Nov2017/"
tycho_directory = "/Users/jaehyeon/Documents/Research/ELG_target_selection/data-repository/"

##############################################################################	
print("Combine all Tractor files and append Tycho-2 stellar mask column.")

trac = combine_tractor_nocut(data_directory)

print("Append Tycho2 stark mask field.")
trac = apply_tycho(trac, tycho_directory+"tycho2.fits", galtype="ELG")
print("Completed.")

##############################################################################	
print("Save the trimmed catalogs.")
cols = fits.ColDefs(trac)
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto("DR5-Tractor-bino-run1.fits", clobber=True)
