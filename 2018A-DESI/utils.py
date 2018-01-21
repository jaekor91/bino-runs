import numpy as np
from astropy.io import ascii, fits
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
from astropy import units as u


import matplotlib as mpl
mpl.rcParams['xtick.major.size'] = 15
mpl.rcParams['xtick.major.width'] = 1.
mpl.rcParams['ytick.major.size'] = 15
mpl.rcParams['ytick.major.width'] = 1.
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15

large_random_constant = -999119283571
deg2arcsec=3600

def default_config(num_config):
	"""
	Return default config to be used.
	"""
	if num_config == 1:
		return [189.2669083333333, 62.21212777777778, 45]
	elif num_config == 2:
		return [189.1928166666667, 62.26021944444445, 45]		
	else:
		assert False

	return 

def load_fits_table(fname):
    """Given the file name, load  the first extension table."""
    return fits.open(fname)[1].data

def save_fits(data, fname):
    """
    Given a rec array and a file name (with "fits" filename), save it.
    """
    cols = fits.ColDefs(np.copy(data)) # This is somehow necessary.
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(fname, clobber=True)
    
    return 


def check_astrometry(ra1,dec1,ra2,dec2,pt_size=0.3):
    """
    Given two sets of ra/dec's return median difference in degrees.
	The first is the reference.
    """
    ra_diff = ra2-ra1
    dec_diff = dec2-dec1
    ra_med_diff = np.median(ra_diff)
    dec_med_diff = np.median(dec_diff)
    return ra_med_diff, dec_med_diff



def crossmatch_cat1_to_cat2(ra1, dec1, ra2, dec2, tol=1./(deg2arcsec+1e-12)):
    """
    Return indices of cat1 (e.g., DR3) and cat2 (e.g., DEE2) cross matched to tolerance. 

    Note: Function used to cross-match DEEP2 and DR3 catalogs in each field 
    and test for any astrometric discrepancies. That is, for every object in 
    DR3, find the nearest object in DEEP2. For each DEEP2 object matched, 
    pick DR3 object that is the closest. The surviving objects after these 
    matching process are the cross-matched set.
    """
    
    # Match cat1 to cat2 using astropy functions.
    idx_cat1_to_cat2, d2d = match_cat1_to_cat2(ra1, dec1, ra2, dec2)
    
    # Indicies of unique cat2 objects that were matched.
    cat2matched = np.unique(idx_cat1_to_cat2)
    
    # For each cat2 object matched, pick cat1 object that is the closest. 
    # Skip if the closest objects more than tol distance away.
    idx1 = [] # Place holder for indices
    idx2 = []
    tag = np.arange(ra1.size,dtype=int)
    for e in cat2matched:
        ibool = (idx_cat1_to_cat2==e)
        candidates = tag[ibool]
        dist2candidates = d2d[ibool]
        # Index of the minimum distance cat1 object
        if dist2candidates.min()<tol:
            idx1.append(candidates[np.argmin(dist2candidates)])
            idx2.append(e)
    
    # Turning list of indices into numpy arrays.
    idx1 = np.asarray(idx1)
    idx2 = np.asarray(idx2)
    
    # Return the indices of cat1 and cat2 of cross-matched objects.
    return idx1, idx2



def match_cat1_to_cat2(ra1, dec1, ra2, dec2):
    """
    "c = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
    catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)  

    idx are indices into catalog that are the closest objects to each of the coordinates in c, d2d are the on-sky distances between them, and d3d are the 3-dimensional distances." -- astropy documentation.  

    Fore more information: http://docs.astropy.org/en/stable/coordinates/matchsep.html#astropy-coordinates-matching 
    """    
    cat1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)  
    cat2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)  
    idx, d2d, d3d = cat1.match_to_catalog_sky(cat2)
    
    return idx, d2d.degree

def closest_idx(arr, val):
    return np.argmin(np.abs(arr-val))   

