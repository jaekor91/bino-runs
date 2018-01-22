import numpy as np
import numpy.lib.recfunctions as rec

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

def flux2mag(flux):
    return 22.5-2.5*np.log10(flux)    

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

def FDR_cut(grz):
    """
    Given a list [g,r,z] magnitudes, apply the cut and return an indexing boolean vector.
    """
    g,r,z=grz; yrz = (r-z); xgr = (g-r)
    ibool = (r<23.4) & (yrz>.3) & (yrz<1.6) & (xgr < (1.15*yrz)-0.15) & (xgr < (1.6-1.2*yrz))
    return ibool


def mag2flux(mag):
    return 10**(0.4*(22.5-mag))
        
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


##############################################################################
# The following is adpated from the URL indicated below.
# """
#     ImagingLSS
#     https://github.com/desihub/imaginglss/blob/master/imaginglss/analysis/tycho_veto.py

#     veto objects based on a star catalogue.
#     The tycho vetos are based on the email discussion at:
#     Date: June 18, 2015 at 3:44:09 PM PDT
#     To: decam-data@desi.lbl.gov
#     Subject: decam-data Digest, Vol 12, Issue 29
#     These objects takes a decals object and calculates the
#     center and rejection radius for the catalogue in degrees.
#     Note : The convention for veto flags is True for 'reject',
#     False for 'preserve'.

#     apply_tycho takes the galaxy catalog and appends a Tychoveto column
#     the code works fine for ELG and LRGs. For other galaxy type, you need to adjust it!
# """

# Import modules
import sys

def BOSS_DR9(tycho):
    bmag = tycho['BMAG']
    # BOSS DR9-11
    b = bmag.clip(6, 11.5)
    R = (0.0802 * b ** 2 - 1.86 * b + 11.625) / 60. #
    return R

def DECAM_LRG(tycho):
    vtmag = tycho['VTMAG']
    R = 10 ** (3.5 - 0.15 * vtmag) / 3600.
    return R

DECAM_ELG = DECAM_LRG

def DECAM_QSO(tycho):
    vtmag = tycho['VTMAG']
    # David Schlegel recommends not applying a bright star mask
    return vtmag - vtmag

def DECAM_BGS(tycho):
    vtmag = tycho['VTMAG']
    R = 10 ** (2.2 - 0.15 * vtmag) / 3600.
    return R

def radec2pos(ra, dec):
    """ converting ra dec to position on a unit sphere.
        ra, dec are in degrees.
    """
    pos = np.empty(len(ra), dtype=('f8', 3))
    ra = ra * (np.pi / 180)
    dec = dec * (np.pi / 180)
    pos[:, 2] = np.sin(dec)
    pos[:, 0] = np.cos(dec) * np.sin(ra)
    pos[:, 1] = np.cos(dec) * np.cos(ra)
    return pos

def tycho(filename):
    """
    read the Tycho-2 catalog and prepare it for the mag-radius relation
    """
    dataf = fits.open(filename)
    data = dataf[1].data
    tycho = np.empty(len(data),
        dtype=[
            ('RA', 'f8'),
            ('DEC', 'f8'),
            ('VTMAG', 'f8'),
            ('VMAG', 'f8'),
            ('BMAG', 'f8'),
            ('BTMAG', 'f8'),
            ('VARFLAG', 'i8'),
            ])
    tycho['RA'] = data['RA']
    tycho['DEC'] = data['DEC']
    tycho['VTMAG'] = data['MAG_VT']
    tycho['BTMAG'] = data['MAG_BT']
    vt = tycho['VTMAG']
    bt = tycho['BTMAG']
    b = vt - 0.09 * (bt - vt)
    v = b - 0.85 * (bt - vt)
    tycho['VMAG']=v
    tycho['BMAG']=b
    dataf.close()
    return tycho


def txts_read(filename):
    obj = np.loadtxt(filename)
    typeobj = np.dtype([
              ('RA','f4'), ('DEC','f4'), ('COMPETENESS','f4'),
              ('rflux','f4'), ('rnoise','f4'), ('gflux','f4'), ('gnoise','f4'),
              ('zflux','f4'), ('znoise','f4'), ('W1flux','f4'), ('W1noise','f4'),
              ('W2flux','f4'), ('W2noise','f4')
              ])
    nobj = obj[:,0].size
    data = np.zeros(nobj, dtype=typeobj)
    data['RA'][:] = obj[:,0]
    data['DEC'][:] = obj[:,1]
    data['COMPETENESS'][:] = obj[:,2]
    data['rflux'][:] = obj[:,3]
    data['rnoise'][:] = obj[:,4]
    data['gflux'][:] = obj[:,5]
    data['gnoise'][:] = obj[:,6]
    data['zflux'][:] = obj[:,7]
    data['znoise'][:] = obj[:,8]
    data['W1flux'][:] = obj[:,9]
    data['W1noise'][:] = obj[:,10]
    data['W2flux'][:] = obj[:,11]
    data['W2noise'][:] = obj[:,12]
    #datas = np.sort(data, order=['RA'])
    return data

def veto(coord, center, R):
    """
        Returns a veto mask for coord. any coordinate within R of center
        is vet.
        Parameters
        ----------
        coord : (RA, DEC)
        center : (RA, DEC)
        R     : degrees
        Returns
        -------
        Vetomask : True for veto, False for keep.
    """
    from sklearn.neighbors import KDTree
    pos_stars = radec2pos(center[0], center[1])
    R = 2 * np.sin(np.radians(R) * 0.5)
    pos_obj = radec2pos(coord[0], coord[1])
    tree = KDTree(pos_obj)
    vetoflag = ~np.zeros(len(pos_obj), dtype='?')
    arg = tree.query_radius(pos_stars, r=R)
    arg = np.concatenate(arg)
    vetoflag[arg] = False
    return vetoflag



def apply_tycho(objgal, tychofn,galtype='LRG'):
    # reading tycho star catalogs
    tychostar = tycho(tychofn)
    #
    # mag-radius relation
    #
    if galtype == 'LRG' or galtype == 'ELG':    # so far the mag-radius relation is the same for LRG and ELG
        radii = DECAM_LRG(tychostar)
    else:
        sys.exit("Check the apply_tycho function for your galaxy type")
    #
    #
    # coordinates of Tycho-2 stars
    center = (tychostar['RA'], tychostar['DEC'])
    #
    #
    # coordinates of objects (galaxies)
    coord = (objgal['ra'], objgal['dec'])
    #
    #
    # a 0.0 / 1.0 array (1.0: means the object is contaminated by a Tycho-2 star, so 0.0s are good)
    tychomask = (~veto(coord, center, radii)).astype('f4')
    objgal = rec.append_fields(objgal, ['TYCHOVETO'], data=[tychomask], dtypes=tychomask.dtype, usemask=False)
    return objgal

def apply_tycho_radec(ra, dec, tychofn,galtype='LRG'):
    """
    Return tycho mask given ra, dec of objects.
    """
    # reading tycho star catalogs
    tychostar = tycho(tychofn)
    #
    # mag-radius relation
    #
    if galtype == 'LRG' or galtype == 'ELG':    # so far the mag-radius relation is the same for LRG and ELG
        radii = DECAM_LRG(tychostar)
    else:
        sys.exit("Check the apply_tycho function for your galaxy type")
    #
    #
    # coordinates of Tycho-2 stars
    center = (tychostar['RA'], tychostar['DEC'])

    # coordinates of objects (galaxies)
    coord = (ra, dec)
    #
    #
    # a 0.0 / 1.0 array (1.0: means the object is contaminated by a Tycho-2 star, so 0.0s are good)
    tychomask = (~veto(coord, center, radii)).astype('f4')
    return tychomask


def apply_tycho_pcat(objgal, tychofn,galtype='LRG'):
    # reading tycho star catalogs
    tychostar = tycho(tychofn)
    #
    # mag-radius relation
    #
    if galtype == 'LRG' or galtype == 'ELG':    # so far the mag-radius relation is the same for LRG and ELG
        radii = DECAM_LRG(tychostar)
    else:
        sys.exit("Check the apply_tycho function for your galaxy type")
    #
    #
    # coordinates of Tycho-2 stars
    center = (tychostar['RA'], tychostar['DEC'])
    #
    #
    # coordinates of objects (galaxies)
    coord = (objgal['RA_DEEP'], objgal['DEC_DEEP'])
    #
    #
    # a 0.0 / 1.0 array (1.0: means the object is contaminated by a Tycho-2 star, so 0.0s are good)
    tychomask = (~veto(coord, center, radii)).astype('f4')
    objgal = rec.append_fields(objgal, ['TYCHOVETO'], data=[tychomask], dtypes=tychomask.dtype, usemask=False)
    return objgal
