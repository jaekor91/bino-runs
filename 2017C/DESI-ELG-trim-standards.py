from utils import *
from os import listdir
from os.path import isfile, join
from astropy.io import ascii

paths = ["./DR5-tractor/0hr-files/", "./DR5-tractor/3hr-files/", "./DR5-tractor/NGC-files/"]


def wrap_ra(ra, new_break_point = 180):
    """
    Given ra [0, 360], wrap it such that the resulting objects have ra [-(360-new_break_point), new_break_point]
    """
    ibool = ra>new_break_point
    ra[ibool] = ra[ibool]-360.

    return ra


# Galaxies
# Extracting columns of interest
BOSS_files = ["boss.r0_45.dm6_p6.txt", "boss.r110_170.dm10_p35.txt"]
ra_gal = []
dec_gal = []
gmag_gal = []
for f in BOSS_files:
    data = ascii.read(f)
    ra_gal.append(data["col1"].data)
    dec_gal.append(data["col2"].data)
    gmag_gal.append(data["col4"].data)
ra_gal = np.concatenate(ra_gal)
dec_gal = np.concatenate(dec_gal)
gmag_gal = np.concatenate(gmag_gal)


# Stars
# Extracting columns of interest
Star_files = ["standard_bino.r0_45.dm6_p6.txt", "standard_bino.r110_170.dm10_p35.txt"]
ra_star = []
dec_star = []
gmag_star = []
for f in Star_files:
    data = ascii.read(f)
    ra_star.append(data["col1"].data)
    dec_star.append(data["col2"].data)
    gmag_star.append(data["col4"].data)
ra_star = np.concatenate(ra_star)
dec_star = np.concatenate(dec_star)
gmag_star = np.concatenate(gmag_star)


# Put galaxies and stars into N x 3 arrays
gal_array = np.array([ra_gal, dec_gal, gmag_gal]).T
star_array = np.array([ra_star, dec_star, gmag_star]).T

star_final = []
gal_final = []

for path in paths:
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    # For each file
    for f in onlyfiles:
        print f
        data = load_fits_table(path+f)
        
        # For each tractor file, compute min/max ra/dec
        ra, dec = data["ra"], data["dec"]
        ra_min, ra_max, dec_min, dec_max = ra.min(), ra.max(), dec.min(), dec.max()
        
        # Select galaxies and stars within the box specified by the boundary above
        iselect_gal = (ra_gal > ra_min) & (ra_gal < ra_max) & (dec_gal > dec_min) & (dec_gal < dec_max)
        iselect_star = (ra_star > ra_min) & (ra_star < ra_max) & (dec_star > dec_min) & (dec_star < dec_max)
        Nstar = iselect_star.sum()  
        Ngal = iselect_gal.sum()
        print "N star/gal: %d/%d" % (Nstar, Ngal)
        print"\n"        
        
        # Make astrometric correction        
        # Star
        if Nstar > 0:
            star_tmp = star_array[iselect_star]        
            star_final.append(star_tmp)
        
        # Gal
        if Ngal>0:
            gal_tmp = gal_array[iselect_gal]                
            gal_final.append(gal_tmp)

# Concatenate all targets
gal_final = np.vstack(gal_final)
np.save("DESI-ELG-boss.npy", gal_final)
star_final = np.vstack(star_final)
np.save("DESI-ELG-standards.npy", star_final)

