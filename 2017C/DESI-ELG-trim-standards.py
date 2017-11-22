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
            ra_tmp, dec_tmp = star_tmp[:, 0], star_tmp[:, 1]
            idx1, idx2 = crossmatch_cat1_to_cat2(ra_tmp, dec_tmp, ra, dec, tol=0.5/(deg2arcsec+1e-12))    
            if idx1.size > 0:
                ra_med_diff = med_x1_minus_x2(ra_tmp[idx1], ra[idx2])
                dec_med_diff = med_x1_minus_x2(dec_tmp[idx1], dec[idx2])
                star_tmp[:, 0] -= ra_med_diff
                star_tmp[:, 1] -= dec_med_diff
                star_final.append(star_tmp)
                print "Stars"
                print "Number of matches: %d" % idx1.size 
                print "Corrections in ra, dec: %.5f, %.5f" % (ra_med_diff*3600, dec_med_diff*3600)
                print"\n"
            else: 
                print "No star matched in this field."
                print "No stars are included from here."
        
        # Gal
        if Ngal>0:
            gal_tmp = gal_array[iselect_gal]                
            ra_tmp, dec_tmp = gal_tmp[:, 0], gal_tmp[:, 1]
            idx1, idx2 = crossmatch_cat1_to_cat2(ra_tmp, dec_tmp, ra, dec, tol=0.5/(deg2arcsec+1e-12))                
            if idx1.size > 0:
                ra_med_diff = med_x1_minus_x2(ra_tmp[idx1], ra[idx2])
                dec_med_diff = med_x1_minus_x2(dec_tmp[idx1], dec[idx2])
                gal_tmp[:, 0] -= ra_med_diff
                gal_tmp[:, 1] -= dec_med_diff
                gal_final.append(gal_tmp)
                print "Gals"
                print "Number of matches: %d" % idx1.size 
                print "Corrections in ra, dec: %.5f, %.5f" % (ra_med_diff*3600, dec_med_diff*3600)
                print"\n"
            else: 
                print "No gal matched in this field."
                print "No gal are included from here."

        print"\n"            

# Concatenate all targets
gal_final = np.vstack(gal_final)
np.save("DESI-ELG-boss.npy", gal_final)
star_final = np.vstack(star_final)
np.save("DESI-ELG-standards.npy", star_final)

