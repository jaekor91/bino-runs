import numpy as np
import matplotlib.pylab as plt
import astropy.io.fits as fits 
from xd_elg_utils import *


def hr2deg(hr):
    return hr * 15

def radec_angular_distance(ra1, dec1, ra2, dec2, indeg = True):
    """
    Assume that angles are provided in degrees to begin with.
    """
    if indeg:
        deg2rad = np.pi/180.        
        ra1, dec1, ra2, dec2 = ra1*deg2rad, dec1*deg2rad, ra2*deg2rad, dec2*deg2rad
    
    return np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1-ra2)) / deg2rad
    

tycho_directory = "/Users/jaehyeon/Documents/Research/ELG_target_selection/data-repository/"




print "Load ccdtosky output file."
# Load in hpix file
data = np.load("/Users/jaehyeon/Documents/Research/DESI-angular-clustering/ccdtosky/outputs/DR5/decals_Nside11/output_arr_chunk0thru50331648.npy")
ra, dec = data["hpix_ra"], data["hpix_dec"]
Ng, Nr, Nz = data["g_Nexp_sum"], data["r_Nexp_sum"], data["z_Nexp_sum"]
N_pix = ra.size# 
iNexp_cut = (Ng >=2) & (Nr >=2) & (Nz >=2)

# Apply further tycho constraint    
print "Generate tycho mask." # (iTycho True is good.)
iTycho = apply_tycho_radec(ra, dec, tycho_directory+"tycho2.fits", galtype="ELG") == 0
print "\n"



print "Spcify regions to target."
# 3 hr. Equatorial.
ifield_cut = (ra < hr2deg(3.5)) & (ra > hr2deg(2.)) & (dec < 5) & (dec > -5)

# For each field, randomly select a pixel until the following criteria area met.
# - hpix pixels that lie within 0.28 x 0.28 brick has a high pass rate.
# - Factor of two multiplication is to give the user more room to choose from
# - Second tolerance factor is used to prevent the fields being too close to each other.
tol = 0.28/2. * 3
tol2 = 1.5

# Controlling random seed
np.random.seed(52)

print "Find candidates for visual inspection."
coverage_threshold = 0.9
tycho_threashold = 0.9
num_candidates = 10 # Number of candidates to choose from.


# Concentrating to the region specified above, and that 
# pass the cuts criteria
ibool = ifield_cut & (iNexp_cut) & (iTycho)
ra_tmp = ra[ibool]
dec_tmp = dec[ibool]

print "Number of centers to try: %d" % ra_tmp.size

list_good = []
list_radec = []
list_coverage_fracs = []
list_junk = []
for j in range(num_candidates):
    print "Fishing for candidate %d" % j
    print "Candidate number: # of Hpix, Nexp>=2 frac, Tycho pass frac"
    
    while True:
        while True: # Pick a pixel that hasn't been considered yet and pass the two pass and tycho requirement.
            idx = np.random.randint(ra_tmp.size)            
            if (idx not in list_junk) & (idx not in list_good):
                # ra/dec of selected bricks
                ra_selected = ra_tmp[idx]
                dec_selected = dec_tmp[idx]                
                
                # If pass above conditions, check whether the selected field is close to any of the already selected field.
                itry = True
                for m in xrange(len(list_good)):
                    # Proper formular should be used here
                    dist = radec_angular_distance(list_radec[m][0], list_radec[m][1], ra_selected, dec_selected)                    
#                     print "dist: %.2f" % dist                                        
                    itry &= (dist > tol2)
                    
                if itry: # If still True, then continue to vet.
                    break
                    
                else:
                    list_junk.append(idx)
    
        # Collect all hpix that lie within the brick centered at the tractor brick center
        # /np.cos(dec_selected*np.pi/180.)
        ibool2 = (ra<(ra_selected+tol)) & (ra>(ra_selected-tol)) &  (dec<(dec_selected+tol)) & (dec>(dec_selected-tol)) 
        
        # Number of pixels inside the region
        Npix_inside = ibool2.sum()
        
        # Tells what fraction of pixels have Nexp_grz >=2
        coverage_fraction = iNexp_cut[ibool2].sum()/float(Npix_inside)

        # Tells what fraction of pixels pass Tycho 2 cut
        tycho_fraction = iTycho[ibool2].sum()/float(Npix_inside)

        print "%d: %d, %.3f, %.3f" % (idx, Npix_inside, coverage_fraction, tycho_fraction)

        if (coverage_fraction > coverage_threshold) and (tycho_fraction > tycho_threashold):
            list_good.append(idx)
            list_radec.append([ra_selected, dec_selected])
            list_coverage_fracs.append([coverage_fraction, tycho_fraction])
            print "Candidate Success"
            break
        else:
            list_junk.append(idx)
            print "Candidate Failed"
    print "\n"
    
    



list_radec= np.array(list_radec)
list_coverage_fracs = np.array(list_coverage_fracs)
ra_selected, dec_selected = list_radec[:, 0], list_radec[:, 1]    

print "Make scatter plots of regions selected"
fig, ax = plt.subplots(1, figsize=(10, 5))
ax.scatter(ra_selected, dec_selected, c="black", edgecolor="none", marker="o", s=50)
ax.axvline(x=hr2deg(3), c= "blue")
ax.axvline(x=hr2deg(2), c= "blue", ls="--")
ax.axvline(x=hr2deg(3.5), c= "blue", ls="--")


n = range(ra_selected.size)
for i, txt in enumerate(n):
    ax.annotate(txt, (ra_selected[i], dec_selected[i]), size=15)

# ax.axvline(x=hr2deg(8+25/60.), c="blue")
ax.axhline(y=5, c="red")
ax.axhline(y=-5, c="red")
ax.set_xlabel("RA")
ax.set_ylabel("DEC")
ax.axis("equal")
# plt.show()
plt.savefig("DESI-ELG-field-candidates-SGC.png", dpi=200, bbox_inches="tight")
plt.close()    


print "RA/DEC (Nexp2 fraction/Tyocho un-masked fracdtion)"
for m in xrange(num_candidates):
    print "%.2f, %.2f (%.2f, %.2f)" % (ra_selected[m], dec_selected[m], list_coverage_fracs[m][0], list_coverage_fracs[m][1])


np.save("DESI-ELG-fields-SGC-radec.npy", list_radec)
np.save("DESI-ELG-fields-SGC-coverage-fracs.npy", list_coverage_fracs)    