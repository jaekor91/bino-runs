# Use this script to move tractor brick files from a designated directory to another.
# bricks-file is not used.

# Loading modules
import numpy as np
from os import listdir
from os.path import isfile, join
import sys
from astropy.io import fits


large_random_constant = -999119283571
deg2arcsec=3600


from_directory = "/global/project/projectdirs/cosmo/data/legacysurvey/dr5/tractor/"
to_directory = "/global/homes/j/jaehyeon/"


##############################################################################
# radec lists
flist = ["DESI-ELG-fields-SGC-0hr-radec.npy", "DESI-ELG-fields-SGC-3hr-radec.npy", "DESI-ELG-fields-NGC-radec.npy"]

ra_list = []
dec_list = []

for f in flist:
    a = np.load(f)
    ra_list.append(a[:,0])
    dec_list.append(a[:,1])
    
# ra = np.concatenate(ra_list)
# dec = np.concatenate(dec_list)


# Using survye-bricks file to get all bricks near the center bricks specified above
data = fits.open("./survey-bricks-dr5.fits")[1].data
ra_bricks = data["RA"]
dec_bricks = data["DEC"]
names_bricks = data["BRICKNAME"]

tol = 0.6
for i in range(3):
	bricks = []
	for j in range(ra_list[i].size):
		ra_c = ra_list[i][j]
		dec_c = dec_list[i][j]
		ibool = (ra_bricks < ra_c+tol) & (ra_bricks > ra_c - tol) & (dec_bricks > dec_c -tol) & (dec_bricks < dec_c+tol)
		bricks.append(list(names_bricks[ibool]))

	bricks = [item for sublist in bricks for item in sublist]
	# print bricks

	postfix = ".fits"
	prefix = "cp "


	f = open("DESI-ELG-tractor-move-binospec-test-%d.sh" % i, "w")
	for brick in bricks:
		if "p" in brick:
			tmp = brick.split("p")
		else:
			tmp = brick.split("m")		
		tractor_directory = tmp[0][-4:-1]
		brick_address = tractor_directory+"/"+"tractor-"+brick
		mv_command = "cp "+from_directory + brick_address + ".fits " + to_directory + "\n"
		f.write(mv_command)
	f.close()
	print("Completed")
