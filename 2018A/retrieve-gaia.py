# scp /Users/jaehyeon/Documents/Research/binospec/bino-runs/2018A/retrieve-gaia.py jaelee@odyssey.rc.fas.harvard.edu:
# scp jaelee@odyssey.rc.fas.harvard.edu:gaia.fits /Users/jaehyeon/Documents/Research/binospec/bino-runs/2018A/data/

import astropy.io.fits as fits
from os import listdir
from os.path import isfile, join
import numpy as np
from functools import reduce

#---- magnitude limits
mlim = [15, 20]
ralim = [187, 191]
declim = [61, 63]

#---- Get all file names
mypath = '/n/fink2/www/dfink/gaia/chunks-gaia_rel1/'
flist = [f for f in listdir(mypath) if isfile(join(mypath, f)) and "chunk" in f]
objs = []

Nf = len(flist)# Number of files

#---- Go through each files and extract objects we want.
nw = 0
for i, fname in enumerate(flist): 
  print "File %d out of %d: %s" % (i, Nf, fname)
  hdr = fits.open(mypath+fname)
  a = hdr[1].data
  ra, dec = a["ra"], a["dec"]  
  Gmag = a["phot_g_mean_mag"]
  w = np.where(reduce(np.logical_and,(Gmag>mlim[0], Gmag<mlim[1], a["phot_g_n_obs"]>10, \
    ra<ralim[1], ra>ralim[0], dec<declim[1], dec>declim[0])))[0]
  nw += w.size 
  print('%d stars found so far' % nw)
  if nw>0:
    objs.append(a[w])
    # Stack all of the fount objects
    objs_all = np.hstack(objs)
    # Save the new array.
    cols = fits.ColDefs(objs_all)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto("/n/home06/jaelee/gaia.fits", clobber=True)    
  hdr.close()



print("Completed.")
