import numpy as np

print "Please select targets based on a *square* region centered at listed ra/dec position and having a length 0.28 deg * 3.\n"
# Printing selected survey fields

print "SGC - 0hr"
a = np.load("DESI-ELG-fields-SGC-0hr-radec.npy")
ra_selected, dec_selected = a[:,0], a[:,1]
a = np.load("DESI-ELG-fields-SGC-0hr-coverage-fracs.npy")
frac1, frac2 = a[:,0], a[:,1]
print "Field Number: RA, DEC (Frac w/ Nexp_grz >= 2, Frac w/ Tycho unamsked)"
for i in range(ra_selected.size):
    print "%d: %.4f, %.4f (%.4f, %.4f)" % (i, ra_selected[i], dec_selected[i], frac1[i], frac2[i])
print "\n"

print "SGC - 3hr"
a = np.load("DESI-ELG-fields-SGC-3hr-radec.npy")
ra_selected, dec_selected = a[:,0], a[:,1]
a = np.load("DESI-ELG-fields-SGC-3hr-coverage-fracs.npy")
frac1, frac2 = a[:,0], a[:,1]
print "Field Number: RA, DEC (Frac w/ Nexp_grz >= 2, Frac w/ Tycho unamsked)"
for i in range(ra_selected.size):
    print "%d: %.4f, %.4f (%.4f, %.4f)" % (i, ra_selected[i], dec_selected[i], frac1[i], frac2[i])
print "\n"


print "NGC - 8h+30"
a = np.load("DESI-ELG-fields-NGC-radec.npy")
ra_selected, dec_selected = a[:,0], a[:,1]
a = np.load("DESI-ELG-fields-NGC-coverage-fracs.npy")
frac1, frac2 = a[:,0], a[:,1]
print "Field Number: RA, DEC (Frac w/ Nexp_grz >= 2, Frac w/ Tycho unamsked)"
for i in range(ra_selected.size):
    print "%d: %.4f, %.4f (%.4f, %.4f)" % (i, ra_selected[i], dec_selected[i], frac1[i], frac2[i])
print "\n"


print "COSMOS"

