import json
from pprint import pprint
import numpy as np

def produce_tally(fname):
    data = json.load(open("./output-catalogs/"+fname))

    name = []
    for i in range(len(data["sidea"])):
        name.append(data["sidea"][i]["OBJECT"])

    for i in range(len(data["sideb"])):
        name.append(data["sideb"][i]["OBJECT"])

    bitmask = []
    nstar = 0
    ngal = 0
    nFDR = 0
    nNDM1 = 0
    nNDM2 = 0
    nNDM3 = 0
    nNDM = 0
    nRF1 = 0
    nRF2 = 0
    nRF3 = 0
    nRF = 0
    nOR = 0
    nAND = 0
    ntot = 0
    for n in name:
        ntot +=1
        if n == "stars":
            nstar +=1
        elif n =="gal":
            ngal += 1
        else:
            bit = int(n)
            if (np.bitwise_and(bit, 2**3) > 0):
                nFDR +=1
            if (np.bitwise_and(bit, 2**4) > 0):
                nNDM1 +=1
            if (np.bitwise_and(bit, 2**6) > 0):
                nNDM2 +=1            
            if (np.bitwise_and(bit, 2**7) > 0):
                nNDM3 +=1                        
            if (np.bitwise_and(bit, 2**4 + 2**6 + 2**7) > 0):
                nNDM +=1
            if (np.bitwise_and(bit, 2**8) > 0):
                nRF1 +=1                                    
            if (np.bitwise_and(bit, 2**9) > 0):
                nRF2 +=1                                                
            if (np.bitwise_and(bit, 2**10) > 0):
                nRF3 +=1                                                
            if (np.bitwise_and(bit, 2**10 + 2**9 + 2**8) > 0):
                nRF +=1
            if ((np.bitwise_and(bit, 2**10 + 2**9 + 2**8) > 0) or (np.bitwise_and(bit, 2**4 + 2**6 + 2**7) > 0)):
                nOR +=1
            if ((np.bitwise_and(bit, 2**10 + 2**9 + 2**8) > 0) and (np.bitwise_and(bit, 2**4 + 2**6 + 2**7) > 0)):
                nAND +=1
    return nstar, ngal, nFDR, nNDM1, nNDM2, nNDM3, nNDM, nRF1, nRF2, nRF3, nRF, nOR, nAND, ntot

print "/---- COSMOS"
print "OR: Union of NDM and RF"
print "AND: Intersection of NDM and RF"
print ("%4s\t" * 14)  % ("star", "gal", "FDR", "NDM1", "NDM2", "NDM3", "NDM", "RF1", "RF2", "RF3", "RF", "OR", "AND", "Tot")
for fname in ["COSMOS-1.json", "COSMOS-2.json"]:
    print (("%4d\t")*14)  % produce_tally(fname)
print "\n"

print "/---- DR6"
print "OR: Union of NDM and RF"
print "AND: Intersection of NDM and RF"
print ("%4s\t" * 14)  % ("star", "gal", "FDR", "NDM1", "NDM2", "NDM3", "NDM", "RF1", "RF2", "RF3", "RF", "OR", "AND", "Tot")
for fname in ["DR6-%d.json" % i for i in range(1, 7)]:
    print (("%4d\t")*14)  % produce_tally(fname)    
print "\n"

print "/---- HSC"
print "OR: Union of NDM and RF"
print "AND: Intersection of NDM and RF"
print ("%4s\t" * 14)  % ("star", "gal", "FDR", "NDM1", "NDM2", "NDM3", "NDM", "RF1", "RF2", "RF3", "RF", "OR", "AND", "Tot")
for fname in ["HSC-%d.json" % i for i in range(1, 7)]:
    print (("%4d\t")*14)  % produce_tally(fname)        