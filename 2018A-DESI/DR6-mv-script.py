from os import listdir
from os.path import isfile, join 

top_dir = "/global/project/projectdirs/cosmo/work/legacysurvey/dr6/tractor/"

print "#---- DR6"
ra_cs = [115.75, 115.75, 117, 118, 119, 119]
dec_cs = [37, 35, 35, 35.5, 35, 33]

# For each center ra/dec
for i in range(6):
    ra_c, dec_c = ra_cs[i], dec_cs[i]
    print "Field", i+1, ra_c, dec_c

    mv_list = []    
    for tractor_dir in ["115/", "116/", "117/", "118/", "119/"]:
        fdir = top_dir + tractor_dir
    
        onlyfiles = [f for f in listdir(fdir) if isfile(join(fdir, f))]
    
        for f in onlyfiles:
            if f.split(".")[-1] =="fits":
                tmp = f[8:-5].split("p")
                ra_tmp = int(tmp[0])/10.
                dec_tmp = int(tmp[1])/10.
                if (ra_tmp > ra_c - 0.5) & (ra_tmp < ra_c + 0.5) & \
                (dec_tmp > dec_c - 0.5) & (dec_tmp < dec_c + 0.5):
                    mv_list.append(fdir+f)
    print len(mv_list)
    new_f = open("DR6-%d-move.sh" % (i+1), "w")
    for f in mv_list:
        new_f.write("cp "+f+" /global/homes/j/jaehyeon/DR6-%d" % (i+1)+"\n")
    new_f.close()
