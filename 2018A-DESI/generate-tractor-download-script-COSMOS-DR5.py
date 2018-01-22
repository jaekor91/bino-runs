from utils import *

# Loading modules
dr_v = "5" # Data release version

bricks = load_fits_table("./data/survey-bricks-dr5.fits")
ra_b, dec_b = bricks["ra"], bricks["dec"]
bricks_COSMOS = bricks["brickname"][(ra_b>149.1) & (ra_b<150.7) & (dec_b>1.3) & (dec_b<2.7)]
print("Generating download scripts for COSMOS")
portal_address = "http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr"+dr_v+"/tractor/"
postfix = ".fits\n"
prefix = "wget "

f = open("tractor-download-COSMOS-DR5.sh","w")
for brick in bricks_COSMOS:
    tractor_directory = brick[:3]
    brick_address = tractor_directory+"/tractor-"+brick+postfix
    download_command = prefix + portal_address + brick_address
    f.write(download_command)
f.close()
print("Completed")