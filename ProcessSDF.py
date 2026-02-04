from os import chdir, listdir
from chimera import runCommand
chdir("/home/niper02/Downloads/vimla/uyt") # change to the SDF file directory
for sdf in listdir("."):
    if not sdf.endswith(".sdf"):
    	continue
    runCommand("open " + sdf)
    runCommand("minimize")
    runCommand("write format mol2 0 " + sdf[:-4] + ".mol2")
