# This case corresponds to: /visu/ScalarMap_On_DeformedShape/F0 case
# Create Scalar Map on Deformed Shape for all fields of the the given MED file

import sys

from paravistest import datadir, pictureext, get_picture_dir
from presentations import CreatePrsForFile, PrsTypeEnum
import paravis


# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "ScalarMap_On_DeformedShape/F0")

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "gro5couches_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.DEFORMEDSHAPESCALARMAP], picturedir, pictureext)


