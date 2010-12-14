# This case corresponds to: /visu/ScalarMap/E9 case
# Create Scalar Map for all fields of the the given MED file

import sys

from paravistest import datadir, pictureext, get_picture_dir
from presentations import CreatePrsForFile, PrsTypeEnum
import paravis


# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "ScalarMap/F4")

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "UO2_250ans_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.SCALARMAP], picturedir, pictureext)


