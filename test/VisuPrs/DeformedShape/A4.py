# This case corresponds to: /visu/DeformedShape/A4 case
# Create Deformed Shape for all fields of the the given MED file

import sys

from paravistest import datadir, pictureext, get_picture_dir
from presentations import CreatePrsForFile, PrsTypeEnum
import paravis


# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "DeformedShape/A4")

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "Hexa8_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.DEFORMEDSHAPE], picturedir, pictureext)

