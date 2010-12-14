# This case corresponds to: /visu/CutPlanes/E0 case
# Create Cut Planes for all fields of the the given MED file

import sys

from paravistest import datadir, pictureext, get_picture_dir
from presentations import CreatePrsForFile, PrsTypeEnum
import paravis


# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "CutPlanes/E0")

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "test_55_solid_concentr_dom_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.CUTPLANES], picturedir, pictureext)


