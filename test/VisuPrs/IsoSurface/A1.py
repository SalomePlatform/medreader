# This case corresponds to: /visu/IsoSurface/A1 case
# Create Iso Surface for all fields of the the given MED file

import sys

from paravistest import datadir, pictureext, get_picture_dir
from presentations import CreatePrsForFile, PrsTypeEnum
import paravis


# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "IsoSurfaces/A1")

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "hexa_28320_ELEM_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.ISOSURFACES], picturedir, pictureext)
