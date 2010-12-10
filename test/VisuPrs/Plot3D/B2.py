#This case corresponds to: /visu/Plot3D/B2 case
#%Create Plot 3D for all fields of the the given MED file%

from paravistest import * 
from presentations import *
import paravis

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "cube_hexa8_quad4_236.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.PLOT3D], picturedir, pictureext)
