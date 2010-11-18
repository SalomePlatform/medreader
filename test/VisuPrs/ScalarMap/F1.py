#This case corresponds to: /visu/ScalarMap/F1 case
#%Create Scalar Map for all fields of the the given MED file%

from paravistest import *
from presentations import *
import paravis

# Create presentations
myParavis = paravis.myParavis

file = datadir +  "ml.med"
print " --------------------------------- "
print "file ", file
print " --------------------------------- "
print "\nCreatePrsForFile..."
CreatePrsForFile(myParavis, file, [PrsTypeEnum.SCALARMAP], picturedir, pictureext)


