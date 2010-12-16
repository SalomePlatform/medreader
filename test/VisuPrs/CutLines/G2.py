# This case corresponds to: /visu/CutLines/G2 case
# Create Cut lines for field of the the given MED file for 10 timestamps%

import sys
from paravistest import datadir, pictureext, get_picture_dir
from presentations import *
from pvsimple import *
import paravis

# Create presentations
myParavis = paravis.myParavis
picturedir = get_picture_dir(sys.argv[1], "CutLines/G2")

theFileName = datadir +  "Bug829_resu_mode_236.med"
print " --------------------------------- "
print "file ", theFileName
print " --------------------------------- "

"""Build presentations of the given types for all fields of the given file."""
#print "Import %s..." % theFileName.split('/')[-1],
result = myParavis.ImportFile(theFileName)
aProxy = GetActiveSource()
if aProxy is None:
	raise RuntimeError, "Error: can't import file."
else: print "OK"
# Get viewScalarMap
aView = GetRenderView()

# Create required presentations for the proxy
# CreatePrsForProxy(aProxy, aView, thePrsTypeList, thePictureDir, thePictureExt, theIsAutoDelete)
aFieldNames = aProxy.PointArrays.GetData()
aNbOnNodes = len(aFieldNames)
aFieldNames.extend(aProxy.CellArrays.GetData())
aTimeStamps = aProxy.TimestepValues.GetData()
aFieldEntity = EntityType.NODE
aFieldName = "MODES_DEPL"

#create Cut Lines presentations for 10 timestamps
for i in range(1,11):
    hide_all(aView, True)
    aPrs = CutLinesOnField(aProxy, aFieldEntity,aFieldName , i)
    if aPrs is None:
        raise RuntimeError, "Presentation is None!!!"
    #display only current scalar map
    aPrs.Visibility=1
    reset_view(aView)
    Render(aView)    
    
    # Add path separator to the end of picture path if necessery
    if not picturedir.endswith(os.sep):
            picturedir += os.sep
    prs_type = PrsTypeEnum.CUTLINES
            
    # Get name of presentation type
    prs_name = PrsTypeEnum.get_name(prs_type)    
    f_prs_type = prs_name.replace(' ', '').upper()
    # Construct image file name
    pic_name = "{folder}{field}_{time}_{type}.{ext}".format(folder=picturedir,
                                                                            field=aFieldName,
                                                                            time=str(i),
                                                                            type=f_prs_type,
                                                                            ext=pictureext)
    
    # Show and record the presentation
    process_prs_for_test(aPrs, aView, pic_name)


