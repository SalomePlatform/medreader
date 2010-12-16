# This case corresponds to: /visu/Vectors/F9 case
# Create Vectors for field of the the given MED file for 10 timestamps%

import sys
from paravistest import datadir, pictureext, get_picture_dir
import paravis
from pvsimple import GetActiveSource, GetRenderView, Render
from presentations import VectorsOnField, hide_all,EntityType,PrsTypeEnum,reset_view,process_prs_for_test

# Create presentations
myParavis = paravis.myParavis

picturedir = get_picture_dir(sys.argv[1], "Vectors/F9")

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
# Get view
aView = GetRenderView()

# Create required presentations for the proxy
# CreatePrsForProxy(aProxy, aView, thePrsTypeList, thePictureDir, thePictureExt, theIsAutoDelete)
aFieldNames = aProxy.PointArrays.GetData()
aNbOnNodes = len(aFieldNames)
aFieldNames.extend(aProxy.CellArrays.GetData())
aTimeStamps = aProxy.TimestepValues.GetData()
aFieldEntity = EntityType.NODE
aFieldName = "MODES_DEPL"

#Creation of a set of non-colored and then colored Vectors presentations, based on time stamps of MODES_DEP field
for colored in [False,True]:
    colored_str = "_non-colored"
    if colored:
        colored_str = "_colored"
    for i in range(1,11):
        hide_all(aView, True)
        aPrs = VectorsOnField(aProxy, aFieldEntity,aFieldName , i,theIsColored=colored)	
        if aPrs is None:
            raise RuntimeError, "Presentation is None!!!"
        # display only current deformed shape
        #display_only(aView,aPrs)
        Prs.Visibility =1	
        reset_view(aView)
        Render(aView)
        # Add path separator to the end of picture path if necessery
        if not picturedir.endswith(os.sep):
                picturedir += os.sep
        prs_type = PrsTypeEnum.VECTORS
                
        # Get name of presentation type
        prs_name = PrsTypeEnum.get_name(prs_type)    
        f_prs_type = prs_name.replace(' ', '').upper()
        # Construct image file name
        pic_name = "{folder}{field}_{time}_{type}.{ext}".format(folder=picturedir,
                                                                                field=aFieldName+colored_str,
                                                                                time=str(i),
                                                                                type=f_prs_type,
                                                                                ext=pictureext)        
        # Show and record the presentation
        process_prs_for_test(aPrs, aView, pic_name)
 
