"""
This module provides auxiliary classes, functions and variables for testing.
"""

from __future__ import print_function
from math import fabs
import os
from datetime import date

import salome

class RepresentationType:
    """
    Types of representation.
    """
    OUTLINE = 0
    POINTS = 1
    WIREFRAME = 2
    SURFACE = 3
    SURFACEEDGES = 4
    VOLUME = 5
    POINTSPRITE = 6
    
    _type2name = {OUTLINE: 'Outline',
                  POINTS: 'Points',
                  WIREFRAME: 'Wireframe',
                  SURFACE:'Surface',
                  SURFACEEDGES:'Surface With Edges',
                  VOLUME:'Volume',
                  POINTSPRITE:'Point Sprite'}
                  
    @classmethod
    def get_name(cls, type):
        """Return paraview representation type by the primitive type."""
        return cls._type2name[type]

        
# Auxiliary classes
class SalomeSession(object):
    def __init__(self):
        import runSalome
        import sys
	sys.argv += ["--killall"]
        sys.argv += ["--show-desktop=1"]
	sys.argv += ["--splash=0"]
        sys.argv += ["--modules=VISU,PARAVIS"]
        clt, d = runSalome.main()
        port = d['port']
        self.port = port
        return
    def __del__(self):
        os.system('killSalomeWithPort.py %s'%(self.port))
        return
    pass

# Auxiliary functions
def test_values(value, et_value, check_error=0):
    """Test values."""
    error = 0
    length = len(value)
    et_length = len(et_value)
    if (length != et_length):
        err_msg = "ERROR!!! There is different number of created {0} and\
        etalon {1} values!!!".format(length, et_length)
        print(err_msg)
	error = error + 1
    else:
        for i in range(et_length):
	    if abs(et_value[i]) > 1:
	        max_val = abs(0.001 * et_value[i])
		if abs(et_value[i] -  value[i]) > max_val:
                    err_msg = "ERROR!!! Got value {0} is not equal\
                    to etalon value {1}!!!".format(value[i],  et_value[i])
		    print(err_msg)
		    error = error + 1
	    else:
	        max_val = 0.001
		if abs(et_value[i] - value[i])> max_val:
		    err_msg = "ERROR!!! Got value {0} is not equal\
                    to etalon value {1}!!!".format(value[i], et_value[i])
		    error = error + 1
    if check_error and error > 0:
        err_msg = "There is(are) some error(s) was(were) found... For more info see ERRORs above..."
        raise RuntimeError(err_msg)
    return error

def get_picture_dir(pic_dir, subdir):
    res_dir = pic_dir
    if not res_dir:
        res_dir = "/tmp/pic"

    # Add current date and subdirectory for the case to the directory path
    cur_date = date.today().strftime("%d%m%Y")
    res_dir += "/test_{date}/{subdir}".format(date=cur_date, subdir=subdir)
    # Create the directory if doesn't exist
    res_dir = os.path.normpath(res_dir)
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    else:
        # Clean the directory
        for root, dirs, files in os.walk(res_dir):
            for f in files:
                os.remove(os.path.join(root, f))

    return res_dir
   
#@VPH: utility function for 3D viewer test for common check of different types of presentation parameters set
def CallAndCheck(prs_name,property_name, value,do_raise = 1, compare_toler = -1.0):
    if isinstance(value,str):
        call_str =  prs_name+"."+property_name+"='" +value + "'" 
    else:
        call_str =  prs_name+"."+property_name+"=" +str(value) 
    exec("from pvsimple import *")    
    exec(prs_name+" = Show()")    
    try:
        exec(call_str)
    except ValueError:
        error_string  ="'"+str(value)+"' is not available for this type of presentations"
    else:
        error_string = None
        
    is_good = (error_string is None)
    if not is_good:
        if do_raise:
            raise RuntimeError, call_str +": "+error_string
        else:
            print(call_str+": "+error_string)
    else:
        # compare just set value and the one got from presentation
        command = "really_set_value = "+prs_name+"."+property_name
        exec(command)
        is_equal = 1
        if compare_toler > 0:
            is_equal = (fabs(really_set_value - value) < compare_toler)
        else:
            is_equal = (really_set_value == value)
        if not is_equal:
            if do_raise:
                raise RuntimeError, call_str + ": " + str(really_set_value) + " has been set instead"
            else:
                print (call_str + ": " + str(really_set_value) + " has been set instead")
                is_good = 0

    return is_good                       

# Auxiliary variables

# Data directory
datadir = os.getenv("DATA_DIR")
if datadir is not None:
    datadir =  os.path.normpath(datadir)
    datadir = datadir+ "/MedFiles/"
    
# Graphica files extension
pictureext = os.getenv("PIC_EXT");
if pictureext == None:
    pictureext="png"

# Run Salome
salome_session = SalomeSession()
salome.salome_init()

# Create new study
print ("Creating new study...", end="")
aStudy = salome.myStudyManager.NewStudy("Study1")
if aStudy is None: 
    raise RuntimeError("Error")
else: 
    print("OK")

salome.myStudy = aStudy
