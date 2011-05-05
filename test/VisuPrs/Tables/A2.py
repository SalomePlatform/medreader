# This case corresponds to: /visu/tables/A2 case
# Import table in Post-Pro specific format from the file;
# create curves.

import sys

from paravistest import tablesdir, get_picture_dir, pictureext
from presentations import *
import paravis
import pvsimple

# Import table from file
print 'Import flux.tab.txt ....',
file_path = tablesdir + "flux.tab.txt"
table_reader = pvsimple.TableReader(FileName=file_path)
if table_reader is None:
    print "FAILED"
else:
    print "OK"

# Get available tables
print 'Get available tables .....'
available_tables = table_reader.GetPropertyValue("AvailableTables")
if (available_tables is None) or (len(available_tables) == 0):
    print "FAILED"
else:
    print available_tables

# Create curves
cur_view = pvsimple.GetRenderView()
if cur_view:
    pvsimple.Delete(cur_view)
xy_view = pvsimple.CreateXYPlotView()

xy_rep = pvsimple.Show(table_reader)
xy_rep.AttributeType = 'Row Data'
xy_rep.UseIndexForXAxis = 0
x_array = xy_rep.GetPropertyValue("SeriesNamesInfo")[0]
xy_rep.XArrayName = x_array
xy_rep.SeriesVisibility = [x_array, '0']

xy_rep.Visibility = 0
pvsimple.Render(xy_view)
xy_rep.Visibility = 1
pvsimple.Render(xy_view)

# Write image

# Directory for saving snapshots
picturedir = get_picture_dir(sys.argv[1], "Tables/A2")
if not picturedir.endswith(os.sep):
    picturedir += os.sep

file_name = picturedir + "flux_tab." + pictureext
pvsimple.WriteImage(file_name, view=xy_view, Magnification=1)


