"""
This module is intended to provide Python API for building presentations
typical for Post-Pro module (Scalar Map, Deformed Shape, Vectors, etc.)
"""


from __future__ import division
from __future__ import print_function

import os
import re
from math import sqrt
from string import upper

import pvsimple as pv

# Constants
EPS = 1E-3
FLT_MIN = 1E-37
VTK_LARGE_FLOAT = 1E+38
GAP_COEFFICIENT = 0.0001


# Globals
_current_bar = None

# Enumerations
class PrsTypeEnum:
    """
    Post-Pro presentation types.
    """
    MESH = 0
    SCALARMAP = 1
    ISOSURFACES = 2
    CUTPLANES = 3
    CUTLINES = 4
    DEFORMEDSHAPE = 5
    DEFORMEDSHAPESCALARMAP = 6
    VECTORS = 7
    PLOT3D = 8
    STREAMLINES = 9
    GAUSSPOINTS = 10

    _type2name = {MESH: 'Mesh',
                  SCALARMAP: 'Scalar Map',
                  ISOSURFACES: 'Iso Surfaces',
                  CUTPLANES: 'Cut Planes',
                  CUTLINES: 'Cut Lines',
                  DEFORMEDSHAPE: 'Deformed Shape',
                  DEFORMEDSHAPESCALARMAP: 'Deformed Shape And Scalar Map',
                  VECTORS: 'Vectors',
                  PLOT3D: 'Plot3D',
                  STREAMLINES: 'Stream Lines',
                  GAUSSPOINTS: 'Gauss Points'}
    
    @classmethod
    def get_name(cls, type):
        """Return presentaion name by its type."""
        return cls._type2name[type]


class EntityType:
    """
    Entity types.
    """
    NODE = 0
    CELL = 1

    _type2name = {NODE: 'OnPoint',
                  CELL: 'OnCell'}
    
    _name2type = {'OnPoint': NODE,
                  'OnCell': CELL}

    _type2pvtype = {NODE: 'POINT_DATA',
                    CELL: 'CELL_DATA'}

    @classmethod
    def get_name(cls, type):
        """Return entity name (used in full group names) by its type."""
        return cls._type2name[type]
    
    @classmethod
    def get_type(cls, name):
        """Return entity type by its name (used in full group names)."""
        return cls._name2type[name]
    
    @classmethod
    def get_pvtype(cls, type):
        """Return entity type from ['CELL_DATA', 'POINT_DATA']"""
        return cls._type2pvtype[type]


class Orientation:
    """
    Orientation types.
    """
    AUTO = 0
    XY = 1
    YZ = 2
    ZX = 3


class GlyphPos:
    """
    Glyph positions.
    """
    CENTER = 0
    TAIL = 1
    HEAD = 2


class GaussType:
    """
    Types of Gauss Points representation.
    """
    SPRITE = 0 	
    POINT = 1 	
    SPHERE = 2

    _type2mode = {SPRITE: 'Texture',
                  POINT: 'SimplePoint',
                  SPHERE: 'Sphere'}

    @classmethod
    def get_mode(cls, type):
        """Return paraview point sprite mode by the primitive type."""
        return cls._type2mode[type]


# Internal classes
class _AuxCounter:
    "Internal class."
    first_render = True


# Auxiliary functions
def process_prs_for_test(prs, view, picture_name, show_bar=True):
    """Show presentation and record image.
    
    Keyword arguments:
    prs -- the presentation to show
    view -- the render view
    picture_name -- the full name of a graphics file to save
    show_bar -- to show scalar bar or not
    
    """
    # Show the presentation only
    display_only(prs, view)

    # Show scalar bar
    if show_bar and _current_bar:
        _current_bar.Visibility = 1
        
    # Reset the view 
    reset_view(view) 
           
    # Create a directory for screenshot if necessary
    file_name = re.sub("\s+","_", picture_name)
    pic_dir = os.path.dirname(picture_name)
    if not os.path.exists(pic_dir):
        os.makedirs(pic_dir)
        
    # Save picture
    pv.WriteImage(file_name, view=view, Magnification=1)


def reset_view(view=None):
    """Reset the view.
    
    If the view is not passed, the active view is used.
    """
    if not view:
        view = pv.GetRenderView()

    if _AuxCounter.first_render:
        # Camera preferences
        view.CameraFocalPoint = [0.0, 0.0, 0.0]
        view.CameraViewUp = [0.0, 0.0, 1.0]
        view.CameraPosition = [738.946, -738.946, 738.946]

        # Turn on the headligth
        view.LightSwitch = 1
        view.LightIntensity = 0.5

        # Use parallel projection
        view.CameraParallelProjection = 1
        
        _AuxCounter.first_render = False

    view.ResetCamera()
    pv.Render(view=view)


def hide_all(view, to_remove=False):
    """Hide all representations in the view."""
    if not view:
        view = pv.GetRenderView()
        
    rep_list = view.Representations
    for rep in rep_list:
        if hasattr(rep, 'Visibility') and rep.Visibility != 0:
            rep.Visibility = 0
        if to_remove:
            view.Representations.remove(rep)
    pv.Render(view=view)


def display_only(prs, view=None):
    """Display only the presentation in the view."""
    hide_all(view)
    if (hasattr(prs, 'Visibility') and prs.Visibility != 1):
        prs.Visibility = 1

    
def _check_vector_mode(vector_mode, nb_components):
    """Check vector mode."""
    if vector_mode not in ('Magnitude', 'X', 'Y', 'Z'):
        raise ValueError("Unexistent vector mode: {0}".format(vector_mode))
    
    if ((nb_components == 1 and (vector_mode == 'Y' or vector_mode == 'Z')) or
        (nb_components == 2 and  vector_mode == 'Z')):
        raise ValueError("Incorrect vector mode {0} for {1}-component field".
                         format(vector_mode, nb_components))


def _get_vector_component(vector_mode):
    vcomponent = -1

    if vector_mode == 'X':
        vcomponent = 0
    elif vector_mode == 'Y':
        vcomponent = 1
    elif vector_mode == 'Z':
        vcomponent = 2

    return vcomponent

def get_data_range(proxy, entity_type, field_name, vector_mode='Magnitude'):
    """Get data range for the field."""
    entity_data_info = None
    if entity_type == EntityType.CELL:
        entity_data_info = proxy.GetCellDataInformation()
    elif entity_type == EntityType.NODE:
        entity_data_info = proxy.GetPointDataInformation()

    data_range = []

    vcomp = _get_vector_component(vector_mode)
    data_range = entity_data_info[field_name].GetComponentRange(vcomp)

    return data_range


def get_bounds(proxy):
    """Get bounds of the proxy."""
    dataInfo = proxy.GetDataInformation()
    bounds_info = dataInfo.GetBounds()
    return bounds_info


def get_x_range(proxy):
    """Get X range."""
    bounds_info = get_bounds(proxy)
    return bounds_info[0:2]


def get_y_range(proxy):
    """Get Y range."""
    bounds_info = get_bounds(proxy)
    return bounds_info[2:4]


def get_z_range(proxy):
    """Get Z range."""
    bounds_info = get_bounds(proxy)
    return bounds_info[4:6]


def is_planar_input(proxy):
    """Check if the given input is planar."""
    bounds_info = get_bounds(proxy)

    if ( abs(bounds_info[0] - bounds_info[1]) <= FLT_MIN or
         abs(bounds_info[2] - bounds_info[3]) <= FLT_MIN or
         abs(bounds_info[4] - bounds_info[5]) <= FLT_MIN ):
        return True

    return False


def is_data_on_cells(proxy, field_name):
    """Check if the field on cells with the given name exists."""
    cell_data_info = proxy.GetCellDataInformation()
    return (field_name in cell_data_info.keys())


def GetOrientation(theProxy):
    """Get the optimum cutting plane orientation for Plot3d."""
    orientation = Orientation.XY

    bounds = get_bounds(theProxy)
    delta = [bounds[1] - bounds[0],
              bounds[3] - bounds[2],
              bounds[5] - bounds[4] ]

    if (delta[0] >= delta[1] and delta[0] >= delta[2]):
        if (delta[1] >= delta[2]):
            orientation = Orientation.XY
        else:
            orientation = Orientation.ZX
    elif (delta[1] >= delta[0] and delta[1] >= delta[2]):
        if (delta[0] >= delta[2]):
            orientation = Orientation.XY
        else:
            orientation = Orientation.YZ
    elif (delta[2] >= delta[0] and delta[2] >= delta[1]):
        if (delta[0] >= delta[1]):
            orientation = Orientation.ZX
        else:
            orientation = Orientation.YZ

    return orientation


def GetNormalByOrientation(theOrientation):
    """Get normal for the plane by its orientation."""
    normal = [0.0, 0.0, 0.0]
    if (theOrientation == Orientation.XY):
        normal[2] = 1.0
    elif (theOrientation == Orientation.ZX):
        normal[1] = 1.0
    elif (theOrientation == Orientation.YZ):
        normal[0] = 1.0
        
    return normal


def GetRangeForOrientation(theProxy, theOrientation):
    """Get source range for cutting plane orientation."""
    val_range = []
    if (theOrientation == Orientation.XY):
        val_range = get_z_range(theProxy)
    elif (theOrientation == Orientation.ZX):
        val_range = get_y_range(theProxy)
    elif (theOrientation == Orientation.YZ):
        val_range = get_x_range(theProxy)

    return val_range


def GetPositions(val_range, theNbPlanes, theDisplacement):
    """Compute plane positions."""
    positions = []
    if (theNbPlanes > 1):
        aDistance = val_range[1] - val_range[0]
        val_range[1] = val_range[0] + (1.0 - EPS)*aDistance
        val_range[0] = val_range[0] + EPS*aDistance
        
        aDistance = val_range[1] - val_range[0]
        step = aDistance/(theNbPlanes-1)
        aDisplacement = step*theDisplacement
        startPos = val_range[0] - 0.5*step + aDisplacement
        
        for i in xrange(theNbPlanes):
            pos = startPos + i*step
            positions.append(pos)
    elif (theNbPlanes == 1):
        pos = val_range[0] + (val_range[1] - val_range[0])*theDisplacement
        positions.append(pos)
        
    return positions

    
def get_nb_components(proxy, entity_type, field_name):
    """Return number of components for the field."""
    entity_data_info = None
    if entity_type == EntityType.CELL:
        entity_data_info = proxy.GetCellDataInformation()
    elif entity_type == EntityType.NODE:
        entity_data_info = proxy.GetPointDataInformation()
        
    nb_components = entity_data_info[field_name].GetNumberOfComponents()
    return nb_components


def get_scale_factor(proxy):
    """Compute scale factor."""
    if not proxy:
        return 0.0

    proxy.UpdatePipeline()
    data_info = proxy.GetDataInformation()

    nb_cells = data_info.GetNumberOfCells()
    nb_points = data_info.GetNumberOfPoints()
    nb_elements = nb_cells if nb_cells > 0  else nb_points
    bounds = get_bounds(proxy)

    volume = 1
    vol = dim = 0
    
    for i in xrange(0, 6, 2):
        vol = abs(bounds[i+1] - bounds[i])
        if (vol > 0):
            dim += 1
            volume *= vol

    if nb_elements == 0 or dim < 1 / VTK_LARGE_FLOAT:
        return 0

    volume /= nb_elements
    
    return pow(volume, 1 / dim)

            
def get_default_scale(prs_type, proxy, entity_type, field_name):
    """Get default scale factor."""
    data_range = get_data_range(proxy, entity_type, field_name)
    
    if prs_type == PrsTypeEnum.DEFORMEDSHAPE:
        EPS = 1.0 / VTK_LARGE_FLOAT
        if ( abs(data_range[1]) > EPS ):
            scale_factor = get_scale_factor(proxy)
            return scale_factor / data_range[1]
    elif prs_type == PrsTypeEnum.PLOT3D:
        bounds = get_bounds(proxy)
        length = sqrt( (bounds[1]-bounds[0])**2 +
                        (bounds[3]-bounds[2])**2 +
                        (bounds[5]-bounds[4])**2 )
        
        EPS = 0.3
        if data_range[1] > 0:
            return length / data_range[1] * EPS
        
    return 0


def get_calc_magnitude(proxy, array_entity, array_name):
    """Compute magnitude for the given vector array via Calculator.
    
    Return the calculator object.
    """
    calculator = None
    
    # Transform vector array to scalar array if possible
    nb_components = get_nb_components(proxy, array_entity, array_name)
    if (nb_components > 1):
        calculator = pv.Calculator(proxy)
        attribute_mode = "point_data"
        if array_entity != EntityType.NODE:
            attribute_mode = "cell_data"
        calculator.AttributeMode = attribute_mode
        if (nb_components == 2):
            # Workaroud: calculator unable to compute magnitude
            # if number of components equal to 2
            calculator.Function = "sqrt({0}_X^2+{0}_Y^2)".format(array_name)
        else:
            calculator.Function = "mag({0})".format(array_name)
        calculator.ResultArrayName = array_name + "_magnitude"
        calculator.UpdatePipeline()
            
    return calculator


def get_add_component_calc(proxy, array_entity, array_name):
    """Creates 3-component array from 2-component.

    The first two components is from the original array. The 3rd component
    is zero.
    If the number of components is not equal to 2 - return original array name.
    
    Return the calculator object.
    """
    calculator = None
    
    nb_components = get_nb_components(proxy, array_entity, array_name)
    if nb_components == 2:
        calculator = pv.Calculator(proxy)
        attribute_mode = "point_data"
        if array_entity != EntityType.NODE:
            attribute_mode = "cell_data"
        calculator.AttributeMode = attribute_mode
        expression = "iHat * {0}_X + jHat * {0}_Y + kHat * 0"
        calculator.Function = expression.format(array_name)
        calculator.ResultArrayName = array_name + "_3c"
        calculator.UpdatePipeline()
            
    return calculator

def select_all_cells(proxy):
    """Select all cell types.
    
    Used in creation of mesh/submesh presentation.
    """
    all_cell_types = proxy.CellTypes.Available
    proxy.CellTypes = all_cell_types
    proxy.UpdatePipeline()


def select_cells_with_data(proxy):
    """Select only cell types which have the data.
    
    Only cell types with data for even one field
    will be selected.
    """
    all_cell_types = proxy.CellTypes.Available
    all_arrays = list(proxy.CellArrays.GetData())
    all_arrays.extend(proxy.PointArrays.GetData())

    if not all_arrays:
        file_name = proxy.FileName.split(os.sep)[-1]
        print("Warning: {0} doesn't contain any data array.".format(file_name))

    cell_types_on = []

    for cell_type in all_cell_types:
        proxy.CellTypes = [cell_type]
        proxy.UpdatePipeline()

        cell_arrays = proxy.GetCellDataInformation().keys()
        point_arrays = proxy.GetPointDataInformation().keys()
        
        in_arrays = lambda array: (array in cell_arrays) or (array in point_arrays)
        if any(in_arrays(array) for array in all_arrays):
            cell_types_on.append(cell_type)

    proxy.CellTypes = cell_types_on
    proxy.UpdatePipeline()


def extract_groups_for_field(proxy, field_name, field_entity):
    """Exctract only groups which have the field.
    
    Return ExtractGroup object, if not all groups have the field.
    Return the initial proxy if no groups had been filtered.
    """
    source = proxy
    
    # Remember the state
    initial_groups = list(proxy.Groups)

    # Get data information for the field entity
    entity_data_info = None
    if field_entity == EntityType.CELL:
        entity_data_info = proxy.GetCellDataInformation()
    elif field_entity == EntityType.NODE:
        entity_data_info = proxy.GetPointDataInformation()

    # Collect groups for extraction
    groups_to_extract = []
    
    for group in initial_groups:
        proxy.Groups = [group]
        proxy.UpdatePipeline()
        if field_name in entity_data_info.keys():
            groups_to_extract.append(group)

    # Restore state
    proxy.Groups = initial_groups
    proxy.UpdatePipeline()

    # Extract groups if necessary
    if len(groups_to_extract) < len(initial_groups):
        extract_group = pv.ExtractGroup(proxy)
        extract_group.Groups = groups_to_extract
        extract_group.UpdatePipeline()
        source = extract_group

    return source

def if_possible(proxy, field_name, entity_type, prs_type):
    """Check if the presentation is possible on the given field."""
    result = True
    if (prs_type == PrsTypeEnum.DEFORMEDSHAPE or
        prs_type == PrsTypeEnum.DEFORMEDSHAPESCALARMAP or
        prs_type == PrsTypeEnum.VECTORS or
        prs_type == PrsTypeEnum.STREAMLINES):
        nb_comp = get_nb_components(proxy, entity_type, field_name)
        result = (nb_comp > 1)
    elif (prs_type == PrsTypeEnum.GAUSSPOINTS):
        result = (entity_type == EntityType.CELL or 
                  field_name in proxy.QuadraturePointArrays.Available)
        #@MZN: quadrature points arrays
    elif (prs_type == PrsTypeEnum.MESH):
        result = len(_get_group_names(proxy, field_name, entity_type)) > 0
                
    return result


def add_scalar_bar(field_name, nb_components,
                   vector_mode, lookup_table, time_value):
    """Add scalar bar."""
    global _current_bar
    
    # Construct bar title
    title = "\n".join([field_name, str(time_value)])
    if nb_components > 1:
        title = "\n".join([title, vector_mode])

    # Create scalar bar
    scalar_bar = pv.CreateScalarBar(Enabled=1, LabelFontSize=10,
                                    TitleFontSize=8, NumberOfLabels=5)
    scalar_bar.Title = title
    scalar_bar.LookupTable = lookup_table

    # Add the scalar bar to the view
    pv.GetRenderView().Representations.append(scalar_bar)

    # Reassign the current bar
    _current_bar = scalar_bar
    
    return scalar_bar


def get_lookup_table(field_name, nb_components, vector_mode='Magnitude'):
    """Get lookup table for the given field."""
    lookup_table = pv.GetLookupTableForArray(field_name, nb_components)

    if vector_mode == 'Magnitude':
        lookup_table.VectorMode = vector_mode
    elif vector_mode == 'X':
        lookup_table.VectorMode = 'Component'
        lookup_table.VectorComponent = 0
    elif vector_mode == 'Y':
        lookup_table.VectorMode = 'Component'
        lookup_table.VectorComponent = 1
    elif vector_mode=='Z':
        lookup_table.VectorMode = 'Component'
        lookup_table.VectorComponent = 2
    else:
        raise ValueError("Incorrect vector mode: {0}".format(vector_mode))
        
    lookup_table.Discretize = 0
    lookup_table.ColorSpace = 'HSV'
    lookup_table.LockScalarRange = 0
    
    return lookup_table


def _get_group_mesh_name(full_group_name):
    """Return mesh name of the group by its full name."""
    group_name = full_group_name.split('/')[1]

    return group_name


def _get_group_entity(full_group_name):
    """Return entity type of the group by its full name."""
    entity = full_group_name.split('/')[2]

    entity_type = EntityType.get_type(entity)
    
    return entity_type


def _get_group_short_name(full_group_name):
    """Return short name of the group by its full name."""
    short_name = full_group_name.split('/')[3]

    return short_name


def _get_mesh_names(proxy):
    """Return all mesh names in the given proxy as a set."""
    groups = proxy.Groups.Available
    mesh_names = set([_get_group_mesh_name(item) for item in groups])
    
    return mesh_names


def _get_group_names(proxy, mesh_name, entity_type, wo_nogroups=False):
    """Return full names of all groups of the given entity type
    from the mesh with the given name as a list.
    """
    groups = proxy.Groups.Available

    condition = lambda item: (_get_group_mesh_name(item) == mesh_name and
                              _get_group_entity(item) == entity_type)
    group_names = [item for item in groups if condition(item)]

    if wo_nogroups:
        # Remove "No_Group" group
        without_no_groups = lambda item: _get_group_short_name(item) != "No_Group"
        group_names = filter(without_no_groups, group_names)
    
    return group_names


def GetTime(theProxy, theTimeStampNumber):
    """Get time value by timestamp number."""
    # Check timestamp number
    timestamps = theProxy.TimestepValues.GetData()
    if ((theTimeStampNumber - 1) not in xrange(len(timestamps))):
        raise ValueError("Timestamp number is out of range: {0}".format(theTimeStampNumber))

    # Return time value
    return timestamps[theTimeStampNumber-1]


# Functions for building Post-Pro presentations
def ScalarMapOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                     theVectorMode='Magnitude'):
    """Creates Scalar Map presentation on the field.

    Keyword arguments:
    theProxy -- the pipeline object, containig data
    theEntityType -- the entity type from PrsTypeEnum
    theFieldName -- the field name
    theTimeStampNumber -- the number of time step (1, 2, ...)
    theVectorMode -- the mode of transformation of vector values into scalar values,
    applicable only if the field contains vector values.
    Possible modes: 'Magnitude' - vector module; 'X', 'Y', 'Z' - vector components.
    
    """
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)
    
    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)

    # Show pipeline object
    scalarmap = pv.GetRepresentation(theProxy)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set properties
    scalarmap.ColorAttributeType = EntityType.get_pvtype(theEntityType)
    scalarmap.ColorArrayName = theFieldName
    scalarmap.LookupTable = lookup_table

    # Add scalar bar
    bar_title = theFieldName + ", " + str(time_value)
    if (nb_components > 1):
        bar_title += "\n" + theVectorMode
    add_scalar_bar(theFieldName, nb_components, theVectorMode,
                   lookup_table, time_value)
        
    return scalarmap


def CutPlanesOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                     theNbPlanes=10, theOrientation = Orientation.YZ, theDisplacement=0.5,
                     theVectorMode='Magnitude'):
    """Creates cut planes on the given field."""
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0
    
    # Create slice filter
    slice = pv.Slice(theProxy)
    slice.SliceType="Plane"

    # Set cut planes normal and get range
    normal = [0.0, 0.0, 0.0]
    val_range = []
    if (theOrientation == Orientation.XY):
        normal[2] = 1.0
        val_range = get_z_range(theProxy)
    elif (theOrientation == Orientation.ZX):
        normal[1] = 1.0
        val_range = get_y_range(theProxy)
    elif (theOrientation == Orientation.YZ):
        normal[0] = 1.0
        val_range = get_x_range(theProxy)
        
    slice.SliceType.Normal = normal

    # Set cut planes positions
    slice.SliceOffsetValues = GetPositions(val_range, theNbPlanes, theDisplacement)
    cut_planes = pv.GetRepresentation(slice)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set properties
    cut_planes.ColorAttributeType = EntityType.get_pvtype(theEntityType)
    cut_planes.ColorArrayName = theFieldName
    cut_planes.LookupTable = lookup_table

    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return cut_planes


def CutLinesOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                    theNbLines=10,
                    theOrientation1 = Orientation.XY, theOrientation2 = Orientation.YZ,
                    theDisplacement1=0.5, theDisplacement2=0.5,
                    theVectorMode='Magnitude'):
    """Creates cut lines on the given field."""
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Create base plane
    base_plane = pv.Slice(theProxy)
    base_plane.SliceType="Plane"

    # Set base plane normal and get position
    base_normal = [0.0, 0.0, 0.0]
    valrange_base = []
    if (theOrientation1 == Orientation.XY):
        base_normal[2] = 1.0
        valrange_base = get_z_range(theProxy)
    elif (theOrientation1 == Orientation.ZX):
        base_normal[1] = 1.0
        valrange_base = get_y_range(theProxy)
    elif (theOrientation1 == Orientation.YZ):
        if (theOrientation2 == theOrientation1):
            theOrientation2 = Orientation.ZX
        base_normal[0] = 1.0
        valrange_base = get_x_range(theProxy)
        
    base_plane.SliceType.Normal = base_normal

    # Set base plane position
    base_plane.SliceOffsetValues = GetPositions(valrange_base, 1, theDisplacement1)

    # Check base plane
    base_plane.UpdatePipeline()
    if (base_plane.GetDataInformation().GetNumberOfCells() == 0):
        base_plane = theProxy
    
    # Create cutting planes
    cut_planes = pv.Slice(base_plane)
    cut_planes.SliceType="Plane"

    # Set cutting planes normal and get positions
    cut_normal = [0.0, 0.0, 0.0]
    valrange_cut = []
    if (theOrientation2 == Orientation.XY):
        cut_normal[2] = 1.0
        valrange_cut = get_z_range(theProxy)
    elif (theOrientation2 == Orientation.ZX):
        cut_normal[1] = 1.0
        valrange_cut = get_y_range(base_plane)
    elif (theOrientation2 == Orientation.YZ):
        cut_normal[0] = 1.0
        valrange_cut = get_x_range(base_plane)
        
    cut_planes.SliceType.Normal = cut_normal

    # Set cutting planes position
    cut_planes.SliceOffsetValues = GetPositions(valrange_cut, theNbLines,
                                                theDisplacement2)
    cut_lines = pv.GetRepresentation(cut_planes)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]
    
    # Set properties
    cut_lines.ColorAttributeType = EntityType.get_pvtype(theEntityType)
    cut_lines.ColorArrayName = theFieldName
    cut_lines.LookupTable = lookup_table

    # Set wireframe represenatation mode
    cut_lines.Representation = 'Wireframe'

    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components,
                   theVectorMode, lookup_table, time_value)
    return cut_lines


def VectorsOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                   theScaleFactor=-1, theGlyphPos=GlyphPos.TAIL,
                   theIsColored=False, theVectorMode='Magnitude'):
    """Creates vectors on the given field."""
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)

    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Extract only groups with data for the field
    new_proxy = extract_groups_for_field(theProxy, theFieldName, theEntityType)
    source = new_proxy
    
    # Cell centers
    if is_data_on_cells(theProxy, theFieldName):
        cell_centers = pv.CellCenters(source)
        cell_centers.VertexCells = 1
        source = cell_centers

    # Glyph
    glyph = pv.Glyph(source)
    glyph.Vectors = theFieldName
    glyph.ScaleMode = 'vector'
    glyph.MaskPoints = 0

    glyph.GlyphType = "2D Glyph"
    glyph.GlyphType.GlyphType = 'Arrow'
    if (theGlyphPos == GlyphPos.TAIL):
        glyph.GlyphType.Center = [0.5, 0.0, 0.0]
    elif (theGlyphPos == GlyphPos.HEAD):
        glyph.GlyphType.Center = [-0.5, 0.0, 0.0]
    elif (theGlyphPos == GlyphPos.CENTER):
        glyph.GlyphType.Center = [0.0, 0.0, 0.0]
        
    if (theScaleFactor > 0):
        glyph.SetScaleFactor = theScaleFactor
    else:
        def_scale = get_default_scale(PrsTypeEnum.DEFORMEDSHAPE, new_proxy, theEntityType, theFieldName)
        #@MZN DEBUG
        print("DEFAULT SCALE FACTOR = ", def_scale)
        glyph.SetScaleFactor = def_scale

    glyph.UpdatePipeline()

    vectors = pv.GetRepresentation(glyph)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set properties
    if (theIsColored):
        vectors.ColorArrayName = 'GlyphVector'
    else:
        vectors.ColorArrayName = ''
    vectors.LookupTable = lookup_table
    
    vectors.LineWidth = 1.0

    # Set wireframe represenatation mode
    vectors.Representation = 'Wireframe'

    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return vectors


def DeformedShapeOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber, theScaleFactor=-1,
                         theIsColored=False, theVectorMode='Magnitude'):
    """Creates Defromed shape on the given field."""
    # We don't need mesh parts with no data on them
    select_cells_with_data(theProxy)
    
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)

    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Extract only groups with data for the field
    new_proxy = extract_groups_for_field(theProxy, theFieldName, theEntityType)

    # Do merge
    source = pv.MergeBlocks(new_proxy)
    
    # Cell data to point data
    if is_data_on_cells(theProxy, theFieldName):
        cell_to_point = pv.CellDatatoPointData()
        cell_to_point.PassCellData = 1
        source = cell_to_point

    vector_array = theFieldName
    # If the given vector array has only 2 components, add the third one
    #@MZN: workaround: paraview doesn't treat 2-component array as vectors
    if nb_components == 2:
        calculator = get_add_component_calc(source, EntityType.NODE, theFieldName)
        vector_array = calculator.ResultArrayName
        source = calculator
        
    # Warp by vector
    warp_vector = pv.WarpByVector(source)
    warp_vector.Vectors = [theFieldName]
    if (theScaleFactor > 0):
        warp_vector.ScaleFactor = theScaleFactor
    else:
        def_scale = get_default_scale(PrsTypeEnum.DEFORMEDSHAPE, theProxy, theEntityType, theFieldName)
        #@MZN DEBUG
        print("DEFAULT SCALE FACTOR = ", def_scale)
        warp_vector.ScaleFactor = def_scale
        
    defshape = pv.GetRepresentation(warp_vector)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set properties
    if (theIsColored):
        defshape.ColorAttributeType = EntityType.get_pvtype(theEntityType)
        defshape.ColorArrayName = theFieldName
    else:
        defshape.ColorArrayName = ''
    defshape.LookupTable = lookup_table

    # Set wireframe represenatation mode
    defshape.Representation = 'Wireframe'

    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return defshape


def DeformedShapeAndScalarMapOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber, theScaleFactor=-1,
                                     theScalarEntityType=None, theScalarFieldName=None,
                                     theVectorMode='Magnitude'):
    """Creates Defromed shape And Scalar Map on the given field."""
    # We don't need mesh parts with no data on them
    select_cells_with_data(theProxy)
    
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
        
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Set scalar field by default
    scalar_entity = theScalarEntityType
    scalar_field = theScalarFieldName
    if (scalar_entity is None) or (scalar_field is None):
        scalar_entity = theEntityType
        scalar_field = theFieldName

    # Extract only groups with data for the field
    new_proxy = extract_groups_for_field(theProxy, theFieldName, theEntityType)

    # Do merge
    source = pv.MergeBlocks(new_proxy)

    # Cell data to point data
    if is_data_on_cells(theProxy, theFieldName):
        cell_to_point = pv.CellDatatoPointData(source)
        cell_to_point.PassCellData = 1
        source = cell_to_point

    vector_array = theFieldName
    # If the given vector array has only 2 components, add the third one
    #@MZN: workaround: paraview doesn't treat 2-component array as vectors
    if nb_components == 2:
        calculator = get_add_component_calc(source, EntityType.NODE, theFieldName)
        vector_array = calculator.ResultArrayName
        source = calculator
        
    # Warp by vector
    warp_vector = pv.WarpByVector(source)
    warp_vector.Vectors = [vector_array]
    if theScaleFactor > 0:
        warp_vector.ScaleFactor = theScaleFactor
    else:
        def_scale = get_default_scale(PrsTypeEnum.DEFORMEDSHAPE, new_proxy, theEntityType, theFieldName)
        #@MZN DEBUG
        print("DEFAULT SCALE FACTOR = ", def_scale)
        warp_vector.ScaleFactor = def_scale
    
    defshapeandmap = pv.GetRepresentation(warp_vector)
   
    # Get lookup table
    lookup_table = get_lookup_table(scalar_field, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, scalar_entity, scalar_field, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]
    
    # Set properties
    defshapeandmap.ColorArrayName = scalar_field
    defshapeandmap.LookupTable = lookup_table
    defshapeandmap.ColorAttributeType = EntityType.get_pvtype(scalar_entity)
     
    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return defshapeandmap


def Plot3DOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                  theOrientation = Orientation.AUTO, thePosition=0.5, theIsRelative=True,
                  theScaleFactor=-1, theIsContourPrs=False, theNbOfContours=32,
                  theVectorMode='Magnitude'):
    """Creates plot 3d on the given field."""
    # We don't need mesh parts with no data on them
    select_cells_with_data(theProxy)
    
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Extract only groups with data for the field
    new_proxy = extract_groups_for_field(theProxy, theFieldName, theEntityType)

    # Do merge
    merge_blocks = pv.MergeBlocks(new_proxy)
    merge_blocks.UpdatePipeline()

    poly_data = None

    # Cutting plane
    
    # Define orientation if necessary (auto mode)
    orientation = theOrientation
    if (orientation == Orientation.AUTO):
        orientation = GetOrientation(theProxy)

    # Get cutting plane normal
    normal = GetNormalByOrientation(orientation)
      
    if (not is_planar_input(theProxy)):
        # Create slice filter
        slice = pv.Slice(merge_blocks)
        slice.SliceType="Plane"

        # Set cutting plane normal
        slice.SliceType.Normal = normal

        # Set cutting plane position
        val_range = GetRangeForOrientation(theProxy, orientation)
        
        if (theIsRelative):
            slice.SliceOffsetValues = GetPositions(val_range, 1, thePosition)
        else:
            slice.SliceOffsetValues = thePosition

        slice.UpdatePipeline()
        poly_data = slice

    use_normal = 0
    # Geometry filter
    if (not poly_data or poly_data.GetDataInformation().GetNumberOfCells() == 0):
        geometry_filter = pv.GeometryFilter(merge_blocks)
        poly_data = geometry_filter
        use_normal = 1 #MZN: temporary workaround
    
    warp_scalar = None
    plot3d = None
    source = poly_data
    
    if is_data_on_cells(poly_data, theFieldName):
        # Cell data to point data
        cell_to_point = pv.CellDatatoPointData(poly_data)
        cell_to_point.PassCellData = 1
        source = cell_to_point
        
    scalars = ['POINTS', theFieldName]
    
    # Transform vector array to scalar array if necessary
    if (nb_components > 1):
        calculator = get_calc_magnitude(source, EntityType.NODE, theFieldName)
        scalars = ['POINTS', calculator.ResultArrayName]
        source = calculator
    
    # Warp by scalar
    warp_scalar = pv.WarpByScalar(source)
    warp_scalar.Scalars = scalars
    warp_scalar.Normal = normal
    warp_scalar.UseNormal = use_normal
    if (theScaleFactor > 0):
        warp_scalar.ScaleFactor = theScaleFactor
    else:
        def_scale = get_default_scale(PrsTypeEnum.PLOT3D, theProxy, theEntityType, theFieldName)
        warp_scalar.ScaleFactor = def_scale
        
    warp_scalar.UpdatePipeline()
    
    if (theIsContourPrs):
        # Contours
        contour = pv.Contour(warp_scalar)
        contour.PointMergeMethod = "Uniform Binning"
        contour.ContourBy = ['POINTS', theFieldName]
        scalar_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
        contour.Isosurfaces = GetPositions(scalar_range, theNbOfContours, 0)
        contour.UpdatePipeline()
        plot3d = pv.GetRepresentation(contour)
    else:
        plot3d = pv.GetRepresentation(warp_scalar)
   
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set properties
    plot3d.ColorAttributeType = EntityType.get_pvtype(theEntityType)
    plot3d.ColorArrayName = theFieldName
    plot3d.LookupTable = lookup_table
  
    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return plot3d


def IsoSurfacesOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                       theRange = None, theNbOfSurfaces=10, theIsColored = True,
                       theColor = None, theVectorMode='Magnitude'):
    """Creates Iso Surfaces on the given field."""
    # We don't need mesh parts with no data on them
    select_cells_with_data(theProxy)
    
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)
    
    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Extract only groups with data for the field
    new_proxy = extract_groups_for_field(theProxy, theFieldName, theEntityType)

    # Do merge
    source = pv.MergeBlocks(new_proxy)
    
    # Transform cell data into point data if necessary
    if is_data_on_cells(theProxy, theFieldName):
        cell_to_point = pv.CellDatatoPointData(source)
        cell_to_point.PassCellData = 1
        source = cell_to_point
        
    contour_by = ['POINTS', theFieldName]

    # Transform vector array to scalar array if necessary
    if (nb_components > 1):
        calculator = get_calc_magnitude(source, EntityType.NODE, theFieldName)
        contour_by = ['POINTS', calculator.ResultArrayName]
        source = calculator

    # Contour filter settings
    contour = pv.Contour(source)
    contour.ComputeScalars = 1
    contour.ContourBy = contour_by

    # Specify the range
    scalar_range = theRange
    if (scalar_range is None):
        scalar_range = get_data_range(theProxy, theEntityType, theFieldName)
        # Cut off the range
        if (scalar_range[0] <= scalar_range[1]):
            delta = abs(scalar_range[1] - scalar_range[0]) * GAP_COEFFICIENT
            scalar_range[0] += delta
            scalar_range[1] -= delta

    # Calculate contour values for the range
    surfaces = []
    for i in xrange(theNbOfSurfaces):
        pos = scalar_range[0] + i*(scalar_range[1] - scalar_range[0])/(theNbOfSurfaces - 1)
        surfaces.append(pos)

    # Set contour values
    contour.Isosurfaces = surfaces
    
    # Get Iso Surfaces representation
    isosurfaces = pv.GetRepresentation(contour)

    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]

    # Set display properties
    if (theIsColored):
        isosurfaces.ColorAttributeType = EntityType.get_pvtype(theEntityType)
        isosurfaces.ColorArrayName = theFieldName
    else:
        isosurfaces.ColorArrayName = ''
        if (theColor is not None):
            isosurfaces.DiffuseColor = theColor
    isosurfaces.LookupTable = lookup_table
  
    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return isosurfaces


def GaussPointsOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber,
                       theIsDeformed=False, theScaleFactor=-1,
                       theIsColored=True, theColor=None,
                       thePrimitiveType=GaussType.SPRITE,
                       theMaxPixelSize=256, theMagnification=1, theVectorMode='Magnitude'):
    """Creates Gauss points on the given field.
    
    If scale factor is defined deformation will be applied.
    """
    # Check vector mode
    nb_components = get_nb_components(theProxy, theEntityType, theFieldName)
    _check_vector_mode(theVectorMode, nb_components)
    
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Quadrature point arrays
    qp_arrays = theProxy.QuadraturePointArrays

    # If no gauss points are defined, use cell centers
    source = None
    if theFieldName in qp_arrays:
        generate_qp = pv.GenerateQuadraturePoints(theProxy)
        generate_qp.SelectSourceArray = ['CELLS', theFieldName] #@MZN
        source = generate_qp
    else:
        cell_centers = pv.CellCenters(theProxy)
        cell_centers.VertexCells = 1
        source = cell_centers

    source.UpdatePipeline()

    # Check if deformation enabled
    if theIsDeformed and nb_components > 1:
        vector_array = theFieldName
        # If the given vector array has only 2 components, add the third one
        #@MZN: workaround: paraview doesn't treat 2-component array as vectors
        if nb_components == 2:
            calculator = get_add_component_calc(source, EntityType.NODE, theFieldName)
            vector_array = calculator.ResultArrayName
            source = calculator
        
        # Warp by vector
        warp_vector = pv.WarpByVector(source)
        warp_vector.Vectors = [vector_array]
        if (theScaleFactor > 0):
            warp_vector.ScaleFactor = theScaleFactor
        else:
            def_scale = get_default_scale(PrsTypeEnum.DEFORMEDSHAPE,
                                          theProxy, theEntityType, theFieldName)
            warp_vector.ScaleFactor = def_scale
        warp_vector.UpdatePipeline()
        source = warp_vector

    # Point sprite settings
    
    #@MZN: Add magnification etc.
    gausspnt = pv.GetRepresentation(source)
    gausspnt.Representation = 'Point Sprite'
    gausspnt.MaxPixelSize = theMaxPixelSize
    gausspnt.RenderMode = GaussType.get_mode(thePrimitiveType)
    #@MZN Debug
    print("Set mode = ", GaussType.get_mode(thePrimitiveType))
    print("Render mode = ", gausspnt.RenderMode)
    if thePrimitiveType == GaussType.SPRITE:
        # Set texture like in visu
        i = 0

    gausspnt.RadiusInitialized = 1
    gausspnt.RadiusTransferFunctionEnabled = 1
    gausspnt.RadiusMode = 'Scalar'
    #gausspnt.RadiusScalarRange = radius_range
    #gausspnt.RadiusArray = [None, theFieldName]
    if nb_components > 0:
        gausspnt.RadiusVectorComponent = _get_vector_component(theVectorMode)
    gausspnt.RadiusUseScalarRange = 1
    gausspnt.RadiusIsProportional = 1
    gausspnt.RadiusProportionalFactor = theMagnification
    
    # Get lookup table
    lookup_table = get_lookup_table(theFieldName, nb_components, theVectorMode)

    # Set field range if necessary
    data_range = get_data_range(theProxy, theEntityType, theFieldName, theVectorMode)
    lookup_table.LockScalarRange = 1
    lookup_table.RGBPoints = [data_range[0], 0, 0, 1, data_range[1], 1, 0, 0]
    
    # Set display properties
    if (theIsColored):
        gausspnt.ColorAttributeType = EntityType.get_pvtype(theEntityType)
        gausspnt.ColorArrayName = theFieldName
    else:
        gausspnt.ColorArrayName = ''
        if theColor:
            gausspnt.DiffuseColor = theColor
            
    gausspnt.LookupTable = lookup_table
    
    # Add scalar bar
    add_scalar_bar(theFieldName, nb_components, theVectorMode, lookup_table, time_value)
    
    return gausspnt


def StreamLinesOnField(theProxy, theEntityType, theFieldName, theTimeStampNumber, theScaleFactor=-1,
                       theIsColored=False, theVectorMode='Magnitude'):
    """Creates Stream Lines on the given field."""
    # Get time value
    time_value = GetTime(theProxy, theTimeStampNumber)

    # Set timestamp
    pv.GetRenderView().ViewTime = time_value
    pv.UpdatePipeline(time_value, theProxy)
    
    # Hide initial object
    rep = pv.GetRepresentation(theProxy)
    rep.Visibility = 0

    # Cell centers
    if is_data_on_cells(theProxy, theFieldName):
        cell_centers = pv.CellCenters(theProxy)
        cell_centers.VertexCells = 1
    
    streamlines = None
    
    return streamlines


def MeshOnEntity(theProxy, theMeshName, theEntityType):
    """Creates submesh of the entity type for the mesh.

    Keyword arguments:
    theProxy -- the pipeline object, containig data
    theMeshName -- the mesh name
    theEntityType -- the entity type
    
    """
    # Select all cell types
    select_all_cells(theProxy)
    
    # Get subset of groups on the given entity
    subset = _get_group_names(theProxy, theMeshName, theEntityType)

    # Select only groups of the given entity type
    theProxy.Groups = subset
    theProxy.UpdatePipeline()

    # Get representation object if the submesh is not empty
    prs = None
    if (theProxy.GetDataInformation().GetNumberOfPoints() or
        theProxy.GetDataInformation().GetNumberOfCells()):
        prs = pv.GetRepresentation(theProxy)
    
    return prs


def MeshOnGroup(theProxy, theGroupName):
    """Creates submesh on the group.

    Keyword arguments:
    theProxy -- the pipeline object, containig data
    theGroupName -- the full group name
    
    """
    # Select all cell types
    select_all_cells(theProxy)
    
    # Select only the group with the given name
    #@MZN: DEBUG
    print ("theProxy.Groups = ", [theGroupName])
    one_group = [theGroupName]
    theProxy.Groups = one_group
    theProxy.UpdatePipeline()
    
    # Get representation object if the submesh is not empty
    prs = None

    #@MZN: DEBUG
    print ("theProxy.GetDataInformation().GetNumberOfPoints() = ", theProxy.GetDataInformation().GetNumberOfPoints())
    print ("theProxy.GetDataInformation().GetNumberOfCells() = ", theProxy.GetDataInformation().GetNumberOfCells())
    print ("theProxy.Groups.GetData() = : ", theProxy.Groups.GetData())

    # Check if the group was set
    if theProxy.Groups.GetData() == one_group:
        group_entity = _get_group_entity(theGroupName)
        # Check if the submesh is not empty
        nb_items = 0
        if group_entity == EntityType.NODE:
            nb_items = theProxy.GetDataInformation().GetNumberOfPoints()
        elif group_entity == EntityType.CELL:
            nb_items = theProxy.GetDataInformation().GetNumberOfCells()
        
        if nb_items:
            prs = pv.GetRepresentation(theProxy)
        
    return prs


def CreatePrsForFile(theParavis, theFileName, thePrsTypeList, thePictureDir,
                     thePictureExt, theIsAutoDelete = 0):
    """Build presentations of the given types for all fields of the given file."""
    print("Import {0}...".format(theFileName.split(os.sep)[-1]), end='')
    theParavis.ImportFile(theFileName)
    proxy = pv.GetActiveSource()
    if proxy is None:
        raise RuntimeError("Error: can't import file.")
    else: print("OK")

    # Get view
    view = pv.GetRenderView()

    # Create required presentations for the proxy
    CreatePrsForProxy(proxy, view, thePrsTypeList, thePictureDir, thePictureExt, theIsAutoDelete)


def _create_prs(prs_type, proxy, field_entity, field_name, timestamp_nb):
    """Internal method. Build presentation of the given type on the given field and
    timestamp number. Set the presentation properties like visu.CreatePrsForResult() do.
    """
    prs = None
    
    if prs_type == PrsTypeEnum.SCALARMAP:
        prs = ScalarMapOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.CUTPLANES:
        prs = CutPlanesOnField(proxy, field_entity, field_name, timestamp_nb,
                               theOrientation = Orientation.ZX)
    elif prs_type == PrsTypeEnum.CUTLINES:
        prs = CutLinesOnField(proxy, field_entity, field_name, timestamp_nb,
                              theOrientation1 = Orientation.XY, theOrientation2 = Orientation.ZX)
    elif prs_type == PrsTypeEnum.DEFORMEDSHAPE:
        prs = DeformedShapeOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.DEFORMEDSHAPESCALARMAP:
        prs = DeformedShapeAndScalarMapOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.VECTORS:
        prs = VectorsOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.PLOT3D:
        prs = Plot3DOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.ISOSURFACES:
        prs = IsoSurfacesOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.GAUSSPOINTS:
        prs = GaussPointsOnField(proxy, field_entity, field_name, timestamp_nb)
    elif prs_type == PrsTypeEnum.STREAMLINES:
        prs = StreamLinesOnField(proxy, field_entity, field_name, timestamp_nb)
    else:
        raise ValueError("Unexistent presentation type.")
        
    return prs

        
def CreatePrsForProxy(theProxy, theView, thePrsTypeList, thePictureDir, thePictureExt, theIsAutoDelete):
    """Build presentations of the given types for all fields of the given proxy.
    Save snapshots in graphics files (type depends on the given extension).
    Stores the files in the given directory.
    
    Keyword arguments:
    theProxy -- the pipeline object, containig data arrays
    theView  -- the render view
    thePrsTypeList  -- the list of presentation types to build
    thePictureDir   -- the directory path for saving snapshots
    thePictureExt   -- graphics files extension (determines file type)
    
    """
    
    field_names = list(theProxy.PointArrays.GetData())
    nb_on_nodes = len(field_names)
    field_names.extend(theProxy.CellArrays.GetData())

    # Add path separator to the end of picture path if necessery
    if not thePictureDir.endswith(os.sep):
        thePictureDir += os.sep
        
    # Mesh Presentation
    if PrsTypeEnum.MESH in thePrsTypeList:
        # Create Mesh presentation. Build all possible submeshes.

        # Remember the current state
        groups = list(theProxy.Groups)

        # Iterate on meshes
        mesh_names = _get_mesh_names(theProxy)
        for mesh_name in mesh_names:
            # Build mesh on nodes and cells
            for entity_type in (EntityType.NODE, EntityType.CELL):
                entity_name = EntityType.get_name(entity_type)
                if if_possible(theProxy, mesh_name, entity_type, PrsTypeEnum.MESH):
                    print("Creating submesh on {0} for '{1}' mesh... ".format(entity_name,
                                                                              mesh_name),
                          end='')
                    prs = MeshOnEntity(theProxy, mesh_name, entity_type)
                    if prs is None :
                        print("FAILED")
                        continue
                    else:
                        print("OK")
                    # Construct image file name
                    pic_name = "{folder}{mesh}_{type}.{ext}".format(folder=thePictureDir,
                                                                    mesh=mesh_name,
                                                                    type=entity_name,
                                                                    ext=thePictureExt)

                    # Show and dump the presentation into a graphics file
                    process_prs_for_test(prs, theView, pic_name, False)
                    
                # Build submesh on all groups of the mesh
                mesh_groups = _get_group_names(theProxy, mesh_name, entity_type, wo_nogroups=True)
                for group in mesh_groups:
                    print("Creating submesh on group {0}... ".format(group), end='')
                    prs = MeshOnGroup(theProxy, group)
                    if prs is None :
                        print("FAILED")
                        continue
                    else:
                        print("OK")
                    # Construct image file name
                    pic_name = "{folder}{group}.{ext}".format(folder=thePictureDir,
                                                              group=group.replace('/','_'),
                                                              ext=thePictureExt)
                
                    # Show and dump the presentation into a graphics file
                    process_prs_for_test(prs, theView, pic_name, False)
                
        # Restore the state
        theProxy.Groups = groups
        theProxy.UpdatePipeline()
 
    # Presentations on fields
    for (i, field_name) in enumerate(field_names):
        # Select only the current field: necessary for getting the right timestamps
        field_entity = None
        if (i >= nb_on_nodes):
            field_entity = EntityType.CELL
            theProxy.PointArrays.DeselectAll()
            theProxy.CellArrays = [field_name]
        else:
            field_entity = EntityType.NODE
            theProxy.CellArrays.DeselectAll()
            theProxy.PointArrays = [field_name]

        # Get timestamps
        theProxy.UpdatePipelineInformation()
        timestamps = theProxy.TimestepValues.GetData()

        for prs_type in thePrsTypeList:
            # Ignore mesh presentation
            if prs_type == PrsTypeEnum.MESH:
                continue
            
            # Get name of presentation type
            prs_name = PrsTypeEnum.get_name(prs_type)

            # Build the presentation if possible
            possible = if_possible(theProxy, field_name, field_entity, prs_type)
            if possible:
                # Presentation type for graphics file name
                f_prs_type = prs_name.replace(' ', '').upper()
                
                for timestamp_nb in xrange(1, len(timestamps)+1):
                    time = timestamps[timestamp_nb-1]
                    print("Creating {0} on {1}, time = {2}... ".format(prs_name, field_name, time),
                          end='')
                    prs = _create_prs(prs_type, theProxy, field_entity, field_name, timestamp_nb)
                    if prs is None :
                        print("FAILED")
                        continue
                    else:
                        print("OK")
            
                    # Construct image file name
                    pic_name = "{folder}{field}_{time}_{type}.{ext}".format(folder=thePictureDir,
                                                                            field=field_name,
                                                                            time=time,
                                                                            type=f_prs_type,
                                                                            ext=thePictureExt)

                    # Show and dump the presentation into a graphics file
                    process_prs_for_test(prs, theView, pic_name)
