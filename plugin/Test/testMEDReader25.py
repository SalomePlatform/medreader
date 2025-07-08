#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2020-2025  CEA, EDF
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : Anthony Geay (EDF R&D)

__doc__ = """
Test of mixed polyhderon and hexahedron in the MEDReader plugin.
"""

import os
import platform
import sys
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
from vtk.util import numpy_support
import numpy as np
import medcoupling as mc
from MEDReaderHelper import WriteInTmpDir,RetriveBaseLine

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")

def generateCase(fname):
    # Inspired by build3DExtrudedUMesh_1 in medcoupling
    coords=[
        0.,0.,0., 1.,1.,0., 1.,1.25,0., 1.,0.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
        3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
        0.,0.,1., 1.,1.,1., 1.,1.25,1., 1.,0.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
        3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
        0.,0.,2., 1.,1.,2., 1.,1.25,2., 1.,0.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
        3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
        0.,0.,3., 1.,1.,3., 1.,1.25,3., 1.,0.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
        3.,2.,3., 0.,1.,3., 1.,3.,3., 2.,2.,3., 2.,3.,3.]

    conn=[
        # 0
        0,11,1,3,15,26,16,18,   1,2,4,7,13,6,-1,1,16,21,6,-1,6,21,28,13,-1,13,7,22,28,-1,7,4,19,22,-1,4,2,17,19,-1,2,1,16,17,-1,16,21,28,22,19,17,
        1,6,5,3,16,21,20,18,   13,10,9,6,28,25,24,21,
        11,8,7,4,2,1,-1,11,26,16,1,-1,1,16,17,2,-1,2,17,19,4,-1,4,19,22,7,-1,7,8,23,22,-1,8,11,26,23,-1,26,16,17,19,22,23,
        7,12,14,13,22,27,29,28,
        # 1
        15,26,16,18,30,41,31,33,   16,17,19,22,28,21,-1,16,31,36,21,-1,21,36,43,28,-1,28,22,37,43,-1,22,19,34,37,-1,19,17,32,34,-1,17,16,31,32,-1,31,36,43,37,34,32,
        16,21,20,18,31,36,35,33,   28,25,24,21,43,40,39,36,
        26,23,22,19,17,16,-1,26,41,31,16,-1,16,31,32,17,-1,17,32,34,19,-1,19,34,37,22,-1,22,23,38,37,-1,23,26,41,38,-1,41,31,32,34,37,38,
        22,27,29,28,37,42,44,43,
        # 2
        30,41,31,33,45,56,46,48,  31,32,34,37,43,36,-1,31,46,51,36,-1,36,51,58,43,-1,43,37,52,58,-1,37,34,49,52,-1,34,32,47,49,-1,32,31,46,47,-1,46,51,58,52,49,47,
        31,36,35,33,46,51,50,48,  43,40,39,36,58,55,54,51,
        41,38,37,34,32,31,-1,41,56,46,31,-1,31,46,47,32,-1,32,47,49,34,-1,34,49,52,37,-1,37,38,53,52,-1,38,41,56,53,-1,56,46,47,49,52,53,
        37,42,44,43,52,57,59,58]
    #

    mesh=mc.MEDCouplingUMesh.New()
    mesh.setName("mesh")
    mesh.setMeshDimension(3)
    mesh.allocateCells(18)
    #
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[0:8])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[8:51])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[51:59])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[59:67])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[67:110])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[110:118])
    #
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[118:126])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[126:169])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[169:177])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[177:185])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[185:228])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[228:236])
    #
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[236:244])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[244:287])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[287:295])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[295:303])
    mesh.insertNextCell(mc.NORM_POLYHED,43,conn[303:346])
    mesh.insertNextCell(mc.NORM_HEXA8,8,conn[346:354])
    #
    mesh.finishInsertingCells()

    myCoords=mc.DataArrayDouble.New()
    myCoords.setValues(coords,60,3)
    mesh.setCoords(myCoords)
    mesh.sortCellsInMEDFileFrmt()

    mm = mc.MEDFileUMesh()
    mm[0] = mesh

    mm.write(fname,2)
    f = mc.MEDCouplingFieldDouble(mc.ON_CELLS)
    f.setMesh(mesh)
    f.setArray(mesh.computeCellCenterOfMass().magnitude())
    f.setName("field")
    mc.WriteFieldUsingAlreadyWrittenMesh(fname,f)
    f2 = mc.MEDCouplingFieldDouble(mc.ON_NODES)
    f2.setMesh(mesh)
    f2.setArray(mesh.getCoords().magnitude())
    f2.setName("field2")
    mc.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
    return f,f2

@WriteInTmpDir
def test(baseline_file):
    fname = "testMEDReader25.med"
    f,f2 = generateCase(fname)
    reader = MEDReader(FileNames=[fname])
    reader.FieldsStatus = ['TS0/mesh/ComSup0/field@@][@@P0','TS0/mesh/ComSup0/field2@@][@@P1','TS0/mesh/ComSup0/mesh@@][@@P0']
    reader.TimesFlagsStatus = ['0000']

    # Check multiblock
    mb_data = servermanager.Fetch(reader)
    MyAssert(mb_data.GetNumberOfBlocks() == 1)

    # Check first block
    ds_data = mb_data.GetBlock(0)
    MyAssert(ds_data.GetNumberOfPoints()) == 60
    MyAssert(ds_data.GetNumberOfCells()) == 18
    if platform.system() == "Windows":
        MyAssert(ds_data.GetPointData().GetNumberOfArrays() == 1) # field2
        MyAssert(ds_data.GetCellData().GetNumberOfArrays() == 3) # familiyIdCell, field, mesh
    else:
        MyAssert(ds_data.GetPointData().GetNumberOfArrays() == 2) # field2, vtkGhostType
        MyAssert(ds_data.GetCellData().GetNumberOfArrays() == 4) # familiyIdCell, field, mesh, vtkGhostType
    # Check array content
    ds_data_ref_conn = np.array([0,11,1,3,15,26,16,18,1,6,5,3,16,21,20,18,13,10,9,6,28,25,24,21,7,12,14,13,22,27,29,28,15,26,16,18,
        30,41,31,33,16,21,20,18,31,36,35,33,28,25,24,21,43,40,39,36,22,27,29,28,37,42,44,43,30,41,31,33,45,56,46,48,31,36,35,33,46,
        51,50,48,43,40,39,36,58,55,54,51,37,42,44,43,52,57,59,58,1,2,4,6,7,13,16,17,19,21,22,28,1,2,4,7,8,11,16,17,19,22,23,26,16,
        17,19,21,22,28,31,32,34,36,37,43,16,17,19,22,23,26,31,32,34,37,38,41,31,32,34,36,37,43,46,47,49,51,52,58,31,32,34,37,38,41,46,47,49,52,53,56],dtype=np.int32)
    MyAssert(np.all(ds_data_ref_conn == numpy_support.vtk_to_numpy(ds_data.GetCells().GetConnectivityArray())))
    MyAssert(mc.DataArrayDouble(numpy_support.vtk_to_numpy(ds_data.GetCellData().GetArray("field") )).isEqual(f.getArray(),1e-12))
    MyAssert(mc.DataArrayDouble(numpy_support.vtk_to_numpy(ds_data.GetPointData().GetArray("field2"))).isEqual(f2.getArray(),1e-12))

    # Rendering test
    if '-D' not in sys.argv:
        # get active view
        renderView1 = GetActiveViewOrCreate('RenderView')

        # show data in view
        repr = Show(reader, renderView1)

        # reset view to fit data
        renderView1.ResetCamera()

        # set scalar coloring
        ColorBy(repr, ('CELLS', 'field'))

        # rescale color and/or opacity maps used to include current data range
        repr.RescaleTransferFunctionToDataRange(True)

        # do not show color bar/color legend
        repr.SetScalarBarVisibility(renderView1, False)

        # render
        Render()

         # compare with baseline image
        import vtk.test.Testing
        from vtk.util.misc import vtkGetTempDir
        vtk.test.Testing.VTK_TEMP_DIR = vtk.util.misc.vtkGetTempDir()
        vtk.test.Testing.compareImage(renderView1.GetRenderWindow(), baseline_file)
        vtk.test.Testing.interact()
        pass


if __name__ == "__main__":
    outImgName="testMEDReader25.png"
    baseline_file = RetriveBaseLine(outImgName)
    test(baseline_file)
