#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2024  CEA, EDF
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
# Author : Anthony Geay

import os
import sys

from medcoupling import *
from paraview.simple import *
from MEDReaderHelper import WriteInTmpDir,RetriveBaseLine

def GenerateCase():
  """ This test is non regression test to check non regression of EDF 8761. ELNO Mesh filter on vector field with 4 comps cut of using GenerateVectors"""

  fname="testMEDReader13.med"
  #

  m=MEDCouplingUMesh("mesh",2)
  m.setCoords(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,1.,1.,0.],4,3))
  m.allocateCells()
  m.insertNextCell(NORM_TRI3,[0,1,3]) ; m.insertNextCell(NORM_TRI3,[1,2,3])
  f=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f.setName("fieldELNO") ; f.setMesh(m)
  arr=DataArrayDouble([0.2,1.1,0.7,0.5,-0.3,0.4])
  f.setArray(DataArrayDouble.Meld(4*[arr]))
  f.checkConsistencyLight()
  WriteField(fname,f,True)
  return fname

@WriteInTmpDir
def test(baseline_file):
  fname = GenerateCase()
  ################### MED write is done -> Go to MEDReader
  testMEDReader13_med = MEDReader( FileNames=[fname] )

  testMEDReader13_med.VectorsProperty = 1
  testMEDReader13_med.FieldsStatus = ['TS0/mesh/ComSup0/fieldELNO@@][@@GSSNE']

  if '-D' not in sys.argv:
    RenderView1 = GetRenderView()
    RenderView1.CameraPosition = [1.0, 0.5, 10000.0]

    RenderView1.CameraPosition = [1.0, 0.5, 4.319751617610021]

    ELNOfieldToSurface3 = ELNOfieldToSurface(Input=testMEDReader13_med)

    DataRepresentation2 = Show()
    #DataRepresentation2.ConstantRadius = 1.9999333620071411
    DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
    #DataRepresentation2.PointGaussianDefaultsInitialized = 1
    DataRepresentation2.SelectionPointFieldDataArrayName = 'fieldELNO'
    DataRepresentation2.SelectionCellFieldDataArrayName = 'FamilyIdCell'
    #DataRepresentation2.SelectInputVectors = ['POINTS', 'fieldELNO_Vector']
    DataRepresentation2.ScalarOpacityUnitDistance = 1.7746382108908556
    DataRepresentation2.Texture = []
    #DataRepresentation2.ExtractedBlockIndex = 1
    #DataRepresentation2.RadiusRange = [6.666666740784422e-05, 1.9999333620071411]
    DataRepresentation2.ScaleFactor = 0.19998666953397334

    #DataRepresentation2.RadiusRange = [6.66667e-05, 1.99993]
    DataRepresentation2.ColorArrayName = ('POINT_DATA', 'fieldELNO_Vector')

    a3_fieldELNO_Vector_PVLookupTable = GetLookupTableForArray( "fieldELNO_Vector", 3, RGBPoints=[0.3464101615137755, 0.23, 0.299, 0.754, 1.1258330249197703, 0.865, 0.865, 0.865, 1.9052558883257653, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

    a3_fieldELNO_Vector_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.3464101615137755, 0.0, 0.5, 0.0, 1.9052558883257653, 1.0, 0.5, 0.0] )

    RenderView1.ViewSize =[300,300]
    Render()

    # compare with baseline image
    import vtk.test.Testing
    from vtk.util.misc import vtkGetTempDir
    vtk.test.Testing.VTK_TEMP_DIR = vtk.util.misc.vtkGetTempDir()
    vtk.test.Testing.compareImage(GetActiveView().GetRenderWindow(), baseline_file,
                                                                threshold=1)
    vtk.test.Testing.interact()

if __name__ == "__main__":
  outImgName="testMEDReader13.png"
  baseline_file = RetriveBaseLine(outImgName)
  test(baseline_file)
  pass
