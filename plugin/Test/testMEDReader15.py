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
# Author : Anthony Geay (EDF R&D)

import os
import sys

from medcoupling import *
from paraview.simple import *
from MEDReaderHelper import WriteInTmpDir,RetriveBaseLine

paraview.simple._DisableFirstRenderCameraReset()

def GenerateCase():
  """ This test is a non regression test that checks the behaviour of MEDReader when a mesh has the same name than a field.
  """

  fname="testMEDReader15.med"

  zeName="zeName"
  c=DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(0.,1.,0.)])
  m=MEDFileUMesh()
  m.setCoords(c)
  m.setName(zeName)
  m.write(fname,2)
  f=MEDCouplingFieldDouble(ON_NODES)
  f.setName(zeName)
  f.setArray(DataArrayDouble([(-1.,1.,0.),(0.,1.,0.),(1.,1.,0.)]))
  tmp=MEDCouplingUMesh.Build0DMeshFromCoords(m.getCoords()) ; tmp.setName(zeName)
  f.setMesh(tmp)
  WriteFieldUsingAlreadyWrittenMesh(fname,f)
  return fname

@WriteInTmpDir
def test(baseline_file):
  fname = GenerateCase()
  reader=MEDReader(FileNames=[fname])
  ExpectedEntries=['TS0/zeName/ComSup0/zeName@@][@@P1','TS0/zeName/ComSup0/MESH@zeName@@][@@P1']
  assert(reader.GetProperty("FieldsTreeInfo")[::2]==ExpectedEntries)

  #
  glyph1=Glyph(Input=reader,GlyphType='Arrow',ScaleArray='FamilyIdNode',OrientationArray='zeName',GlyphMode='All Points',ScaleFactor=0.1,GlyphTransform='Transform2')

  if '-D' not in sys.argv:
    renderView1=GetActiveViewOrCreate('RenderView')
    renderView1.InteractionMode='3D'
    zeNameLUT = GetColorTransferFunction('zeName')
    zeNameLUT.RGBPoints = [1.0, 0.231373, 0.298039, 0.752941, 1.2071067811865475, 0.865003, 0.865003, 0.865003, 1.4142135623730951, 0.705882, 0.0156863, 0.14902]
    zeNameLUT.ScalarRangeInitialized = 1.
    zeNameLUT.VectorMode = 'Component'

    glyph1Display=Show(glyph1,renderView1)
    glyph1Display.ColorArrayName = ['POINTS', 'FamilyIdNode']
    glyph1Display.LookupTable = zeNameLUT
    # set scalar coloring
    ColorBy(glyph1Display, ('POINTS', 'zeName'))
    # rescale color and/or opacity maps used to include current data range
    glyph1Display.RescaleTransferFunctionToDataRange(True)
    # do not show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, False)
    #
    renderView1.ViewSize =[300,300]
    renderView1.GetRenderWindow().DoubleBufferOff()
    Render()

    # compare with baseline image
    import vtk.test.Testing
    from vtk.util.misc import vtkGetTempDir
    vtk.test.Testing.VTK_TEMP_DIR = vtk.util.misc.vtkGetTempDir()
    vtk.test.Testing.compareImage(GetActiveView().GetRenderWindow(), baseline_file,
                                                                threshold=1)
    vtk.test.Testing.interact()

if __name__ == "__main__":
  outImgName="testMEDReader15.png"
  baseline_file = RetriveBaseLine(outImgName)
  test(baseline_file)
  pass
