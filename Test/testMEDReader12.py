#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2014  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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

from MEDLoader import *

""" This test is non regression test to check that ExtractGroup filter works well on non unstructured meshes."""

fname="testMEDReader12.med"
outImgName="testMEDReader12.png"
#
p=DataArrayDouble([0.6075913659107736,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.5814307903048253,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.6075913659107737,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.4734426322438601,0.4734426322438599,0.5535261763251467,0.4734426322438599,0.4734426322438599,0.5535261763251467,0.4734426322438599,0.4734426322438599,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.5814307903048253,0.553526176325147,0.4734426322438601,0.4734426322438602,0.553526176325147,0.4734426322438603,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.5814307903048255,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.6075913659107737,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048255,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.6075913659107739,0.6075913659107737,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.5814307903048254,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.581430790304825,0.6075913659107737,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224386013,0.4734426322438599,0.5535261763251467,0.4734426322438599,0.4734426322438599,0.5535261763251467,0.4734426322438599,0.4734426322438599,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.5814307903048254,0.553526176325147,0.47344263224386013,0.4734426322438602,0.553526176325147,0.47344263224386035,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.4734426322438602,0.4734426322438602,0.553526176325147,0.5814307903048255,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.4734426322438599,0.47344263224385996,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.4734426322438602,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.47344263224385996,0.47344263224385996,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.581430790304825,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.553526176325147,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5535261763251467,0.5814307903048253,0.6075913659107737,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048255,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.5814307903048253,0.6075913659107739])
p.setInfoOnComponents(["- [-]"])
#
famIds=DataArrayInt([-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-2,-2,-3,-3,-4,-4,-4,-4,-3,-4,-4,-3,-4,-4,-4,-4,-3,-3,-2,-2,-3,-3,-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-4,-4,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-4,-4,-3,-3,-2,-2,-3,-3,-4,-4,-4,-4,-3,-4,-4,-3,-4,-4,-4,-4,-3,-3,-2,-2,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-2,-2,-3,-3,-4,-4,-4,-4,-3,-4,-4,-3,-4,-4,-4,-4,-3,-3,-2,-2,-3,-3,-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-4,-4,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-2,-2,-3,-3,-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-4,-4,-3,-3,-2,-2,-3,-3,-4,-4,-4,-4,-3,-4,-4,-3,-4,-4,-4,-4,-3,-3,-2,-2,-3,-3,-3,-3,-4,-4,-3,-4,-4,-3,-4,-4,-3,-3,-3,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-1,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1])
#
arrX=DataArrayDouble([0.,0.00672,0.01932,0.03192,0.04452,0.05712,0.06972,0.08232,0.09492,0.10752,0.12012,0.13272,0.14532,0.15792,0.17052,0.18312,0.19572,0.20832,0.21504])
arrY=DataArrayDouble([0.,0.00672,0.01932,0.03192,0.04452,0.05712,0.06972,0.08232,0.09492,0.10752,0.12012,0.13272,0.14532,0.15792,0.17052,0.18312,0.19572,0.20832,0.21504])
arrZ=DataArrayDouble([4.6025,4.695,4.805])
#
cm=MEDCouplingCMesh("Maillage_THYC") ; cm.setCoords(arrX,arrY,arrZ)
mm=MEDFileCMesh() ; mm.setMesh(cm)
mm.setFamilyFieldArr(0,famIds)
mm.setFamilyId("FAMILLE -1",-1)
mm.setFamilyId("FAMILLE -2",-2)
mm.setFamilyId("FAMILLE -3",-3)
mm.setFamilyId("FAMILLE -4",-4)
mm.write(fname,2)
f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(cm) ; f.setTime(0.,-1,0)
f.setName("POROSITE") ; f.setArray(p)
f1ts=MEDFileField1TS()
f1ts.setFieldNoProfileSBT(f)
f1ts.write(fname,0)
################### MED write is done -> Go to MEDReader
from paraview.simple import *

testMEDReader12_med = MEDReader( FileName='/export/home/geay/Salome7/V7_3_0/BUG_7972/TRY4/testMEDReader12.med' )
testMEDReader12_med.AllArrays = ['TS0/Maillage_THYC/ComSup0/POROSITE@@][@@P0']

RenderView1 = GetRenderView()
RenderView1.CenterOfRotation = [0.10751999914646149, 0.10751999914646149, 4.703749895095825]

DataRepresentation1 = Show()
DataRepresentation1.Representation = 'Outline'
DataRepresentation1.ScaleFactor = 0.021504000000000002
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation1.SelectionCellFieldDataArrayName = 'FamilyIdCell'

RenderView1.CameraPosition = [0.10751999914646149, 0.10751999914646149, 5.409578268564619]
RenderView1.CameraFocalPoint = [0.10751999914646149, 0.10751999914646149, 4.703749895095825]
RenderView1.CameraClippingRange = [0.4972827225809383, 0.9700468818440346]
RenderView1.CameraParallelScale = 0.18268182562745858

ExtractGroup1 = ExtractGroup(Input=testMEDReader12_med)
ExtractGroup1.UpdatePipelineInformation()
ExtractGroup1.AllGroups = ['FAM_FAMILLE -2@@][@@-2', 'FAM_FAMILLE -4@@][@@-4']

DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionCellFieldDataArrayName = 'FamilyIdCell'
DataRepresentation2.ScalarOpacityUnitDistance = 0.053416489858379865
DataRepresentation2.ExtractedBlockIndex = 1
DataRepresentation2.ScaleFactor = 0.021503999829292297

DataRepresentation1.Visibility = 0

a1_POROSITE_PVLookupTable = GetLookupTableForArray( "POROSITE", 1, RGBPoints=[0.4734426322438599, 0.23, 0.299, 0.754, 0.5274367112743427, 0.865, 0.865, 0.865, 0.5814307903048255, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

a1_POROSITE_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.4734426322438599, 0.0, 0.5, 0.0, 0.5814307903048255, 1.0, 0.5, 0.0] )

DataRepresentation2.ScalarOpacityFunction = a1_POROSITE_PiecewiseFunction
DataRepresentation2.ColorArrayName = ('CELL_DATA', 'POROSITE')
DataRepresentation2.LookupTable = a1_POROSITE_PVLookupTable
DataRepresentation2.ColorAttributeType = 'CELL_DATA'

Render()
RenderView1.ViewSize =[300,300]
WriteImage(outImgName)