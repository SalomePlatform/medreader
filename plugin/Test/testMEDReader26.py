#  -*- coding: utf-8 -*-
# Copyright (C) 2020-2026  CEA, EDF
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

__doc__ = """
EDF34655 : test management of pipeline medreader* -> extractgroup in multifiles context
"""

import medcoupling as mc
from MEDReaderHelper import WriteInTmpDir,RetriveBaseLine

from math import pi

Radius = 200.0
NbOfCellsOrthoRadPerProc = 10
NbCellsRad = 10
NbCellsAx = 10
NbOfProcs = 10 # must be even for good testing
# to be in [0., 1.]
RadiusThresForGroup = 0.3
grpName = "MyGroup"

def MyAssert(clue):
    if not clue:
        raise RuntimeError("Assertion failed !")
    
def CreateMEDFiles():

    def ManipulateGroups( mm, iProc ):
        """
        This method modifies intentionally families to be sure that extract group works correctly part by part.
        """
        global grpName
        MyAssert( len( mm.getFamiliesOnGroup( grpName ) ) == 1 )
        famName = mm.getFamiliesOnGroup( grpName )[0]
        famId = mm.getFamilyId( famName )
        tupleId = mm.getFamilyFieldAtLevel(0).findIdsEqual( famId )
        mm.setFamilyId(famName,famId - iProc)
        mm.getFamilyFieldAtLevel(0)[ tupleId ] -= iProc
        pass

    arrX = mc.DataArrayDouble( NbCellsRad + 1 )
    arrX.iota()
    arrX *= Radius / NbCellsRad
    arrY = mc.DataArrayDouble( NbOfCellsOrthoRadPerProc * NbOfProcs + 1 )
    arrY.iota()
    arrY = ( arrY * 2 * pi ) / ( NbOfCellsOrthoRadPerProc * NbOfProcs )
    m = mc.MEDCouplingCMesh()
    m.setCoords( arrX, arrY )
    m = m.buildUnstructured()
    m.setCoords( m.getCoords().fromPolarToCart() )
    m.changeSpaceDimension(3,0.)

    ax = mc.DataArrayDouble( NbCellsAx + 1 ) ; ax.iota()
    ax *= Radius / NbCellsAx
    ax = ax.changeNbOfComponents(3,0.)
    ax = ax[:,[1,2,0]]
    ax = mc.MEDCouplingUMesh.Build1DMeshFromCoords( ax )
    m = m.buildExtrudedMesh(ax,0)
    cc = m.computeIsoBarycenterOfNodesPerCell()
    cc = cc.fromCartToCyl()
    cc[cc[:,1].findIdsLowerThan(0)] += 2*pi

    ret = [ ]

    for iProc in range( NbOfProcs ):
        cids = cc[:,1].findIdsInRange( iProc * 2 * pi / NbOfProcs , ( iProc + 1) * 2 *pi / NbOfProcs )
        part = m[ cids ]
        part.setName("mesh")
        part.zipCoords()
        #
        mm = mc.MEDFileUMesh()
        mm[0] = part

        if iProc%2 == 0: # here important to force different groups on procs especially for group #NbOfProcs -1 which is the representant
            cellIdsInGroup = part.computeIsoBarycenterOfNodesPerCell().fromCartToCyl()[:,0].findIdsLowerThan( Radius * RadiusThresForGroup )
            cellIdsInGroup.setName( grpName )
            mm.addGroup(0,cellIdsInGroup)
            #
            ManipulateGroups( mm, iProc )
            #
        fname = f"mesh_{iProc}.med"
        ret.append( fname )
        mm.write( fname, 2 )
    
    return ret

@WriteInTmpDir
def test26():
    fileNames = CreateMEDFiles( )

    from MEDReaderVTK import MEDReaderIO
    import vtk
    file_names = vtk.vtkStringArray()

    procZero = 0
    fileName = fileNames[procZero]
    medreader = MEDReaderIO.vtkMEDReader()
    medreader.SetFileName( fileName )
    medreader.UpdateInformation()
    arrName = [ medreader.GetFieldsTreeArrayName(i) for i in range( medreader.GetNumberOfFieldsTreeArrays() ) ]
    medreader.SetFieldsStatus( arrName[0], True)
    medreader.Update()
    data = medreader.GetOutputDataObject(0)
    data_ug = data.GetBlock(0)
    MyAssert( data_ug.GetNumberOfCells() == NbCellsRad * NbOfCellsOrthoRadPerProc * NbCellsAx )

    medreader2 = MEDReaderIO.vtkFileSeriesGroupReader()
    medreaderReal = MEDReaderIO.vtkMEDReader()
    medreaderReal.SetFileName( fileName )
    medreader2.SetReader( medreaderReal )
    medreader2.SetFileNameMethod("SetFileName")
    for fileName in fileNames:
        file_names.InsertNextValue( fileName )
        medreader2.AddFileName( fileName )
        medreader2.CanReadFile( fileName )
    medreader2.Update()
    resu = medreader2.GetOutputDataObject(0)
    MyAssert( resu.GetNumberOfBlocks() == NbOfProcs )
    MyAssert( resu.GetNumberOfCells() == NbCellsRad * NbOfCellsOrthoRadPerProc * NbCellsAx * NbOfProcs )

    extractGroup = MEDReaderIO.vtkExtractGroup()
    extractGroup.SetInputConnection( medreader2.GetOutputPort(0) )
    extractGroup.UpdateInformation()
    groupsNames = [ extractGroup.GetGroupsFlagsArrayName(i) for i in range( extractGroup.GetNumberOfGroupsFlagsArrays() ) ]
    MyAssert( groupsNames == ["GRP_MyGroup"] )# <<< important test is here
    extractGroup.SetGroupsFlagsStatus(f"GRP_{grpName}",True)
    extractGroup.Update()
    resuExtractGrp = extractGroup.GetOutputDataObject(0)
    MyAssert( resuExtractGrp.GetNumberOfBlocks() == NbOfProcs // 2 )
    nbCellsExpected = int( RadiusThresForGroup * NbCellsRad * NbOfCellsOrthoRadPerProc * NbCellsAx ) * ( NbOfProcs // 2 )
    MyAssert( resuExtractGrp.GetNumberOfCells() == nbCellsExpected )

if __name__ == "__main__":
    test26()

