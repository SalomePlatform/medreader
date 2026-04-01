// Copyright (C) 2020-2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDReaderCommon.hxx"

#include <vtkMutableDirectedGraph.h>
#include <vtkStringArray.h>

namespace
{
  template< class FamGroupFeeder >
  void FeedSILForFamsAndGrps(FamGroupFeeder feeder, vtkMutableDirectedGraph* sil, vtkIdType root, vtkStringArray *names)
  {
    vtkIdType typeId0(sil->AddChild(root));
    names->InsertNextValue( feeder.getName() );
    //
    vtkIdType typeId1(sil->AddChild(typeId0));
    names->InsertNextValue(std::string(ROOT_OF_GRPS_IN_TREE));
    std::vector<std::string> grps(feeder.getGroupsNames());
    for( std::string it0 : grps )
      {
        vtkIdType typeId2(sil->AddChild(typeId1));
        names->InsertNextValue( it0 );
        std::vector<std::string> famsOnGrp( feeder.getFamiliesOnGroup( it0 ));
        for( std::string it1 : famsOnGrp )
          {
            sil->AddChild(typeId2);
            names->InsertNextValue( it1.c_str());
          }
      }
    //
    vtkIdType typeId11(sil->AddChild(typeId0));
    names->InsertNextValue(std::string(ROOT_OF_FAM_IDS_IN_TREE));
    std::vector<std::string> fams( feeder.getFamiliesNames() );
    for( std::string it00 : fams )
      {
        sil->AddChild(typeId11);
        int famId( feeder.getFamilyId( it00.c_str()) );
        names->InsertNextValue( ExtractGroupFam::FromFamInfoToEntry( it00, famId ) );
      }
  }
}
