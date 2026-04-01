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
// Author : Anthony Geay (EDF R&D)

#include "ExtractGroupHelper.h"
#include "MEDFileFieldRepresentationTree.hxx"
#include "MEDReaderCommon.hxx"

#include "vtkInformation.h"
#include "vtkInformationDataObjectMetaDataKey.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkDataSetAttributes.h"
#include "vtkStringArray.h"
#include "vtkTree.h"

#include <cstring>
#include <algorithm>

const char ExtractGroupGrp::START[]="GRP_";

const char ExtractGroupFam::START[]="FAM_";

ExtractGroupStatus::ExtractGroupStatus(const char *name):_status(false),_name(name)
{
}

void ExtractGroupStatus::printMySelf(std::ostream& os) const
{
  os << "      -" << _ze_key_name << "(";
  if(_status)
    os << "X";
  else
    os << " ";
  os << ")" << std::endl;
}

bool ExtractGroupStatus::isSameAs(const ExtractGroupStatus& other) const
{
  return _name==other._name && _ze_key_name==other._ze_key_name;
}

bool ExtractGroupGrp::isSameAs(const ExtractGroupGrp& other) const
{
  bool ret(ExtractGroupStatus::isSameAs(other));
  if(ret)
    return _fams==other._fams;
  else
    return false;
}

ExtractGroupFam::ExtractGroupFam(const char *name):ExtractGroupStatus(name),_id(0)
{
  std::size_t pos(_name.find(MEDFileFieldRepresentationLeavesArrays::ZE_SEP));
  std::string name0(_name.substr(0,pos)),name1(_name.substr(pos+strlen(MEDFileFieldRepresentationLeavesArrays::ZE_SEP)));
  std::istringstream iss(name1);
  iss >> _id;
  std::ostringstream oss; oss << START << name; _ze_key_name=oss.str(); _name=name0;
}

std::string ExtractGroupFam::FromFamInfoToEntry( const std::string& famName, int famId )
{
  std::ostringstream oss;
  oss << famName << MEDFileFieldRepresentationLeavesArrays::ZE_SEP << famId;
  return oss.str();
}

bool ExtractGroupFam::isSameAs(const ExtractGroupFam& other) const
{
  bool ret(ExtractGroupStatus::isSameAs(other));
  if(ret)
    return _id==other._id;
  else
    return false;
}

void ExtractGroupFam::printMySelf(std::ostream& os) const
{
  os << "      -" << _ze_key_name << " famName : \"" << _name << "\" id : " << _id << " (";
  if(_status)
    os << "X";
  else
    os << " ";
  os << ")" << std::endl;
}

void ExtractGroupFam::fillIdsToKeep(std::set<int>& s) const
{
  s.insert(_id);
}
///////////////////

bool ExtractGroupInternal::IndependantIsInformationOK(vtkInformationDataObjectMetaDataKey *medReaderMetaData, vtkInformation *info)
{
  // Check the information contain meta data key
  if(!info->Has(medReaderMetaData))
    return false;

  // Recover Meta Data
  vtkMutableDirectedGraph *sil(vtkMutableDirectedGraph::SafeDownCast(info->Get(medReaderMetaData)));
  if(!sil)
    return false;
  int idNames(0);
  vtkStringArray *verticesNames( vtkStringArray::SafeDownCast( sil->GetVertexData()->GetAbstractArray("Names",idNames) ) );
  if(!verticesNames)
    return false;
  for(int i=0;i<verticesNames->GetNumberOfValues();i++)
    {
      vtkStdString &st(verticesNames->GetValue(i));
      if( st == MESHES_FAMS_GROUPS )
        return true;
    }
  return false;
}

/*!
 * Take care of the returned pointer. Must be steady pointer.
 */
const char *ExtractGroupInternal::getMeshNameYellowSubmarine() const
{
  return this->_data.getMeshNameYellowSubmarine();
}

vtkIdType ExtractGroupData::LocateFamGroupsTreeTopVertexInSIL(vtkMutableDirectedGraph *sil, vtkIdType *idMeshesFamsGrpsMulti)
{
  int nameId( -1 );
  vtkStringArray *verticesNames( vtkStringArray::SafeDownCast(sil->GetVertexData()->GetAbstractArray("Names",nameId)) );
  vtkIdType idMeshesFamsGrps( -1 );
  bool found(false);
  for( int i = 0; i<verticesNames->GetNumberOfValues() ; i++ )
    {
      vtkStdString &st(verticesNames->GetValue(i));
      if( st == MESHES_FAMS_GROUPS )
        {
          idMeshesFamsGrps = i;
          found=true;
        }
      if( st == MESHES_FAMS_GROUPS_MULTI )
      {
        if( idMeshesFamsGrpsMulti )
        {
          *idMeshesFamsGrpsMulti = i;
        }
      }
    }
  if(!found)
    throw INTERP_KERNEL::Exception("There is an internal error ! The tree on server side has not the expected look !");
  return idMeshesFamsGrps;
}

void ExtractGroupData::fuseForAggregated( const std::vector< ExtractGroupData >& egds )
{
  if( egds.empty() )
  {
    throw INTERP_KERNEL::Exception("FuseForAggregated : empty input vector !");
  }
  const ExtractGroupData& ref( egds[0] );
  _mesh_name = ref.getMeshName();
  std::set<std::string> grpsOnAllInstances;
  for( const ExtractGroupData& egd : egds )
  {
    for( std::pair<std::string,std::vector<int> > grpData : egd.getAllGroups() )
    {
      grpsOnAllInstances.insert( grpData.first );
    }
  }
  _groups.clear();
  _fams.clear();
  for( std::string grp : grpsOnAllInstances )
  {
    _groups.emplace_back( ExtractGroupGrp( grp.c_str() ) );
  }
}

/*!
 * Load tree data in META_DATA information key under MeshesFamsGrps.
 */
void ExtractGroupData::loadFromSingle(vtkMutableDirectedGraph *sil, vtkStringArray *verticesNames, vtkIdType vertexIdOfMeshesFamsGrps)
{
  std::vector<ExtractGroupGrp> oldGrps(_groups); _groups.clear();
  std::vector<ExtractGroupFam> oldFams(_fams); _fams.clear();
  vtkNew<vtkAdjacentVertexIterator> it0;
  sil->GetAdjacentVertices(vertexIdOfMeshesFamsGrps,it0);
  while(it0->HasNext())
    {
      vtkIdType id1(it0->Next());
      std::string meshName(verticesNames->GetValue(id1).c_str());
      this->_mesh_name=meshName;
      vtkNew<vtkAdjacentVertexIterator> it1;
      sil->GetAdjacentVertices(id1,it1);
      vtkIdType idZeGrps(it1->Next());//zeGroups
      vtkNew<vtkAdjacentVertexIterator> itGrps;
      sil->GetAdjacentVertices(idZeGrps,itGrps);
      while(itGrps->HasNext())
        {
          vtkIdType idg(itGrps->Next());
          ExtractGroupGrp grp(verticesNames->GetValue(idg).c_str());
          vtkNew<vtkAdjacentVertexIterator> itGrps2;
          sil->GetAdjacentVertices(idg,itGrps2);
          std::vector<std::string> famsOnGroup;
          while(itGrps2->HasNext())
            {
              vtkIdType idgf(itGrps2->Next());
              famsOnGroup.push_back(std::string(verticesNames->GetValue(idgf)).c_str());
            }
          grp.setFamilies(famsOnGroup);
          _groups.push_back(grp);
        }
      vtkIdType idZeFams(it1->Next());//zeFams
      vtkNew<vtkAdjacentVertexIterator> itFams;
      sil->GetAdjacentVertices(idZeFams,itFams);
      while(itFams->HasNext())
        {
          vtkIdType idf(itFams->Next());
          ExtractGroupFam fam(verticesNames->GetValue(idf).c_str());
          _fams.push_back(fam);
        }
    }
  //
  std::size_t szg(_groups.size()),szf(_fams.size());
  if(szg==oldGrps.size() && szf==oldFams.size())
    {
      bool isSame(true);
      for(std::size_t i=0;i<szg && isSame;i++)
        isSame=_groups[i].isSameAs(oldGrps[i]);
      for(std::size_t i=0;i<szf && isSame;i++)
        isSame=_fams[i].isSameAs(oldFams[i]);
      if(isSame)
        {
          for(std::size_t i=0;i<szg;i++)
            _groups[i].cpyStatusFrom(oldGrps[i]);
          for(std::size_t i=0;i<szf;i++)
            _fams[i].cpyStatusFrom(oldFams[i]);
        }
    }
}

/*!
 * Method exposable to ParaView.
 */
const char *ExtractGroupData::getMeshNameYellowSubmarine() const
{
  return _mesh_name.c_str();
}

int ExtractGroupData::getNumberOfEntries() const
{
  std::size_t sz0(_groups.size()),sz1(_fams.size());
  return (int)(sz0+sz1);
}

const char *ExtractGroupData::getKeyOfEntry(int i) const
{
  int sz0((int)_groups.size());
  if(i>=0 && i<sz0)
    return _groups[i].getKeyOfEntry();
  else
    return _fams[i-sz0].getKeyOfEntry();
}

const ExtractGroupStatus& ExtractGroupData::getEntry(const char *entry) const
{
  std::string entryCpp(entry);
  for(std::vector<ExtractGroupGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  for(std::vector<ExtractGroupFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  std::ostringstream oss; oss << "ExtractGroupData::getEntry : no such entry \"" << entry << "\"!";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

const ExtractGroupStatus *ExtractGroupData::getEntryNoThrow(const char *entry) const
{
  std::string entryCpp(entry);
  for(std::vector<ExtractGroupGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return &(*it0);
  for(std::vector<ExtractGroupFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return &(*it0);
  return nullptr;
}

ExtractGroupStatus& ExtractGroupData::getEntry(const char *entry)
{
  std::string entryCpp(entry);
  for(std::vector<ExtractGroupGrp>::iterator it0=_groups.begin();it0!=_groups.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  for(std::vector<ExtractGroupFam>::iterator it0=_fams.begin();it0!=_fams.end();it0++)
    if(entryCpp==(*it0).getKeyOfEntry())
      return *it0;
  std::ostringstream oss; oss << "ExtractGroupData::getEntry : no such entry \"" << entry << "\"!";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

void ExtractGroupData::printMySelf(std::ostream& os) const
{
  os << "Groups :" << std::endl;
  for(std::vector<ExtractGroupGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    (*it0).printMySelf(os);
  os << "Families :" << std::endl;
  for(std::vector<ExtractGroupFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    (*it0).printMySelf(os);
}

int ExtractGroupData::getIdOfFamily(const std::string& famName) const
{
  for(std::vector<ExtractGroupFam>::const_iterator it=_fams.begin();it!=_fams.end();it++)
    {
      if((*it).getName()==famName)
        return (*it).getId();
    }
  return std::numeric_limits<int>::max();
}

std::map<std::string,int> ExtractGroupData::computeFamStrIdMap() const
{
  std::map<std::string,int> ret;
  for(std::vector<ExtractGroupFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    ret[(*it0).getName()]=(*it0).getId();
  return ret;
}

std::set<int> ExtractGroupData::getIdsHighlighted() const
{
  std::map<std::string,int> m(this->computeFamStrIdMap());
  std::set<int> s;
  for(std::vector<ExtractGroupGrp>::const_iterator it0=_groups.begin();it0!=_groups.end();it0++)
    {
      if((*it0).getStatus())
        {
          const std::vector<std::string>& fams((*it0).getFamiliesLyingOn());
          for(std::vector<std::string>::const_iterator it1=fams.begin();it1!=fams.end();it1++)
            {
              std::map<std::string,int>::iterator it2(m.find((*it1)));
              if(it2!=m.end())
                s.insert((*it2).second);
            }
        }
     }
  for(std::vector<ExtractGroupFam>::const_iterator it0=_fams.begin();it0!=_fams.end();it0++)
    if((*it0).getStatus())
      (*it0).fillIdsToKeep(s);
  return s;
}

void ExtractGroupData::clearSelection() const
{
  for(auto it : _groups)
    it.resetStatus();
  for(auto it : _fams)
    it.resetStatus();
}

// see reference : https://en.cppreference.com/w/cpp/iterator/iterator
class FamilyIterator : public std::iterator< std::input_iterator_tag, long, long, int*, int >
{
  long _num = 0;
  const ExtractGroupData *_egi = nullptr;
  const std::vector<std::string> *_fams = nullptr;
public:
  explicit FamilyIterator(long num , const ExtractGroupData *egi, const std::vector<std::string>& fams) : _num(num),_egi(egi),_fams(&fams) {}
  FamilyIterator& operator++() { ++_num; return *this;}
  bool operator==(const FamilyIterator& other) const {return _num == other._num;}
  bool operator!=(const FamilyIterator& other) const {return !(*this == other);}
  reference operator*() const {return _egi->getIdOfFamily((*_fams)[_num]);}
};

std::vector< std::pair<std::string,std::vector<int> > > ExtractGroupData::getAllGroups() const
{
  std::vector< std::pair<std::string,std::vector<int> > > ret;
  for(const auto&  grp : _groups)
  {
      const std::vector<std::string>& fams(grp.getFamiliesLyingOn());
      std::vector<int> famIds(FamilyIterator(0,this,fams),FamilyIterator(fams.size(),this,fams));
      std::pair<std::string,std::vector<int> > elt(grp.getName(),std::move(famIds));
      ret.emplace_back(std::move(elt));
  }
  return ret;
}

std::vector<std::string> ExtractGroupData::getGroupsNames() const
{
  std::vector<std::string> ret;
  std::for_each( _groups.cbegin(), _groups.cend(), [&ret]( const ExtractGroupGrp& grp )
  {
    ret.push_back( grp.getName() );
  }
   );
  return ret;
}

std::vector<std::string> ExtractGroupData::getFamiliesNames() const
{
  std::vector<std::string> ret;
  std::for_each( _fams.cbegin(), _fams.cend(), [&ret]( const ExtractGroupFam& fam )
  {
    ret.push_back( fam.getName() );
  }
   );
  return ret;
}

std::vector<std::string> ExtractGroupData::getFamiliesOnGroup(const std::string& grpName) const
{
  for( const ExtractGroupGrp& grp : _groups )
  {
    if( grp.getName() == grpName )
    {
      return grp.getFamiliesLyingOn();
    }
  }
  THROW_IK_EXCEPTION( "ExtractGroupData::getFamiliesOnGroup : Impossible to locate group \"" << grpName << "\"" );
}

int ExtractGroupData::getFamilyId( const std::string& famName ) const
{
  for( const ExtractGroupFam& fam : _fams )
  {
    if( fam.getName() == famName )
    {
      return fam.getId();
    }
  }
  THROW_IK_EXCEPTION( "ExtractGroupData::getFamiliesOnGroup : Impossible to locate family \"" << famName << "\"" );
}

void ExtractGroupInternal::loadFromMultiple(vtkMutableDirectedGraph *sil, vtkStringArray *verticesNames, vtkIdType idMeshesFamsGrpsMulti)
{
  vtkNew<vtkTree> tree;
  tree->CheckedShallowCopy( sil );
  _data_multi.resize( tree->GetNumberOfChildren( idMeshesFamsGrpsMulti ) );
  for (vtkIdType iFile = 0; iFile < tree->GetNumberOfChildren( idMeshesFamsGrpsMulti ) ; iFile++ )
  {
    vtkIdType iFileEntry( tree->GetChild(idMeshesFamsGrpsMulti, iFile) );
    int iFileExpected( std::stoi( verticesNames->GetValue(iFileEntry) ) );
    if( iFile != iFileExpected )
    {
      throw INTERP_KERNEL::Exception("There is an internal error inside SIL subtree of aggregated files.");
    }
    _data_multi[ iFile ].loadFromSingle( sil, verticesNames, iFileEntry);
  }
}

void ExtractGroupInternal::loadFrom(vtkMutableDirectedGraph *sil)
{
  int nameId( -1 );
  vtkStringArray *verticesNames( vtkStringArray::SafeDownCast(sil->GetVertexData()->GetAbstractArray("Names",nameId)) );
  vtkIdType idMeshesFamsGrpsMulti( -1 );
  vtkIdType idMeshesFamsGrps( ExtractGroupData::LocateFamGroupsTreeTopVertexInSIL( sil, &idMeshesFamsGrpsMulti ) );
  _data.loadFromSingle( sil, verticesNames, idMeshesFamsGrps );
  if( idMeshesFamsGrpsMulti != -1 )
  {
    loadFromMultiple( sil, verticesNames, idMeshesFamsGrpsMulti );
  }
}

bool ExtractGroupInternal::getStatusOfEntryStr(const char *entry) const
{
  const ExtractGroupStatus& elt(_data.getEntry(entry));
  return elt.getStatus();
}

void ExtractGroupInternal::setStatusOfEntryStr(const char *entry, bool status)
{
  _selection.emplace_back(entry,status);
}

std::set<int> ExtractGroupInternal::getIdsToKeep() const
{
  for(auto it: _selection)
    {
      const ExtractGroupStatus& elt(_data.getEntry(it.first.c_str()));
      elt.setStatus(it.second);
    }
  return _data.getIdsHighlighted();
}

std::set<int> ExtractGroupInternal::getIdsToKeepMulti( unsigned int iBlock ) const
{
  if( iBlock >= ( unsigned int )_data_multi.size() )
  {
    THROW_IK_EXCEPTION( "vtkExtractGroup::RequestData : in multi mode invalid block " << iBlock );
  }
  const ExtractGroupData& cur( _data_multi[iBlock] );
  for(auto it: _selection)
    {
      if( it.second )
      {
        const ExtractGroupStatus *elt( cur.getEntryNoThrow( it.first.c_str() ) );
        if( elt )
        {
          elt->setStatus( it.second );
        }
      }
    }
  return _data_multi[iBlock].getIdsHighlighted();
}

bool ExtractGroupInternal::isMultiMode( unsigned int nbBlock ) const
{
  if( _data_multi.empty() )
    return false;
  return (unsigned int) _data_multi.size() == nbBlock;
}

void ExtractGroupInternal::clearSelection() const
{
  _selection.clear();
  _data.clearSelection();
}
