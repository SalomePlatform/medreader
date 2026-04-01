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

#pragma once

#include <vtkType.h>

#include <sstream>
#include <vector>
#include <set>
#include <map>

#include "MEDLoaderForPV.h"

class MEDLOADERFORPV_EXPORT ExtractGroupStatus
{
public:
  ExtractGroupStatus():_status(false) { }
  ExtractGroupStatus(const char *name);
  bool getStatus() const { return _status; }
  void setStatus(bool status) const { _status=status; }
  void cpyStatusFrom(const ExtractGroupStatus& other) { _status=other._status; }
  std::string getName() const { return _name; }
  void resetStatus() const { _status=false; }
  const char *getKeyOfEntry() const { return _ze_key_name.c_str(); }
  virtual void printMySelf(std::ostream& os) const;
  virtual bool isSameAs(const ExtractGroupStatus& other) const;
protected:
mutable bool _status;
std::string _name;
std::string _ze_key_name;
};

class MEDLOADERFORPV_EXPORT ExtractGroupGrp : public ExtractGroupStatus
{
public:
  ExtractGroupGrp(const char *name):ExtractGroupStatus(name) { std::ostringstream oss; oss << START << name; _ze_key_name=oss.str(); }
  void setFamilies(const std::vector<std::string>& fams) { _fams=fams; }
  const std::vector<std::string>& getFamiliesLyingOn() const { return _fams; }
  bool isSameAs(const ExtractGroupGrp& other) const;
public:
  static const char START[];
  std::vector<std::string> _fams;
};

class MEDLOADERFORPV_EXPORT ExtractGroupFam : public ExtractGroupStatus
{
public:
  ExtractGroupFam(const char *name);
  static std::string FromFamInfoToEntry( const std::string& famName, int famId );
  void printMySelf(std::ostream& os) const;
  void fillIdsToKeep(std::set<int>& s) const;
  int getId() const { return _id; }
  bool isSameAs(const ExtractGroupFam& other) const;
public:
  static const char START[];
private:
  int _id;
};

class vtkMutableDirectedGraph;
class vtkStringArray;

class MEDLOADERFORPV_EXPORT ExtractGroupData
{
public:
  static vtkIdType LocateFamGroupsTreeTopVertexInSIL(vtkMutableDirectedGraph *sil, vtkIdType *idMeshesFamsGrpsMulti = nullptr);
  void fuseForAggregated( const std::vector< ExtractGroupData >& egds );
  void loadFromSingle(vtkMutableDirectedGraph *sil, vtkStringArray *verticesNames, vtkIdType vertexIdOfMeshesFamsGrps);
  const char *getMeshNameYellowSubmarine() const;
  std::string getMeshName() const { return _mesh_name; }
  int getNumberOfEntries() const;
  const char *getKeyOfEntry(int i) const;
  const ExtractGroupStatus *getEntryNoThrow(const char *entry) const;
  const ExtractGroupStatus& getEntry(const char *entry) const;
  ExtractGroupStatus& getEntry(const char *entry);
  void printMySelf(std::ostream& os) const;
  int getIdOfFamily(const std::string& famName) const;
  std::set<int> getIdsHighlighted() const;
  std::vector< std::pair<std::string,std::vector<int> > > getAllGroups() const;
  std::vector<std::string> getGroupsNames() const;
  std::vector<std::string> getFamiliesNames() const;
  std::vector<std::string> getFamiliesOnGroup(const std::string& grpName) const;
  int getFamilyId( const std::string& famName ) const;
  void clearSelection() const;
private:
  std::map<std::string,int> computeFamStrIdMap() const;
private:
  std::vector<ExtractGroupGrp> _groups;
  std::vector<ExtractGroupFam> _fams;
  std::string _mesh_name;
};

class vtkInformationDataObjectMetaDataKey;
class vtkInformation;

class MEDLOADERFORPV_EXPORT ExtractGroupInternal
{
public:
  void loadFrom(vtkMutableDirectedGraph *sil);
  void loadFromSingle(vtkMutableDirectedGraph *sil, vtkStringArray *verticesNames, vtkIdType vertexIdOfMeshesFamsGrps);
  void loadFromMultiple(vtkMutableDirectedGraph *sil, vtkStringArray *verticesNames, vtkIdType idMeshesFamsGrpsMulti);
  int getNumberOfEntries() const { return _data.getNumberOfEntries(); }
  //! Warning : supposed to be const char * steady for paraview
  const char *getMeshNameYellowSubmarine() const;
  const char *getKeyOfEntry(int i) const { return _data.getKeyOfEntry(i); }
  bool getStatusOfEntryStr(const char *entry) const;
  void setStatusOfEntryStr(const char *entry, bool status);
  void printMySelf(std::ostream& os) const { _data.printMySelf(os); }
  bool isMultiMode( unsigned int nbBlock ) const;
  std::set<int> getIdsToKeep() const;
  std::set<int> getIdsToKeepMulti( unsigned int iBlock ) const;
  void clearSelection() const;
  int getIdOfFamily(const std::string& famName) const { return _data.getIdOfFamily( famName ); }
  std::vector< std::pair<std::string,std::vector<int> > > getAllGroups() const { return _data.getAllGroups(); }
  static bool IndependantIsInformationOK(vtkInformationDataObjectMetaDataKey *medReaderMetaData, vtkInformation *info);
private:
  ExtractGroupData _data;
  std::vector< ExtractGroupData > _data_multi;
  mutable std::vector< std::pair<std::string,bool> > _selection;
};
