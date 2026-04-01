// Copyright (C) 2022-2026  CEA, EDF
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

#include "vtkFileSeriesGroupReader.h"

#include <vtkCollection.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkTree.h>

#include <vtkTreeBFSIterator.h>
#include <vtkStringArray.h>
#include <vtkDataSetAttributes.h>

#include "vtkMEDReader.h"
#include "MEDReaderAndFileSeriesGateway.hxx"
#include "ExtractGroupHelper.h"
#include "ExtractGroupHelper.txx"

#include <vector>
#include <string>
#include <memory>

vtkStandardNewMacro(vtkFileSeriesGroupReader);

//=============================================================================
struct vtkFileSeriesGroupReaderInternals
{
  std::vector<std::string> FileNames;
  vtkNew<vtkCollection> ReaderCollection;
};

//=============================================================================
vtkFileSeriesGroupReader::vtkFileSeriesGroupReader()
    : Internals(new vtkFileSeriesGroupReaderInternals())
{
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//-----------------------------------------------------------------------------
vtkFileSeriesGroupReader::~vtkFileSeriesGroupReader() = default;

//----------------------------------------------------------------------------
void vtkFileSeriesGroupReader::AddFileName(const char *name)
{
  // Make sure the reader always has a filename set
  this->ReaderSetFileName(name);

  this->AddFileNameInternal(name);
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkFileSeriesGroupReader::RemoveAllFileNames()
{
  this->RemoveAllFileNamesInternal();
  this->Modified();
}

//----------------------------------------------------------------------------
void vtkFileSeriesGroupReader::RemoveAllFileNamesInternal()
{
  this->Internals->FileNames.clear();
}

//----------------------------------------------------------------------------
void vtkFileSeriesGroupReader::AddFileNameInternal(const char *name)
{
  this->Internals->FileNames.emplace_back(name);
}

//----------------------------------------------------------------------------
unsigned int vtkFileSeriesGroupReader::GetNumberOfFileNames()
{
  return static_cast<unsigned int>(this->Internals->FileNames.size());
}

//----------------------------------------------------------------------------
const char *vtkFileSeriesGroupReader::GetFileName(unsigned int idx)
{
  if (idx >= this->Internals->FileNames.size())
  {
    return nullptr;
  }
  return this->Internals->FileNames[idx].c_str();
}

//----------------------------------------------------------------------------
int vtkFileSeriesGroupReader::CanReadFile(const char *filename)
{
  if (!this->Reader)
  {
    return 0;
  }

  return this->ReaderCanReadFile(filename);
}

namespace
{
  void RetrieveProcInfo(unsigned int &iProc, unsigned int &nProc)
  {
    vtkMultiProcessController *vmpc(vtkMultiProcessController::GetGlobalController());
    iProc = vmpc ? vmpc->GetLocalProcessId() : 0;
    nProc = vmpc ? vmpc->GetNumberOfProcesses() : 1;
  }

  enum class FileSerieMode
  {
    OneFile_Sequential_And_Parallel,
    MultiFile_No_Aggregation,
    MultiFile_Aggregation
  };

  FileSerieMode GetMode(unsigned int nbOfFiles, unsigned int iProc, unsigned int nProc)
  {
    if (nbOfFiles == 1)
    {
      return FileSerieMode::OneFile_Sequential_And_Parallel;
    }
    else if (nbOfFiles <= nProc && iProc < nProc)
    {
      return FileSerieMode::MultiFile_No_Aggregation;
    }
    else
    {
      return FileSerieMode::MultiFile_Aggregation;
    }
  }

  void ComputeInCaseOfMultiFileWithAggregation(
    unsigned int nbOfFiles, unsigned int iProc, unsigned int nProc,
    unsigned int& nFiles, unsigned int& offFile, unsigned int& supFiles
  )
  {
    nFiles = nbOfFiles / nProc;
    offFile = iProc * nFiles;
    supFiles = nbOfFiles % nProc;
  }

  void DeepCopy( vtkMutableDirectedGraph *destSIL, vtkIdType distRootWhereToInsertLocalHierarchy, vtkTree *localTree, vtkStringArray *localVertexData, vtkIdType localVertexToCopyFrom )
  {
    vtkStringArray *destVertexData( vtkStringArray::SafeDownCast( destSIL->GetVertexData()->GetAbstractArray(0) ) );
    vtkIdType destVertex( destSIL->AddChild( distRootWhereToInsertLocalHierarchy ) );
    destVertexData->InsertNextValue( localVertexData->GetValue(localVertexToCopyFrom) );
    for (vtkIdType iter = 0; iter < localTree->GetNumberOfChildren(localVertexToCopyFrom) ; iter++ )
    {
      vtkIdType localVertexIter( localTree->GetChild(localVertexToCopyFrom, iter) );
      DeepCopy(destSIL, destVertex, localTree, localVertexData, localVertexIter);
    }
  }

  void FeedGrpsFamsInfoInSILOfLocal( vtkMutableDirectedGraph *destSIL, vtkIdType distRootWhereToInsertLocalHierarchy, vtkMutableDirectedGraph *localSIL )
  {
    vtkStringArray *destVertexData( vtkStringArray::SafeDownCast( destSIL->GetVertexData()->GetAbstractArray(0) ) );
    vtkStringArray *localVertexData( vtkStringArray::SafeDownCast( localSIL->GetVertexData()->GetAbstractArray(0) ) );
    vtkNew<vtkTree> localTree;
    localTree->CheckedShallowCopy( localSIL );
    vtkIdType localRoot( localTree->GetRoot() );
    vtkIdType mfg( localTree->GetChild(localRoot, 1) );// MeshesFamsGrps
    for (vtkIdType iterMesh = 0; iterMesh < localTree->GetNumberOfChildren(mfg) ; iterMesh++ )// iterate over meshes in MeshesFamsGrps
    {
      vtkIdType localMesh( localTree->GetChild(mfg, iterMesh) );// mesh
      DeepCopy( destSIL, distRootWhereToInsertLocalHierarchy, localTree, localVertexData, localMesh );
    }
  }

  void PrintSIL( vtkMutableDirectedGraph *sil, vtkStringArray *vertexData )
  {
    std::cout << "Printing............" << std::endl;
    vtkNew<vtkTree> tree;
    tree->CheckedShallowCopy( sil );
    vtkNew<vtkTreeBFSIterator> bfsIterator;
    bfsIterator->SetTree(tree);
    bfsIterator->SetStartVertex(tree->GetRoot());
    while (bfsIterator->HasNext())
    {
        vtkIdType vertex = bfsIterator->Next();
        vtkIdType lev( tree->GetLevel(vertex) );
        vtkStdString val( vertexData->GetValue( vertex ) );
        std::cout << " level = " << lev << " vertex = " << vertex << " " << val << std::endl; 
    }
    std::cout << "End printing............" << std::endl;
  }
  
  class InternalFamGroupFeeder
  {
  public:
    InternalFamGroupFeeder( ExtractGroupData&& egd ):_egd( std::move( egd ) )
    {
    }

    std::string getName() const
    {
      return _egd.getMeshName();
    }
    
    std::vector<std::string> getGroupsNames() const
    {
      return _egd.getGroupsNames();
    }
    
    std::vector<std::string> getFamiliesNames() const
    {
      return _egd.getFamiliesNames();
    }
    
    std::vector<std::string> getFamiliesOnGroup(const std::string& grp) const
    {
      return _egd.getFamiliesOnGroup( grp );
    }
    
    int getFamilyId( const std::string &famName ) const
    {
      return _egd.getFamilyId( famName );
    }
  private:
    ExtractGroupData _egd;
  };

  /*!
   * EDF34655 :  Method in charge to compute the fusion of groups / families over a collection of vtkMEDReader.
   */
  void FuseFamGrpsOverParts( vtkMutableDirectedGraph *destSIL, vtkCollection *collection )
  {
    vtkNew<vtkTree> destTree;
    destTree->CheckedShallowCopy( destSIL );
    vtkIdType destRoot( destTree->GetRoot() );
    vtkStringArray *distVertexData( vtkStringArray::SafeDownCast( destSIL->GetVertexData()->GetAbstractArray(0) ) );
    // Warning destSIL is the origin in this part of code
    vtkNew< vtkMutableDirectedGraph > newDestSIL;
    newDestSIL->Initialize();
    // CrossEdge is an edge linking hierarchies.
    vtkNew<vtkStringArray> crossEdgesArray;
    newDestSIL->GetEdgeData()->AddArray(crossEdgesArray);
    vtkNew<vtkStringArray> namesArray;
    namesArray->SetName("Names");
    newDestSIL->GetVertexData()->AddArray(namesArray);
    // Go Feeding
    // First compute the MESHES_FAMS_GROUPS tree coming from the fusion of all collection of MEDReaders
    std::vector< ExtractGroupData > egds;
    collection->InitTraversal();
    vtkObject *obj( nullptr );
    while ((obj = collection->GetNextItemAsObject()) != nullptr)
    {
      vtkMEDReader *localReader( vtkMEDReader::SafeDownCast( obj ) );
      vtkMutableDirectedGraph *localSIL( localReader->GetSIL() );
      vtkStringArray *localVertexData( vtkStringArray::SafeDownCast( localSIL->GetVertexData()->GetAbstractArray(0) ) );
      vtkIdType idMeshesFamsGrps( ExtractGroupData::LocateFamGroupsTreeTopVertexInSIL( localSIL ) );
      ExtractGroupData egd;
      egd.loadFromSingle( localSIL, localVertexData, idMeshesFamsGrps );
      egds.emplace_back( std::move( egd ) );
    }
    ExtractGroupData egd;
    egd.fuseForAggregated( egds );
    // Fusion of collection of fam/grps MEDReaders info done -> now feed newDestSIL
    {
      vtkIdType newDestRootId( newDestSIL->AddVertex() );
      namesArray->InsertNextValue( distVertexData->GetValue( destRoot ) );
      for( vtkIdType i = 0 ; i < destTree->GetNumberOfChildren( destRoot ) ; ++i )
      {
        if( i != MESHES_FAMS_GROUPS_ID ) // Copy all from destSIL except MESHES_FAMS_GROUPS part
        {
          vtkIdType fieldsStatusTreeDest( destTree->GetChild(destRoot, i) );
          DeepCopy( newDestSIL, newDestRootId, destTree, distVertexData, fieldsStatusTreeDest );
        }
        else
        {
          vtkIdType meshesFamsGrpsRoot(newDestSIL->AddChild(newDestRootId));
          namesArray->InsertNextValue( MESHES_FAMS_GROUPS );
          InternalFamGroupFeeder feeder( std::move( egd ) );
          FeedSILForFamsAndGrps<InternalFamGroupFeeder>( feeder, newDestSIL, meshesFamsGrpsRoot, namesArray);
        }
      }
    }
    // End of feeding -> patch the original
    destSIL->ShallowCopy( newDestSIL );
  }

  void FillGrpsFamsInfoInSIL( vtkMutableDirectedGraph *destSIL, vtkCollection *collection )
  {
    FuseFamGrpsOverParts(destSIL, collection);
    ///////
    vtkIdType fgmV( -1 );
    {
      vtkNew<vtkTree> destTree;
      destTree->CheckedShallowCopy( destSIL );
      vtkIdType destRoot( destTree->GetRoot() );
      fgmV = destTree->GetChild(destRoot, MESHES_FAMS_GROUPS_MULTI_ID);
    }
    //////
    collection->InitTraversal();
    vtkObject *obj( nullptr );
    int i(0);
    vtkStringArray *vertexData( vtkStringArray::SafeDownCast( destSIL->GetVertexData()->GetAbstractArray(0) ) );
    while ((obj = collection->GetNextItemAsObject()) != nullptr)
    {
      vtkMEDReader *localReader( vtkMEDReader::SafeDownCast( obj ) );
      vtkMutableDirectedGraph *localSIL( localReader->GetSIL() );
      vtkIdType rootOfLocalMeshGrpFams( destSIL->AddChild( fgmV) );
      std::string s( std::to_string( i ) );
      vertexData->InsertNextValue( s );
      FeedGrpsFamsInfoInSILOfLocal( destSIL, rootOfLocalMeshGrpFams, localSIL );
      i++;
    }
    //PrintSIL( destSIL, vertexData );
  }
}

//----------------------------------------------------------------------------
int vtkFileSeriesGroupReader::RequestInformation(
    vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  unsigned int iProc, nProc;
  RetrieveProcInfo(iProc, nProc);
  unsigned int nbFiles( this->GetNumberOfFileNames() );
  CallbackForSILManager cbmanager;
  if( GetMode(nbFiles, iProc, nProc) == FileSerieMode::MultiFile_Aggregation )
  {
    unsigned int nFiles, offFile, supFiles;
    ComputeInCaseOfMultiFileWithAggregation(
      nbFiles, iProc, nProc,
      nFiles, offFile, supFiles
    );
    vtkMEDReader *exposedReader( vtkMEDReader::SafeDownCast(this->Reader) );
    cbmanager.setReader( exposedReader );
    this->Internals->ReaderCollection->RemoveAllItems();
    for(unsigned int i = 0; i < nFiles; i++)
    {
      // Create as many MEDReader as we need to avoid deep copy
      vtkNew<vtkMEDReader> localReader;
      this->Internals->ReaderCollection->AddItem(localReader.Get());

      for (int iField = 0; iField < exposedReader->GetNumberOfFieldsTreeArrays(); iField++)
      {
        const char *name = exposedReader->GetFieldsTreeArrayName(iField);
        localReader->SetFieldsStatus(name, exposedReader->GetFieldsTreeArrayStatus(name));
      }
      for (int iTimes = 0; iTimes < exposedReader->GetNumberOfTimesFlagsArrays(); iTimes++)
      {
        const char *name = exposedReader->GetTimesFlagsArrayName(iTimes);
        localReader->SetTimesFlagsStatus(name, exposedReader->GetTimesFlagsArrayStatus(name));
      }
      localReader->GenerateVectors(exposedReader->GetGenerateVect());
      localReader->ChangeMode(exposedReader->GetIsStdOrMode());
      localReader->GhostCellGeneratorCallForPara(exposedReader->GetGCGCP());
      localReader->GetRidOffDebugArrays(exposedReader->GetRemoveDebugArrays());

      // Configure the localReader for usage with the files
      localReader->SetFileName(this->GetFileName(i + offFile));
      localReader->SetDistributeWithMPI(false);
      localReader->UpdateInformation();
    }
    cbmanager.setCallbackInfo( FillGrpsFamsInfoInSIL, this->Internals->ReaderCollection );
  }
  return this->Reader->ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int vtkFileSeriesGroupReader::RequestData(vtkInformation *request,
                                          vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
  auto output = vtkMultiBlockDataSet::GetData(outputVector, 0);
  unsigned int nbFiles( this->GetNumberOfFileNames() );
  
  vtkInformation *info = outputVector->GetInformationObject(0);
  double time = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  unsigned int iProc, nProc;
  RetrieveProcInfo(iProc, nProc);

  // Simple case, one file/bloc for n proc
  switch (GetMode(nbFiles, iProc, nProc))
  {
  case FileSerieMode::OneFile_Sequential_And_Parallel:
  {
    // Make sure the information is up to date
    this->ReaderSetFileName(this->GetFileName(0));

    // Update the reader
    this->Reader->ProcessRequest(request, inputVector, outputVector);

    // Copy the GAUSS_DATA info key
    vtkInformation *mInfo = this->Reader->GetOutputInformation(0);
    if (mInfo->Has(vtkMEDReader::GAUSS_DATA()))
    {
      info->CopyEntry(mInfo, vtkMEDReader::GAUSS_DATA());
    }
    break;
  }
  case FileSerieMode::MultiFile_No_Aggregation:
  { // N file/block read by m proc, with n <= m, means 0/1 file/block per proc
    // Distribute in MEDReader only when reading a single file in a single block
    vtkMEDReader::SafeDownCast(this->Reader)->SetDistributeWithMPI(false);

    // Needed as the MEDReader do not support changing its filename
    // without reloading everything.
    this->ReaderSetFileName(this->GetFileName(iProc));
    vtkMEDReader::SafeDownCast(this->Reader)->ReloadInternals();

    // Creating new output information.
    vtkNew<vtkMultiBlockDataSet> mb;
    vtkNew<vtkInformationVector> newOutInfoVec;
    newOutInfoVec->Copy(outputVector, true);
    vtkInformation *newOutInfo = newOutInfoVec->GetInformationObject(0);
    newOutInfo->Set(vtkDataObject::DATA_OBJECT(), mb);

    // Update the reader
    this->Reader->ProcessRequest(request, inputVector, newOutInfoVec);

    // Recover the output
    output->SetNumberOfBlocks( nProc );
    output->SetBlock(iProc, mb->GetBlock(0));

    // Copy the GAUSS_DATA info key
    vtkInformation *mInfo = this->Reader->GetOutputInformation(0);
    if (mInfo->Has(vtkMEDReader::GAUSS_DATA()))
    {
      info->CopyEntry(mInfo, vtkMEDReader::GAUSS_DATA());
    }
    break;
  }
  case FileSerieMode::MultiFile_Aggregation:
  {
    output->SetNumberOfBlocks( nbFiles );
    unsigned int nFiles, offFile, supFiles;
    ComputeInCaseOfMultiFileWithAggregation(
      nbFiles, iProc, nProc,
      nFiles, offFile, supFiles
    );

    // Last proc handle remaining files/block
    if (iProc + 1 == nProc)
    {
      nFiles += supFiles;
    }

    vtkMEDReader *exposedReader = vtkMEDReader::SafeDownCast(this->Reader);
    this->Internals->ReaderCollection->InitTraversal();
    vtkObject *obj( nullptr );
    int i(0);
    while ((obj = this->Internals->ReaderCollection->GetNextItemAsObject()) != nullptr)
    {
      vtkMEDReader *localReader( vtkMEDReader::SafeDownCast( obj ) );
      localReader->UpdateTimeStep(time);
      vtkDataObject *outputReader = vtkMultiBlockDataSet::SafeDownCast(localReader->GetOutputDataObject(0))->GetBlock(0);
      output->SetBlock(i + offFile, outputReader);

      if (i == 0)
      {
        // Copy the GAUSS_DATA info key of the first filename
        vtkInformation *mInfo = localReader->GetOutputInformation(0);
        if (mInfo->Has(vtkMEDReader::GAUSS_DATA()))
        {
          info->CopyEntry(mInfo, vtkMEDReader::GAUSS_DATA());
        }
      }
      i++;
    }
    break;
  }
  default:
  {
    vtkErrorMacro("vtkFileSeriesGroupReader::RequestData : fail to find a corresponding mode.");
    return 0;
  }
  }
  return 1;
}

//------------------------------------------------------------------------------
int vtkFileSeriesGroupReader::FillOutputPortInformation(
    int vtkNotUsed(port), vtkInformation *info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

//-----------------------------------------------------------------------------
void vtkFileSeriesGroupReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "MetaFileName: " << (this->_MetaFileName ? this->_MetaFileName : "(none)")
     << endl;
}
