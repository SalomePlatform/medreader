#include "vtkTableTo3D.h"

#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkVariantArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkTable.h"
#include "vtkInformation.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkWarpScalar.h"
#include "vtkContourFilter.h"

vtkStandardNewMacro(vtkTableTo3D);
vtkCxxRevisionMacro(vtkTableTo3D, "$Revision$");


vtkTableTo3D::vtkTableTo3D()
{
  this->ScaleFactor = 1.0;
  this->UseOptimusScale = true;
  this->PresentationType = TABLETO3D_SURFACE;
  this->NumberOfContours = 32;
}

vtkTableTo3D::~vtkTableTo3D()
{
}

int vtkTableTo3D::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkTable");
  return 1;
}

int vtkTableTo3D::RequestData(vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  vtkTable* input = vtkTable::GetData(inputVector[0], 0);
  vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

  if (input->GetNumberOfRows() == 0 ||input->GetNumberOfColumns() < 2)
    {
      return 1;
    }
  
  vtkIdType xSize = input->GetNumberOfRows();
  vtkIdType ySize = input->GetNumberOfColumns() - 1; 
  vtkIdType nbPoints = xSize * ySize;

  vtkDataArray* xAxis = vtkDataArray::SafeDownCast(input->GetColumn(0));
  if (!xAxis)
    {
      vtkErrorMacro("The first column is not numeric.");
      return 1;
    }
    
  double xRange = xAxis->GetTuple1(xSize - 1) - xAxis->GetTuple1(0);
  double yDelta = xRange / ySize;
  
  vtkSmartPointer<vtkDoubleArray> yAxis = 
    vtkSmartPointer<vtkDoubleArray>::New();
  yAxis->SetNumberOfValues(ySize);
  for (vtkIdType i = 0; i < ySize; i++ )
    {
      yAxis->SetValue(i, i*yDelta);
    }

  vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nbPoints);

  vtkSmartPointer<vtkIntArray> pointsIdMapper = 
    vtkSmartPointer<vtkIntArray>::New();
  pointsIdMapper->SetName("POINTS_ID_MAPPER");
  pointsIdMapper->SetNumberOfComponents(2);
  pointsIdMapper->SetNumberOfTuples(nbPoints);
  int *pointsIdMapperPtr = pointsIdMapper->GetPointer(0);

  for (vtkIdType i = 0, pntId = 0; i < ySize; i++) 
    {
      for (vtkIdType j = 0; j < xSize; j++, pntId++) 
	{
	  points->SetPoint(pntId, xAxis->GetTuple1(j), 
			          yAxis->GetValue(i), 
                                  0.0);

	  *pointsIdMapperPtr++ = pntId;
	  *pointsIdMapperPtr++ = 0;
	}
    }

  vtkSmartPointer<vtkDoubleArray> scalars = 
    vtkSmartPointer<vtkDoubleArray>::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(nbPoints);
  double *scalarsPtr = scalars->GetPointer(0);
  for (vtkIdType i = 0; i < ySize; i++) 
    {
      vtkDataArray* col = 
	vtkDataArray::SafeDownCast(input->GetColumn(i + 1));
      
      if (!col)
	{
	  vtkErrorMacro("Column "<< i <<"is not numeric.");
	  return 1;
	}
      
      for ( vtkIdType j = 0; j < xSize; j++ ) 
	{
	  double value = col->GetTuple1(j);
	  *scalarsPtr++ = value;
	}
    }

  vtkSmartPointer<vtkStructuredGrid> structuredGrid = 
    vtkSmartPointer<vtkStructuredGrid>::New();
  structuredGrid->SetPoints(points);

  structuredGrid->SetDimensions(xSize, ySize, 1);

  // structuredGrid->GetPointData()->AddArray(pointsIdMapper);
  if (input->GetInformation()->Has(vtkDataObject::FIELD_NAME()))
    {
      scalars->SetName(input->GetInformation()->Get(vtkDataObject::FIELD_NAME()));
    }
  else
    {
      scalars->SetName("Table");
    }
  structuredGrid->GetPointData()->SetScalars(scalars);

  vtkSmartPointer<vtkStructuredGridGeometryFilter> geomFilter = 
    vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
  geomFilter->SetInput(structuredGrid);
  geomFilter->Update();
  
  vtkSmartPointer<vtkWarpScalar> warpScalar = 
    vtkSmartPointer<vtkWarpScalar>::New();
  
  double scaleFactor = this->ScaleFactor;
  if (this->UseOptimusScale)
    {
      double range[2];
      geomFilter->GetOutput()->GetScalarRange(range);
      double length = geomFilter->GetOutput()->GetLength();
      if (range[1] > 0)
	{
	  scaleFactor = length / range[1] * 0.3;
	}
      else
	{
	  scaleFactor = 0;
	}
    }

  if (this->PresentationType == TABLETO3D_SURFACE)
    {
      warpScalar->SetInput(geomFilter->GetOutput());
      warpScalar->SetScaleFactor(scaleFactor);
    }
  else
    {
      vtkSmartPointer<vtkContourFilter> contourFilter = 
	vtkSmartPointer<vtkContourFilter>::New();
      contourFilter->SetInput(geomFilter->GetOutput());
      contourFilter->GenerateValues(this->NumberOfContours, 
				    geomFilter->GetOutput()->GetScalarRange());
      warpScalar->SetInput(contourFilter->GetOutput());
      warpScalar->SetScaleFactor(scaleFactor);
    }

  warpScalar->Update();
  output->ShallowCopy(warpScalar->GetPolyDataOutput());
  
  return 1;
}

//----------------------------------------------------------------------------
void vtkTableTo3D::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "ScaleFactor: " << this->ScaleFactor << endl;
  os << indent << "UseOptimusScale: "
     << (this->UseOptimusScale? "true" : "false") << endl;
  os << indent << "PresentationType: " 
     << ((this->PresentationType == TABLETO3D_SURFACE)? "Surface" : "Contour")
     << endl;
  os << indent << "NumberOfContours: " << this->NumberOfContours << endl;
}