#ifndef __vtkMedInterpolation_h_
#define __vtkMedInterpolation_h_

#include "vtkObject.h"
#include "vtkMed.h"
#include "vtkMedSetGet.h"

class vtkMedString;
class vtkMedFraction;

class VTK_EXPORT vtkMedInterpolation : public vtkObject
{
public:
  static vtkMedInterpolation* New();
  vtkTypeRevisionMacro(vtkMedInterpolation, vtkObject)
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // name of the interpolation function
  vtkGetObjectMacro(Name, vtkMedString);

  // Description:
  // This is the iterator that should be used to read this interpolation
  // in the med file
  vtkSetMacro(MedIterator, med_int);
  vtkGetMacro(MedIterator, med_int);

	// Description:
	// Type geometrique des mailles
	vtkSetMacro(GeometryType, med_geometry_type);
	vtkGetMacro(GeometryType, med_geometry_type);

	// Description:
	// 1 if the basis functions are relative to the vertices of the cell.
	vtkSetMacro(IsCellNode, int);
	vtkGetMacro(IsCellNode, int);

  // Description:
  // Maximum degree of any coefficient of any basis function
  vtkSetMacro(MaximumDegree, int);
  vtkGetMacro(MaximumDegree, int);

  // Description:
  // Maximum number of coefficients for any basis function
  vtkSetMacro(MaximumNumberOfCoefficient, int);
  vtkGetMacro(MaximumNumberOfCoefficient, int);

  // Description:
  // Maximum number of coefficients for any basis function
  vtkSetMacro(NumberOfVariable, int);
  vtkGetMacro(NumberOfVariable, int);

  // Description:
  // The basis functions
  vtkGetObjectVectorMacro(BasisFunction, vtkMedFraction);
  vtkSetObjectVectorMacro(BasisFunction, vtkMedFraction);

protected :
	vtkMedInterpolation();
	~vtkMedInterpolation();

	med_int MedIterator;
	med_geometry_type GeometryType;
	int IsCellNode;
	int MaximumNumberOfCoefficient;
	int MaximumDegree;
	int NumberOfVariable;
	vtkMedString* Name;
	vtkObjectVector<vtkMedFraction>* BasisFunction;
};

#endif //__vtkMedInterpolation_h_
