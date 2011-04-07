#ifndef __vtkMedField_h_
#define __vtkMedField_h_

#include "vtkObject.h"
#include "vtkMedSetGet.h"
#include "vtkMed.h"

#include "vtkSmartPointer.h"

#include <set>

class vtkMedInterpolation;
class vtkMedFieldOverEntity;
class vtkMedString;
class vtkMedFieldStep;
class vtkMedComputeStep;
template <class T>
class vtkMedComputeStepMap;
class vtkMedFile;

class VTK_EXPORT vtkMedField: public vtkObject
{
public:
	static vtkMedField* New();
	vtkTypeRevisionMacro(vtkMedField, vtkObject);
	void PrintSelf(ostream& os, vtkIndent indent);

	// Description:
	// The number of component of this field
	virtual void	SetNumberOfComponent(int);
	vtkGetMacro(NumberOfComponent, int);

	// Description:
	// The type of data stored in this field
	vtkSetMacro(DataType, med_field_type);
	vtkGetMacro(DataType, med_field_type);

	// Description:
	// The name of this field
	vtkGetObjectMacro(Name, vtkMedString);

	// Description:
	// The name of this mesh this field is on
	vtkGetObjectMacro(MeshName, vtkMedString);

	// Description:
	// The name of this mesh this field is on
	vtkGetObjectMacro(TimeUnit, vtkMedString);

	// Description:
	// The units of each component of this field
	vtkGetObjectVectorMacro(Unit, vtkMedString);
	vtkSetObjectVectorMacro(Unit, vtkMedString);

	// Description:
	// The name of each component of this field
	vtkGetObjectVectorMacro(ComponentName, vtkMedString);
	vtkSetObjectVectorMacro(ComponentName, vtkMedString);

	// Description:
	// add a cell type as support to this field
	void	AddFieldStep(vtkMedFieldStep*);
	void	ClearFieldStep();
	vtkMedFieldStep* GetFieldStep(const vtkMedComputeStep&);
	vtkMedFieldStep* FindFieldStep(const vtkMedComputeStep&, int);
	med_int GetNumberOfFieldStep();
	vtkMedFieldStep* GetFieldStep(med_int);
	void  GatherFieldTimes(std::set<med_float>&);
	void  GatherFieldIterations(med_float,std::set<med_int>&);

	// Description:
	// returns if the field is on point, cell, quadrature point or elno
	//BTX
	enum {
		UnknownFieldType = 0x00,
		PointField = 0x01,
		CellField = 0x02,
		QuadratureField = 0x04,
		ElnoField = 0x08};
	//ETX
	//Description:
	// returns the type of field this is. The returned code is and OR between
	// the different possible types.
	vtkGetMacro(FieldType, int);

	// This computes the FieldType
	// (currently, it does it by looking only at the first compute step)
	virtual void	ComputeFieldType();

	// Description:
	// This returns true if the FieldType is composed of several types
	virtual int HasManyFieldTypes();

	// Description:
	// returns the first support type this field is on.
	virtual int GetFirstType();

	// Description:
	// This methods extracts from the other field all the fields that are
	// on the given support type and add them to the current field.
	// It also updates the other FieldType ivar.
	virtual void	ExtractFieldType(vtkMedField* otherfield, int type);

	// Description:
	// The index of this field in the med file
	vtkSetMacro(MedIterator, med_int);
	vtkGetMacro(MedIterator, med_int);

	// Description:
	// if the mesh is local or not.
	vtkSetMacro(Local, med_int);
	vtkGetMacro(Local, med_int);

	// Description:
	// The interpolation functions associated with this field
	vtkGetObjectVectorMacro(Interpolation, vtkMedInterpolation);
	vtkSetObjectVectorMacro(Interpolation, vtkMedInterpolation);

	// Description:
	// This stores the file this field is stored on.
	virtual void	SetParentFile(vtkMedFile*);
	vtkGetObjectMacro(ParentFile, vtkMedFile);

protected:
	vtkMedField();
	virtual ~vtkMedField();

	vtkSetMacro(FieldType, int);

	int NumberOfComponent;
	med_field_type DataType;
	med_int MedIterator;
	med_int Local;
	vtkMedString* Name;
	vtkMedString* MeshName;
	vtkMedString* TimeUnit;
	int FieldType;
	vtkMedFile* ParentFile;

	//BTX
	vtkMedComputeStepMap<vtkMedFieldStep>* FieldStep;

	vtkObjectVector<vtkMedString>* Unit;
	vtkObjectVector<vtkMedString>* ComponentName;
	vtkObjectVector<vtkMedInterpolation>* Interpolation;
	//ETX

private:
	vtkMedField(const vtkMedField&); // Not implemented.
	void operator=(const vtkMedField&); // Not implemented.

};

#endif //__vtkMedField_h_
