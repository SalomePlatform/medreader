// Copyright (C) 2010-2024  CEA, EDF
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
// Author : Anthony Geay

#ifndef vtkGenerateVectors_h__
#define vtkGenerateVectors_h__

#include "vtkDataSetAlgorithm.h"

class vtkFieldData;
class vtkDoubleArray;

class VTK_EXPORT vtkGenerateVectors
{
public:
  static void ChangeMeshTimeToUpdateCache(vtkDataSet *dataset);
  static void Operate(vtkFieldData *fd);
  static vtkDoubleArray *Operate2Compo(vtkDoubleArray *oldArr);
  static vtkDoubleArray *OperateMoreThan3Compo(vtkDoubleArray *oldArr);
  static std::string SuffixFieldName(const std::string& name);
public:
  static const char VECTOR_SUFFIX[];
protected:
  static void UpdateInformationOfArray(vtkDoubleArray *oldArr, vtkDoubleArray *arr);
};

#endif
