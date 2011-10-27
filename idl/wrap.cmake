# Copyright (C) 2010-2011  CEA/DEN, EDF R&D
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

ADD_CUSTOM_COMMAND(
 OUTPUT vtkWrapIDL.h
 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/create_header.py ${CMAKE_BINARY_DIR}
 DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt
)

ADD_CUSTOM_COMMAND(
 OUTPUT hints
 COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/create_hints.py ${PARAVIEW_LIBRARY_DIRS}
 DEPENDS ${PARAVIEW_LIBRARY_DIRS}/hints ${CMAKE_CURRENT_SOURCE_DIR}/hints_paravis
)

SET(WRAP_IDL)
SET(WRAP_SK_FILES)

IF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)
 EXECUTE_PROCESS(
  COMMAND ${PYTHON_EXECUTABLE} -c "f = open('${CMAKE_BINARY_DIR}/wrapfiles.txt') ; print f.read(), ; f.close()"
  OUTPUT_VARIABLE WRAP_LIST_FULL
 )

 STRING(REGEX  MATCHALL "[^\n]+" WRAP_LIST_REG ${WRAP_LIST_FULL})
 FOREACH(STR ${WRAP_LIST_REG})

  SEPARATE_ARGUMENTS(STR)
  LIST(LENGTH STR WRAP_LEN)
  SET(DEP)
 
  LIST(GET STR 0 VAL)

  IF(WRAP_LEN GREATER 1)
   MATH(EXPR WRAP_LEN1 "${WRAP_LEN} - 1" )

   FOREACH(IND RANGE 1 ${WRAP_LEN1})
    LIST(GET STR ${IND} DEP_VAL)
    SET(DEP ${DEP} PARAVIS_Gen_${DEP_VAL}.idl)
   ENDFOREACH(IND RANGE 1 ${WRAP_LEN1})

  ENDIF(WRAP_LEN GREATER 1)

  SET(WRAP_IDL ${WRAP_IDL} PARAVIS_Gen_${VAL}.idl)
  SET(WRAP_SK_FILES ${WRAP_SK_FILES} PARAVIS_Gen_${VAL}SK.cc)
  SET(vtkWrapIDL_EXEFILE ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL)
  IF(WINDOWS)
    IF(CMAKE_BUILD_TOOL STREQUAL nmake)
      SET(vtkWrapIDL_EXEFILE ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL.exe)
    ELSE(CMAKE_BUILD_TOOL STREQUAL nmake)
      SET(vtkWrapIDL_EXEFILE ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_BUILD_TYPE}/vtkWrapIDL.exe)
    ENDIF(CMAKE_BUILD_TOOL STREQUAL nmake)
  ENDIF(WINDOWS)
  ADD_CUSTOM_COMMAND(
   OUTPUT PARAVIS_Gen_${VAL}.idl
   COMMAND ${vtkWrapIDL_EXEFILE} ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h hints 0 PARAVIS_Gen_${VAL}.idl
   DEPENDS vtkWrapIDL ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h hints ${DEP}
  )

 ENDFOREACH(STR ${WRAP_LIST_REG})
ENDIF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)

ADD_CUSTOM_TARGET(generate_txt DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt vtkWrapIDL.h hints)
ADD_CUSTOM_TARGET(generate_idl ALL DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt vtkWrapIDL.h hints ${WRAP_IDL})
