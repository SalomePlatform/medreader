ADD_CUSTOM_COMMAND(
 OUTPUT vtkWrapIDL.h
 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_header.sh ${CMAKE_BINARY_DIR}
 DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt
)

ADD_CUSTOM_COMMAND(
 OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/vtkWrapIDL.c
 DEPENDS vtkWrapIDL.h
)

ADD_CUSTOM_COMMAND(
 OUTPUT hints
 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_hints.sh ${PARAVIEW_LIBRARY_DIRS}
 DEPENDS ${PARAVIEW_LIBRARY_DIRS}/hints ${CMAKE_CURRENT_SOURCE_DIR}/hints_paravis
)

SET(WRAP_IDL)
SET(WRAP_SK_FILES)

IF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)
 EXECUTE_PROCESS(
  COMMAND cat ${CMAKE_BINARY_DIR}/wrapfiles.txt
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
  ADD_CUSTOM_COMMAND(
   OUTPUT PARAVIS_Gen_${VAL}.idl
   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h hints 0 PARAVIS_Gen_${VAL}.idl
   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h hints ${DEP}
  )

 ENDFOREACH(STR ${WRAP_LIST_REG})
ENDIF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)

ADD_CUSTOM_TARGET(generate_idl ALL DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt vtkWrapIDL.h ${CMAKE_CURRENT_SOURCE_DIR}/vtkWrapIDL.c hints ${WRAP_IDL})
