SET(WRAP_IDL_I_HH)
SET(WRAP_IDL_I_CC)

IF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)
 EXECUTE_PROCESS(
  COMMAND cat ${CMAKE_BINARY_DIR}/wrapfiles.txt
  OUTPUT_VARIABLE WRAP_LIST_FULL
 )

 STRING(REGEX  MATCHALL "[^\n]+" WRAP_LIST_REG ${WRAP_LIST_FULL})
 FOREACH(STR ${WRAP_LIST_REG})

  SEPARATE_ARGUMENTS(STR)
  LIST(LENGTH STR WRAP_LEN)
  SET(DEP_HH)
  SET(DEP_CC)
 
  LIST(GET STR 0 VAL)

  IF(WRAP_LEN GREATER 1)
   MATH(EXPR WRAP_LEN1 "${WRAP_LEN} - 1" )

   FOREACH(IND RANGE 1 ${WRAP_LEN1})
    LIST(GET STR ${IND} DEP_VAL)
    SET(DEP_HH ${DEP_HH} PARAVIS_Gen_${DEP_VAL}_i.hh)
    SET(DEP_CC ${DEP_CC} PARAVIS_Gen_${DEP_VAL}_i.cc)
   ENDFOREACH(IND RANGE 1 ${WRAP_LEN1})

  ENDIF(WRAP_LEN GREATER 1)

  SET(WRAP_IDL_I_HH ${WRAP_IDL_I_HH} PARAVIS_Gen_${VAL}_i.hh)
  SET(WRAP_IDL_I_CC ${WRAP_IDL_I_CC} PARAVIS_Gen_${VAL}_i.cc)

  ADD_CUSTOM_COMMAND(
   OUTPUT PARAVIS_Gen_${VAL}_i.hh
   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_HH_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h ${CMAKE_BINARY_DIR}/idl/hints 0 PARAVIS_Gen_${VAL}_i.hh
   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_HH_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h ${CMAKE_BINARY_DIR}/idl/hints ${DEP_HH}
  ) 

  ADD_CUSTOM_COMMAND(
   OUTPUT PARAVIS_Gen_${VAL}_i.cc
   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_CC_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h ${CMAKE_BINARY_DIR}/idl/hints 0 PARAVIS_Gen_${VAL}_i.cc
   DEPENDS PARAVIS_Gen_${VAL}_i.hh ${CMAKE_CURRENT_BINARY_DIR}/vtkWrapIDL_CC_exe ${PARAVIEW_INCLUDE_DIRS}/${VAL}.h ${CMAKE_BINARY_DIR}/idl/hints ${DEP_CC}
  )

 ENDFOREACH(STR ${WRAP_LIST_REG})
ENDIF(EXISTS ${CMAKE_BINARY_DIR}/wrapfiles.txt)

ADD_CUSTOM_COMMAND(
 OUTPUT PARAVIS_CreateClass.cxx
 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/create_class.sh ${CMAKE_SOURCE_DIR}
 DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt ${WRAP_IDL_I_HH}
)
ADD_CUSTOM_TARGET(generate_pvgui ALL DEPENDS ${CMAKE_BINARY_DIR}/wrapfiles.txt PARAVIS_CreateClass.cxx ${WRAP_IDL_I_HH} ${WRAP_IDL_I_CC})
