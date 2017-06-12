SET(VTK_LIBS vtkCommonExecutionModel vtkParallelCore)

vtk_module(vtkMEDReader
  DEPENDS
    ${VTK_LIBS}
  TEST_DEPENDS
    vtkRenderingCore
    vtkTestingRendering
    vtkInteractionStyle
  EXCLUDE_FROM_WRAP_HIERARCHY
)
