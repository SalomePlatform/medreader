# - Find MED library
# Find the MED includes and library
# This module defines
#  MED2HOME, a directory where MED was installed. This directory is used to help find trhe other values.
#  MED_INCLUDE_DIR, where to find med.h
#  MED_INCLUDE_DIRS, where to find med.h file, concatenated with other include dirs from HDF5 and MPI (if parallel)
#  MED_LIBRARIES, libraries to link against to use MED. (including HDF5 and MPI if parallel)
#  MED_FOUND, If false, do not try to use MED.
# also defined, but not for general use are
#  MED_LIBRARY, the med library.
#  MEDC_LIBRARY, the medC library

SET(MED2HOME $ENV{MED2HOME} CACHE PATH "Path to the med install dir")

FIND_PATH(MED_INCLUDE_DIR med.h
  HINTS
  ${MED2HOME}/include
  PATHS
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(MED_LIBRARY med
  HINTS
  ${MED_INCLUDE_DIR}/../lib
  ${MED2HOME}/lib
  PATHS
  /usr/local/lib
  /usr/lib
)

get_filename_component(MED_LIBRARY_DIR ${MED_LIBRARY} PATH)

FIND_LIBRARY(MEDC_LIBRARY medC
  HINTS
  ${MED_LIBRARY_DIR}
  ${MED2HOME}/lib
  PATHS
  /usr/local/lib
  /usr/lib
)

IF(MED_INCLUDE_DIR)
  IF(MED_LIBRARY)
    IF(MEDC_LIBRARY)
      SET(MED_LIBRARIES ${MED_LIBRARY} ${MEDC_LIBRARY} )
      SET( MED_FOUND "YES" )
    ENDIF(MEDC_LIBRARY)
  ENDIF(MED_LIBRARY)
ENDIF(MED_INCLUDE_DIR)

IF(${MED_FOUND})
  FIND_PACKAGE(HDF5 REQUIRED)
  SET(MED_LIBRARIES ${MED_LIBRARIES} ${HDF5_LIBRARIES})
  SET(MED_INCLUDE_DIRS ${MED_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS})
  IF(${HDF5_IS_PARALLEL})
    FIND_PACKAGE(MPI REQUIRED)
    SET(MED_LIBRARIES ${MED_LIBRARIES} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    SET(MED_INCLUDE_DIRS ${MED_INCLUDE_DIRS} ${MPI_INCLUDE_PATH})
  ENDIF(${HDF5_IS_PARALLEL})
ENDIF(${MED_FOUND})

SET(MED_INCLUDE_DIR ${MED_INCLUDE_DIRS})
