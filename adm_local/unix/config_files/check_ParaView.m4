dnl  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
dnl  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS 
dnl 
dnl  This library is free software; you can redistribute it and/or 
dnl  modify it under the terms of the GNU Lesser General Public 
dnl  License as published by the Free Software Foundation; either 
dnl  version 2.1 of the License. 
dnl 
dnl  This library is distributed in the hope that it will be useful, 
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of 
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
dnl  Lesser General Public License for more details. 
dnl 
dnl  You should have received a copy of the GNU Lesser General Public 
dnl  License along with this library; if not, write to the Free Software 
dnl  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA 
dnl 
dnl  See http://www.opencascade.org/SALOME/ or email : webmaster.salome@opencascade.org 
dnl
dnl
dnl

dnl
dnl  OPTIONS_PARAVIEW
dnl  ------------------------------------------------------------------------
dnl  Adds the --with-paraview=path 
dnl
AC_DEFUN([OPTIONS_PARAVIEW], [
  AC_ARG_WITH([paraview],
              [AC_HELP_STRING([--with-paraview], [The prefix where ParaView is located (default "" means taking from PVHOME environment variable)])],
              [with_paraview=$withval], [with_paraview=""])
])

dnl
dnl  CHECK_PARAVIEW
dnl  ------------------------------------------------------------------------
dnl  a simplified check that assumes VTK is a part of ParaView, separate location 
dnl  for VTK is not supported
dnl
AC_DEFUN([CHECK_PARAVIEW], [
AC_REQUIRE([AC_PROG_CC])dnl
AC_REQUIRE([AC_PROG_CXX])dnl
AC_REQUIRE([AC_PROG_CPP])dnl
AC_REQUIRE([AC_PROG_CXXCPP])dnl
AC_REQUIRE([AC_LINKER_OPTIONS])dnl

AC_REQUIRE([OPTIONS_PARAVIEW])dnl

AC_CHECKING(for ParaView)

AC_LANG_SAVE
AC_LANG_CPLUSPLUS

AC_SUBST(PARAVIEW_INCLUDES)
AC_SUBST(PARAVIEW_LIBS)
dnl AC_SUBST(VTKPY_MODULES)

PARAVIEW_INCLUDES=""
PARAVIEW_LIBS=""
dnl VTKPY_MODULES=""

paraview_ok=no

if test "x$OpenGL_ok" != "xyes" ; then
   AC_MSG_WARN(ParaView needs OpenGL correct configuration, check configure output)
fi


if test "x$qt_ok" != "xyes" ; then
   AC_MSG_WARN(ParaView needs Qt correct configuration, check configure output)
fi

if test "x$hdf5_ok" != "xyes" ; then
   AC_MSG_WARN(ParaView needs correct HDF5 configuration, check configure output)
fi

AC_PATH_X
if test "x$x_libraries" != "x"
then
   LXLIB="-L$x_libraries"
else
   LXLIB=""
fi

if test "x$x_libraries" = "x/usr/lib"
then
   LXLIB=""
fi

if test "x$x_libraries" = "x/usr/lib${LIB_LOCATION_SUFFIX}"
then
   LXLIB=""
fi

LOCAL_INCLUDES="$OGL_INCLUDES"
LOCAL_LIBS="-lvtksys -lvtkzlib -lvtkpng -lvtkjpeg -lvtktiff -lvtkexpat -lvtksqlite -lvtkmetaio -lvtkverdict -lvtkNetCDF -lvtkDICOMParser -lvtkfreetype -lvtkftgl -lvtkexoIIc -lvtklibxml2 -lvtkCommon -lvtkGraphics -lvtkImaging -lvtkFiltering -lvtkIO -lvtkRendering -lvtkHybrid -lvtkParallel -lvtkWidgets -lvtkClientServer -lvtkCommonCS -lvtkFilteringCS -lvtkIOCS -lvtkImagingCS -lvtkGraphicsCS -lvtkGenericFiltering -lvtkGenericFilteringCS -lvtkRenderingCS -lvtkHybridCS -lvtkParallelCS -lvtkPVCommandOptions -lKWCommon -lvtkInfovis -lvtkInfovisCS -lvtkPVServerCommon -lvtkPVServerCommonCS -lvtkCommonPythonD -lvtkFilteringPythonD -lvtkIOPythonD -lXdmf -lvtkXdmf -lvtkXdmfCS -lvtkWidgetsCS -lvtkVolumeRendering -lvtkVolumeRenderingCS -lvtkPVFilters -lvtkPVFiltersCS -lvtkViews -lvtkPVServerManagerPythonD -lvtkPVServerCommonPythonD -lvtkPVPythonInterpretor -lvtkPVServerManager -lQVTK -lpqWidgets -lQtTesting -lvtkQtChart -lpqCore -lQtPython -lpqApplicationComponents -lpqComponents -lvtkalglib $LXLIB -lX11 -lXt"

dnl VTK install dir
if test -z $with_paraview ; then
  with_paraview=""
fi
if test "x$with_paraview" = "xyes" ; then
  dnl in case user wrote --with-paraview=yes
  with_paraview=""
fi
if test "x$with_paraview" = "xno" ; then
  dnl in case user wrote --with-paraview=no
  with_paraview=""
  AC_MSG_WARN(Value "no", specified for option --with-paraview, is not supported)
fi

if test "x$with_paraview" != "x" ; then
  PVHOME="$with_paraview"
fi

  if test -z $PVHOME ; then
    AC_MSG_WARN(undefined PVHOME and PVSRCHOME variables which specify where ParaView is located)
  fi  
  
LOCAL_INCLUDES="-I$PVHOME/include/paraview-$PVVERSION $LOCAL_INCLUDES"

LOCAL_LIBS="-L$PVHOME/lib/paraview-$PVVERSION $LOCAL_LIBS"

dnl vtk headers
CPPFLAGS_old="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $LOCAL_INCLUDES $QT_INCLUDES"

AC_CHECK_HEADER(pqApplicationCore.h,paraview_ok="yes",paraview_ok="no")

CPPFLAGS="$CPPFLAGS_old"

if test "x$paraview_ok" = "xyes"
then
   PARAVIEW_INCLUDES="$LOCAL_INCLUDES"

   dnl vtk libraries

   AC_MSG_CHECKING(linking ParaView component library)

   LIBS_old="$LIBS"
   LIBS="$LIBS $QT_LIBS $LOCAL_LIBS"
   CPPFLAGS_old="$CPPFLAGS"
   CPPFLAGS="$CPPFLAGS $QT_INCLUDES $PARAVIEW_INCLUDES"

   dnl VTKPY_MODULES="$VTKHOME/python"

   AC_CACHE_VAL(salome_cv_lib_paraview,[
     AC_TRY_LINK([#include "pqApplicationCore.h"
                 ],
		 [pqApplicationCore::instance()],
		 [salome_cv_lib_paraview=yes],
		 [salome_cv_lib_paraview=no])
   ])
   paraview_ok="$salome_cv_lib_paraview"
   LIBS="$LIBS_old"
   CPPFLAGS="$CPPFLAGS_old"
fi

if  test "x$paraview_ok" = "xno"
then
  AC_MSG_RESULT(no)
  AC_MSG_WARN(unable to link with ParaView component library)
else
  AC_MSG_RESULT(yes)
  PARAVIEW_LIBS="$LOCAL_LIBS"
fi

AC_MSG_RESULT(for ParaView: $paraview_ok)

AC_LANG_RESTORE

# Save cache
AC_CACHE_SAVE

])dnl
