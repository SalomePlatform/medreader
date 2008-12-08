// Copyright (C) 2005  OPEN CASCADE, CEA/DEN, EDF R&D, PRINCIPIA R&D
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either 
// version 2.1 of the License.
// 
// This library is distributed in the hope that it will be useful 
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
// Plot2d_ViewModel.cxx: implementation of the Plot2d_ViewModel class.

#include "PVGUI_ViewModel.h"
#include "PVGUI_ViewWindow.h"

/*!
  Constructor
*/
PVGUI_Viewer::PVGUI_Viewer()
:SUIT_ViewModel() 
{
}

/*!
  Destructor
*/
PVGUI_Viewer::~PVGUI_Viewer()
{
}

/*!
  Create new instance of view window on desktop \a theDesktop.
  \retval SUIT_ViewWindow* - created view window pointer.
*/
SUIT_ViewWindow* PVGUI_Viewer::createView(SUIT_Desktop* theDesktop)
{
  PVGUI_ViewWindow* aPVView = new PVGUI_ViewWindow(theDesktop, this);
  return aPVView;
}