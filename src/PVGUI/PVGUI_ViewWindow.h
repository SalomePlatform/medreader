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
// File   : Plot2d_ViewWindow.h
// Author : Vadim SANDLER, Open CASCADE S.A.S. (vadim.sandler@opencascade.com)
//

#ifndef PVGUI_VIEWWINDOW_H
#define PVGUI_VIEWWINDOW_H

#include <SUIT_ViewWindow.h>
#include <QMap>

class SUIT_Desktop;
class PVGUI_Viewer;
class pqViewManager;

class PVGUI_ViewWindow : public SUIT_ViewWindow  
{
  Q_OBJECT

public:
  PVGUI_ViewWindow( SUIT_Desktop*, PVGUI_Viewer* );
  virtual ~PVGUI_ViewWindow();

  virtual QString   getVisualParameters();
  virtual void      setVisualParameters( const QString& );
  
  void              setMultiViewManager( pqViewManager* );
  pqViewManager*    getMultiViewManager() const;

private:
  PVGUI_Viewer*     myModel;
  pqViewManager*    myPVMgr;
};

#endif // PLOT2D_VIEWWINDOW_H