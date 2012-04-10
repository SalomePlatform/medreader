// PARAVIS : ParaView wrapper SALOME module
//
// Copyright (C) 2010-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
// File   : PVGUI_Module_MenuActions.cxx
// Author : Margarita KARPUNINA
//

#include "PVGUI_Module.h"

#include <QtxActionToolMgr.h>
#include <LightApp_Application.h>
#include <SUIT_Desktop.h>

#include <QAction>
#include <QDockWidget>
#include <QToolBar>
#include <QStatusBar>
#include <QShortcut>
#include <QScrollArea>
#include <QVBoxLayout>

#include <pqAnimationViewWidget.h> 

#include <pqApplicationCore.h>
#include <pqComparativeVisPanel.h>
#include <pqObjectInspectorWidget.h>
#include <pqPipelineBrowserWidget.h>
//#include <pqProxyTabWidget.h>
#include <pqObjectInspectorWidget.h>
#include <pqProxyInformationWidget.h>
#include <pqDisplayProxyEditorWidget.h>
#include <pqSettings.h>
#include <pqDataInformationWidget.h>
#include <pqPVAnimationWidget.h>
#include <pqSelectionInspectorWidget.h>
#include <pqProgressWidget.h>

#include <pqAlwaysConnectedBehavior.h>
#include <pqApplicationCore.h>
#include <pqAutoLoadPluginXMLBehavior.h>
#include <pqCommandLineOptionsBehavior.h>
#include <pqCrashRecoveryBehavior.h>
#include <pqDataTimeStepBehavior.h>
#include <pqDefaultViewBehavior.h>
#include <pqDeleteBehavior.h>
#include <pqPersistentMainWindowStateBehavior.h>
#include <pqPluginActionGroupBehavior.h>
#include <pqPluginDockWidgetsBehavior.h>
#include <pqPluginManager.h>
#include <pqPVNewSourceBehavior.h>
#include <pqSpreadSheetVisibilityBehavior.h>
#include <pqStandardViewModules.h>
#include <pqUndoRedoBehavior.h>
#include <pqViewFrameActionsBehavior.h>
#include <pqParaViewMenuBuilders.h>
#include <pqCollaborationPanel.h>
#include <pqMemoryInspector.h>

/*!
  \brief Create dock widgets for ParaView widgets such as object inspector, pipeline browser, etc.
  ParaView pqMainWIndowCore class is fully responsible for these dock widgets' contents.
*/
void PVGUI_Module::setupDockWidgets()
{
  SUIT_Desktop* desk = application()->desktop();
 
  desk->setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
  desk->setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);

  // Pipeline
  QDockWidget* pipelineBrowserDock = new QDockWidget( tr( "TTL_PIPELINE_BROWSER" ), desk );
  pipelineBrowserDock->setObjectName("pipelineBrowserDock");
  pipelineBrowserDock->setAllowedAreas( Qt::LeftDockWidgetArea|Qt::NoDockWidgetArea|Qt::RightDockWidgetArea );
  desk->addDockWidget( Qt::LeftDockWidgetArea, pipelineBrowserDock );
  pqPipelineBrowserWidget* browser = new pqPipelineBrowserWidget(pipelineBrowserDock);
  pqParaViewMenuBuilders::buildPipelineBrowserContextMenu(*browser);
  pipelineBrowserDock->setWidget(browser);
  myDockWidgets[pipelineBrowserDock] = true;


  

  //Object inspector
  QDockWidget* objectInspectorDock = new QDockWidget( tr( "TTL_OBJECT_INSPECTOR" ), desk );
  objectInspectorDock->setObjectName("objectInspectorDock");
  objectInspectorDock->setAllowedAreas( Qt::LeftDockWidgetArea|Qt::NoDockWidgetArea|Qt::RightDockWidgetArea );
  desk->addDockWidget( Qt::LeftDockWidgetArea, objectInspectorDock );

  pqObjectInspectorWidget* objectInspectorWidget = new pqObjectInspectorWidget(objectInspectorDock);
  objectInspectorDock->setObjectName("objectInspectorWidget");
  objectInspectorWidget->setShowOnAccept(true);
  objectInspectorDock->setWidget(objectInspectorWidget);
  connect( objectInspectorWidget, SIGNAL( helpRequested(const QString&, const QString&) ),  this, SLOT( showHelpForProxy(const QString&, const QString&) ) );
  myDockWidgets[objectInspectorDock] = true;

  //Display Dock
  QDockWidget* displayDock = new QDockWidget( tr( "TTL_DISPLAY" ), desk );
  displayDock->setObjectName("displayDock");
  QWidget* displayWidgetFrame = new QWidget(displayDock);
  displayWidgetFrame->setObjectName("displayWidgetFrame");
  displayDock->setWidget(displayWidgetFrame);

  QScrollArea* displayScrollArea = new QScrollArea(displayWidgetFrame);
  displayScrollArea->setObjectName("displayScrollArea");
  displayScrollArea->setWidgetResizable(true);

  QVBoxLayout* verticalLayout = new QVBoxLayout(displayWidgetFrame);
  verticalLayout->setSpacing(0);
  verticalLayout->setContentsMargins(0, 0, 0, 0);

  pqDisplayProxyEditorWidget* displayWidget = new pqDisplayProxyEditorWidget(displayDock);
  displayWidget->setObjectName("displayWidget");
  displayScrollArea->setWidget(displayWidget);
  verticalLayout->addWidget(displayScrollArea);

  myDockWidgets[displayDock] = true;

  // information dock
  QDockWidget* informationDock = new QDockWidget(tr( "TTL_INFORMATION" ), desk);
  informationDock->setObjectName("informationDock");

  QWidget* informationWidgetFrame = new QWidget(informationDock);
  informationWidgetFrame->setObjectName("informationWidgetFrame");
  
  QVBoxLayout* verticalLayout_2 = new QVBoxLayout(informationWidgetFrame);
  verticalLayout_2->setSpacing(0);
  verticalLayout_2->setContentsMargins(0, 0, 0, 0);

  QScrollArea* informationScrollArea = new QScrollArea(informationWidgetFrame);
  informationScrollArea->setObjectName("informationScrollArea") ;
  informationScrollArea->setWidgetResizable(true);

  pqProxyInformationWidget* informationWidget = new pqProxyInformationWidget();
  informationWidget->setObjectName("informationWidget");
  informationWidget->setGeometry(QRect(0, 0, 77, 214));
  informationScrollArea->setWidget(informationWidget);

  verticalLayout_2->addWidget(informationScrollArea);
  informationDock->setWidget(informationWidgetFrame);

  myDockWidgets[informationDock] = true;

  desk->setTabPosition(Qt::LeftDockWidgetArea, QTabWidget::North);
  desk->tabifyDockWidget(objectInspectorDock, displayDock);
  desk->tabifyDockWidget(objectInspectorDock, informationDock);
  objectInspectorDock->raise();

  // Statistic View
  QDockWidget* statisticsViewDock  = new QDockWidget( tr( "TTL_STATISTICS_VIEW" ), desk );
  statisticsViewDock->setObjectName("statisticsViewDock");
  statisticsViewDock->setAllowedAreas(Qt::BottomDockWidgetArea|Qt::LeftDockWidgetArea|
                                      Qt::NoDockWidgetArea|Qt::RightDockWidgetArea );
  desk->addDockWidget( Qt::BottomDockWidgetArea, statisticsViewDock );
  pqDataInformationWidget* aStatWidget = new pqDataInformationWidget(statisticsViewDock);
  statisticsViewDock->setWidget(aStatWidget);
  myDockWidgets[statisticsViewDock] = false; // hidden by default

  //Animation view
  QDockWidget* animationViewDock     = new QDockWidget( tr( "TTL_ANIMATION_VIEW" ), desk );
  animationViewDock->setObjectName("animationViewDock");
  desk->addDockWidget( Qt::BottomDockWidgetArea, animationViewDock );
  pqPVAnimationWidget* animation_panel = new pqPVAnimationWidget(animationViewDock);
  animationViewDock->setWidget(animation_panel);
  myDockWidgets[animationViewDock] = false; // hidden by default

  desk->tabifyDockWidget(animationViewDock,  statisticsViewDock);

  // Selection view
  QDockWidget* selectionInspectorDock = new QDockWidget( tr( "TTL_SELECTION_INSPECTOR" ), desk );
  selectionInspectorDock->setObjectName("selectionInspectorDock");
  selectionInspectorDock->setAllowedAreas( Qt::AllDockWidgetAreas );
  desk->addDockWidget( Qt::LeftDockWidgetArea, selectionInspectorDock );
  pqSelectionInspectorPanel* aSelInspector = new pqSelectionInspectorWidget(selectionInspectorDock);
  selectionInspectorDock->setWidget(aSelInspector);
  myDockWidgets[selectionInspectorDock] = false; // hidden by default

  // Comparative View
  QDockWidget* comparativePanelDock  = new QDockWidget( tr( "TTL_COMPARATIVE_VIEW_INSPECTOR" ), desk );
  comparativePanelDock->setObjectName("comparativePanelDock");
  desk->addDockWidget( Qt::LeftDockWidgetArea, comparativePanelDock );
  pqComparativeVisPanel* cv_panel = new pqComparativeVisPanel( comparativePanelDock );
  comparativePanelDock->setWidget(cv_panel);
  myDockWidgets[comparativePanelDock] = false; // hidden by default

  // Collaboration view
  QDockWidget* collaborationPanelDock = new QDockWidget(tr( "TTL_COLLABORATIVE_DOCK" ), desk);
  collaborationPanelDock->setObjectName("collaborationPanelDock");
  pqCollaborationPanel* collaborationPanel = new pqCollaborationPanel();
  collaborationPanel->setObjectName("collaborationPanel");
  collaborationPanelDock->setWidget(collaborationPanel);
  desk->addDockWidget(Qt::RightDockWidgetArea, collaborationPanelDock);
  myDockWidgets[collaborationPanelDock] = false; // hidden by default
  
  // Memory inspector dock
  QDockWidget* memoryInspectorDock = new QDockWidget(tr( "TTL_MEMORY_INSPECTOR" ), desk);
  memoryInspectorDock->setObjectName("memoryInspectorDock");
  pqMemoryInspector* dockWidgetContents = new pqMemoryInspector();
  dockWidgetContents->setObjectName("dockWidgetContents");
  memoryInspectorDock->setWidget(dockWidgetContents);
  desk->addDockWidget(Qt::RightDockWidgetArea, memoryInspectorDock);
  myDockWidgets[memoryInspectorDock] = false; // hidden by default


  // Setup the statusbar ...
  pqProgressWidget* aProgress = new pqProgressWidget(desk->statusBar());
  desk->statusBar()->addPermanentWidget(aProgress);
  
  // Set up the dock window corners to give the vertical docks more room.
  desk->setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
  desk->setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);
  
  // Setup the default dock configuration ...
  statisticsViewDock->hide();
  comparativePanelDock->hide();
  animationViewDock->hide();
  selectionInspectorDock->hide();
  collaborationPanelDock->hide();
  memoryInspectorDock->hide();


  // Setup quick-launch shortcuts.
  QShortcut *ctrlSpace = new QShortcut(Qt::CTRL + Qt::Key_Space, desk);
  QObject::connect(ctrlSpace, SIGNAL(activated()), pqApplicationCore::instance(), SLOT(quickLaunch()));
  QShortcut *altSpace = new QShortcut(Qt::ALT + Qt::Key_Space, desk);
  QObject::connect(altSpace, SIGNAL(activated()), pqApplicationCore::instance(), SLOT(quickLaunch()));

}

/*!
  \brief Save states of dockable ParaView widgets.
*/
void PVGUI_Module::saveDockWidgetsState()
{
  SUIT_Desktop* desk = application()->desktop();

  // VSR: 19/06/2011: do not use Paraview's methods, since it conflicts with SALOME GUI architecture
  // ... the following code is commented...
  // Save the state of the window ...
  // pqApplicationCore::instance()->settings()->saveState(*desk, "MainWindow");
  //
  //for (int i = 0; i < myDockWidgets.size(); ++i)
  //  myDockWidgets.at(i)->setParent(0);
  // ... and replaced - manually hide dock windows

  // store dock widgets visibility state and hide'em all
  QMapIterator<QWidget*, bool> it1( myDockWidgets );
  while( it1.hasNext() ) {
    it1.next();
    QDockWidget* dw = qobject_cast<QDockWidget*>( it1.key() );
    myDockWidgets[dw] = dw->isVisible();
    dw->setVisible( false );
    dw->toggleViewAction()->setVisible( false );
  }
  // store toolbar breaks state and remove all tollbar breaks 
  QMapIterator<QWidget*, bool> it2( myToolbarBreaks );
  while( it2.hasNext() ) {
    it2.next();
    QToolBar* tb = qobject_cast<QToolBar*>( it2.key() );
    myToolbarBreaks[tb] = desk->toolBarBreak( tb );
    if ( myToolbarBreaks[tb] )
      desk->removeToolBarBreak( tb );
  }
  // store toolbars visibility state and hide'em all
  QMapIterator<QWidget*, bool> it3( myToolbars );
  while( it3.hasNext() ) {
    it3.next();
    QToolBar* tb = qobject_cast<QToolBar*>( it3.key() );
    myToolbars[tb] = tb->isVisible();
    tb->setVisible( false );
    tb->toggleViewAction()->setVisible( false );
  }
}

/*!
  \brief Restore states of dockable ParaView widgets.
*/
void PVGUI_Module::restoreDockWidgetsState()
{
  SUIT_Desktop* desk = application()->desktop();

  // VSR: 19/06/2011: do not use Paraview's methods, since it conflicts with SALOME GUI architecture
  // ... the following code is commented...
  //for (int i = 0; i < myDockWidgets.size(); ++i)
  //  myDockWidgets.at(i)->setParent(desk);
  //
  // Restore the state of the window ...
  //pqApplicationCore::instance()->settings()->restoreState("MainWindow", *desk);
  // ... and replaced - manually hide dock windows

  // restore dock widgets visibility state
  QMapIterator<QWidget*, bool> it1( myDockWidgets );
  while( it1.hasNext() ) {
    it1.next();
    QDockWidget* dw = qobject_cast<QDockWidget*>( it1.key() );
    dw->setVisible( it1.value() );
    dw->toggleViewAction()->setVisible( true );
  }
  // restore toolbar breaks state
  QMapIterator<QWidget*, bool> it2( myToolbarBreaks );
  while( it2.hasNext() ) {
    it2.next();
    QToolBar* tb = qobject_cast<QToolBar*>( it2.key() );
    if ( myToolbarBreaks[tb] )
      desk->insertToolBarBreak( tb );
  }
  // restore toolbar visibility state
  QMapIterator<QWidget*, bool> it3( myToolbars );
  while( it3.hasNext() ) {
    it3.next();
    QToolBar* tb = qobject_cast<QToolBar*>( it3.key() );
    tb->setVisible( it3.value() );
    tb->toggleViewAction()->setVisible( true );
  }
}
