#include "pqSetModeStarter.h"

#include <pqServerManagerModel.h>

#include <pqApplicationCore.h>
#include <pqRenderView.h>
#include <pqSettings.h>


//-----------------------------------------------------------------------------
pqSetModeStarter::pqSetModeStarter(QObject* p/*=0*/)
  : QObject(p)
{
}

//-----------------------------------------------------------------------------
pqSetModeStarter::~pqSetModeStarter()
{
}


//-----------------------------------------------------------------------------
void pqSetModeStarter::onStartup()
{
  this->setStandardMode();
}

//-----------------------------------------------------------------------------
void pqSetModeStarter::setStandardMode()
{
  pqSettings* settings = pqApplicationCore::instance()->settings();
  settings->beginGroup("renderModule");
  if (!settings->contains("InteractorStyle/CameraManipulators")) {
    // Set Post-Pro-like settings
    QStringList strs;
    pqRenderView::ManipulatorType manips[9];
    const pqRenderView::ManipulatorType* default3DManips = pqRenderView::getDefaultManipulatorTypes();

    // Copy default settings, make changes for Ctrl+MB and MB modes
    for(int i=0; i<9; i++)
      {
	manips[i] = default3DManips[i];

	// Ctrl+MB
	if (manips[i].Shift == 0 && manips[i].Control == 1) {
	  if (manips[i].Mouse == 1)
	    manips[i].Name = QByteArray("Zoom");
	  else  if (manips[i].Mouse == 2)
	    manips[i].Name = QByteArray("Pan");
	  else  if (manips[i].Mouse == 3)
	    manips[i].Name = QByteArray("Rotate");
	}

	// MB only
	if (manips[i].Shift == 0 && manips[i].Control == 0) {
	  if (manips[i].Mouse == 1)
	    manips[i].Name = QByteArray("Rotate");
	  else  if (manips[i].Mouse == 2)
	    manips[i].Name = QByteArray("Pan");
	  else  if (manips[i].Mouse == 3)
	    manips[i].Name = QByteArray("Zoom");
	}
      }

    // Save settings
    for(int i=0; i<9; i++)
      {
	strs << QString("Manipulator%1Mouse%2Shift%3Control%4Name%5")
	  .arg(i+1)
	  .arg(manips[i].Mouse)
	  .arg(manips[i].Shift)
	  .arg(manips[i].Control)
	  .arg(QString(manips[i].Name));
      }
      
    settings->setValue("InteractorStyle/CameraManipulators", strs);
  }
  settings->endGroup();

  // Loop through render views and apply the settings
  QList<pqRenderViewBase*> views =
    pqApplicationCore::instance()->getServerManagerModel()->
    findItems<pqRenderViewBase*>();

  foreach(pqRenderViewBase* view, views) {
    view->restoreSettings(true);
  }
}