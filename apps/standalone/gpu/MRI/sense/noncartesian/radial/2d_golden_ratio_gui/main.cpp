////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

#include "radialSenseAppMainWidget.h"

#include <stdlib.h>

#include <QtGui/QApplication>

int
main( int argc, char** argv) 
{
  QApplication app(argc, argv);
  radialSenseAppMainWindow window;
  window.show();
  
  return app.exec();
}
