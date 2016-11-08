#include <iostream>
#include <QApplication>

#include "gui/base_visualization.h"
#include "gui/gl_widget.h"
#include "gui/main_window.h"

#include "example_visualization.h"

int main(int _argc, char* _argv[]) {
  // Create application and main window.
  QApplication app(_argc, _argv);
  main_window window;

  example_visualization sim;
  window.gl()->visualization(&sim);

  // Show window and execute app.
  window.show();
  app.exec();

  return 0;
}
