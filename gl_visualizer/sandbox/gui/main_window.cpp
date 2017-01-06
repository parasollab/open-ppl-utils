#include "main_window.h"

#include <iostream>

#include "gl_widget.h"
#include "tool_bar.h"

/*-------------------------------- Construction ------------------------------*/

main_window::
main_window()
  : QMainWindow(nullptr)
{
  // Set window properties.
  setMinimumSize(800, 600);
  setWindowTitle("OpenGL Sandbox");

  // Create the status bar.
  m_status_bar = new QStatusBar(this);
  m_status_bar->setSizeGripEnabled(false);
  m_status_bar->showMessage("Ready.");
  this->setStatusBar(m_status_bar);

  // Create stack layout.
  m_stack = new QStackedWidget(this);
  this->setCentralWidget(m_stack);

  // Create the gl window.
  m_gl_widget = new gl_widget(this);
  m_stack->addWidget(m_gl_widget);

  // Create the tool bar.
  m_tool_bar = new tool_bar(this);

  // Connect components.

  // Connect the tool bar to the gl widget.
  connect(m_tool_bar, SIGNAL(start()),
          m_gl_widget,  SLOT(start()));
  connect(m_tool_bar, SIGNAL(reset()),
          m_gl_widget,  SLOT(reset()));

  // Connect the gl widget to the status bar.
  connect(m_gl_widget, SIGNAL(status_message(const QString&, int)),
          m_status_bar,  SLOT(showMessage(const QString&, int)));
}

/*------------------------------- Accessors ----------------------------------*/

gl_widget*
main_window::
gl() const
{
  return m_gl_widget;
}

/*--------------------------- Visualization Interface ------------------------*/

void
main_window::
visualization(base_visualization* const _v)
{
  m_gl_widget->visualization(_v);
}


base_visualization*
main_window::
visualization()
{
  return m_gl_widget->visualization();
}

/*---------------------------- Window Controls -------------------------------*/

void
main_window::
show_window(QWidget* const _w)
{
  int i = m_stack->indexOf(_w);
  if(m_stack->widget(i))
    m_stack->setCurrentIndex(i);
  else
    std::cerr << "main_window::show_window error: requested display of non-"
              << "existent stack widget." << std::endl;
}


void
main_window::
add_window(QWidget* const _w)
{
  int i = m_stack->indexOf(_w);
  if(!m_stack->widget(i))
    m_stack->addWidget(_w);
  else
    std::cerr << "main_window::add_window error: requested re-addition of a "
              << "widget that is already on the stack." << std::endl;
}


void
main_window::
delete_window(QWidget* const _w)
{
  int i = m_stack->indexOf(_w);
  if(m_stack->widget(i)) {
    m_stack->removeWidget(_w);
    delete _w;
  }
  else
    std::cerr << "main_window::delete_window error: requested deletion of a "
              << "widget that isn't on the stack." << std::endl;
}

/*----------------------------------------------------------------------------*/
