#include "tool_bar.h"

/*------------------------------ Construction --------------------------------*/

tool_bar::
tool_bar(QWidget* _parent)
  : QToolBar(_parent)
{
  // Define buttons.
  m_actions["play"] = this->addAction("Play", this, SIGNAL(start()));
  m_actions["reset"] = this->addAction("Reset", this, SIGNAL(reset()));

  // Make toolbar a permanent feature of the parent widget.
  this->setContextMenuPolicy(Qt::PreventContextMenu);
  this->setFloatable(false);
  this->setMovable(false);
  this->setToolButtonStyle(Qt::ToolButtonTextOnly);
  this->adjustSize();
}

/*----------------------------------------------------------------------------*/
