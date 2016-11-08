#ifndef TOOL_BAR_H_
#define TOOL_BAR_H_

#include <map>
#include <string>

#include <QtGui>


////////////////////////////////////////////////////////////////////////////////
/// Provides the application toolbar.
////////////////////////////////////////////////////////////////////////////////
class tool_bar
  : public QToolBar
{

  Q_OBJECT

  ///@name Internal State
  ///@{

  std::map<std::string, QAction*> m_actions;  ///< Action list.

  ///@}

  public:

    ///@name Construction
    ///@{

    tool_bar(QWidget* _parent);

    ///@}

  signals:

    ///@name Simulation State Signals
    ///@{

    void start(); ///< Signal the simulation to start.
    void reset(); ///< Signal the simulation to reset.

    ///@}
};

#endif
