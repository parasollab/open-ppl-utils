#ifndef MAIN_WINDOW_H_
#define MAIN_WINDOW_H_

#include <QtGui>

#include <vector>

class gl_widget;
class tool_bar;


////////////////////////////////////////////////////////////////////////////////
/// Owns and manages all the application's widgets.
////////////////////////////////////////////////////////////////////////////////
class main_window
  : public QMainWindow
{

  Q_OBJECT

  ///@name Layout Components
  ///@{

  QStackedWidget* m_stack;      ///< Handles window switching.
  gl_widget*      m_gl_widget;  ///< Manages and displays the OpenGL scene.
  tool_bar*       m_tool_bar;   ///< The tool bar at the top of the window.
  QStatusBar*     m_status_bar; ///< Status bar at the bottom of the window.

  ///@}

  public:

    ///@name Construction
    ///@{

    main_window();
    main_window(const main_window&) = delete;
    virtual ~main_window() = default;

    ///@}
    ///@name Accessors
    ///@{

    gl_widget* gl() const;

    ///@}

  public slots:

    ///@name Window Controls
    ///@{

    /// Display a widget from the window stack.
    /// @param[in] _w The widget to display.
    void show_window(QWidget* const _w);

    /// Add a widget to the window stack.
    /// @param[in] _w The widget to add.
    void add_window(QWidget* const _w);

    /// Delete a widget from the window stack and de-allocate its memory.
    /// @param[in] _w The widget to delete.
    void delete_window(QWidget* const _w);

    ///@}
};

#endif
