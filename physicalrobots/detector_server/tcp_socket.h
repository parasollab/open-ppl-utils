#ifndef NONSTD_TCP_SOCKET_H_
#define NONSTD_TCP_SOCKET_H_

#include <atomic>
#include <functional>
#include <string>
#include <thread>
#include <sys/socket.h>

#include "runtime.h"

namespace nonstd {

  //////////////////////////////////////////////////////////////////////////////
  /// Managed socket for TCP connections, providing both client and server
  /// functionality.
  ///
  /// Client: use the client constructor, or call connect.
  /// Server: define a handler function, then call listen.
  //////////////////////////////////////////////////////////////////////////////
  class tcp_socket final
  {

    public:

      ///@name Local Types
      ///@{

      /// This kind of function will be called when a listening socket receives
      /// incoming connections. The socket reference is used to send and receive
      /// data.
      typedef std::function<void(tcp_socket&&)> handler_function;

      ///@}

    private:

      ///@name Universal State
      ///@{

      static constexpr bool m_debug{true};   ///< Show debugging messages?

      int m_id{-1};                          ///< The socket descriptor.
      std::string m_server;                  ///< The connected server.
      std::string m_port;                    ///< The connection port.

      ///@}
      ///@name Server-Specific State
      ///@{

      std::atomic<bool> m_listening{false};  ///< Listening for new connections?
      std::thread m_listening_thread;        ///< Separate thread for listening.
      handler_function m_handler;            ///< The handler function to use.

      ///@}
      ///@name Construction
      ///@{

      /// Construct a socket from an existing connection. Only for internal use
      /// when receiving new connections.
      /// @param[in] _socket An open TCP socket.
      /// @param[in] _server The connected host.
      /// @param[in] _port The port in use.
      tcp_socket(const int _socket, const std::string& _server,
          const std::string& _port);

    public:

      /// Construct a socket without making a connection.
      tcp_socket() = default;

      /// Construct a new socket and connect (as a client) to the designated
      /// server.
      /// @param[in] _server The host address to connect to.
      /// @param[in] _port The port to use.
      tcp_socket(const std::string& _server, const std::string& _port);

      /// No copy, only move allowed.
      tcp_socket(const tcp_socket&) = delete;

      /// The owned socket is closed on destruction.
      ~tcp_socket();

      ///@}
      ///@name Transmission Interface
      ///@{

      /// Send a transmission to the server.
      /// @tparam T The transmission object type.
      /// @param[in] _t The object to send.
      /// @warning This does not perform any serialization, so the host and
      ///          client machines must use the same representation for this
      ///          type.
      template <typename T>
      tcp_socket& operator<<(const T& _t);

      /// Receive a transmission from the server.
      /// @tparam T The transmission object type.
      /// @param[in] _t A writable location for the received object.
      /// @warning This does not perform any serialization, so the host and
      ///          client machines must use the same representation for this
      ///          type.
      template <typename T>
      tcp_socket& operator>>(T& _t);

      ///@}
      ///@name Connection Interface
      ///@{

      /// Close the current connection.
      void disconnect();

      /// Connect to a listening server as a client.
      /// @param[in] _server The server IP address.
      /// @param[in] _port The port to use.
      /// @return True if the connection was successful.
      bool connect(const std::string& _server, const std::string& _port);

      /// Listen for incoming client connections on the designated port.
      /// @param[in] _port The port to listen on.
      /// @param[in] _backlog The number of backlogged connections to allow.
      /// @param[in] _concurrent Use separate thread? Blocking listen if false.
      /// @return A bool indicating success or failure.
      bool listen(const std::string& _port, int _backlog = 12,
          bool _concurrent = true);

      /// Set the handler function to use on incoming connections.
      /// @param[in] _f The handler function to use. If the socket will be
      ///               listening concurrently, this must be thread-safe.
      void set_handler(const handler_function& _f) {m_handler = _f;}

    private:

      /// Listen for new connections. When one is received, call m_handler with
      /// the corresponding tcp_socket. This function is executed in a separate
      /// thread when using concurrent listening.
      static void listener(tcp_socket* _t);

      /// Get the description of a socket error code (stored in errno).
      /// @return The descriptive meaning of the error.
      static std::string print_error();

      ///@}
  };

  /*---------------------- Transmission Interface ----------------------------*/

  template <typename T>
  tcp_socket&
  tcp_socket::
  operator<<(const T& _t)
  {
    assert_msg(m_id != -1,
        "nonstd::tcp_socket error: tried to write to unopened connection.");

    assert_msg(send(m_id, (void*)&_t, sizeof(_t), 0) != -1,
        "nonstd::tcp_socket error: couldn't write to socket.");

    return *this;
  }


  template <typename T>
  tcp_socket&
  tcp_socket::
  operator>>(T& _t)
  {
    assert_msg(m_id != -1,
        "nonstd::tcp_socket error: tried to read from unopened connection.");

    assert_msg(recv(m_id, (void*)&_t, sizeof(_t), MSG_WAITALL) != -1,
        "nonstd::tcp_socket error: couldn't read from socket.");

    return *this;
  }

  /*--------------------------------------------------------------------------*/
}

#endif
