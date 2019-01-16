#ifndef PACKET_H_
#define PACKET_H_


////////////////////////////////////////////////////////////////////////////////
/// A packet of information about the camera's position in the local frame of a
/// detected aruco marker.
///
/// The packet is used to send data over the network in integer format. The
/// floating point values (distance, heading) will be multiplied by s_factor so
/// that they can be transported as integers.
///
/// @warning We do not currently worry about byte-ordering as we are always
///          going between intel (little endian) machines. We will need to
///          address this in the future if we wish to transport packets across
///          machines with different byte ordering.
////////////////////////////////////////////////////////////////////////////////
struct packet {

  static constexpr double s_factor = 10000;

  int32_t id{0}, x{0}, y{0}, t{0};

  packet() = default;

  packet(int _id, int _x, int _y, int _t)
      : id(_id), x(_x), y(_y), t(_t) { }

};

#endif
