#include <chrono>
#include <ctgmath>

#include "nonstd/random.h"

using namespace std;

namespace nonstd {

  mt19937
  random_generator(chrono::system_clock::now().time_since_epoch().count());

  double
  drand()
  {
    static uniform_real_distribution<double> distribution(0., 1.);
    return distribution(random_generator);
  }


  double
  grand()
  {
    double x = 2. * drand48() - 1.;
    double y = 2. * drand48() - 1.;
    return x / sqrt(x * x + y * y);
  }

}
