#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <complex>
#include <string>
#include <cstddef>

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>

#include <iomanip>
#include <ostream>
#include <sstream>

#include <array>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>

using std::size_t;

/* Forward declarations */
namespace Config {
  extern unsigned nBit;
  extern size_t arSize;
  extern size_t popSize;
  extern size_t maxGen;
  extern double selectBias;
  extern double expLengthIni;
  extern double expMutationCount;
  extern double expSliceLength;
  extern const double pControl;
  extern const double dAlpha;
  extern const size_t circLineLength;
}

/* Useful constants and typedefs */
namespace QGA {
  namespace Const {
    const std::complex<double> i{0, 1};
    const double pi = std::acos(-1);
    const double v12 = 1/std::sqrt(2);
  }
}

#endif // !defined QGA_COMMONS_HPP
