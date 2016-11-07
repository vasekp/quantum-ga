#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <atomic>
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
#include <regex>

using std::size_t;

/* Forward declarations */
namespace Config {
  extern const unsigned nBit;
  extern const size_t arSize;
  extern const size_t popSize;
  extern const double selectBias;
  extern const double heurFactor;
  extern const double expLengthIni;
  extern const double expMutationCount;
  extern const double pCrossUniform;
  extern const double pControl;
  extern const double dAlpha;
}

/* Useful constants and typedefs */
namespace QGA {
  namespace Const {
    const std::complex<double> i{0, 1};
    const double pi = std::acos(-1);
    const double v12 = 1/std::sqrt(2);
  }
}

#ifdef USE_QPP
  #include "QGA_Backend/QPP.hpp"
#elif defined USE_QICLIB
  #include "QGA_Backend/QIClib.hpp"
#else
  #error Either USE_QPP or USE_QICLIB needed.
#endif

#include "QGA_bits/CircuitPrinter.hpp"
#include "QGA_bits/Fitness.hpp"
#include "QGA_bits/GateBase.hpp" // uses Fitness.hpp and CircuitPrinter.hpp
#include "QGA_bits/Gene.hpp"     // uses GateBase.hpp
#include "QGA_bits/Tools.hpp"
#include "QGA_bits/Gates.hpp"    // uses Tools.hpp and GateBase.hpp
#include "QGA_bits/CandidateCounter.hpp"
#include "QGA_bits/CandidateBase.hpp"
#include "QGA_bits/CandidateFactory.hpp"

#endif // !defined QGA_COMMONS_HPP
