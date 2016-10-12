#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <atomic>
#include <string>
#include <iomanip>
#include <cstddef>
#include <complex>
#include <string>
#include <cmath>
#include <sstream>

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
  extern const double pChoiceUniform;
  extern const double pCrossUniform;
  extern const double pControl;
}

/* Useful constants and typedefs */
namespace QGA {
  namespace Const {
    const std::complex<double> i{0, 1};
    const double pi = std::acos(-1);
    const double v12 = 1/std::sqrt(2);
  }
}

#include "GateBase.hpp"
#include "Gene.hpp"
#include "Fitness.hpp"
#include "CandidateCounter.hpp"
#include "CandidateBase.hpp"
#include "CandidateFactory.hpp"

#include "Colours.hpp"
#include "Tools.hpp"

#endif // !defined QGA_COMMONS_HPP
