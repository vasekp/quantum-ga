#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <atomic>
#include <string>
#include <iomanip>
#include <cstddef>

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

#include "Fitness.hpp"
#include "CandidateCounter.hpp"
#include "CandidateBase.hpp"
#include "CandidateFactory.hpp"

#include "colours.hpp"

#endif // !defined QGA_COMMONS_HPP
