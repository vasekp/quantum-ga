#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <atomic>
#include <string>
#include <iomanip>

/* Forward declarations */
namespace Config {
  extern const unsigned nBit;
  extern const size_t arSize;
  extern const size_t popSize;
  extern const float selectBias;
  extern const float heurFactor;
  extern const float expLengthIni;
  extern const float expLengthAdd;
  extern const float pDeleteUniform;
  extern const float pControl;
}

#include "GeneTools.hpp"
#include "Fitness.hpp"
#include "CandidateCounter.hpp"
#include "CandidateBase.hpp"
#include "CandidateFactory.hpp"

#include "colours.hpp"

#endif // !defined QGA_COMMONS_HPP
