#ifndef QGA_COMMONS_HPP
#define QGA_COMMONS_HPP

#include <atomic>
#include <string>
#include <iomanip>

/* Forward declarations */
namespace Config {
  extern const unsigned nBit;
  extern const size_t popSize;
  extern const float selectBias;
  extern const float heurFactor;
  extern const float expLengthIni;
  extern const float expLengthAdd;
  extern const float pDeleteUniform;
  extern const float pControl;
}

namespace Wrapper {
  std::string gate_name(unsigned);
  extern const unsigned gate_cnt;
}
/* End forward declarations */

#include "GeneTools.hpp"
#include "Fitness.hpp"
#include "CandidateCounter.hpp"
#include "CandidateBase.hpp"
#include "CandidateFactory.hpp"

#endif // !defined QGA_COMMONS_HPP
