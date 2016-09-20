// allow only one wrapper
#ifndef QGA_WRAPPER_HPP
#define QGA_WRAPPER_HPP

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include "qpp.h"
#include <string>
#include <cmath>
#include <sstream>

namespace Wrapper {

namespace internal {

struct Gate {
  qpp::cmat op;
  std::string name;
};

std::vector<Gate> gates {
  { qpp::gt.H, "H" },
/*{ qpp::gt.X, "X" },
  { qpp::gt.Y, "Y" },
  { qpp::gt.Z, "Z" },*/
  { qpp::gt.T, "T" }
};

qpp::ket out{};

} // namespace internal


class GeneFactory;


class Gene {

  unsigned op;
  unsigned tgt;
  unsigned hw;
  std::vector<qpp::idx> ixv{};

public:

  using Factory = GeneFactory;

  NOINLINE Gene(unsigned op_, unsigned tgt_, unsigned control_enc):
      op(op_), tgt(tgt_), hw(0) {
    ixv.reserve(Config::nBit);
    unsigned ctrl = QGA::GeneTools::ctrlBitString(control_enc, tgt);
    for(unsigned i = 0; i < Config::nBit; i++) {
      if(ctrl & 1) {
        ixv.push_back(i);
        hw++;
      }
      ctrl >>= 1;
    }
  }

  const std::vector<qpp::idx>& ix_vector() const {
    return ixv;
  }

  unsigned target() const {
    return tgt;
  }

  unsigned gate() const {
    return op;
  }

  unsigned weight() const {
    return hw;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << internal::gates[g.gate()].name << g.target() + 1;
    if(g.ixv.size()) {
      os << '[';
      for(auto ctrl : g.ixv)
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

}; // class Gene


class GeneFactory {

  // distribution of possible gates
  std::uniform_int_distribution<unsigned> dOp{0,
    (unsigned)internal::gates.size() - 1};
  // distribution of targets
  std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
  // distribution of controls
  std::uniform_int_distribution<unsigned> dCtrl{};

public:

  GeneFactory() { }

  Gene getNew() {
    return {dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
  }

}; // class GeneFactory


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    return 1 - std::abs(sim().dot(internal::out));
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
        " ", " "); // row, col separators
    os << sim().format(fmt) << '\n';
    return os.str();
  }

private:

  qpp::ket sim() const {
    std::vector<qpp::idx> dims(Config::nBit, 2);
    qpp::ket psi = qpp::mket(qpp::n2multiidx(0, dims));
    for(const Gene& g : gt) {
      /* control-gate (QIClib) */
      psi = qpp::applyCTRL(
          psi,                          // state
          internal::gates[g.gate()].op, // operator
          g.ix_vector(),                // vector<qpp::idx> of control systems
          {g.target()});                // vector<qpp:idx> of target systems
    }
    return psi;
  }

}; // class Candidate


const unsigned gate_cnt = internal::gates.size();


std::string gate_name(unsigned ix) {
  return internal::gates[ix].name;
}


void init() {
  std::vector<qpp::idx> dims(Config::nBit, 2);
  internal::out = qpp::mket(qpp::n2multiidx(3, dims));
}

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
