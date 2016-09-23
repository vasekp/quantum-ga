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
  int inv;
  int sq;
};

std::vector<Gate> gates {
  { qpp::gt.Id2, "I", 0, 0 },
  { qpp::gt.H, "H", 0, -1 },
/*{ qpp::gt.X, "X", 0, -2 },
  { qpp::gt.Y, "Y", 0, -3 },
  { qpp::gt.Z, "Z", 0, -4 },*/
  { qpp::gt.T, "T", +1, 0/*+2*/ },
  { qpp::gt.T.conjugate(), "Ti", -1, 0/*+2*/ },
/*{ qpp::gt.S, "S", +1, -3 },
  { qpp::gt.S.conjugate(), "Si", -1, -4 }*/
};

qpp::ket out{};

} // namespace internal


class Gene {

  unsigned op;
  unsigned tgt;
  unsigned hw;
  std::vector<qpp::idx> ixv{};

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<unsigned> dOp{1,
      (unsigned)internal::gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    return {dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
  }

  const std::vector<qpp::idx>& ix_vector() const {
    return ixv;
  }

  unsigned target() const {
    return tgt;
  }

  const internal::Gate& gate() const {
    return internal::gates[op];
  }

  unsigned weight() const {
    return hw;
  }

  bool invert() {
    int dIx = gate().inv;
    op += dIx;
    return dIx ? true : false;
  }

  bool merge(const Gene& g) {
    if(op == 0) {
      // Identity * G = G
      *this = g;
      return true;
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixv == ixv
        && internal::gates[op].sq != 0) {
      // G * G = square(G) if also among our operations
      op += internal::gates[op].sq;
      return true;
    } else return false;
  }

  bool mutate() {
    /* no-op */
    return false;
  }

  bool simplify() {
    /* no-op */
    return false;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << g.gate().name << g.target() + 1;
    if(g.ixv.size()) {
      os << '[';
      for(auto ctrl : g.ixv)
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

private:

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

}; // class Gene


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
          psi,            // state
          g.gate().op,    // operator
          g.ix_vector(),  // vector<qpp::idx> of control systems
          {g.target()});  // vector<qpp:idx> of target systems
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
