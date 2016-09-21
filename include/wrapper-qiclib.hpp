// allow only one wrapper
#ifndef QGA_WRAPPER_HPP
#define QGA_WRAPPER_HPP

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include "QIClib"
#include <string>
#include <cmath>
#include <sstream>

namespace Wrapper {

namespace internal {

arma::cx_mat22 H {
  1/std::sqrt(2),  1/std::sqrt(2),
  1/std::sqrt(2), -1/std::sqrt(2)
};

arma::cx_mat22 X {
  0, 1, 1, 0
};

arma::cx_mat22 Y {
  0, {0,-1}, {0,1}, 0
};

arma::cx_mat22 Z {
  1, 0, 0, -1
};

arma::cx_mat22 T {
  1, 0, 0, {1/std::sqrt(2), 1/std::sqrt(2)}
};

arma::cx_mat22 Ti {
  1, 0, 0, {1/std::sqrt(2), -1/std::sqrt(2)}
};

struct Gate {
  arma::cx_mat22 op;
  std::string name;
  int inv;
};

std::vector<Gate> gates {
  { H, "H", 0 },
/*{ X, "X", 0 },
  { Y, "Y", 0 },
  { Z, "Z", 0 },*/
  { T, "T", +1 },
  { Ti, "Ti", -1 }
};

arma::cx_vec out{};

} // namespace internal


class GeneFactory;


class Gene {

  unsigned op;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;

  friend class GeneFactory;

public:

  using Factory = GeneFactory;

  NOINLINE Gene(unsigned op_, unsigned tgt_, unsigned control_enc):
      op(op_), tgt(tgt_), hw(0) {
    std::vector<arma::uword> ixv;
    ixv.reserve(Config::nBit);
    unsigned ctrl = QGA::GeneTools::ctrlBitString(control_enc, tgt - 1);
    for(unsigned i = 0; i < Config::nBit; i++) {
      if(ctrl & 1) {
        ixv.push_back(i + 1);
        hw++;
      }
      ctrl >>= 1;
    }
    ixs = ixv;
  }

  const arma::uvec& ix_vector() const {
    return ixs;
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

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << g.gate().name << g.target();
    if(g.ixs.size()) {
      os << '[';
      for(auto ctrl : g.ixs)
        os << ctrl;
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
  std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
  // distribution of controls
  std::uniform_int_distribution<unsigned> dCtrl{};

public:

  GeneFactory() { }

  Gene getNew() {
    return {dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
  }

  Gene invert(const Gene& g) {
    Gene ret = g;
    ret.op += g.gate().inv;
    return ret;
  }

  Gene&& invert(Gene&& g) {
    g.op += g.gate().inv;
    return std::move(g);
  }

}; // class GeneFactory


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    return 1 - std::abs(arma::cdot(internal::out, sim()));
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    sim().st().raw_print(os);
    return os.str();
  }

private:

  arma::cx_vec sim() const {
    arma::cx_vec psi = qic::mket({0}, {arma::uword(1) << Config::nBit});
    for(const Gene& g : gt) {
      /* control-gate (QIClib) */
      psi = qic::apply_ctrl(
          psi,            // state
          g.gate().op,    // operator
          g.ix_vector(),  // arma::uvec of control systems
          {g.target()});  // arma::uvec of target systems
    }
    return psi;
  }

}; // class Candidate


void init() {
  internal::out = qic::mket({3}, {arma::uword(1) << Config::nBit});
}

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
