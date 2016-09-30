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

using cxd = arma::cx_double;

const double pi = std::acos(-1);

const cxd i{0,1};

const double v12 = 1/std::sqrt(2);

arma::cx_mat22 I {
  1, 0, 0, 1
};

arma::cx_mat22 H {
  v12, v12, v12, -v12
};

arma::cx_mat22 X {
  0, 1, 1, 0
};

arma::cx_mat22 Y {
  0, -i, i, 0
};

arma::cx_mat22 Z {
  1, 0, 0, -1
};

arma::cx_mat22 T {
  1, 0, 0, std::exp(i*pi/4.)
};

arma::cx_mat22 Ti {
  1, 0, 0, std::exp(-i*pi/4.)
};

arma::cx_mat22 S {
  1, 0, 0, i
};

arma::cx_mat22 Si {
  1, 0, 0, -i
};

struct Gate {
  arma::cx_mat22 op;
  std::string name;
  int inv;
  int sq;
};

std::vector<Gate> gates {
  { I, "I", 0, 0 },
  { H, "H", 0, -1 },
/*{ X, "X", 0, -2 },
  { Y, "Y", 0, -3 },
  { Z, "Z", 0, -4 },*/
  { T, "T", +1, 0/*+2*/ },
  { Ti, "Ti", -1, 0/*+2*/ },
/*{ S, "S", +1, -3 },
  { Si, "Si", -1, -4 }*/
};

arma::cx_vec out{};

} // namespace internal


class Gene {

  size_t op;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{1, internal::gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    return {dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
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
    } else if(g.op == 0) {
      // G * Identity = G
      return true;
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)
        && internal::gates[op].sq != 0) {
      // G * G = square(G) if also among our operations
      op += internal::gates[op].sq;
      return true;
    } else return false;
  }

  void mutate() {
    /* no-op */
  }

  bool simplify() {
    /* no-op */
    return false;
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

private:

  NOINLINE Gene(size_t op_, unsigned tgt_, unsigned control_enc):
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

}; // class Gene


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
    for(const auto& g : gt) {
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
