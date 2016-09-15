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
  1, 0, 0, {0,1}
};

struct Gate {
  arma::cx_mat22 op;
  std::string name;
};

std::vector<Gate> gates {
  { H, "H" },
  { X, "X" },
  { Y, "Y" },
  { Z, "Z" },
  { T, "T" }
};

arma::cx_vec out{};

} // namespace internal


class Candidate : public ::CBase<Candidate> {

  using Base = ::CBase<Candidate>;

public:

  using Base::Base;

  double error() const {
    return 1 - std::abs(arma::cdot(internal::out, sim()));
  }

  std::string dump() const {
    std::ostringstream os{};
    sim().t().raw_print(os);
    return os.str();
  }

private:

  arma::cx_vec sim() const {
    arma::cx_vec psi = qic::mket({0}, {arma::uword(1) << Config::nBit});
    for(const Gene& g : gt) {
      /* convert std::vector<unsigned> to arma::uvec, adding 1 */
      auto& ixv = g.ix_vector();
      arma::uvec ixs(ixv.size());
      for(size_t i = 0; i < ixv.size(); i++)
        ixs.at(i) = ixv[i] + 1;
      /* control-gate (QIClib) */
      psi = qic::apply_ctrl(
          psi,                          // state
          internal::gates[g.gate()].op, // operator
          ixs,                          // arma::uvec of control systems
          {1 + g.target()});            // arma::uvec of target systems
    }
    return psi;
  }

}; // class CBase


const unsigned gate_cnt = internal::gates.size();


std::string gate_name(unsigned ix) {
  return internal::gates[ix].name;
}


void init() {
  internal::out = qic::mket({3}, {arma::uword(1) << Config::nBit});
}

} // namespace Wrapper
