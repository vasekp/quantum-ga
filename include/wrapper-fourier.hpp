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

arma::cx_mat22 phi(double a) {
  return { std::exp(i*a), 0, 0, std::exp(i*a) };
}

arma::cx_mat22 xrot(double a) {
  return { std::cos(a), i*std::sin(a), i*std::sin(a), std::cos(a) };
}

arma::cx_mat22 yrot(double a) {
  return { std::cos(a), std::sin(a), -std::sin(a), std::cos(a) };
}

arma::cx_mat22 zrot(double a) {
  return { std::exp(i*a), 0, 0, std::exp(-i*a) };
}

} // namespace internal


class Gene {

  char names[4] = { 'P', 'X', 'Y', 'Z' };

  unsigned op;
  double angle;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;
  arma::cx_mat22 mat;

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<unsigned> dOp{0, 3};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{0, 2.0*internal::pi};
    return {dOp(gen::rng), dAng(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
  }

  const arma::uvec& ix_vector() const {
    return ixs;
  }

  unsigned target() const {
    return tgt;
  }

  arma::cx_mat22 gate() const {
    switch(op) {
      case 0:
        return internal::phi(angle);
      case 1:
        return internal::xrot(angle);
      case 2:
        return internal::yrot(angle);
      case 3:
        return internal::zrot(angle);
      default:
        throw std::logic_error("gate must be between 0 and 3");
    }
  }

  unsigned weight() const {
    return hw;
  }

  bool invert() {
    angle = -angle;
    return true;
  }

  bool mutate() {
    std::normal_distribution<> dAng{0.0, 0.1};
    angle += dAng(gen::rng);
    return true;
  }

  bool merge(const Gene& g) {
    if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)) {
      angle += g.angle;
      return true;
    } else return false;
  }

  bool simplify() {
    /* TODO: parametrize */
    constexpr int N = 5;
    int coeffs[N];
    double a = std::fmod(angle, 2*internal::pi) / internal::pi;
    for(int i = 0; i < N; i++) {
      coeffs[i] = std::floor(a);
      a = 1/(a - coeffs[i]);
    }
    int t;
    for(t = 1; t < N; t++)
      if(std::abs(coeffs[t]) > 20)
        break;
    // 1 ≤ t ≤ N
    if(t == N) // no simplification
      return false;
    // 1 ≤ t < N
    angle = coeffs[--t];
    // 0 ≤ t < N-1
#pragma GCC diagnostic ignored "-Warray-bounds"
    while(t--) // t = 0 breaks loop
      // 0 ≤ t
      angle = coeffs[t] + 1/angle;
#pragma GCC diagnostic pop
    angle *= internal::pi;
    return true;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << g.names[g.op] << g.target()
      << '(' << g.angle / internal::pi << "π)";
    if(g.ixs.size()) {
      os << '[';
      for(auto ctrl : g.ixs)
        os << ctrl;
      os << ']';
    }
    return os;
  }

private:

  NOINLINE Gene(unsigned op_, double angle_, unsigned tgt_, unsigned control_enc):
      op(op_), angle(angle_), tgt(tgt_), hw(0) {
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
    double error{0};
    size_t dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim), in, out;
    for(size_t i = 0; i < dim; i++) {
      psi.fill(0);
      psi[i] = 1;
      in = psi;
      out = arma::fft(psi) / sqrt(dim);
      //error += 1 - std::pow(std::abs(arma::cdot(out, sim(psi))), 2);
      error += 1 - std::real(arma::cdot(out, sim(psi)));
    }
    return error;
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    arma::cx_vec psi(arma::uword(1) << Config::nBit);
    os << '\n';
    size_t dim = arma::uword(1) << Config::nBit;
    for(size_t i = 0; i < dim; i++) {
      psi.fill(0);
      psi[i] = 1;
      //sim(psi).st().raw_print(os);
      for(auto& p : sim(psi))
        os << std::abs(p)*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(p)/3.14159 << "π " << std::noshowpos;
      os << '\n';
    }
    os << '\n';
    return os.str();
  }

private:

  arma::cx_vec sim(arma::cx_vec& psi) const {
    for(const Gene& g : gt) {
      /* control-gate (QIClib) */
      psi = qic::apply_ctrl(
          psi,            // state
          g.gate(),       // operator
          g.ix_vector(),  // arma::uvec of control systems
          {g.target()});  // arma::uvec of target systems
    }
    return psi;
  }

}; // class Candidate


void init() {
}

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
