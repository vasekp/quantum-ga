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

  char names[3] = { 'X', 'Y', 'Z' };

  unsigned op;
  double angle;
  double gphase;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<unsigned> dOp{0, 2};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{0, 2.0*internal::pi};
    return {dOp(gen::rng), dAng(gen::rng), dAng(gen::rng),
      dTgt(gen::rng), dCtrl(gen::rng)};
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
        return internal::xrot(angle) * std::exp(gphase * internal::i);
      case 1:
        return internal::yrot(angle) * std::exp(gphase * internal::i);
      case 2:
        return internal::zrot(angle) * std::exp(gphase * internal::i);
      default:
        throw std::logic_error("gate must be between 0 and 2");
    }
  }

  unsigned weight() const {
    return hw;
  }

  double phase() const {
    return gphase;
  }

  bool invert() {
    angle = -angle;
    return true;
  }

  bool mutate() {
    std::normal_distribution<> dAng{0.0, 0.1};
    std::bernoulli_distribution dWhich{};
    (dWhich(gen::rng) ? angle : gphase) += dAng(gen::rng);
    return true;
  }

  bool merge(const Gene& g) {
    if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)) {
      angle += g.angle;
      gphase += g.gphase;
      return true;
    } else return false;
  }

  double rationalize(double x) {
    double a = x;
    /* TODO: parametrize */
    constexpr int N = 10;
    int coeffs[N];
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
      return x;
    // 1 ≤ t < N
    a = coeffs[--t];
    // 0 ≤ t < N-1
#pragma GCC diagnostic ignored "-Warray-bounds"
    while(t--) // t = 0 breaks loop
      // 0 ≤ t
      a = coeffs[t] + 1/a;
#pragma GCC diagnostic pop
    return a;
  }

  bool simplify() {
    angle = rationalize(std::fmod(angle, 2*internal::pi) / internal::pi)
      * internal::pi;
    gphase = rationalize(std::fmod(gphase, 2*internal::pi) / internal::pi)
      * internal::pi;
    return true;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << g.names[g.op] << g.target();
    if(g.ixs.size()) {
      os << '[';
      for(auto ctrl : g.ixs)
        os << ctrl;
      os << ']';
    }
    os << '(' << g.angle / internal::pi << "π)";
    return os;
  }

private:

  NOINLINE Gene(unsigned op_, double angle_, double phase_,
    unsigned tgt_, unsigned control_enc):
      op(op_), angle(angle_), gphase(phase_), tgt(tgt_), hw(0) {
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
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    size_t dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim), in, out;
    for(size_t i = 0; i < dim; i++) {
      psi.fill(0);
      psi[i] = 1;
      in = psi;
      out = arma::fft(psi) / sqrt(dim);
      error += std::max(1 - std::real(arma::cdot(out, sim(psi))), 0.0);
    }
    return error < 1E-6 ? 0 : error;
  }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    os << (Base&)c;
    double phase = 0;
    for(auto& g : c.gt)
      phase += g.phase();
    os << "φ " << phase / internal::pi << "π";
    return os;
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
      for(auto& p : sim(psi))
        os << std::abs(p)*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(p)/internal::pi << "π " << std::noshowpos;
      os << '\n';
    }
    return os.str();
  }

private:

  arma::cx_vec sim(arma::cx_vec& psi) const {
    for(auto& g : gt) {
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
