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


struct Gate {
  arma::cx_mat22(*fn)(double);
  std::string name;
  bool ctrl;
};

Gate gates[] = {
  {xrot, "X", true},
  {yrot, "Y", true},
  {zrot, "Z", true}
};

constexpr size_t gate_count = std::extent<decltype(gates)>::value;

} // namespace internal


class Gene {

  unsigned op;
  double angle;
  double gphase;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;
  arma::cx_mat22 mat;

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates
    std::uniform_int_distribution<unsigned> dOp{1, internal::gate_count};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*internal::pi, 0.5*internal::pi};
    return {dOp(gen::rng), dAng(gen::rng), dAng(gen::rng),
      dTgt(gen::rng), dCtrl(gen::rng)};
  }

  const arma::uvec& ix_vector() const {
    return ixs;
  }

  unsigned target() const {
    return tgt;
  }

  const arma::cx_mat22& gate() const {
    return mat;
  }

  unsigned weight() const {
    return hw;
  }

  double phase() const {
    return gphase;
  }

  bool invert() {
    angle = -angle;
    update();
    return true;
  }

  void mutate() {
    std::normal_distribution<> dAng{0.0, 0.1};
    std::bernoulli_distribution dWhich{};
    (dWhich(gen::rng) ? angle : gphase) += dAng(gen::rng);
    update();
  }

  bool merge(const Gene& g) {
    if(angle == 0) {
      // op1 = (phase*)identity
      op = g.op;
      tgt = g.tgt;
      ixs = g.ixs;
      hw = g.hw;
      angle = g.angle;
      gphase += g.gphase;
      update();
      return true;
    } else if(g.angle == 0) {
      // op2 = (phase*)identity
      gphase += g.gphase;
      update();
      return true;
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)) {
      angle += g.angle;
      gphase += g.gphase;
      update();
      return true;
    } else
      return false;
  }

  double rationalize(double x) {
    double a = std::abs(x);
    constexpr int N = 8;
    double coeffs[N];
    int t;
    for(t = 0; t < N; t++) {
      coeffs[t] = std::floor(a);
      if(coeffs[t] > 100) {
        coeffs[t++] = 100;
        break;
      }
      a = 1/(a - coeffs[t]);
    }
    std::discrete_distribution<> dStop(&coeffs[1], &coeffs[t]);
    int cut = dStop(gen::rng) + 1;
    if(cut == t)
      return x;
    a = coeffs[--cut];
    while(cut > 0)
      a = coeffs[--cut] + 1/a;
    return x < 0 ? -a : a;
  }

  bool simplify() {
    angle = rationalize(std::fmod(angle / internal::pi, 2.0)) * internal::pi;
    gphase = rationalize(std::fmod(gphase / internal::pi, 2.0)) * internal::pi;
    update();
    return true;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << internal::gates[g.op - 1].name << g.target();
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
    if(internal::gates[op - 1].ctrl) {
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
    } else
      ixs.clear();
    update();
  }

  void update() {
    if(op < 1 || op > internal::gate_count)
      throw std::logic_error("gate must be between 0 and 2");
    mat = std::exp(gphase * internal::i) * internal::gates[op - 1].fn(angle);
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
    return error < 1E-8 ? 0 : error;
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
