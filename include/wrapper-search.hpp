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

// An assymetric version of zrot
arma::cx_mat22 phase(double a) {
  return { 1, 0, 0, std::exp(i*a) };
}


struct Gate {
  arma::cx_mat22(*fn)(double);
  std::string name;
  bool ctrl;
};

Gate gates[] = {
  {xrot, "X", false},
//{yrot, "Y", true},
//{zrot, "Z", true}
  {phase, "P", true}
};

constexpr size_t gate_count = std::extent<decltype(gates)>::value;

} // namespace internal


class Gene {

  bool is_oracle;
  size_t op;
  double angle;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;
  arma::cx_mat22 mat;

public:

  static Gene getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of being an oracle
    std::bernoulli_distribution dOracle{0.1};
    if(dOracle(gen::rng))
      return {true};
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{0, internal::gate_count - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*internal::pi, 0.5*internal::pi};
    return {dOp(gen::rng), dAng(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)};
  }

  arma::cx_vec apply(const arma::cx_vec& psi, unsigned mark) const {
    if(is_oracle) {
      // oracle
      arma::cx_vec ret{psi};
      ret[mark] = -ret[mark];
      return ret;
    } else {
      // control-gate (QIClib)
      return qic::apply_ctrl(
          psi,    // state
          mat,    // operator
          ixs,    // arma::uvec of control systems
          {tgt}); // arma::uvec of target systems
    }
  }

  bool isOracle() const {
    return is_oracle;
  }

  unsigned weight() const {
    return hw;
  }

  bool invert() {
    if(is_oracle)
      return false;
    angle = -angle;
    update();
    return true;
  }

  void mutate() {
    if(is_oracle)
      return;
    std::normal_distribution<> dAng{0.0, 0.1};
    angle += dAng(gen::rng);
    update();
  }

  bool merge(const Gene& g) {
    if(is_oracle || g.is_oracle)
      return false;
    if(angle == 0) {
      // op1 = identity
      op = g.op;
      tgt = g.tgt;
      ixs = g.ixs;
      hw = g.hw;
      angle = g.angle;
      update();
      return true;
    } else if(g.angle == 0) {
      // op2 = identity
      return true;
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)) {
      angle += g.angle;
      update();
      return true;
    } else
      return false;
  }

  bool simplify() {
    if(is_oracle)
      return false;
    angle = rationalize(std::fmod(angle / internal::pi, 2.0)) * internal::pi;
    update();
    return true;
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    if(g.is_oracle)
      return os << "Oracle";
    os << internal::gates[g.op].name << g.tgt;
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

  NOINLINE Gene(size_t op_, double angle_, unsigned tgt_, unsigned control_enc):
      is_oracle(false), op(op_), angle(angle_), tgt(tgt_), hw(0) {
    if(internal::gates[op].ctrl) {
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

  Gene(bool oracle): is_oracle(oracle), op(0), angle(0), tgt(0), hw(0) {
    if(!oracle)
      throw std::logic_error("Gene(bool) called with false");
    ixs.clear();
  }

  double rationalize(double x) {
    double a = std::abs(x);
    constexpr unsigned N = 8;
    double coeffs[N];
    unsigned t;
    for(t = 0; t < N; t++) {
      coeffs[t] = std::floor(a);
      if(coeffs[t] > 100) {
        coeffs[t++] = 100;
        break;
      }
      a = 1/(a - coeffs[t]);
    }
    std::discrete_distribution<unsigned> dStop(&coeffs[1], &coeffs[t]);
    unsigned cut = dStop(gen::rng) + 1;
    if(cut == t)
      return x;
    a = coeffs[--cut];
    while(cut > 0)
      a = coeffs[--cut] + 1/a;
    return x < 0 ? -a : a;
  }

  void update() {
    if(op > internal::gate_count - 1)
      throw std::logic_error("gate must be between 0 and 2");
    mat = internal::gates[op].fn(angle);
  }

}; // class Gene


struct Fitness {

  double error;
  size_t length;
  size_t ocalls;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{'
       << f.error << ','
       << f.length << ','
       << f.ocalls << '}';
  }

  friend NOINLINE bool operator< (const Fitness& a, const Fitness& b) {
    return a.error < b.error || (a.error == b.error && a.ocalls < b.ocalls);
  }

  friend NOINLINE bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error
        && a.length <= b.length
        && a.ocalls <= b.ocalls
        && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error
        && a.length == b.length
        && a.ocalls == b.ocalls;
  }

}; // struct Fitness


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  NOINLINE Fitness fitness() const {
    size_t ocalls = 0;
    for(const auto& g : gt)
      if(g.isOracle())
        ocalls++;
    QGA::counter.hit();
    return {error(), gt.size(), ocalls};
  }

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    arma::uword dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim), in, out(dim);
    for(unsigned mark = 0; mark < dim; mark++) {
      psi.fill(0);
      psi[0] = 1;
      in = psi;
      out.fill(0);
      out[mark] = 1;
      error += std::max(1 -
          std::pow(std::abs(arma::cdot(out, sim(psi, mark))), 2), 0.0);
    }
    error /= dim;
    return (unsigned long)(error * (1UL<<24)) / (double)(1UL<<24);
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    arma::uword dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim);
    for(arma::uword mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      psi.fill(0);
      psi[0] = 1;
      for(auto& p : sim(psi, mark))
        os << p << ' ';
      /*for(auto& g : gt) {
        os << psi << '\n';
        psi = g.apply(psi, mark);
      }
      for(auto& p : psi)
        os << p << ' ';*/
      os << '\n';
    }
    return os.str();
  }

private:

  arma::cx_vec sim(arma::cx_vec& psi, unsigned mark) const {
    for(auto& g : gt)
      psi = g.apply(psi, mark);
    return psi;
  }

}; // class Candidate


void init() {
}

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
