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

using State = arma::cx_vec;


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
  {xrot, "X", true},
  {yrot, "Y", true},
  {zrot, "Z", true}
};

constexpr size_t gate_count = std::extent<decltype(gates)>::value;

} // namespace internal


template<class GeneBase>
class XYZGene : public GeneBase {

  size_t op;
  double angle;
  double gphase;
  unsigned tgt;
  unsigned hw;
  arma::uvec ixs;
  arma::cx_mat22 mat;

  using SP = std::shared_ptr<GeneBase>;

public:

  static SP getNew() {
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{0, internal::gate_count - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{1, Config::nBit};
    // distribution of controls
    std::uniform_int_distribution<unsigned> dCtrl{};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*internal::pi, 0.5*internal::pi};
    return std::make_shared<XYZGene>(
        dOp(gen::rng), dAng(gen::rng), dAng(gen::rng),
        dTgt(gen::rng), dCtrl(gen::rng));
  }

  State apply(const State& psi) const override {
    // control-gate (QIClib)
    return qic::apply_ctrl(
        psi,    // state
        mat,    // operator
        ixs,    // arma::uvec of control systems
        {tgt}); // arma::uvec of target systems
  }

  unsigned complexity() const override {
    return hw * hw;
  }

  /*double phase() const {
    return gphase;
  }*/

  SP invert(const SP&) override {
    return std::make_shared<XYZGene>(op, -angle, gphase, hw, tgt, ixs);
  }

  SP mutate(const SP&) override {
    std::normal_distribution<> dAng{0.0, 0.1};
    std::bernoulli_distribution dWhich{};
    bool modAngle = dWhich(gen::rng);
    return std::make_shared<XYZGene>(
        op,
        modAngle ? angle + dAng(gen::rng) : angle,
        modAngle ? gphase : gphase + dAng(gen::rng),
        hw, tgt, ixs);
  }

  SP simplify(const SP&) override {
    return std::make_shared<XYZGene>(op,
        rationalize(std::fmod(angle / internal::pi, 2.0)) * internal::pi,
        rationalize(std::fmod(gphase / internal::pi, 2.0)) * internal::pi,
        hw, tgt, ixs);
  }

  SP invite(const SP& g) const override {
    return g.get()->visit(g, *this);
  }

  SP visit(const SP& self, const XYZGene& g) override {
    if(angle == 0) {
      // op1 = (phase*)identity
      return std::make_shared<XYZGene>(
          g.op, g.angle, g.gphase + gphase, g.hw, g.tgt, g.ixs);
    } else if(g.angle == 0) {
      // op2 = (phase*)identity
      return std::make_shared<XYZGene>(
          op, angle, gphase + g.gphase, hw, tgt, ixs);
    } else if(g.op == op
        && g.tgt == tgt
        && g.ixs.size() == ixs.size()
        && arma::all(g.ixs == ixs)) {
      return std::make_shared<XYZGene>(
          op, angle + g.angle, gphase + g.gphase, hw, tgt, ixs);
    } else
      return self;
  }

  std::ostream& write(std::ostream& os) const override {
    os << internal::gates[op].name << tgt;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs)
        os << ctrl;
      os << ']';
    }
    os << '(' << angle / internal::pi << "π)";
    return os;
  }

// Should be private, but constructors are needed my std::make_shared().

  NOINLINE XYZGene(size_t op_, double angle_, double phase_,
    unsigned tgt_, unsigned control_enc):
      op(op_), angle(angle_), gphase(phase_), tgt(tgt_), hw(0) {
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
    computeMat();
  }

  NOINLINE XYZGene(size_t op_, double angle_, double phase_,
    unsigned hw_, unsigned tgt_, arma::uvec ixs_):
      op(op_), angle(angle_), gphase(phase_), tgt(tgt_), hw(hw_), ixs(ixs_) {
    computeMat();
  }

private:

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

  void computeMat() {
    if(op >= internal::gate_count)
      throw std::logic_error("gate must be between 0 and 2");
    mat = std::exp(gphase * internal::i) * internal::gates[op].fn(angle);
  }

}; // class XYZGene

} // namespace Wrapper

#endif // !defined QGA_WRAPPER_HPP
