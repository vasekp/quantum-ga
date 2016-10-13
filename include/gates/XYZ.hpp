#ifndef GATE_XYZ_HPP
#define GATE_XYZ_HPP

namespace QGA {

struct gate_struct {
  Backend::Gate(*fn)(double);
  char name;
};

static const std::vector<gate_struct> gates {
  {Backend::xrot, 'X'},
  {Backend::yrot, 'Y'},
  {Backend::zrot, 'Z'}
};

using Tools::Controls;


template<class GateBase, const std::vector<gate_struct>* gates, Controls cc>
class Param : public GateBase {

  size_t op;
  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{0, gates->size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    Tools::controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<Param>(gates->size() == 1 ? 0 : dOp(gen::rng),
        tgt, dAng(gen::rng), dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle == 0;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  void invert(Pointer& self) const override {
    self = std::make_shared<Param>(op, tgt, -angle, ixs);
  }

  void mutate(Pointer& self) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<Param>(op, tgt, angle + dAng(gen::rng), ixs);
  }

  void simplify(Pointer& self) const override {
    self = std::make_shared<Param>(op, tgt,
        Tools::rationalize_angle(angle), ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const Param& g) const override {
    if(g.op == op && g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<Param>(op, tgt, angle + g.angle, ixs);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << (*gates)[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    os << '(' << angle / Const::pi << "π)";
    return os;
  }

  Param(size_t op_, unsigned tgt_, double angle_,
      std::vector<bool> ctrl):
    op(op_), tgt(tgt_), angle(angle_), ixs(ctrl), mat((*gates)[op].fn(angle))
  { }

  Param(size_t op_, unsigned tgt_, double angle_,
      const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), angle(angle_), ixs(ixs_), mat((*gates)[op].fn(angle))
  { }

}; // class Param


template<class GateBase>
using CnXYZ = Param<GateBase, &gates, Controls::ANY>;

} // namespace QGA

#endif // !defined GATE_XYZ_HPP
