namespace QGA {

namespace Gates {

using Tools::Controls;

struct gate_struct_p {
  Backend::Gate(*fn)(double);
  char name;
};


namespace internal {

static const std::vector<gate_struct_p> gates_param {
  {Backend::xrot, 'X'},
  {Backend::yrot, 'Y'},
  {Backend::zrot, 'Z'}
};


template<class GateBase, const std::vector<gate_struct_p>* gates, Controls cc>
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

  Pointer invert(const Pointer&) const override {
    return std::make_shared<Param>(op, tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    return std::make_shared<Param>(op, tgt, angle + dAng(gen::rng), ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<Param>(op, tgt,
        Tools::rationalize_angle(angle), ixs);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Param& g) const override {
    if(g.op == op && g.tgt == tgt && g.ixs == ixs)
      return std::make_shared<Param>(op, tgt, angle + g.angle, ixs);
    else
      return {};
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
      const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), angle(angle_), ixs(ixs_), mat((*gates)[op].fn(angle))
  { }

}; // class Param

} // namespace internal


template<Controls cc = Controls::NONE,
  const std::vector<gate_struct_p>* gates = &internal::gates_param>
struct XYZ {

  template<class GateBase>
  using Template = internal::Param<GateBase, gates, cc>;

}; // struct XYZ

} // namespace Gates

} // namespace QGA
