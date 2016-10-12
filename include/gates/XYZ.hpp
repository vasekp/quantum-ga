namespace QGA {

struct gate_struct {
  Backend::Gate(*fn)(double);
  char name;
};

std::vector<gate_struct> gates {
  {Backend::xrot, 'X'},
  {Backend::yrot, 'Y'},
  {Backend::zrot, 'Z'}
};


template<class GateBase>
class XYZ : public GateBase {

  size_t op;
  double angle;
  unsigned tgt;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    // distribution of possible gates
    std::uniform_int_distribution<size_t> dOp{0, gates.size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt_ = dTgt(gen::rng);
    Tools::controls_distribution dCtrl{Config::nBit, Config::pControl, tgt_};
    // distribution of angle
    std::uniform_real_distribution<> dAng{-0.5*Const::pi, 0.5*Const::pi};
    return std::make_shared<XYZ>(
        dOp(gen::rng), dAng(gen::rng), tgt_, dCtrl(gen::rng));
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
    self = std::make_shared<XYZ>(op, -angle, tgt, ixs);
  }

  void mutate(Pointer& self) const override {
    std::normal_distribution<> dAng{0.0, 0.1};
    self = std::make_shared<XYZ>(op, angle + dAng(gen::rng), tgt, ixs);
  }

  void simplify(Pointer& self) const override {
    self = std::make_shared<XYZ>(op, Tools::rationalize_angle(angle), tgt, ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const XYZ& g) const override {
    if(g.op == op && g.tgt == tgt && g.ixs == ixs) {
      first = std::make_shared<XYZ>(op, angle + g.angle, tgt, ixs);
      return true;
    } else
      return false;
  }

  std::ostream& write(std::ostream& os) const override {
    os << gates[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    os << '(' << angle / Const::pi << "π)";
    return os;
  }

  NOINLINE XYZ(size_t op_, double angle_, unsigned tgt_,
      std::vector<bool> ctrl):
      op(op_), angle(angle_), tgt(tgt_), ixs(ctrl) {
    mat = gates[op].fn(angle);
  }

  NOINLINE XYZ(size_t op_, double angle_, unsigned tgt_,
      const Backend::Controls& ixs_):
      op(op_), angle(angle_), tgt(tgt_), ixs(ixs_) {
    mat = gates[op].fn(angle);
  }

}; // class XYZ

} // namespace QGA
