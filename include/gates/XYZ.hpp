namespace QGA {

namespace Gates {

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

static const std::vector<gate_struct_p> gates_x {
  {Backend::xrot, 'X'},
};

static const std::vector<gate_struct_p> gates_y {
  {Backend::yrot, 'Y'},
};

static const std::vector<gate_struct_p> gates_z {
  {Backend::zrot, 'Z'},
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
    controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
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
    return std::make_shared<Param>(op, tgt, rationalize_angle(angle), ixs);
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

  static Pointer read(const std::string& s) {
    std::string reS{};
    for(const gate_struct_p& g : *gates)
      reS = reS + "|(" + g.name + ")";
    std::regex re{"(?:" + reS.substr(1) + ")" +
      "(\\d)(\\[(\\d+)\\])?\\((-?[0-9.]+)(?:π)?\\)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    size_t num = gates->size();
    size_t op;
    for(op = 0; op < num; op++)
      if(m[op + 1].matched)
        break;
    // TODO overflow
    unsigned tgt = m[num + 1].str()[0] - '1';
    std::vector<bool> ctrl(Config::nBit, false);
    if(m[num + 2].matched)
      for(const char& c : m[num + 3].str())
        ctrl[c - '1'] = true;
    double angle = std::stod(m[num + 4].str()) * Const::pi;
    return std::make_shared<Param>(op, tgt, angle, Backend::Controls{ctrl});
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
};

template<Controls cc = Controls::NONE>
struct X {
  template<class GateBase>
  using Template = internal::Param<GateBase, &internal::gates_x, cc>;
};

template<Controls cc = Controls::NONE>
struct Y {
  template<class GateBase>
  using Template = internal::Param<GateBase, &internal::gates_y, cc>;
};

template<Controls cc = Controls::NONE>
struct Z {
  template<class GateBase>
  using Template = internal::Param<GateBase, &internal::gates_z, cc>;
};

} // namespace Gates

} // namespace QGA
