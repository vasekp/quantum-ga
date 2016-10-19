namespace QGA {

namespace Gates {

struct gate_struct_p {
  Backend::Gate(*fn)(double);
  char name;
};


namespace internal {

static const std::vector<gate_struct_p> gates_param {
  {func::xrot, 'X'},
  {func::yrot, 'Y'},
  {func::zrot, 'Z'}
};

static const std::vector<gate_struct_p> gates_x {
  {func::xrot, 'X'},
};

static const std::vector<gate_struct_p> gates_y {
  {func::yrot, 'Y'},
};

static const std::vector<gate_struct_p> gates_z {
  {func::zrot, 'Z'},
};

static const std::vector<gate_struct_p> gates_r {
  {func::rrot, 'R'},
};


template<class GateBase, const std::vector<gate_struct_p>* gates, Controls cc>
class Param : public GateBase {

  size_t op;
  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

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
    angle_distribution<> dAng{};
    return std::make_shared<Param>(gates->size() == 1 ? 0 : dOp(gen::rng),
        tgt, dAng(gen::rng), dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return psi.apply_ctrl(mat, ixs, tgt);
  }

  bool isTrivial() const override {
    return angle == 0;
  }

  unsigned controls() const override {
    return ixs.size();
  }

  Pointer invert(const Pointer&) const override {
    return std::make_shared<Param>(op, tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    angle_distribution<true> dAng{};
    return std::make_shared<Param>(op, tgt, angle + dAng(gen::rng), ixs);
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<Param>(op, tgt, rationalize_angle(angle), ixs);
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
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
    unsigned tgt = m[num + 1].str()[0] - '1';
    if(tgt < 0 || tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(m[num + 2].matched)
      for(const char& c : m[num + 3].str()) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[pos] = true;
      }
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

template<Controls cc = Controls::NONE>
struct R {
  template<class GateBase>
  using Template = internal::Param<GateBase, &internal::gates_r, cc>;
};

} // namespace Gates

} // namespace QGA
