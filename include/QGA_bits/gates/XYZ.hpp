namespace QGA {

namespace Gates {

struct gate_struct_p {
  Backend::Gate(*fn)(double);
  std::string name;
};


namespace internal {

static const std::vector<gate_struct_p> gates_param_xyz {
  {func::xrot, "X"},
  {func::yrot, "Y"},
  {func::zrot, "Z"}
};

static const std::vector<gate_struct_p> gates_param_x {
  {func::xrot, "X"},
};

static const std::vector<gate_struct_p> gates_param_y {
  {func::yrot, "Y"},
};

static const std::vector<gate_struct_p> gates_param_z {
  {func::zrot, "Z"},
};

/* Do not use: does not represent a 1-parametric group!
 * This interferes with the isTrivial and merge mechanisms of Param.
 * An equivalent of R(φ) is P(π) Y(φ). */
/*static const std::vector<gate_struct_p> gates_param_r {
  {func::rrot, "R"},
};*/


template<Controls cc, const std::vector<gate_struct_p>* gates>
struct Param {

template<class GateBase>
class ParamTemp : public GateBase {

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  ParamTemp():
    op(std::uniform_int_distribution<size_t>{0, gates->size() - 1}(gen::rng)),
    tgt(std::uniform_int_distribution<unsigned>{0, Config::nBit - 1}
        (gen::rng)),
    angle(angle_distribution<>{}(gen::rng)),
    ixs(controls_distribution<cc>{Config::nBit, tgt, Config::pControl}
        (gen::rng)),
    mat((*gates)[op].fn(angle))
  { }

  // construct using parameters
  ParamTemp(size_t op_, unsigned tgt_, double angle_,
      const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), angle(angle_), ixs(ixs_), mat((*gates)[op].fn(angle))
  { }

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
    return std::make_shared<ParamTemp>(op, tgt, -angle, ixs);
  }

  Pointer mutate(const Pointer&) const override {
    std::bernoulli_distribution dCont{};
    if(dCont(gen::rng)) {
      // Continuous
      angle_distribution<true> dAng{};
      return std::make_shared<ParamTemp>(op, tgt, angle + dAng(gen::rng), ixs);
    } else
      // Discrete
      return std::make_shared<ParamTemp>();
  }

  Pointer simplify(const Pointer&) const override {
    return std::make_shared<ParamTemp>(op, tgt, rationalize_angle(angle), ixs);
  }

  Pointer swapQubits(const Pointer&, unsigned s1, unsigned s2) const override {
    return std::make_shared<ParamTemp>(
        op,
        tgt == s1 ? s2 : tgt == s2 ? s1 : tgt,
        angle,
        Backend::Controls::swapQubits(ixs, s1, s2));
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const ParamTemp& g) const override {
    if(g.op == op && g.tgt == tgt && g.ixs == ixs)
      return std::make_shared<ParamTemp>(op, tgt, angle + g.angle, ixs);
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

  void printOn(QGA::CircuitPrinter& p) const override {
    p.addControlledGate((*gates)[op].name, tgt, ixs.as_vector());
  }

  static Pointer read(const std::string& s) {
    std::string reS{};
    for(const gate_struct_p& g : *gates)
      reS = reS + "|(" + g.name + ")";
    regex::regex re{"(?:" + reS.substr(1) + ")" +
      "(\\d)(\\[(\\d+)\\])?\\((-?[0-9.]+)(?:π)?\\)"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    size_t num = gates->size();
    size_t op;
    for(op = 0; op < num; op++)
      if(ms.matched(op + 1))
        break;
    unsigned tgt = ms.match(num + 1)[0] - '1';
    if(tgt < 0 || tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(ms.matched(num + 2))
      for(const char& c : ms.match(num + 3)) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[pos] = true;
      }
    double angle = std::stod(ms.match(num + 4)) * Const::pi;
    return std::make_shared<ParamTemp>(op, tgt, angle, Backend::Controls{ctrl});
  }

private:

  size_t op;
  unsigned tgt;
  double angle;
  Backend::Controls ixs;
  Backend::Gate mat;

}; // class Param<Controls, Gates>::ParamTemp<GateBase>

template<class GateBase>
using Template = ParamTemp<GateBase>;

template<Controls cc_>
using WithControls = Param<cc_, gates>;

template<const std::vector<gate_struct_p>* gates_>
using WithGates = Param<cc, gates_>;

}; // struct Param<Controls, Gates>

} // namespace internal

using XYZ = internal::Param<Controls::NONE, &internal::gates_param_xyz>;
using X = internal::Param<Controls::NONE, &internal::gates_param_x>;
using Y = internal::Param<Controls::NONE, &internal::gates_param_y>;
using Z = internal::Param<Controls::NONE, &internal::gates_param_z>;

} // namespace Gates

} // namespace QGA
