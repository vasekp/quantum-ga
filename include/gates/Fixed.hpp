namespace QGA {

namespace Gates {

struct gate_struct_f {
  Backend::Gate op;
  std::string name;
  int inv;
  int sq;
};


namespace internal {

static const std::vector<gate_struct_f> gates_fixed {
  { Backend::I, "I", 0, 0 },
  { Backend::H, "H", 0, -1 },
  { Backend::X, "X", 0, -2 },
  { Backend::Y, "Y", 0, -3 },
  { Backend::Z, "Z", 0, -4 },
  { Backend::T, "T", +1, +2 },
  { Backend::Ti, "Ti", -1, +2 },
  { Backend::S, "S", +1, -3 },
  { Backend::Si, "Si", -1, -4 }
};

template<class GateBase, const std::vector<gate_struct_f>* gates, Controls cc>
class Fixed : public GateBase {

  size_t op;
  unsigned tgt;
  Backend::Controls ixs;

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    /* Distributions: cheap and safer in MT environment this way */
    // distribution of possible gates (except of identity)
    std::uniform_int_distribution<size_t> dOp{1, gates->size() - 1};
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    return std::make_shared<Fixed>(
        dOp(gen::rng), tgt, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return psi.apply_ctrl((*gates)[op].op, ixs, tgt);
  }

  bool isTrivial() const override {
    return op == 0;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  Pointer invert(const Pointer& self) const override {
    int dIx = (*gates)[op].inv;
    if(dIx != 0)
      return std::make_shared<Fixed>(op + dIx, tgt, ixs);
    else
      return self;
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Fixed& g) const override {
    // G * G = square(G) if also among our operations
    if(g.op == op && g.tgt == tgt && g.ixs == ixs && (*gates)[op].sq != 0)
      return std::make_shared<Fixed>(op + (*gates)[op].sq, tgt, ixs);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    os << (*gates)[op].name << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto& ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  static Pointer read(const std::string& s) {
    std::string reS{};
    for(const gate_struct_f& g : *gates)
      reS = reS + "|(" + g.name + ")";
    std::regex re{"(?:" + reS.substr(1) + ")(\\d)(\\[(\\d+)\\])?"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    size_t num = gates->size();
    size_t op;
    for(op = 0; op < num; op++)
      if(m[op + 1].matched)
        break;
    unsigned tgt = m[num + 1].str()[0] - '1';
    if(tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(m[num + 2].matched)
      for(const char& c : m[num + 3].str()) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[pos] = true;
      }
    return std::make_shared<Fixed>(op, tgt, Backend::Controls{ctrl});
  }

  Fixed(size_t op_, unsigned tgt_, const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), ixs(ixs_) { }

}; // class Fixed

} // namespace internal


template<Controls cc = Controls::NONE,
  const std::vector<gate_struct_f>* gates = internal::gates_fixed>
struct Fixed {

  template<class GateBase>
  using Template = internal::Fixed<GateBase, gates, cc>;

}; // struct Fixed

} // namespace Gates

} // namespace QGA
