namespace QGA {

namespace Gates {

using Tools::Controls;

namespace internal {

struct gate_struct_f {
  Backend::Gate op;
  std::string name;
  int inv;
  int sq;
};

static const std::vector<gate_struct_f> fixed_full {
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

static const std::vector<gate_struct_f> fixed_reduced {
  { Backend::I, "I", 0, 0 },
  { Backend::H, "H", 0, -1 },
  { Backend::T, "T", +1, 0 },
  { Backend::Ti, "Ti", -1, 0 },
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
    Tools::controls_distribution<Controls::ANY>
      dCtrl{Config::nBit, tgt, Config::pControl};
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

  void invert(Pointer& self) const override {
    int dIx = (*gates)[op].inv;
    if(dIx != 0)
      self = std::make_shared<Fixed>(op + dIx, tgt, ixs);
  }

  bool invite(Pointer& first, Pointer& second) const override {
    return first->merge(first, second, *this);
  }

  bool merge(Pointer& first, Pointer&, const Fixed& g) const override {
    // G * G = square(G) if also among our operations
    if(g.op == op && g.tgt == tgt && g.ixs == ixs && (*gates)[op].sq != 0) {
      first = std::make_shared<Fixed>(op + (*gates)[op].sq, tgt, ixs);
      return true;
    } else
      return false;
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

  Fixed(size_t op_, unsigned tgt_, std::vector<bool> ctrl):
      op(op_), tgt(tgt_), ixs(ctrl) { }

  Fixed(size_t op_, unsigned tgt_, const Backend::Controls& ixs_):
    op(op_), tgt(tgt_), ixs(ixs_) { }

}; // class Fixed

} // namespace internal


template<class GateBase>
using Fixed = internal::Fixed<GateBase,
        &internal::fixed_full, Controls::NONE>;

template<class GateBase>
using FixedRed = internal::Fixed<GateBase,
        &internal::fixed_reduced, Controls::NONE>;

} // namespace Gates

} // namespace QGA
