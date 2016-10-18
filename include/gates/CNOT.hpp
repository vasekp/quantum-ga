namespace QGA {

namespace Gates {

namespace internal {

template<class GateBase, Controls cc>
class CNOT : public GateBase {

  unsigned tgt;
  Backend::Controls ixs;
  bool odd;  // parity of the power

  using typename GateBase::Pointer;

public:

  static Pointer getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    // distribution of controls
    unsigned tgt = dTgt(gen::rng);
    controls_distribution<cc> dCtrl{Config::nBit, tgt, Config::pControl};
    return std::make_shared<CNOT>(tgt, dCtrl(gen::rng));
  }

  Backend::State applyTo(const Backend::State& psi) const override {
    return odd ? psi.apply_ctrl(Backend::X, ixs, tgt) : psi;
  }

  bool isTrivial() const override {
    // CNOT^(2k) = CNOT^0 = identity
    return !odd;
  }

  unsigned complexity() const override {
    return ixs.size() * ixs.size();
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const CNOT& g) const override {
    if(g.tgt == tgt && g.ixs == ixs)
      return std::make_shared<CNOT>(tgt, ixs, odd ^ g.odd);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    if(!odd)
      return os << "[Id]";
    os << "NOT" << tgt + 1;
    if(ixs.size()) {
      os << '[';
      for(auto ctrl : ixs.as_vector())
        os << ctrl + 1;
      os << ']';
    }
    return os;
  }

  static Pointer read(const std::string& s) {
    std::regex re{"(\\[Id\\])|NOT(\\d)(\\[(\\d+)\\])?"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    if(m[1].matched)
      return std::make_shared<CNOT>(0, Backend::Controls{}, false);
    unsigned tgt = m[2].str()[0] - '1';
    if(tgt >= Config::nBit)
      return {};
    std::vector<bool> ctrl(Config::nBit, false);
    if(m[3].matched)
      for(const char& c : m[4].str()) {
        size_t pos = c - '1';
        if(pos >= 0 && pos < Config::nBit && pos != tgt)
          ctrl[c - '1'] = true;
      }
    return std::make_shared<CNOT>(tgt, Backend::Controls{ctrl});
  }

  CNOT(unsigned tgt_, const Backend::Controls& ixs_, bool odd_ = true):
      tgt(tgt_), ixs(ixs_), odd(odd_) { }

}; // class CNOT

} // namespace internal


template<Controls cc = Controls::ONE>
struct CNOT {

  template<class GateBase>
  using Template = internal::CNOT<GateBase, cc>;

}; // struct CNOT

} // namespace Gates

} // namespace QGA
