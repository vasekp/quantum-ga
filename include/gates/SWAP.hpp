namespace QGA {

namespace Gates {

struct SWAP {

template<class GateBase>
class Inner : public GateBase {

  unsigned s1, s2;
  Backend::Controls ixs;
  bool odd;  // parity of the power

  using typename GateBase::Pointer;
  using typename GateBase::Counter;
  using typename GateBase::Context;

public:

  static Pointer getNew() {
    // distribution of targets
    std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 2};
    // distribution of controls
    unsigned s1 = dTgt(gen::rng),
             s2 = dTgt(gen::rng);
    // ensure that the two subsystem indices are different
    s2 += s2 >= s1;
    return std::make_shared<Inner>(s1, s2);
  }

  Backend::State applyTo(const Backend::State& psi, const Context*) const
  override {
    return odd ? psi.swap(ixs) : psi;
  }

  bool isTrivial() const override {
    // SWAP^(2k) = SWAP^0 = identity
    return !odd;
  }

  unsigned complexity() const override {
    return 0;
  }

  void hit(Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Inner& g) const override {
    if(g.s1 == s1 && g.s2 == s2)
      return std::make_shared<Inner>(s1, s2, ixs, odd ^ g.odd);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    if(!odd)
      return os << "[Id]";
    else
      return os << "SWAP" << s1 + 1 << s2 + 1;
  }

  static Pointer read(const std::string& s) {
    std::regex re{"(\\[Id\\])|SWAP(\\d)(\\d)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    if(m[1].matched)
      return std::make_shared<Inner>(false);
    unsigned s1 = m[2].str()[0] - '1',
             s2 = m[3].str()[0] - '1';
    if(s1 >= Config::nBit || s2 >= Config::nBit || s2 == s1)
      return {};
    return std::make_shared<Inner>(s1, s2);
  }

  Inner(unsigned s1_, unsigned s2_, const Backend::Controls& ixs_, bool odd_):
      s1(s1_), s2(s2_), ixs(ixs_), odd(odd_) { }

  Inner(unsigned s1_, unsigned s2_): s1(s1_), s2(s2_),
    ixs(std::move(Backend::Controls::swap(s1, s2))), odd(true) { }

  Inner(bool): s1(), s2(), ixs(), odd(false) { }

}; // class Inner

template<class GateBase>
using Template = Inner<GateBase>;

}; // struct SWAP

} // namespace Gates

} // namespace QGA
