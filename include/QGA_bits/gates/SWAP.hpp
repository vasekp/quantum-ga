namespace QGA {

namespace Gates {

struct SWAP {

template<class GateBase>
class SWAPTemp : public GateBase {

  using typename GateBase::Pointer;
  using Ctx = typename GateBase::Context;

public:

  // construct a random gate
  SWAPTemp():
    s1(std::uniform_int_distribution<unsigned>{0, Config::nBit - 2}
        (gen::rng)), // NB the -2 is important!
    s2(std::uniform_int_distribution<unsigned>{0, Config::nBit - 2}
        (gen::rng)),
    ixs(), odd(true)
  {
    // ensure that the two systems are different and in ascending order
    if(s2 < s1)
      std::swap(s1, s2);
    s2 += s2 >= s1;
    ixs = Backend::Controls::swap(s1, s2);
  }

  // construct using parameters
  SWAPTemp(unsigned s1_, unsigned s2_, const Backend::Controls& ixs_, bool odd_):
      s1(s1_), s2(s2_), ixs(ixs_), odd(odd_) { }

  SWAPTemp(unsigned s1_, unsigned s2_): s1(s1_), s2(s2_),
    ixs(std::move(Backend::Controls::swap(s1, s2))), odd(true) { }

  SWAPTemp(bool): s1(), s2(), ixs(), odd(false) { }

  Backend::State applyTo(const Backend::State& psi, const Ctx*) const override {
    return odd ? psi.swap(ixs) : psi;
  }

  bool isTrivial() const override {
    // SWAP^(2k) = SWAP^0 = identity
    return !odd;
  }

  Pointer mutate(const Pointer&) const override {
    return std::make_shared<SWAPTemp>();
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const SWAPTemp& g) const override {
    if(g.s1 == s1 && g.s2 == s2)
      return std::make_shared<SWAPTemp>(s1, s2, ixs, odd ^ g.odd);
    else
      return {};
  }

  std::ostream& write(std::ostream& os) const override {
    if(!odd)
      return os << "[Id]";
    else
      return os << "SWAP" << s1 + 1 << s2 + 1;
  }

  void printOn(QGA::internal::CircuitPrinter& p) const override {
    p.addGates({
        {{s1}, "X"},
        {{s2}, "X"},
    });
  }

  static Pointer read(const std::string& s) {
    std::regex re{"(\\[Id\\])|SWAP(\\d)(\\d)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    if(m[1].matched)
      return std::make_shared<SWAPTemp>(false);
    unsigned s1 = m[2].str()[0] - '1',
             s2 = m[3].str()[0] - '1';
    if(s1 >= Config::nBit || s2 >= Config::nBit || s2 == s1)
      return {};
    return std::make_shared<SWAPTemp>(s1, s2);
  }

private:

  unsigned s1, s2;
  Backend::Controls ixs;
  bool odd;  // parity of the power

}; // class SWAP::SWAPTemp<GateBase>

template<class GateBase>
using Template = SWAPTemp<GateBase>;

}; // struct SWAP

} // namespace Gates

} // namespace QGA
