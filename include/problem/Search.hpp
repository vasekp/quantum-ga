// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;


/* The Context to be used in Gene below, holding a mark for the Oracle */
struct Context {
  unsigned mark;
};


/* The oracle gate. */

struct Oracle {

template<class GateBase>
class Inner : public GateBase {

  using typename GateBase::Pointer;
  using typename GateBase::Context;

  bool odd;  // parity of the power

public:

  Inner(bool odd_ = true): odd(odd_) { }

  State applyTo(const State& psi, const Context* pMark) const override {
    unsigned mark = pMark->mark;
    State ret{psi};
    if(odd)
      ret[mark] = -ret[mark];
    return ret;
  }

  bool isTrivial() const override {
    // oracle^(2k) = oracle^0 = identity
    return !odd;
  }

  void hit(typename GateBase::Counter& c) const {
    c.hit(this);
  }

  Pointer invite(const Pointer& first) const override {
    return first->merge(*this);
  }

  Pointer merge(const Inner& g) const override {
    // oracle * oracle = oracle^2 → true ^ true = false
    return std::make_shared<Inner>(odd ^ g.odd);
  }

  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  static Pointer read(const std::string& s) {
    std::regex re{"\\[Id\\]|(Oracle)"};
    std::smatch m{};
    if(!std::regex_match(s, m, re))
      return {};
    return std::make_shared<Inner>(m[1].matched);
  }

}; // class Oracle::Inner

template<class GateBase>
using Template = Inner<GateBase>;

}; // struct Oracle


using Gene = typename QGA::Gene<Oracle, QGA::Gates::X<>, QGA::Gates::CPhase>
                ::WithContext<Context>;

class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      State out{mark};
      error += std::max(1 -
          std::pow(std::abs(State::overlap(out, sim(psi, mark))), 2), 0.0);
    }
    return error / dim;
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      os << sim(psi, mark);
    }
    return os.str();
  }

private:

  State sim(const State& psi, unsigned mark) const {
    State ret{psi};
    Context c{mark};
    for(const auto& g : gt)
      ret = ret.apply(g, &c);
    return ret;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
