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
class OracleTemp : public GateBase {

  using typename GateBase::Pointer;
  using typename GateBase::Context;

  bool odd;  // parity of the power

public:

  OracleTemp(bool odd_ = true): odd(odd_) { }

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

  Pointer swapQubits(const Pointer& self, unsigned, unsigned) const override {
    return self;
  }

  unsigned type() const override {
    return GateBase::Indexer::index(this);
  }

  const OracleTemp* cast(const OracleTemp*) const override {
    return this;
  }

  bool sameType(const GateBase& other) const override {
    const OracleTemp* c = other.cast(this);
    return c != nullptr;
  }

  Pointer getAnother() const override {
    return std::make_shared<OracleTemp>();
  }

  Pointer merge(const GateBase& other) const override {
    if(!sameType(other))
      return {};
    // oracle * oracle = oracle^2 → true ^ true = false
    const OracleTemp* c = other.cast(this);
    return std::make_shared<OracleTemp>(odd ^ c->odd);
  }

  std::ostream& write(std::ostream& os) const override {
    return os << (odd ? "Oracle" : "[Id]");
  }

  void printOn(QGA::CircuitPrinter& p) const override {
    if(odd)
      p.addBarrierGate("U_f");
  }

  static Pointer read(const std::string& s) {
    regex::regex re{"\\[Id\\]|(Oracle)"};
    regex::matches ms{};
    if(!re.match(s, ms))
      return {};
    return std::make_shared<OracleTemp>(ms.matched(1));
  }

}; // class Oracle::OracleTemp<GateBase>

template<class GateBase>
using Template = OracleTemp<GateBase>;

}; // struct Oracle


using Gene = typename QGA::Gene<
                Oracle,
                QGA::Gates::X,
                QGA::Gates::CPhase
              >::WithContext<Context>;


class Candidate : public QGA::CandidateBase<Candidate, Gene, double, unsigned, unsigned> {

  using Base = QGA::CandidateBase<Candidate, Gene, double, unsigned, unsigned>;

public:

  using Base::Base;

  Base::Fitness fitness() const {
    if(genotype().size() > 1000)
      return {};
    double errMax = 0;
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    for(unsigned mark = 0; mark < dim; mark++) {
      State out{mark};
      double error = std::max(1 -
          std::pow(std::abs(State::overlap(out, sim(psi, mark))), 2), 0.0);
      if(error > errMax)
        errMax = error;
    }
    unsigned oracles = 0;
    for(const auto& g : genotype())
      if(g->type() == Gene::gateType<Oracle>())
        oracles++;
    return {
      trimError(errMax),
      genotype().size(),
      oracles
    };
  }

  std::ostream& print_full(std::ostream& os) const {
    unsigned dim = 1 << Config::nBit;
    State psi{0};
    os << '\n';
    for(unsigned mark = 0; mark < dim; mark++) {
      os << mark << ": ";
      os << sim(psi, mark);
    }
    return os;
  }

private:

  State sim(const State& psi, unsigned mark) const {
    State ret{psi};
    Context c{mark};
    for(const auto& g : genotype())
      ret = g->applyTo(ret, &c);
    return ret;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
