// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;

static const std::vector<QGA::Gates::gate_struct_f> reduced_set {
  { &QGA::Backend::I, "I", 0, 0 },
  { &QGA::Backend::H, "H", 0, -1 },
  { &QGA::Backend::T, "T", +1, 0 },
  { &QGA::Backend::Ti, "Ti", -1, 0 },
};

using Gene = QGA::Gene<
               QGA::Gates::Fixed
                 ::WithControls<QGA::Controls::ANY>
                 ::WithGates<&reduced_set>
             >;

const State out{3};


class Candidate : public QGA::CandidateBase<Candidate, Gene, double, unsigned>
{

  using Base = QGA::CandidateBase<Candidate, Gene, double, unsigned>;

public:

  using Base::Base;

  Base::FitnessMain fitness_main() const {
    return {
      this->trimError(1 - std::abs(State::overlap(out, sim()))), // error
      this->controls() // total number of control qubits
    };
  }

  std::ostream& print_full(std::ostream& os) const {
    return os << sim();
  }

private:

  State sim() const {
    State psi{0};
    for(const auto& g : genotype())
      psi = g->applyTo(psi);
    return psi;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
