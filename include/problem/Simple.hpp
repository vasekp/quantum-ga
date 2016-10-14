// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;

static const std::vector<QGA::Gates::gate_struct_f> reduced_set {
  { QGA::Backend::I, "I", 0, 0 },
  { QGA::Backend::H, "H", 0, -1 },
  { QGA::Backend::T, "T", +1, 0 },
  { QGA::Backend::Ti, "Ti", -1, 0 },
};

using Gene = QGA::Gene<QGA::Gates::Fixed<QGA::Tools::Controls::ANY, &reduced_set>>;

const State out{3};


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    return 1 - std::abs(State::overlap(out, sim()));
  }

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << sim();
    return os.str();
  }

private:

  State sim() const {
    State psi{0};
    for(const auto& g : gt)
      psi = psi.apply(g);
    return psi;
  }

}; // class Candidate

} // anonymous namespace

#endif // !defined QGA_PROBLEM_HPP
