// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

using QGA::Backend::State;

using Gene = QGA::Gene<QGA::Gates::FixedRed>;

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

#endif // !defined QGA_PROBLEM_HPP
