// allow only one problem
#ifndef QGA_PROBLEM_HPP
#define QGA_PROBLEM_HPP

namespace {

using QGA::Backend::State;

using Gene = QGA::Gene<QGA::Gates::Y, QGA::Gates::CPhase, QGA::Gates::SWAP>;


class Candidate : public QGA::CandidateBase<Candidate, Gene, double, size_t> {

  using Base = QGA::CandidateBase<Candidate, Gene, double, size_t>;

public:

  using Base::Base;

  Base::Fitness fitness() const {
    if(genotype().size() > 1000)
      return {};
    using cxd = std::complex<double>;
    cxd overlapTotal{0};
    unsigned dim = 1 << Config::nBit;
    State psi{};
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = State::fourier(psi);
      cxd overlap = State::overlap(out, sim(psi));
      overlapTotal += overlap;
    }
    double errorAvg = std::max(1.0 - std::abs(overlapTotal / cxd(dim)), 0.0);
    return {
      trimError(errorAvg),
      genotype().size()
    };
  }

  std::ostream& print_full(std::ostream& os) const {
    unsigned dim = 1 << Config::nBit;
    State psi{};
    os << '\n';
    for(unsigned i = 0; i < dim; i++) {
      psi.reset(i);
      State out = sim(psi);
      for(unsigned j = 0; j < dim; j++)
        os << std::abs(out[j])*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(out[j]) / QGA::Const::pi << "π "
          << std::noshowpos;
      os << '\n';
    }
    return os;
  }

private:

  State sim(const State& psi) const {
    State ret{psi};
    for(const auto& g : genotype())
      ret = g->applyTo(ret);
    return ret;
  }

}; // class Candidate

} // anonymous namespace
#endif // !defined QGA_PROBLEM_HPP
