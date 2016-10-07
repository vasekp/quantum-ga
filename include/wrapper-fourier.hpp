#include "wrapper-common.hpp"
#include "GeneBase.hpp"

namespace Wrapper {

class Gene : public QGA::GeneBase<Gene, XYZGene> {

public:

  static std::shared_ptr<Gene> getNew() {
    return Wrapper::XYZGene<Gene>::getNew();
  }

}; // class Gene


class Candidate : public QGA::CandidateBase<Candidate, Gene> {

  using Base = QGA::CandidateBase<Candidate, Gene>;

public:

  using Base::Base;

  double error() const {
    if(gt.size() > 1000)
      return INFINITY;
    double error{0};
    arma::uword dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim), in, out;
    for(arma::uword i = 0; i < dim; i++) {
      psi.fill(0);
      psi[i] = 1;
      in = psi;
      out = arma::fft(psi) / sqrt(dim);
      error += std::max(1 - std::real(arma::cdot(out, sim(psi))), 0.0);
    }
    error /= dim;
    return error < 1E-8 ? 0 : error;
  }

  /*friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    os << (Base&)c;
    double phase = 0;
    for(auto& g : c.gt)
      phase += g.phase();
    os << "φ " << phase / internal::pi << "π";
    return os;
  }*/

  std::string dump(const std::ostream& ex) const {
    std::ostringstream os{};
    os.flags(ex.flags());
    os.precision(ex.precision());
    os << '\n';
    arma::uword dim = arma::uword(1) << Config::nBit;
    arma::cx_vec psi(dim);
    for(arma::uword i = 0; i < dim; i++) {
      psi.fill(0);
      psi[i] = 1;
      for(auto& p : sim(psi))
        os << std::abs(p)*std::sqrt(dim) << "/√" << dim << "∠"
          << std::showpos << std::arg(p)/internal::pi << "π " << std::noshowpos;
      os << '\n';
    }
    return os.str();
  }

private:

  arma::cx_vec sim(arma::cx_vec& psi) const {
    for(auto& g : gt)
      psi = g->apply(psi);
    return psi;
  }

}; // class Candidate


void init() {
}

} // namespace Wrapper
