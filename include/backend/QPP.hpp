// allow only one backend
#ifndef QGA_BACKEND_HPP
#define QGA_BACKEND_HPP

#include "qpp.h"
//#include <unsupported/Eigen/FFT>

namespace QGA {

namespace Backend {

class Gate : public Eigen::Matrix2cd {

  using Base = Eigen::Matrix2cd;
  using cxd = std::complex<double>;

public:

  Gate(const Base& mat) : Base(mat) { }

  Gate(cxd u11, cxd u12, cxd u21, cxd u22) : Base(2, 2) {
    *this << u11, u12, u21, u22;
  }

}; // class Gate


/* Fixed gates */

const Gate I { qpp::gt.Id2 };
const Gate H { qpp::gt.H };
const Gate X { qpp::gt.X };
const Gate Y { qpp::gt.Y };
const Gate Z { qpp::gt.Z };
const Gate T { qpp::gt.T };
const Gate Ti { qpp::gt.T.conjugate() };
const Gate S { qpp::gt.S };
const Gate Si { qpp::gt.S.conjugate() };


class Controls : public std::vector<qpp::idx> {

  using Base = std::vector<qpp::idx>;

public:

  Controls() = default;

  Controls(const std::vector<bool>& bits): Base() {
    for(unsigned i = 0; i < Config::nBit; i++)
      if(bits[i])
        Base::push_back(i);
  }

  const Base& as_vector() const {
    return static_cast<const Base&>(*this);
  }

}; // class Controls


class State : public qpp::ket {

  using Base = qpp::ket;

public:

  State(const Base& base): Base(base) { }

  State(Base&& base): Base(std::move(base)) { }

  // initializes in a basis state
  State(size_t index = 0):
    Base(qpp::mket(qpp::n2multiidx(index, dims()))) { }

  // resets in a basis state
  void reset(size_t index) {
    *this = qpp::mket(qpp::n2multiidx(index, dims()));
  }

  /* BROKEN in Eigen 3.2.9: fft expects DenseCoeffsBase::operator[] to return
   * by reference but it returns by value. Won't work in any modification!

  static State fourier(const State& in) {
    Eigen::VectorXcd ret{in.size()};
    Eigen::FFT<double> fft;
    fft.fwd(ret, in.col(0));
    return {ret};
  }*/

  static std::complex<double> overlap(const State& lhs, const State& rhs) {
    return rhs.dot(lhs);
  }

  template<class Gene, class... Args>
  State apply(const Gene& g, Args... args) {
    return g->applyTo(*this, args...);
  }

  State apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) const {
    return {qpp::applyCTRL(*this, mat, ixs, {tgt})};
  }

  friend std::ostream& operator<< (std::ostream& os, const State& state) {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
        " ", " "); // row, col separators
    return os << state.format(fmt) << '\n';
  }

private:

  std::vector<qpp::idx> dims() {
    return std::vector<qpp::idx>(Config::nBit, 2);
  }

}; // class State

} // namespace Backend

} // namespace QGA

#endif // !defined QGA_BACKEND_HPP
