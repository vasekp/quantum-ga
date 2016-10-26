// allow only one backend
#ifndef QGA_BACKEND_HPP
#define QGA_BACKEND_HPP

#include "qpp.h"
//#include <unsupported/Eigen/FFT>

namespace QGA {

namespace Backend {

class Gate : Eigen::Matrix2cd {

  using Base = Eigen::Matrix2cd;
  using cxd = std::complex<double>;

public:

  Gate(const Base& mat) : Base(mat) { }

  Gate(cxd u11, cxd u12, cxd u21, cxd u22) : Base(2, 2) {
    *this << u11, u12, u21, u22;
  }

  const Base& rep() const {
    return static_cast<const Base&>(*this);
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


class Controls : std::vector<qpp::idx> {

  using Base = std::vector<qpp::idx>;

public:

  Controls() = default;

  Controls(const std::vector<bool>& bits): Base() {
    for(unsigned i = 0; i < Config::nBit; i++)
      if(bits[i])
        Base::push_back(i);
  }

  friend bool operator==(const Controls& lhs, const Controls& rhs) {
    return lhs.rep() == rhs.rep();
  }

  const Base& as_vector() const {
    return static_cast<const Base&>(*this);
  }

  static Controls swap(unsigned s1, unsigned s2) {
    Base ret(Config::nBit);
    for(unsigned i = 0; i < Config::nBit; i++)
      ret[i] = i;
    ret[s1] = s2;
    ret[s2] = s1;
    return {std::move(ret)};
  }

  const Base& rep() const {
    return static_cast<const Base&>(*this);
  }

  using Base::size;

private:

  Controls(Base&& vec): Base(std::move(vec)) { }

}; // class Controls


class State : qpp::ket {

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

  template<class Gene, class Context = void>
  State apply(const Gene& g, const Context* c = nullptr) {
    return g->applyTo(*this, c);
  }

  State apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) const {
    return {qpp::applyCTRL(rep(), mat.rep(), ixs.rep(), {tgt})};
  }

  State swap(const Controls& ixs) const {
    return {qpp::syspermute(rep(), ixs.rep())};
  }

  friend std::ostream& operator<< (std::ostream& os, const State& state) {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
        " ", " "); // row, col separators
    return os << state.format(fmt) << '\n';
  }

  using Base::operator[];

private:

  const Base& rep() const {
    return static_cast<const Base&>(*this);
  }

  std::vector<qpp::idx> dims() {
    return std::vector<qpp::idx>(Config::nBit, 2);
  }

}; // class State

} // namespace Backend

} // namespace QGA

#endif // !defined QGA_BACKEND_HPP
