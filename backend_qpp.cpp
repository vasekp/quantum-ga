#include "QGA_commons.hpp"
#include "QGA_bits/Backend.hpp"
#include "make_unique.hpp"

#include "qpp.h"
//#include <unsupported/Eigen/FFT>

namespace QGA {

namespace Backend {


/* Not defined by the standard until C++14 */
template<class T, typename... Args>
static std::unique_ptr<T> make_unique(Args... args) {
  return std::unique_ptr<T>{new T(std::forward<Args>(args)...)};
}


// class Gate

class Gate::GateImpl : public Eigen::Matrix2cd {
  
  public:
    using Eigen::Matrix2cd::Matrix2cd;

  Eigen::Matrix2cd& rep() {
    return static_cast<Eigen::Matrix2cd&>(*this);
  }

};

Gate::Gate(const Gate& other):
  pImpl(make_unique<GateImpl>(other.impl()))
{ }

Gate::Gate(GateImpl&& otherImpl):
  pImpl(make_unique<GateImpl>(std::move(otherImpl)))
{ }

Gate::Gate(cxd u11, cxd u12, cxd u21, cxd u22):
  pImpl(make_unique<GateImpl>(2, 2))
{
  pImpl->rep() << u11, u12, u21, u22;
}

Gate::~Gate() { }

const Gate::GateImpl& Gate::impl() const {
  return *pImpl;
}

Gate operator*(const Gate& lhs, const Gate& rhs) {
  return {lhs.impl() * rhs.impl()};
}


// Gate constants

using QGA::Const::i;
using QGA::Const::pi;
using QGA::Const::v12;

const Gate I { qpp::gt.Id2 };
const Gate H { qpp::gt.H };
const Gate X { qpp::gt.X };
const Gate Y { qpp::gt.Y };
const Gate Z { qpp::gt.Z };
const Gate T { qpp::gt.T };
const Gate Ti { qpp::gt.T.conjugate() };
const Gate S { qpp::gt.S };
const Gate Si { qpp::gt.S.conjugate() };


// class Controls

class Controls::ControlsImpl : public std::vector<qpp::idx> {

  public:
    using std::vector<qpp::idx>::vector;

  ControlsImpl() = default;

  ControlsImpl(std::vector<qpp::idx>&& other):
    std::vector<qpp::idx>(std::move(other))
  { }

  const std::vector<qpp::idx>& rep() const {
    return static_cast<const std::vector<qpp::idx>&>(*this);
  }

};

Controls::Controls():
  pImpl(make_unique<ControlsImpl>())
{ }

Controls::Controls(const Controls& other):
  pImpl(make_unique<ControlsImpl>(other.impl()))
{ }

Controls::Controls(ControlsImpl&& otherImpl):
  pImpl(make_unique<ControlsImpl>(std::move(otherImpl)))
{ }

Controls::Controls(const std::vector<bool>& bits):
  pImpl(make_unique<ControlsImpl>())
{
  for(unsigned i = 0; i < Config::nBit; i++)
    if(bits[i])
      pImpl->push_back(i);
}

Controls::~Controls() { }

const Controls::ControlsImpl& Controls::impl() const {
  return *pImpl;
}

Controls& Controls::operator=(const Controls& other) {
  pImpl = make_unique<ControlsImpl>(other.impl());
  return *this;
}

Controls& Controls::operator=(Controls&& other) {
  pImpl = make_unique<ControlsImpl>(std::move(other.impl()));
  return *this;
}

bool operator==(const Controls& lhs, const Controls& rhs) {
  return lhs.impl() == rhs.impl();
}

size_t Controls::size() const {
  return impl().size();
}

std::vector<unsigned> Controls::as_vector() const {
  std::vector<unsigned> ret{};
  ret.reserve(size());
  for(auto ix : impl())
    ret.push_back(ix);
  return ret;
}

Controls Controls::swapGate(unsigned s1, unsigned s2) {
  std::vector<qpp::idx> ret(Config::nBit);
  for(unsigned i = 0; i < Config::nBit; i++)
    ret[i] = i;
  ret[s1] = s2;
  ret[s2] = s1;
  return {std::move(ret)};
}

Controls Controls::swapQubits(const Controls& orig, unsigned s1, unsigned s2) {
  std::vector<qpp::idx> vector{orig.impl().rep()};
  for(auto&& v : vector)
    v = v == s1 ? s2 : v == s2 ? s1 : v;
  std::sort(vector.begin(), vector.end());
  return {std::move(vector)};
}


// class State

class State::StateImpl : public qpp::ket {

  public:
    using qpp::ket::ket;

    qpp::ket& rep() {
    return static_cast<qpp::ket&>(*this);
  }

  const qpp::ket& rep() const {
    return static_cast<const qpp::ket&>(*this);
  }

};

static std::vector<qpp::idx> dims() {
  return std::vector<qpp::idx>(Config::nBit, 2);
}

State::State(const State& other):
  pImpl(make_unique<StateImpl>(other.impl()))
{ }

State::State(State&& other) = default;

State::State(StateImpl&& otherImpl):
  pImpl(make_unique<StateImpl>(std::move(otherImpl)))
{ }

State::State(size_t index):
  pImpl(make_unique<StateImpl>(qpp::mket(qpp::n2multiidx(index, dims()))))
{ }

State::~State() { }

State& State::operator=(const State& other) {
  pImpl = make_unique<StateImpl>(other.impl());
  return *this;
}

State& State::operator=(State&& other) {
  pImpl = make_unique<StateImpl>(std::move(other.impl()));
  return *this;
}

const State::StateImpl& State::impl() const {
  return *pImpl;
}

State::StateImpl& State::impl() {
  return *pImpl;
}

void State::reset(size_t index) {
  impl() = qpp::mket(qpp::n2multiidx(index, dims()));
}

/* BROKEN in Eigen 3.2.9: fft expects DenseCoeffsBase::operator[] to return
 * by reference but it returns by value. Won't work in any modification!

State State::fourier(const State& in) {
  Eigen::VectorXcd ret{in.impl().size()};
  Eigen::FFT<double> fft;
  fft.fwd(ret, in.impl().col(0));
  return {ret};
}*/

cxd State::overlap(const State& lhs, const State& rhs) {
  return rhs.impl().dot(lhs.impl());
}

State State::apply_ctrl(
    const Gate& mat,
    const Controls& ixs,
    unsigned tgt) const
{
  return {qpp::applyCTRL(impl(), mat.impl(), ixs.impl(), {tgt})};
}

State State::swapQubits(const Controls& ixs) const {
  return {qpp::syspermute(impl(), ixs.impl())};
}

std::ostream& operator<< (std::ostream& os, const State& state) {
  Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
      " ", " "); // row, col separators
  return os << state.impl().format(fmt) << '\n';
}

cxd& State::operator[](size_t index) {
  return impl()[index];
}

} // namespace Backend

} // namespace QGA
