#include "QGA_commons.hpp"
#include "QGA_bits/Backend.hpp"
#include "make_unique.hpp"

#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER
#include "QIClib"

namespace QGA {

namespace Backend {


// class Gate

class Gate::GateImpl : public arma::cx_mat22 {
  public:
    using arma::cx_mat22::cx_mat22;
};

Gate::Gate(const Gate& other):
  pImpl(make_unique<GateImpl>(other.impl()))
{ }

Gate::Gate(GateImpl&& otherImpl):
  pImpl(make_unique<GateImpl>(std::move(otherImpl)))
{ }

Gate::Gate(cxd u11, cxd u12, cxd u21, cxd u22):
  pImpl(make_unique<GateImpl>(
        std::initializer_list<cxd>{u11, u12, u21, u22}))
{ }

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

const Gate I { 1, 0, 0, 1 };
const Gate H { v12, v12, v12, -v12 };
const Gate X { 0, 1, 1, 0 };
const Gate Y { 0, -i, i, 0 };
const Gate Z { 1, 0, 0, -1 };
const Gate T { 1, 0, 0, std::exp(i*pi/4.) };
const Gate Ti { 1, 0, 0, std::exp(-i*pi/4.) };
const Gate S { 1, 0, 0, i };
const Gate Si { 1, 0, 0, -i };


// class Controls

class Controls::ControlsImpl : public arma::uvec {

  public:
    using arma::uvec::uvec;

  const arma::uvec& rep() const {
    return static_cast<const arma::uvec&>(*this);
  }

};

static std::vector<arma::uword> ix_vector(const std::vector<bool>& bits) {
  std::vector<arma::uword> ret{};
  for(unsigned i = 0; i < Config::nBit; i++)
    if(bits[i])
      ret.push_back(i + 1);
  return ret;
}

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
  pImpl(make_unique<ControlsImpl>(ix_vector(bits)))
{ }

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
  return lhs.size() == rhs.size() &&
    arma::all(lhs.impl().rep() == rhs.impl().rep());
}

size_t Controls::size() const {
  return impl().size();
}

std::vector<unsigned> Controls::as_vector() const {
  std::vector<unsigned> ret{};
  ret.reserve(size());
  for(auto ix : impl())
    ret.push_back(ix - 1);
  return ret;
}

Controls Controls::swapGate(unsigned s1, unsigned s2) {
  arma::uvec ret(Config::nBit);
  for(unsigned i = 0; i < Config::nBit; i++)
    ret[i] = i + 1;
  ret[s1] = s2 + 1;
  ret[s2] = s1 + 1;
  return {std::move(ret)};
}

Controls Controls::swapQubits(const Controls& orig, unsigned s1, unsigned s2) {
  arma::uvec vector{orig.impl().rep()};
  for(auto&& v : vector)
    v = v == s1 + 1 ? s2 + 1 : v == s2 + 1 ? s1 + 1 : v;
  std::sort(vector.begin(), vector.end());
  return {std::move(vector)};
}


// class State

class State::StateImpl : public arma::cx_vec {

  public:
    using arma::cx_vec::cx_vec;

  arma::cx_vec& rep() {
    return static_cast<arma::cx_vec&>(*this);
  }

  const arma::cx_vec& rep() const {
    return static_cast<const arma::cx_vec&>(*this);
  }

};

static arma::uword dim() {
  return arma::uword(1) << Config::nBit;
}

State::State(const State& other):
  pImpl(make_unique<StateImpl>(other.impl()))
{ }

State::State(State&& other) = default;

State::State(StateImpl&& otherImpl):
  pImpl(make_unique<StateImpl>(std::move(otherImpl)))
{ }

State::State(size_t index):
  pImpl(make_unique<StateImpl>(dim(), arma::fill::zeros))
{
  reset(index);
}

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
  impl().fill(0);
  impl()[index] = 1;
}

State State::fourier(const State& in) {
  return {arma::fft(in.impl().rep()) / sqrt(dim())};
}

cxd State::overlap(const State& lhs, const State& rhs) {
  return arma::cdot(lhs.impl().rep(), rhs.impl().rep());
}

State State::apply_ctrl(
    const Gate& gate,
    const Controls& ixs,
    unsigned tgt) const
{
  return {qic::apply_ctrl(
      impl().rep(),
      gate.impl(),
      ixs.impl(),
      {tgt + 1}
    )};
}

State State::swapQubits(const Controls& ixs) const {
  return {qic::sysperm(impl().rep(), ixs.impl())};
}

std::ostream& operator<< (std::ostream& os, const State& state) {
  state.impl().st().raw_print(os);
  return os;
}

cxd& State::operator[](size_t index) {
  return impl()[index];
}

} // namespace Backend

} // namespace QGA
