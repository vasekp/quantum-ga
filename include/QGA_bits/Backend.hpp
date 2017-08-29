namespace QGA {

namespace Backend {

using cxd = std::complex<double>;


class Gate {

  class GateImpl;
  std::unique_ptr<GateImpl> pImpl;

public:

  Gate(const Gate&);
  Gate(GateImpl&&);
  Gate(cxd u11, cxd u12, cxd u21, cxd u22);
  ~Gate();

  friend Gate operator*(const Gate& lhs, const Gate& rhs);

  cxd operator() (size_t rx, size_t cx);

private:

  const GateImpl& impl() const;

  friend class State;

}; // class Gate

extern const Gate I;
extern const Gate H;
extern const Gate X;
extern const Gate Y;
extern const Gate Z;
extern const Gate T;
extern const Gate Ti;
extern const Gate S;
extern const Gate Si;


class Controls {

  class ControlsImpl;
  std::unique_ptr<ControlsImpl> pImpl;

public:

  Controls();
  Controls(const Controls&);
  Controls(ControlsImpl&&);
  Controls(const std::vector<bool>& bits);
  ~Controls();

  Controls& operator=(const Controls&);
  Controls& operator=(Controls&&);

  friend bool operator==(const Controls& lhs, const Controls& rhs);

  size_t size() const;

  static Controls swapGate(unsigned s1, unsigned s2);
  static Controls swapQubits(const Controls& orig, unsigned s1, unsigned s2);

  std::vector<unsigned> as_vector() const;

private:

  const ControlsImpl& impl() const;

  friend class State;

}; // class Controls


class State {

  class StateImpl;
  std::unique_ptr<StateImpl> pImpl;

public:

  State(const State&);
  State(State&&);
  State(StateImpl&&);
  State(size_t index = 0);
  ~State();

  State& operator=(const State&);
  State& operator=(State&&);

  void reset(size_t index);

  State apply_ctrl(const Gate& mat, const Controls& ixs, unsigned tgt) const;
  State swapQubits(const Controls& ixs) const;

  static State fourier(const State& in);
  static cxd overlap(const State& lhs, const State& rhs);

  friend std::ostream& operator<< (std::ostream& os, const State& state);

  cxd& operator[](size_t);

private:

  const StateImpl& impl() const;
  StateImpl& impl();

}; // class State

} // namespace Backend

} // namespace QGA
