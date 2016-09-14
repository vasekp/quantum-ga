namespace Config {

  const float selectBias = 1.0;
  const size_t popSize = 10;
  const size_t popSize2 = 100;
  const int nGen = 50;

  const float expLengthIni = 30;      // expected length of circuits in 0th generation
  const float expLengthAdd = 1.5;     // expected length of gates inserted in mutation
  const float pDeleteUniform = 0.10;  // probability of single gate deletion

  const float heurFactor = 0.15;      // how much prior success of genetic ops should influence future choices

  const float pControl = 0.25;        // how much each bit is likely to be a control bit at gate creation

  const int nBit = 3;

} // namespace Config


class Gene {

  unsigned op;      // which gate to use (see Glogals::gates)
  unsigned tgt;     // target qubit
  std::vector<unsigned> ixs;   // list of control bits
  unsigned ctrlEnc; // 0 through UINT_MAX
  uint32_t hw;      // Hamming weight of ctrl

public:

  NOINLINE Gene(unsigned op_, unsigned target_, unsigned control_):
      op(op_), tgt(target_), ixs{}, ctrlEnc(control_), hw(0) {
    size_t ctrl = 0;
    double c = (double)control_ / std::numeric_limits<unsigned>::max();
    /* Convert an unsigned between 0 and UINT_MAX to a bit string where the
     * probability of 1 in each position is given by Config::pControl. A value
     * less than 0.5 means that plain NOTs and C-NOTs will be generated more
     * often than CC-NOTs and higher. */
    for(int i = 0; i < Config::nBit - 1; i++) {
      ctrl <<= 1;
      if(c < Config::pControl) {
        ctrl |= 1;
        hw++;
        c /= Config::pControl;
      } else {
        c = (c - Config::pControl)/(1 - Config::pControl);
      }
    }
    /* At this point ctrl has nBit-1 bits. We use this to guarantee that
     * 1<<tgt is left unoccupied. */
    ctrl =
      ((ctrl >> tgt) << (tgt+1))  // shift bits left of tgt to the left
        |
      (ctrl & ((1 << tgt) - 1));  // keep bits right of tgt
    ixs.reserve(Config::nBit);
    for(int i = 0; i < Config::nBit; i++) {
      if(ctrl & 1)
        ixs.push_back(i);
      ctrl >>= 1;
    }
  }

  unsigned gate() const {
    return op;
  }

  unsigned target() const {
    return tgt;
  }

  unsigned control() const {
    return ctrlEnc;
  }

  const std::vector<unsigned>& ix_vector() const {
    return ixs;
  }

  unsigned weight() const {
    return hw;
  }

}; // class Gene


struct Fitness {

  double error;
  size_t length;
  size_t cplx;

  friend std::ostream& operator<< (std::ostream& os, const Fitness& f) {
    return os << '{' << f.error << ',' << f.length << ',' << f.cplx << '}';
  }

  friend NOINLINE bool operator<< (const Fitness& a, const Fitness& b) {
    return a.error <= b.error && a.length <= b.length && a.cplx <= b.cplx && !(a == b);
  }

  friend bool operator== (const Fitness& a, const Fitness& b) {
    return a.error == b.error && a.length == b.length && a.cplx == b.cplx;
  }

}; // struct Fitness
