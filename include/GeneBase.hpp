class GeneBase {

  unsigned op;      // which gate to use (see Glogals::gates)
  unsigned tgt;     // target qubit
  unsigned ctrlEnc; // 0 through UINT_MAX
  unsigned ctrlDec; // binary control string with 0 at target
  uint32_t hw;      // Hamming weight of ctrl

public:

  NOINLINE GeneBase(unsigned op_, unsigned target_, unsigned control_):
      op(op_), tgt(target_), ctrlEnc(control_), hw(0) {
    size_t ctrl = 0;
    double c = (double)control_ / std::numeric_limits<unsigned>::max();
    /* Convert an unsigned between 0 and UINT_MAX to a bit string where the
     * probability of 1 in each position is given by Config::pControl. A value
     * less than 0.5 means that plain NOTs and C-NOTs will be generated more
     * often than CC-NOTs and higher. */
    for(unsigned i = 0; i < Config::nBit - 1; i++) {
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
    ctrlDec =
      ((ctrl >> tgt) << (tgt+1))  // shift bits left of tgt to the left
        |
      (ctrl & ((1 << tgt) - 1));  // keep bits right of tgt
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

  unsigned controlDec() const {
    return ctrlDec;
  }

  unsigned weight() const {
    return hw;
  }

}; // class GeneBase
