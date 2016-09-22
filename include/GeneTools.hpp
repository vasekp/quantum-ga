namespace QGA {

namespace GeneTools {

size_t ctrlBitString(unsigned rand, unsigned target) {
  size_t ctrl = 0;
  double c = (double)rand / std::numeric_limits<unsigned>::max();
  /* Convert an unsigned between 0 and UINT_MAX to a bit string where the
   * probability of 1 in each position is given by Config::pControl. A value
   * less than 0.5 means that plain NOTs and C-NOTs will be generated more
   * often than CC-NOTs and higher. */
  for(unsigned i = 0; i < Config::nBit - 1; i++) {
    ctrl <<= 1;
    if(c < Config::pControl) {
      ctrl |= 1;
      c /= Config::pControl;
    } else {
      c = (c - Config::pControl)/(1 - Config::pControl);
    }
  }
  /* At this point ctrl has nBit - 1 bits. We use this to guarantee that
   * 1 << target is left unoccupied. */
  return
    ((ctrl >> target) << (target + 1))  // shift bits left of tgt to the left
      |
    (ctrl & ((1 << target) - 1));  // keep bits right of tgt
}

} // namespace GeneTools

} // namespace QGA
