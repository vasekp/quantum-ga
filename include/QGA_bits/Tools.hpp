namespace QGA {

double rationalize(double x);
double rationalize_angle(double a);


/* An enum for possible settings for controls_distribution below */

enum class Controls {
  NONE,
  ONE,
  LEAST1,
  ANY
};

/* A distribution generating bit strings of length nBit where the probability
 * of 1 in each position is given by pTrue. The bit at position nSkip is left
 * off. */

template<Controls cc>
class controls_distribution {

public:

  controls_distribution(unsigned nBit_, unsigned iSkip_, double pTrue_):
    nBit(nBit_), pTrue(pTrue_), iSkip(iSkip_) {
#ifdef DEBUG
      if((cc == Controls::ONE || cc == Controls::LEAST1) && nBit <= 1)
        throw std::logic_error("nBit < 2 in controls_distribution<≥1>!");
#endif
    }

  template<class URNG>
  std::vector<bool> operator() (URNG& rng) {
    std::vector<bool> bits(nBit, false);

    if(cc == Controls::ANY || cc == Controls::LEAST1) {
      std::bernoulli_distribution dist(pTrue);
      for(unsigned i = 0; i < nBit; i++) {
        if(i == iSkip)
          continue;
        bits[i] = dist(rng);
      }
    }

    if(cc == Controls::ONE || cc == Controls::LEAST1) {
      std::uniform_int_distribution<> dist(0, nBit - 2);
      unsigned res = dist(rng);
      bits[res + (res >= iSkip)] = true;
    }

    return bits;
  }

private:

  const unsigned nBit;
  const double pTrue;
  const unsigned iSkip;

}; // class controls_distribution<Controls>


/* Two related pre-initialized distributions: for generating initial values of
 * angle for parametric gates and for generating angle deviations for
 * continuous gate mutation. */

template<bool diff = false>
struct angle_distribution : std::uniform_real_distribution<> {

  angle_distribution():
    std::uniform_real_distribution<>(-Const::pi, Const::pi) { }

}; // class angle_distribution<false>

template<>
struct angle_distribution<true> : std::normal_distribution<> {

  angle_distribution():
    std::normal_distribution<>(0.0, Config::dAlpha) { }

}; // class angle_distribution<true>

} // namespace QGA
