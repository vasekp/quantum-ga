#include <atomic>
#include <string>

/* Forward declarations */
namespace Config {
  extern const unsigned nBit;
  extern const float pControl;
}

namespace Wrapper {
  std::string gate_name(unsigned);
}
/* End forward declarations */


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
    ctrl =
      ((ctrl >> tgt) << (tgt+1))  // shift bits left of tgt to the left
        |
      (ctrl & ((1 << tgt) - 1));  // keep bits right of tgt
    ixs.reserve(Config::nBit);
    for(unsigned i = 0; i < Config::nBit; i++) {
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


template<class Derived>
class CBase {

protected:

  std::vector<Gene> gt{};

private:

  int origin = -1;

  /* Total fitness() evaluation counter */
  static std::atomic_ulong count;

  const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

public:

  CBase(std::vector<Gene>&) = delete;

  CBase(std::vector<Gene>&& gt_): gt(std::move(gt_)) { }

  NOINLINE Fitness fitness() const {
    /* Complexity = square sum of numbers of control bits per gate */
    size_t cplx = 0;
    for(const Gene& g : gt) {
      unsigned h = g.weight();
      cplx += h*h;
    }
    count++;
    return {derived().error(), gt.size(), cplx};
  }

  friend std::ostream& operator<< (std::ostream& os, const CBase& c) {
    auto first = c.gt.begin(), last = c.gt.end();
    for(auto it = first; it != last; it++) {
      if(it != first) os << ' ';
      os << Wrapper::gate_name(it->gate()) << it->target()+1;
      auto& ixv = it->ix_vector();
      if(ixv.size()) {
        os << '[';
        for(auto ix : ixv)
          os << ix;
        os << ']';
      }
    }
    return os;
  }

  const std::vector<Gene>& genotype() const {
    return gt;
  }

  void setOrigin(int _origin) {
    origin = _origin;
  }

  int getOrigin() const {
    return origin;
  }

  static unsigned long totalCount() {
    return count;
  }

  friend class CandidateFactory;

}; // class Candidate
