#include <atomic>
#include <string>
#include <iomanip>

/* Forward declarations */
namespace Config {
  extern const unsigned nBit;
  extern const size_t popSize;
  extern const float selectBias;
  extern const float heurFactor;
  extern const float expLengthIni;
  extern const float expLengthAdd;
  extern const float pDeleteUniform;
  extern const float pControl;
}

namespace Wrapper {
  std::string gate_name(unsigned);
  extern const unsigned gate_cnt;
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

  template<class, class>
  friend class CandidateFactory;

}; // class Candidate


template<class Candidate, class Population = gen::NSGAPopulation<Candidate>>
class CandidateFactory {

  using GenOp = Candidate (CandidateFactory::*)();

  std::uniform_real_distribution<> dUni{0, 1};
  std::uniform_int_distribution<unsigned> dOp{0, Wrapper::gate_cnt - 1};
  std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
  std::uniform_int_distribution<unsigned> dCtrl{};

  static const std::vector<std::pair<GenOp, std::string>> func;
  static std::vector<unsigned> weights;
  std::discrete_distribution<> dFun{};

  Population& pop;

public:

  CandidateFactory(Population& _pop): pop(_pop) {
    if(weights.size() == 0) {
      weights = std::vector<unsigned>(func.size(), 1);
      normalizeWeights();
    }
    applyWeights();
  }

  static NOINLINE Candidate genInit() {
    const static double probTerm = 1/Config::expLengthIni;  // probability of termination; expLength = expected number of genes
    static thread_local std::uniform_real_distribution<> dUni{0, 1};
    static thread_local std::uniform_int_distribution<unsigned> dOp{0, Wrapper::gate_cnt - 1};
    static thread_local std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    static thread_local std::uniform_int_distribution<unsigned> dCtrl{};
    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do {
      gt.emplace_back(dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gt)};
  }

  NOINLINE Candidate getNew() {
    int index = dFun(gen::rng);
    Candidate c = (this->*func[index].first)();
    c.setOrigin(index);
    return c;
  }

  static void hit(int ix) {
    if(ix >= 0)
      weights[ix]++;
  }

  static void normalizeWeights() {
    unsigned total = std::accumulate(weights.begin(), weights.end(), 0);
    float factor = 1/Config::heurFactor * (float)func.size()*Config::popSize / total;
    for(auto& w : weights)
      w *= factor;
  }

  void applyWeights() {
    dFun = std::discrete_distribution<>(weights.begin(), weights.end());
  }

  static void dumpWeights(std::ostream& os) {
    float total = std::accumulate(weights.begin(), weights.end(), 0);
    int sz = func.size();
    /* Find the longest GenOp name */
    using cmp = typename decltype(func)::value_type;
    auto max = std::max_element(func.begin(), func.end(),
        [](const cmp& v1, const cmp& v2) {
          return v1.second.length() < v2.second.length();
        });
    auto maxw = max->second.length();
    auto _flags = os.flags(std::ios_base::left | std::ios_base::fixed);
    auto _precision = os.precision(4);
    for(int i = 0; i < sz; i++)
      os << std::setw(maxw+3) << func[i].second + ':' << weights[i] / total << '\n';
    os.flags(_flags);
    os.precision(_precision);
  }

  private:
  const Candidate& get() {
    return pop.NSGASelect(Config::selectBias);
  }

  Candidate mAlterGate() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    unsigned pos = gen::rng() % sz;
    gm[pos] = Gene(dOp(gen::rng), gm[pos].target(), gm[pos].control());
    return Candidate{std::move(gm)};
  }

  Candidate mAlterTarget() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    unsigned pos = gen::rng() % sz;
    gm[pos] = Gene(gm[pos].gate(), dTgt(gen::rng), gm[pos].control());
    return Candidate{std::move(gm)};
  }

  Candidate mAlterControl() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    unsigned pos = gen::rng() % sz;
    gm[pos] = Gene(gm[pos].gate(), gm[pos].target(),
        gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate{std::move(gm)};
  }

  Candidate mAlterSingle() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    unsigned pos = gen::rng() % sz;
    gm[pos] = Gene(dOp(gen::rng), dTgt(gen::rng), gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate{std::move(gm)};
  }

  Candidate mAddSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    unsigned pos = gen::rng() % (sz + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(sz + ins.size());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), gt.begin() + pos, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mAddPairs() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins;
    ins.reserve(2*Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(sz + 2*ins.size());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), gt.begin() + pos1, gt.begin() + pos2);
    gm.insert(gm.end(), std::make_move_iterator(ins.rbegin()), std::make_move_iterator(ins.rend()));
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(sz - (pos2 - pos1));
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteSliceShort() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1);
    /* Integer with the same distribution in mAddSlice */
    int len = 1 + floor(log(dUni(gen::rng)) / log(1 - 1/Config::expLengthAdd));
    unsigned pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gm;
    gm.reserve(sz - (pos2 - pos1));
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteUniform() {
    auto &p = get();
    auto gt = p.genotype();
    std::vector<Gene> gm;
    gm.reserve(gt.size());
    for(auto& g : gt)
      if(dUni(gen::rng) >= Config::pDeleteUniform)
        gm.push_back(g);
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap2() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos = gen::rng() % (sz + 1);
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin() + pos, gt.end());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos);
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap4() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1),
             pos3 = gen::rng() % (sz + 1);
    if(pos2 < pos1) std::swap(pos1, pos2);
    if(pos3 < pos1) std::swap(pos1, pos3);
    if(pos3 < pos2) std::swap(pos2, pos3);
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.begin() + pos3);
    gm.insert(gm.end(), gt.begin() + pos1, gt.begin() + pos2);
    gm.insert(gm.end(), gt.begin() + pos3, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap5() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    std::array<unsigned, 4> pos;
    for(int i = 0; i < 4; i++)
      pos[i] = gen::rng() % (sz + 1);
    std::sort(pos.begin(), pos.end());
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos[0]);
    gm.insert(gm.end(), gt.begin() + pos[2], gt.begin() + pos[3]);
    gm.insert(gm.end(), gt.begin() + pos[1], gt.begin() + pos[2]);
    gm.insert(gm.end(), gt.begin() + pos[0], gt.begin() + pos[1]);
    gm.insert(gm.end(), gt.begin() + pos[3], gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mReverseSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), gt.rbegin() + sz - pos2, gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossover1() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype();
    unsigned pos1 = gen::rng() % (gt1.size() + 1),
             pos2 = gen::rng() % (gt2.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(pos1 + (gt2.size() - pos2));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossover2() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype();
    unsigned pos1l = gen::rng() % (gt1.size() + 1),
             pos1r = gen::rng() % (gt1.size() + 1),
             pos2l = gen::rng() % (gt2.size() + 1),
             pos2r = gen::rng() % (gt2.size() + 1);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm;
    gm.reserve(gt1.size() - (pos1r - pos1l) + (pos2r - pos2r));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate{std::move(gm)};
  }

}; // class CandidateFactory
