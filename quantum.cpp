#define QICLIB_DONT_USE_NLOPT
#define ARMA_DONT_USE_WRAPPER

#include <iostream>
#include <cmath>
#include <atomic>

#include "QIClib"
#include "genetic.hpp"

namespace Config {

  const float selectBias = 1.0;
  const size_t popSize = 100;
  const size_t popSize2 = 1000;
  const int nGen = 50;

  const float expLengthIni = 30;      // expected length of circuits in 0th generation
  const float expLengthAdd = 1.5;     // expected length of gates inserted in mutation
  const float pDeleteUniform = 0.10;  // probability of single gate deletion

  const float heurFactor = 0.15;      // how much prior success of genetic ops should influence future choices

  const float pControl = 0.25;        // how much each bit is likely to be a control bit at gate creation

  const int nBit = 3;

} // namespace Config


namespace Globals {

  arma::cx_mat22 H {
    1/std::sqrt(2),  1/std::sqrt(2),
    1/std::sqrt(2), -1/std::sqrt(2)
  };

  arma::cx_mat22 X {
    0, 1, 1, 0
  };

  arma::cx_mat22 Y {
    0, {0,-1}, {0,1}, 0
  };

  arma::cx_mat22 Z {
    1, 0, 0, -1
  };

  arma::cx_mat22 T {
    1, 0, 0, {0,1}
  };

  struct Gate {
    arma::cx_mat22 op;
    char name;
  };

  std::vector<Gate> gates {
    { H, 'H' },
    /*{ X, 'X' },
    { Y, 'Y' },
    { Z, 'Z' },*/
    { T, 'T' }
  };

  arma::cx_vec out{};

} // namespace Globals


class Gene {

  unsigned op;      // which gate to use (see Glogals::gates)
  unsigned tgt;     // target qubit
  arma::uvec ixs;   // list of control bits
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
    std::vector<arma::uword> ixv{};
    for(int i = 0; i < Config::nBit; i++) {
      if(ctrl & 1)
        ixv.push_back(i+1);
      ctrl >>= 1;
    }
    ixs = ixv;
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

  unsigned weight() const {
    return hw;
  }

  arma::cx_vec apply(arma::cx_vec psi) const {
    return qic::apply_ctrl(psi, Globals::gates[op].op, ixs, {1+tgt});
  }

  friend std::ostream& operator<< (std::ostream& os, const Gene& g) {
    os << Globals::gates[g.op].name << g.tgt+1;
    if(g.ixs.size()) {
      os << '[';
      for(auto ix : g.ixs)
        os << ix;
      os << ']';
    }
    return os;
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


class Candidate {

  std::vector<Gene> gt{};
  int origin = -1;
  static std::atomic_ulong count;

public:

  Candidate(std::vector<Gene>& _gt) = delete;

  Candidate(std::vector<Gene>&& _gt): gt(std::move(_gt)) { }

  NOINLINE Fitness fitness() const {
    arma::cx_vec psi = qic::mket({0}, {1 << Config::nBit});
    for(const Gene& g : gt)
      psi = g.apply(psi);
    size_t cplx = 0;
    for(const Gene& g : gt) {
      unsigned h = g.weight();
      cplx += h*h;
    }
    count++;
    return {1 - std::abs(cdot(psi, Globals::out)), gt.size(), cplx};
  }

  friend std::ostream& operator<< (std::ostream& os, const Candidate& c) {
    auto first = c.gt.begin(), last = c.gt.end();
    for(auto it = first; it != last; it++)
      os << (it == first ? "" : " ") << *it;
    return os;
  }

  void setOrigin(int _origin) {
    origin = _origin;
  }

  int getOrigin() const {
    return origin;
  }

  void dump(std::ostream& os) const {
    arma::cx_vec psi = qic::mket({0}, {1 << Config::nBit});
    for(const Gene& g : gt)
      psi = g.apply(psi);
    psi.t().raw_print(os);
  }

  static unsigned long totalCount() {
    return count;
  }

  friend class CandidateFactory;

};


typedef gen::NSGAPopulation<Candidate> Population;
typedef gen::Candidate<Candidate> GenCandidate;


class CandidateFactory {

  typedef Candidate (CandidateFactory::*GenOp)();

  std::uniform_real_distribution<> dUni{0, 1};
  std::uniform_int_distribution<unsigned> dOp{0, Globals::gates.size() - 1};
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
    static thread_local std::uniform_int_distribution<unsigned> dOp{0, Globals::gates.size() - 1};
    static thread_local std::uniform_int_distribution<unsigned> dTgt{0, Config::nBit - 1};
    static thread_local std::uniform_int_distribution<unsigned> dCtrl{};
    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do {
      gt.push_back({dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng)});
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
    typedef decltype(func)::value_type cmp;
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
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(dOp(gen::rng), gm[pos].target(), gm[pos].control());
    return Candidate{std::move(gm)};
  }

  Candidate mAlterTarget() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(gm[pos].gate(), dTgt(gen::rng), gm[pos].control());
    return Candidate{std::move(gm)};
  }

  Candidate mAlterControl() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(gm[pos].gate(), gm[pos].target(),
        gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate{std::move(gm)};
  }

  Candidate mAlterSingle() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    auto gm = p.gt;
    unsigned pos = gen::rng() % p.gt.size();
    gm[pos] = Gene(dOp(gen::rng), dTgt(gen::rng), gm[pos].control() ^ dCtrl(gen::rng));
    return Candidate{std::move(gm)};
  }

  Candidate mAddSlice() {
    auto &p = get();
    unsigned pos = gen::rng() % (p.gt.size() + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() + ins.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mAddPairs() {
    auto &p = get();
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins;
    ins.reserve(2*Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do {
      ins.emplace_back(dOp(gen::rng), dTgt(gen::rng), dCtrl(gen::rng));
    } while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() + 2*ins.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.insert(gm.end(), std::make_move_iterator(ins.rbegin()), std::make_move_iterator(ins.rend()));
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteSlice() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() - (pos2 - pos1));
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteSliceShort() {
    auto &p = get();
    auto sz = p.gt.size();
    if(sz == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1);
    /* Integer with the same distribution in mAddSlice */
    int len = 1 + floor(log(dUni(gen::rng)) / log(1 - 1/Config::expLengthAdd));
    unsigned pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gm;
    gm.reserve(p.gt.size() - (pos2 - pos1));
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteUniform() {
    auto &p = get();
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    for(auto& g : p.gt)
      if(dUni(gen::rng) >= Config::pDeleteUniform)
        gm.push_back(g);
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap2() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos = gen::rng() % (p.gt.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin() + pos, p.gt.end());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos);
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap4() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    unsigned pos1 = gen::rng() % (p.gt.size() + 1),
             pos2 = gen::rng() % (p.gt.size() + 1),
             pos3 = gen::rng() % (p.gt.size() + 1);
    if(pos2 < pos1) std::swap(pos1, pos2);
    if(pos3 < pos1) std::swap(pos1, pos3);
    if(pos3 < pos2) std::swap(pos2, pos3);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.begin() + pos3);
    gm.insert(gm.end(), p.gt.begin() + pos1, p.gt.begin() + pos2);
    gm.insert(gm.end(), p.gt.begin() + pos3, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mSplitSwap5() {
    auto &p = get();
    if(p.gt.size() == 0)
      return p;
    std::vector<unsigned> pos;
    for(int i = 0; i < 4; i++)
      pos.push_back(gen::rng() % (p.gt.size() + 1));
    std::sort(pos.begin(), pos.end());
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos[0]);
    gm.insert(gm.end(), p.gt.begin() + pos[2], p.gt.begin() + pos[3]);
    gm.insert(gm.end(), p.gt.begin() + pos[1], p.gt.begin() + pos[2]);
    gm.insert(gm.end(), p.gt.begin() + pos[0], p.gt.begin() + pos[1]);
    gm.insert(gm.end(), p.gt.begin() + pos[3], p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mReverseSlice() {
    auto &p = get();
    auto sz = p.gt.size();
    if(sz == 0)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(p.gt.size());
    gm.insert(gm.end(), p.gt.begin(), p.gt.begin() + pos1);
    gm.insert(gm.end(), p.gt.rbegin() + sz - pos2, p.gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), p.gt.begin() + pos2, p.gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossover1() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
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
    auto &gt1 = p1.gt,
         &gt2 = p2.gt;
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


std::atomic_ulong Candidate::count { 0 };

std::vector<unsigned> CandidateFactory::weights;

const std::vector<std::pair<CandidateFactory::GenOp, std::string>> CandidateFactory::func {
    { &CandidateFactory::mAlterGate,         "MGate" },
    { &CandidateFactory::mAlterTarget,       "MTarget" },
    { &CandidateFactory::mAlterControl,      "MControl" },
    //{ &CandidateFactory::mAlterSingle,       "MSingle" },
    { &CandidateFactory::mAddSlice,          "AddSlice" },
    { &CandidateFactory::mAddPairs,          "AddPairs" },
    { &CandidateFactory::mDeleteSlice,       "DelSlice" },
    { &CandidateFactory::mDeleteSliceShort,  "DelShort" },
    //{ &CandidateFactory::mDeleteUniform,     "DelUnif" },
    //{ &CandidateFactory::mSplitSwap2,        "SpltSwp2"  },
    { &CandidateFactory::mSplitSwap4,        "SpltSwp4"  },
    { &CandidateFactory::mSplitSwap5,        "SpltSwp5"  },
    //{ &CandidateFactory::mReverseSlice,      "InvSlice" },
    { &CandidateFactory::crossover1,         "C/Over1" },
    { &CandidateFactory::crossover2,         "C/Over2" },
  };


int main() {
#ifdef BENCH
  gen::rng.seed(1);
  omp_set_num_threads(1);
#endif
  Globals::out = qic::mket({3}, {1 << Config::nBit});

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  Population pop{Config::popSize, [&] { return CandidateFactory::genInit(); }};

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do popSize */
    auto nondom = pop.front();
    //nondom.randomTrim(Config::popSize); // not necessary: this is usually about 10
    size_t nd = nondom.size();

    /* Top up to popSize2 candidates, precomputing fitnesses */
    Population pop2{Config::popSize2};
    pop.precompute();
    CandidateFactory cf{pop};
    pop2.add(Config::popSize2 - nd, [&]() -> const Candidate { return cf.getNew(); }, true);

    /* Merge the nondominated subset of the previous population */
    pop2.add(nondom);
    pop = std::move(pop2);

    for(auto& c : pop.front().randomSelect(Config::popSize))
      CandidateFactory::hit(c.getOrigin());
    pop.prune([](const GenCandidate& a, const GenCandidate& b) -> bool {
      return a.fitness() == b.fitness();
      });

    /* Summarize */
    nondom = pop.front();
    std::cout << "Gen " << gen << ": " << nondom.size() << " nondominated";
    if(nondom.size() > 0) {
      auto& e = nondom.randomSelect();
      std::cout << ", e.g. " << e.fitness() << ' ' << e;
    }
    std::cout << std::endl;

    /* Make older generations matter less in the choice of gen. op. */
    CandidateFactory::normalizeWeights();

  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre;
  std::cout << std::endl << "Run took " << dur.count() << " s (" << dur.count()/Config::nGen << " s/gen avg), " <<
    Candidate::totalCount() << " candidates tested, best of run:" << std::endl;

  /* Dump the heuristic distribution */
  std::cout << "\nGenetic operator distribution:\n";
  CandidateFactory::dumpWeights(std::cout);

  /* List results */
  auto nondom = pop.front();
  std::cout << '\n' << nondom.size() << " nondominated candidates, ";
  nondom.prune([](const GenCandidate& a, const GenCandidate& b) -> bool {
      return a.fitness() == b.fitness();
    });
  std::cout << nondom.size() << " unique fitnesses:\n";
  for(auto& c : nondom) {
    std::cout << c.fitness() << ' ' << c << ": ";
    c.dump(std::cout);
  }
}
