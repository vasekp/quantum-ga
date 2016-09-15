#include <iostream>

#include "genetic.hpp"
#include "commons.hpp"
#include "wrapper-qiclib.hpp"

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

  const unsigned nBit = 3;

} // namespace Config

using Candidate = Wrapper::Candidate;
using Population = gen::NSGAPopulation<Candidate>;
using GenCandidate = gen::Candidate<Candidate>;


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
    using cmp = decltype(func)::value_type;
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


template<class Candidate>
std::atomic_ulong CBase<Candidate>::count{0};

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

  Wrapper::init();

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
    std::cout << c.fitness() << ' ' << c << ": " << c.dump();
  }
}
