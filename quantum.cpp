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
using CF = CandidateFactory<Candidate>;


template<>
std::atomic_ulong CBase<Candidate>::count{0};

template<>
std::vector<unsigned> CF::weights{};

template<>
const std::vector<std::pair<CF::GenOp, std::string>> CF::func {
    { &CF::mAlterGate,         "MGate" },
    { &CF::mAlterTarget,       "MTarget" },
    { &CF::mAlterControl,      "MControl" },
    //{ &CF::mAlterSingle,       "MSingle" },
    { &CF::mAddSlice,          "AddSlice" },
    { &CF::mAddPairs,          "AddPairs" },
    { &CF::mDeleteSlice,       "DelSlice" },
    { &CF::mDeleteSliceShort,  "DelShort" },
    //{ &CF::mDeleteUniform,     "DelUnif" },
    //{ &CF::mSplitSwap2,        "SpltSwp2"  },
    { &CF::mSplitSwap4,        "SpltSwp4"  },
    { &CF::mSplitSwap5,        "SpltSwp5"  },
    //{ &CF::mReverseSlice,      "InvSlice" },
    { &CF::crossover1,         "C/Over1" },
    { &CF::crossover2,         "C/Over2" },
  };


int main() {
#ifdef BENCH
  gen::rng.seed(1);
  omp_set_num_threads(1);
#endif

  Wrapper::init();

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  Population pop{Config::popSize, [&] { return CF::genInit(); }};

  for(int gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do popSize */
    auto nondom = pop.front();
    //nondom.randomTrim(Config::popSize); // not necessary: this is usually about 10
    size_t nd = nondom.size();

    /* Top up to popSize2 candidates, precomputing fitnesses */
    Population pop2{Config::popSize2};
    pop.precompute();
    CF cf{pop};
    pop2.add(Config::popSize2 - nd, [&]() -> const Candidate { return cf.getNew(); }, true);

    /* Merge the nondominated subset of the previous population */
    pop2.add(nondom);
    pop = std::move(pop2);

    for(auto& c : pop.front().randomSelect(Config::popSize))
      CF::hit(c.getOrigin());
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
    CF::normalizeWeights();

  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre;
  std::cout << std::endl << "Run took " << dur.count() << " s (" << dur.count()/Config::nGen << " s/gen avg), " <<
    Candidate::totalCount() << " candidates tested, best of run:" << std::endl;

  /* Dump the heuristic distribution */
  std::cout << "\nGenetic operator distribution:\n";
  CF::dumpWeights(std::cout);

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
