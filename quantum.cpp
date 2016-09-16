#include <iostream>

#include "genetic.hpp"

#include "include/commons.hpp"

#ifdef USE_QPP
#include "include/wrapper-qpp.hpp"
#elif defined USE_QICLIB
#include "include/wrapper-qiclib.hpp"
#else
#error Either USE_QPP or USE_QICLIB needed.
#endif

namespace Config {

  // Circuit width (constant)
  const unsigned nBit = 3;

  // strength parameter of NSGA selection
  const float selectBias = 1.0;

  // Archive (external population) size
  const size_t popSize = 10;

  // Internal population size
  const size_t popSize2 = 500;

  // Number of generations (constant)
  const int nGen = 100;

  // Expected curcuit depth in 0th generation
  const float expLengthIni = 30;

  // Expected number of gates inserted in mutation
  const float expLengthAdd = 1.5;

  // Probability of single gate deletion
  const float pDeleteUniform = 0.10;

  // How much prior success of genetic ops should influence future choices
  const float heurFactor = 0.15;

  // How much each bit is likely to be a control bit at gate creation
  const float pControl = 0.25;

} // namespace Config


using Candidate = Wrapper::Candidate;
using Population = gen::NSGAPopulation<Candidate>;
using GenCandidate = gen::Candidate<Candidate>;
using CandidateFactory = QGA::CandidateFactory<Candidate>;


// The following initializations are needed in the .cpp due to language
// restrictions (ODR)

// Initialize the candidate counter
std::atomic_ulong QGA::CandidateCounter::count{0};

// Define the weights field
template<>
std::vector<unsigned> CandidateFactory::weights{};

// Choose the genetic operators
template<>
const std::vector<
  std::pair<CandidateFactory::GenOp,
  std::string>>
CandidateFactory::func {
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
    // trimming not necessary: this is usually about 10
    //nondom.randomTrim(Config::popSize);
    size_t nd = nondom.size();

    /* Top up to popSize2 candidates in parallel */
    Population pop2{Config::popSize2};
    pop.precompute();
    CandidateFactory cf{pop};
    pop2.add(Config::popSize2 - nd,
            [&]() -> const Candidate { return cf.getNew(); });

    /* Merge the nondominated subset of the previous population */
    pop2.add(nondom);
    pop = std::move(pop2);

    /* Take a record which GenOps were successful in making good candidates */
    for(auto& c : pop.front().randomSelect(Config::popSize))
      CandidateFactory::hit(c.getOrigin());

    /* Leave only one representative of each fitness */
    pop.prune([](const GenCandidate& a, const GenCandidate& b) -> bool {
      return a.fitness() == b.fitness();
      });

    /* Summarize */
    nondom = pop.front();
    std::cout << "Gen " << gen << ": "
      << pop.size() << " unique fitnesses, "
      << nondom.size() << " nondominated";
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
  std::cout << std::endl << "Run took " << dur.count() << " s"
    << " (" << dur.count()/Config::nGen << " s/gen avg), "
    << QGA::CandidateCounter::total() << " candidates tested, "
    << "best of run:" << std::endl;

  /* Dump the heuristic distribution */
  std::cout << "\nGenetic operator distribution:\n";
  CandidateFactory::dumpWeights(std::cout);

  /* List results */
  auto nondom = pop.front();
  std::vector<GenCandidate> vec{};

  /* Sort by error, from high to low */
  vec.reserve(nondom.size());
  for(auto& c : nondom)
    vec.push_back(c);
  std::sort(vec.begin(), vec.end(),
      [](const GenCandidate& a, const GenCandidate&b ) -> bool {
        return a.fitness().error > b.fitness().error;
      }
  );
  std::cout << '\n' << vec.size() << " nondominated candidates:\n";
  for(auto& c : vec) {
    std::cout << c.fitness() << ' ' << c << ": " << c.dump();
  }
}
