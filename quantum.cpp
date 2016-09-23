#include <iostream>
#include <csignal>

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
  const size_t arSize = 20;

  // Internal population size
  const size_t popSize = 500;

  // Number of generations (constant)
  const int nGen = 100;

  // Expected curcuit depth in 0th generation
  const float expLengthIni = 30;

  // Expected number of gates inserted in mutation
  const float expLengthAdd = 1.5;

  // Probability of single gate deletion
  const float pDeleteUniform = 0.10;

  // How much prior success of genetic ops should influence future choices
  const float heurFactor = 0.10;

  // How much each bit is likely to be a control bit at gate creation
  const float pControl = 0.25;

} // namespace Config


namespace SigComm {

  enum StopState {
    RUNNING,
    HANDLER,
    STOPPING
  };

  volatile sig_atomic_t stopped = RUNNING;

  std::atomic<std::chrono::duration<double>> timeOut{};

} // namespace SigComm


using Candidate = Wrapper::Candidate;
using Population = gen::NSGAPopulation<Candidate>;
using GenCandidate = gen::Candidate<Candidate>;
using CandidateFactory = QGA::CandidateFactory<Candidate>;


// Initialize the candidate counter
// Needs to appear in the .cpp
QGA::CandidateCounter QGA::counter{};


void int_handler(int);


int main() {
#ifdef BENCH
  gen::rng.seed(1);
  omp_set_num_threads(1);
#endif

  Colours::use = isatty(1);
  Wrapper::init();

  std::signal(SIGINT, int_handler);

  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();

  Population pop{Config::popSize, [&] { return CandidateFactory::genInit(); }};

  std::cout << std::fixed << std::setprecision(4);

  CandidateFactory::Selector sel = CandidateFactory::getInitSelector();

  int gen;
  for(gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do arSize */
    auto nondom = pop.front();
    // trimming not necessary: this is usually about 10
    //nondom.randomTrim(Config::arSize);
    size_t nd = nondom.size();

    /* Top up to popSize candidates in parallel */
    Population pop2{Config::popSize};
    pop.precompute();
    CandidateFactory cf{pop, sel};
    pop2.add(Config::popSize - nd,
            [&]() -> const Candidate { return cf.getNew(); });

    /* Merge the nondominated subset of the previous population */
    pop2.add(nondom);
    pop = std::move(pop2);

    /* Take a record which GenOps were successful in making good candidates */
    for(auto& c : pop.front())
      sel.hit(c.getOrigin());

    /* Leave only one representative of each fitness */
    pop.prune([](const GenCandidate& a, const GenCandidate& b) -> bool {
      return a.fitness() == b.fitness();
      });

    /* Summarize */
    nondom = pop.front();
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset()
      << Colours::yellow() << pop.size() << Colours::reset()
      << " unique fitnesses, "
      << Colours::yellow() << nondom.size() << Colours::reset()
      << " nondominated";
    if(nondom.size() > 0) {
      auto& e = nondom.randomSelect();
      std::cout << ", e.g. " << e.fitness() << ' ' << e;
    }
    std::cout << std::endl;

    if(SigComm::stopped == SigComm::STOPPING)
      break;
  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre - SigComm::timeOut.load();
  std::cout << std::endl << "Run took " << dur.count() << " s ("
    << Colours::blue() << dur.count()/gen
    << " s/gen " << Colours::reset() << "avg), "
    << Colours::blue() << QGA::counter.total() << Colours::reset()
    << " candidates tested" << std::endl;

  /* Dump the heuristic distribution */
  std::cout << "\nGenetic operator distribution:\n";
  sel.dump(std::cout);

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
  std::cout << '\n'
    << Colours::yellow() << vec.size() << Colours::reset()
    << " nondominated candidates:\n";
  for(auto& c : vec)
    std::cout << Colours::green() << c.fitness() << Colours::reset()
      << ' ' << c << ": " << c.dump(std::cout);
}


void int_handler(int) {
  if(SigComm::stopped == SigComm::HANDLER)
    return;
  else if(SigComm::stopped == SigComm::STOPPING)
    // we got stuck during a stop request (e.g., popSize too large)
    std::_Exit(1);
  SigComm::stopped = SigComm::HANDLER;
  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();
  std::cerr << "\n\nComputation stopped. Choose action:\n"
    << Colours::blue() << "a: " << Colours::reset() << "abort,\n"
    << Colours::blue() << "c: " << Colours::reset() << "continue,\n"
    << Colours::blue() << "s: " << Colours::reset() << "stop now and output results.\n";
  char c;
  do {
    std::cerr << "\nYour choice: ";
    std::cin >> c;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  } while(c != 'a' && c != 'c' && c != 's');
  switch(c) {
    case 'a':
      std::_Exit(1);
    case 'c':
      SigComm::stopped = SigComm::RUNNING;
      break;
    case 's':
      SigComm::stopped = SigComm::STOPPING;
      break;
  }
  post = std::chrono::steady_clock::now();
  SigComm::timeOut.store(SigComm::timeOut.load() + (post - pre));
}
