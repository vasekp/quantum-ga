#include <iostream>
#include <csignal>

#ifdef USE_QPP
  #include "include/wrapper-qpp.hpp"
#elif defined USE_QICLIB
  #include "include/wrapper-qiclib.hpp"
#else
  #error Either USE_QPP or USE_QICLIB needed.
#endif

#include "genetic.hpp"
#include "include/commons.hpp"

#ifdef FOURIER
  #include "include/problem-fourier.hpp"
#elif defined(SEARCH)
  #include "include/problem-search.hpp"
#else
  #include "include/problem-simple.hpp"
#endif

namespace Config {

  // Circuit width (constant)
  const unsigned nBit = 3;

  // strength parameter of NSGA selection
  const double selectBias = 0.3;

  // Archive (external population) size
  const size_t arSize = 100;

  // Internal population size
  const size_t popSize = 2000;

  // Number of candidates to keep from parent generation
  const size_t popKeep = 300;

  // Number of generations (constant)
  const unsigned long nGen = std::numeric_limits<unsigned long>::max();

  // Expected curcuit depth in 0th generation
  const double expLengthIni = 30;

  // Expected number of gates inserted / modified / removed in mutation
  const double expMutationCount = 4.0;

  // Probability of removing a single gate in uniform deletion
  const double pChoiceUniform = 0.1;

  // Probability of a crossover at any given point
  const double pCrossUniform = 0.1;

  // How much prior success of genetic ops should influence future choices
  const double heurFactor = 1.0 / nGen;

  // How much each bit is likely to be a control bit at gate creation
  const double pControl = 0.25;

  // Size of random subset of candidates to list on demand at interrupt
  const size_t nIntList = 20;

} // namespace Config


namespace SigComm {

  enum StopState {
    RUNNING,
    INTERRUPTED,
    STOPPING
  };

  enum Response {
    CONTINUE,
    DUMP,
    RESTART,
    LIST,
    STOP
  };

  volatile sig_atomic_t state = RUNNING;

  std::chrono::duration<double> timeOut{};

} // namespace SigComm


using Candidate = Wrapper::Candidate;
using Population = gen::NSGAPopulation<Candidate>;
using GenCandidate = gen::Candidate<Candidate>;
using CandidateFactory = QGA::CandidateFactory<Candidate>;


// Initialize the candidate counter
// Needs to appear in the .cpp
QGA::CandidateCounter QGA::counter{};


void int_handler(int);
int int_response();
void dumpResults(Population&, CandidateFactory::Selector&,
    std::chrono::time_point<std::chrono::steady_clock>, unsigned long);
void listRandom(Population&);


int main() {
#ifdef BENCH
  gen::rng.seed(1);
  omp_set_num_threads(1);
#endif

  Colours::use = isatty(1);
  Wrapper::init();

  std::signal(SIGINT, int_handler);

  std::chrono::time_point<std::chrono::steady_clock>
    start{std::chrono::steady_clock::now()};

  Population pop{Config::popSize,
    [&] { return CandidateFactory::genInit().setGen(0); }};

  std::cout << std::fixed << std::setprecision(4);

  CandidateFactory::Selector sel = CandidateFactory::getInitSelector();

  unsigned long gen;
  for(gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do arSize */
    Population pop2 = pop.front();
    pop2.rankTrim(Config::arSize);

    /* Randomly select popKeep candidates for survival without modification */
    pop2.reserve(Config::popSize);
    pop2.add(pop.randomSelect(Config::popKeep));

    /* Top up to popSize candidates in parallel */
    CandidateFactory cf{pop, sel};
    pop.precompute();
    pop2.add(Config::popSize - pop2.size(),
        [&] { return cf.getNew().setGen(gen); });

    /* We don't need the original population anymore */
    pop = std::move(pop2);

    /* Leave only one representative of each fitness */
    pop.prune([](const GenCandidate& a, const GenCandidate& b) -> bool {
        return a.fitness() == b.fitness();
      }, 0, false);

    /* Take a record which GenOps were successful in making good candidates */
    auto nondom = pop.front();
    for(auto& c : nondom)
      if(c.getGen() == gen)
        sel.hit(c.getOrigin());

    /* Summarize */
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset()
      << Colours::yellow() << pop.size() << Colours::reset()
      << " unique fitnesses, lowest error "
      << Colours::green() << pop.best().fitness() << Colours::blue()
      << " [" << pop.best().getGen() << ']' << Colours::reset() << ", "
      << Colours::yellow() << nondom.size() << Colours::reset()
      << " nondominated, newest: ";
    auto& newest = *std::min_element(nondom.begin(), nondom.end(),
        [](const GenCandidate& c1, const GenCandidate& c2) {
          return c1.getGen() > c2.getGen();
        });
    std::cout
      << Colours::green() << newest.fitness() << Colours::blue()
      << " [" << newest.getGen() << ']' << Colours::reset() << std::endl;

    while(SigComm::state == SigComm::INTERRUPTED) {
      int res = int_response();
      switch(res) {
        case SigComm::DUMP:
          dumpResults(pop, sel, start, gen);
          break;
        case SigComm::LIST:
          listRandom(pop);
          break;
        case SigComm::RESTART:
          pop = Population{Config::popSize,
            [&] { return CandidateFactory::genInit().setGen(0); }};
          sel = CandidateFactory::getInitSelector();
          start = std::chrono::steady_clock::now();
          SigComm::timeOut = std::chrono::duration<double>(0);
          gen = 0;
          break;
      }
    }
    if(SigComm::state == SigComm::STOPPING)
      break;
  }

  dumpResults(pop, sel, start, gen);
}


void dumpResults(Population& pop, CandidateFactory::Selector& sel,
    std::chrono::time_point<std::chrono::steady_clock> start,
    unsigned long gen) {
  /* Timing information */
  std::chrono::time_point<std::chrono::steady_clock>
    now{std::chrono::steady_clock::now()};
  std::chrono::duration<double> dur = now - start - SigComm::timeOut;
  std::cout << "\nRun took " << dur.count() << " s ("
    << Colours::blue() << dur.count()/gen
    << " s/gen " << Colours::reset() << "avg), "
    << Colours::blue() << QGA::counter.total() << Colours::reset()
    << " candidates tested\n";

  /* List results */
  auto nondom = pop.front();
  nondom.sort();
  std::cout << '\n'
    << Colours::yellow() << nondom.size() << Colours::reset()
    << " nondominated candidates:\n";
  for(auto& c : nondom.reverse()) {
    std::cout << Colours::green() << c.fitness() << Colours::reset()
      << " [" << Colours::blue() << 'g' << c.getGen() << Colours::reset()
      << "] " << c;
    if(c.fitness().error < 0.01)
      std::cout << ": " << c.dump(std::cout);
    else
      std::cout << '\n';
  }

  /* Dump the heuristic distribution */
  std::cout << "\nGenetic operator distribution:\n";
  sel.dump(std::cout);
}


void listRandom(Population& pop) {
  auto sel = pop.randomSelect(Config::nIntList);
  sel.sort();
  for(auto& c : sel.reverse())
    std::cout << Colours::green() << c.fitness() << Colours::reset()
      << " [" << Colours::blue() << 'g' << c.getGen() << Colours::reset()
      << "] " << c << '\n';
}


void int_handler(int) {
  if(SigComm::state != SigComm::RUNNING)
    // we got stuck while processing another signal (e.g., popSize too large
    // or a deadlock)
    std::_Exit(1);
  SigComm::state = SigComm::INTERRUPTED;
}


int int_response() {
  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();
  std::cerr << "\nComputation stopped. Choose action:\n"
    << Colours::blue() << "a: " << Colours::reset()
      << "abort,\n"
    << Colours::blue() << "c: " << Colours::reset()
      << "continue,\n"
    << Colours::blue() << "d: " << Colours::reset()
      << "diagnose / list current results,\n"
    << Colours::blue() << "l: " << Colours::reset()
      << "list " << Config::nIntList << " random candidates,\n"
    << Colours::blue() << "r: " << Colours::reset()
      << "restart,\n"
    << Colours::blue() << "q: " << Colours::reset()
      << "quit after this generation.\n";
  int ret = -1;
  do {
    char c;
    std::cerr << "\nYour choice: ";
    std::cin >> c;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(std::cin.eof())
      c = 'a';
    switch(c) {
      case 'a':
        std::_Exit(1);
      case 'c':
        SigComm::state = SigComm::RUNNING;
        ret = SigComm::CONTINUE;
        break;
      case 'd':
        ret = SigComm::DUMP;
        break;
      case 'l':
        ret = SigComm::LIST;
        break;
      case 'r':
        SigComm::state = SigComm::RUNNING;
        ret = SigComm::RESTART;
        break;
      case 'q':
        SigComm::state = SigComm::STOPPING;
        ret = SigComm::STOP;
        break;
    }
  } while(ret < 0);
  post = std::chrono::steady_clock::now();
  SigComm::timeOut += post - pre;
  return ret;
}
