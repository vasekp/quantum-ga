#include <iostream>
#include <csignal>

#include "genetic.hpp"

#include "include/commons.hpp"

#ifdef FOURIER
  #ifdef USE_QICLIB
    #include "include/wrapper-fourier.hpp"
  #else
    #error Fourier problem is only implemented using QIClib.
  #endif
#else
  #ifdef USE_QPP
    #include "include/wrapper-qpp.hpp"
  #elif defined USE_QICLIB
    #include "include/wrapper-qiclib.hpp"
  #else
    #error Either USE_QPP or USE_QICLIB needed.
  #endif
#endif

namespace Config {

  // Circuit width (constant)
  const unsigned nBit = 3;

  // strength parameter of NSGA selection
  const float selectBias = 2.0;

  // Archive (external population) size
  const size_t arSize = 100;

  // Internal population size
  const size_t popSize = 2000;

  // Number of generations (constant)
  const size_t nGen = std::numeric_limits<size_t>::max();

  // Expected curcuit depth in 0th generation
  const float expLengthIni = 30;

  // Expected number of gates inserted in mutation
  const float expLengthAdd = 1.5;

  // Probability of altering / deleting a single gate in uniform operators
  const float pChoiceUniform = 0.03;

  // Probability of a crossover at any given point
  const float pCrossUniform = 0.1;

  // How much prior success of genetic ops should influence future choices
  const float heurFactor = 1.0 / nGen;

  // How much each bit is likely to be a control bit at gate creation
  const float pControl = 0.25;

} // namespace Config


namespace SigComm {

  enum StopState {
    RUNNING,
    INTERRUPTED
  };

  enum Response {
    CONTINUE,
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
SigComm::Response int_response();
void list(Population&, CandidateFactory::Selector&);


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

  Population pop{Config::popSize, [&] { return CandidateFactory::genInit().setGen(0); }};

  std::cout << std::fixed << std::setprecision(4);

  CandidateFactory::Selector sel = CandidateFactory::getInitSelector();

  size_t gen;
  for(gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do arSize */
    Population pop2 = pop.front();
    pop2.rankTrim(Config::arSize);
    size_t nd = pop2.size();

    /* Top up to popSize candidates in parallel */
    pop2.reserve(Config::popSize);
    pop.precompute();
    CandidateFactory cf{pop, sel};
    pop2.add(Config::popSize - nd,
            [&]() -> const Candidate { return cf.getNew().setGen(gen); });

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
      << " nondominated";
    if(nondom.size() > 0) {
      auto& e = nondom.randomSelect();
      std::cout << ", e.g. " << e.fitness() << ' ' << e;
    }
    std::cout << std::endl;

    if(SigComm::state != SigComm::RUNNING) {
      SigComm::Response res = int_response();
      while(res == SigComm::LIST) {
        list(pop, sel);
        res = int_response();
      }
      if(res == SigComm::STOP)
        break;
    }
  }

  post = std::chrono::steady_clock::now();
  std::chrono::duration<double> dur = post - pre - SigComm::timeOut;
  std::cout << std::endl << "Run took " << dur.count() << " s ("
    << Colours::blue() << dur.count()/gen
    << " s/gen " << Colours::reset() << "avg), "
    << Colours::blue() << QGA::counter.total() << Colours::reset()
    << " candidates tested" << std::endl;

  list(pop, sel);
}


void list(Population& pop, CandidateFactory::Selector& sel) {
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


void int_handler(int) {
  if(SigComm::state != SigComm::RUNNING)
    // we got stuck while processing another signal (e.g., popSize too large
    // or a deadlock)
    std::_Exit(1);
  SigComm::state = SigComm::INTERRUPTED;
}


SigComm::Response int_response() {
  std::chrono::time_point<std::chrono::steady_clock> pre, post;
  pre = std::chrono::steady_clock::now();
  std::cerr << "\nComputation stopped. Choose action:\n"
    << Colours::blue() << "a: " << Colours::reset() << "abort,\n"
    << Colours::blue() << "c: " << Colours::reset() << "continue,\n"
    << Colours::blue() << "d: " << Colours::reset() << "diagnose / list current results,\n"
    << Colours::blue() << "q: " << Colours::reset() << "quit after this generation.\n";
  char c;
  do {
    std::cerr << "\nYour choice: ";
    std::cin >> c;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    if(std::cin.eof()) {
      c = 'a';
      break;
    }
  } while(c != 'a' && c != 'c' && c != 'd' && c != 'q');
  post = std::chrono::steady_clock::now();
  SigComm::timeOut += post - pre;
  switch(c) {
    case 'a':
      std::_Exit(1);
    case 'c':
      SigComm::state = SigComm::RUNNING;
      return SigComm::CONTINUE;
      break;
    case 'd':
      return SigComm::LIST;
      break;
    case 'q':
    default:
      return SigComm::STOP;
      break;
  }
}
