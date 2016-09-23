#include <iostream>
#include <csignal>

#include "genetic.hpp"

#include "include/commons.hpp"

/*#ifdef USE_QPP
#include "include/wrapper-qpp.hpp"
#elif defined USE_QICLIB*/
#include "include/wrapper-fourier.hpp"
/*#else
#error Either USE_QPP or USE_QICLIB needed.
#endif*/

namespace Config {

  // Circuit width (constant)
  const unsigned nBit = 2;

  // strength parameter of NSGA selection
  const float selectBias = 1.7;

  // Archive (external population) size
  const size_t arSize = 100;

  // Internal population size
  const size_t popSize = 2000;

  // Number of generations (constant)
  const int nGen = 500;

  // Expected curcuit depth in 0th generation
  const float expLengthIni = 30;

  // Expected number of gates inserted in mutation
  const float expLengthAdd = 1.5;

  // Probability of altering / deleting a single gate in uniform operators
  const float pChoiceUniform = 0.03;

  // Probability of a crossover at any given point
  const float pCrossUniform = 0.1;

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

  Population pop{Config::popSize, [&] { return CandidateFactory::genInit().setGen(0); }};

  std::cout << std::fixed << std::setprecision(4);

  CandidateFactory::Selector sel = CandidateFactory::getInitSelector();

  int gen;
  for(gen = 0; gen < Config::nGen; gen++) {

    /* Find the nondominated subset and trim down do arSize */
    Population pop2 = pop.front();
    pop2.randomTrim(Config::arSize);
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
        return a.fitness() == b.fitness() ||
          std::abs(a.fitness().error - b.fitness().error) < 0.0001;
      }, 0, false);

    /* Take a record which GenOps were successful in making good candidates */
    auto nondom = pop.front();
    for(auto& c : nondom)
      if(c.getGen() == gen)
        sel.hit(c.getOrigin());

    /* Summarize */
    auto minerr = std::min_element(nondom.begin(), nondom.end(),
        [](const GenCandidate& a, const GenCandidate& b) -> bool {
          return a.fitness().error + 0.05*a.fitness().length
               < b.fitness().error + 0.05*b.fitness().length; });
    std::cout << Colours::bold() << "Gen " << gen << ": " << Colours::reset()
      << Colours::yellow() << pop.size() << Colours::reset()
      << " unique fitnesses, lowest err+len "
      << Colours::green() << minerr->fitness() << Colours::reset() << ", "
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
  for(auto& c : vec) {
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
