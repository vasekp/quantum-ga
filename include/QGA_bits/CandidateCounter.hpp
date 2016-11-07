namespace QGA {

class CandidateCounter {

public:

  CandidateCounter(): count{0} { };

  void hit() {
    count++;
  }

  unsigned long total() {
    return count;
  }

private:

  /* Total fitness() evaluation counter */
  std::atomic_ulong count;

}; // class CandidateCounter


namespace {

  // It's not the purest technique but one global variable here is acceptable
  CandidateCounter counter{};

} // anonymous namespace

} // namespace QGA
