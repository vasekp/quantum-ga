namespace QGA {

class CandidateCounter {

  /* Total fitness() evaluation counter */
  std::atomic_ulong count;

public:

  CandidateCounter(): count{0} { };

  void hit() {
    count++;
  }

  unsigned long total() {
    return count;
  }

}; // class CandidateCounter


extern CandidateCounter counter;

} // namespace QGA
