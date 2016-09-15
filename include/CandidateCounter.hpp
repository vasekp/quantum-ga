class CandidateCounter {

  /* Total fitness() evaluation counter */
  static std::atomic_ulong count;

public:

  static void hit() {
    count++;
  }

  static unsigned long total() {
    return count;
  }

}; // class CandidateCounter
