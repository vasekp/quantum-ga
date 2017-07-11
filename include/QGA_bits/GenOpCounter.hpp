namespace QGA {

template<class CandidateFactory>
class GenOpCounter {

  using GenOp = typename CandidateFactory::GenOp;
  // private members declared at bottom

public:

  void hit(size_t ix) {
    if(ix < hits.size())
      hits[ix]++;
  }

  void reset() {
    std::fill(hits.begin(), hits.end(), 0);
  }

  friend std::ostream& operator<< (std::ostream& os, const GenOpCounter& trk) {
    /* Find the longest GenOp name */
    auto max = std::max_element(ops.begin(), ops.end(),
        [](const GenOp& a, const GenOp& b) {
          return a.name.length() < b.name.length();
        });
    auto maxw = max->name.length();

    /* Preserve settings of os */
    auto flags_ = os.flags(std::ios_base::left);

    /* List all op names and probabilities */
    for(size_t ix = 0; ix < ops.size(); ix++)
      os << ops[ix].name
         << std::setw(maxw + 3 - ops[ix].name.length()) << ':'
         << trk.hits[ix] << '\n';

    os.flags(flags_);
    return os;
  }

private:

  static constexpr auto& ops = CandidateFactory::ops;
  std::array<size_t, ops.size()> hits{};

}; // class GenOpCounter

} // namespace QGA
