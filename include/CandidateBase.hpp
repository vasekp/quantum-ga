template<class Derived, class Gene>
class CandidateBase {

protected:

  std::vector<Gene> gt{};

private:

  int origin = -1;

  const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

public:

  using GeneType = Gene;

  CandidateBase(std::vector<Gene>&) = delete;

  CandidateBase(std::vector<Gene>&& gt_): gt(std::move(gt_)) { }

  NOINLINE Fitness fitness() const {
    /* Complexity = square sum of numbers of control bits per gate */
    size_t cplx = 0;
    for(const Gene& g : gt) {
      unsigned h = g.weight();
      cplx += h*h;
    }
    CandidateCounter::hit();
    return {derived().error(), gt.size(), cplx};
  }

  friend std::ostream& operator<< (std::ostream& os, const CandidateBase& c) {
    auto first = c.gt.begin(), last = c.gt.end();
    for(auto it = first; it != last; it++) {
      if(it != first) os << ' ';
      os << Wrapper::gate_name(it->gate()) << it->target()+1;
      auto& ixv = it->ix_vector();
      if(ixv.size()) {
        os << '[';
        for(auto ix : ixv)
          os << ix;
        os << ']';
      }
    }
    return os;
  }

  const std::vector<Gene>& genotype() const {
    return gt;
  }

  void setOrigin(int _origin) {
    origin = _origin;
  }

  int getOrigin() const {
    return origin;
  }

  template<class, class, class>
  friend class CandidateFactory;

}; // class CandidateBase
