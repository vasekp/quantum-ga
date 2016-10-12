namespace QGA {

template<class Derived, class Gene>
class CandidateBase {

protected:

  std::vector<Gene> gt{};

private:

  size_t origin = (size_t)(~0);
  unsigned long gen = (unsigned long)(~0);

  const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

public:

  using GeneType = Gene;

  CandidateBase(std::vector<Gene>&& gt_):
      gt(std::move(gt_)) {
    if(gt.size() == 0)
      return;
    auto end = gt.end(), last = gt.begin();
    for(auto cur = last + 1; cur != end; cur++) {
      // Can be merged: done, go to next cur
      // Can not: put *cur after *last and increase both
      bool consumed = (*last).merge(*cur);
      if(!consumed)
        std::swap(*++last, *cur);
    }
    gt.erase(++last, gt.end());
  }

  NOINLINE Fitness fitness() const {
    unsigned cplx{0};
    for(const auto& g : gt)
      cplx += g->complexity();
    counter.hit();
    return {trimError(derived().error()), gt.size(), cplx};
  }

  friend std::ostream& operator<< (std::ostream& os, const CandidateBase& c) {
    for(const auto& g : c.gt)
      os << g << ' ';
    return os;
  }

  const std::vector<Gene>& genotype() const {
    return gt;
  }

  Derived& setOrigin(size_t origin_) {
    if(origin == (size_t)(~0))
      origin = origin_;
    return static_cast<Derived&>(*this);
  }

  size_t getOrigin() const {
    return origin;
  }

  Derived& setGen(unsigned long gen_) {
    if(gen == (unsigned long)(~0))
      gen = gen_;
    return static_cast<Derived&>(*this);
  }

  unsigned long getGen() const {
    return gen;
  }

  template<class, class>
  friend class CandidateFactory;

protected:

  static double trimError(double error) {
    // Ignore deviations of roughly 10^-7
    return (unsigned long)(error * (1UL<<24)) / (double)(1UL<<24);
  }

}; // class CandidateBase

} // namespace QGA
