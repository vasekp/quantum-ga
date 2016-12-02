namespace QGA {

namespace internal {
  // Defined below
  template<class>
  class FullPrinter;
}

template<class Derived, class Gene, typename... Elements>
class CandidateBase {

protected:

  using FitnessMain = Fitness<Elements...>;
  using FitnessFull = Fitness<typename Gene::Counter, Elements...>;

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

  FitnessFull fitness() const {
    typename Gene::Counter cc{};
    for(const auto& g : gt)
      g->hit(cc);
    counter.hit();
    return {derived().fitness_main(), cc};
  }

  friend bool sameCirc(const CandidateBase& lhs, const CandidateBase& rhs) {
    auto& gt1 = lhs.gt;
    auto& gt2 = rhs.gt;
    if(gt1.size() != gt2.size())
      return false;
    for(size_t i = 0; i < gt1.size(); i++)
      if(!sameType(gt1[i], gt2[i]))
        return false;
    return true;
  }

  friend std::ostream& operator<< (std::ostream& os, const CandidateBase& c) {
    for(const auto& g : c.gt)
      os << g << ' ';
    return os;
  }

  internal::FullPrinter<Derived> full() const {
    return {derived()};
  }

  internal::CircuitPrinter circuit() const {
    internal::CircuitPrinter printer{Config::nBit};
    for(const auto& g : gt)
      printer.print(g);
    return printer;
  }

  static Derived read(const std::string str) {
    std::istringstream is{str};
    std::vector<Gene> gt{};
    Gene gene{};
    while(is >> gene)
      gt.push_back(std::move(gene));
    return {std::move(gt)};
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

protected:

  unsigned controls() const {
    unsigned controls = 0;
    for(const auto& g : gt)
      controls += g->controls();
    return controls;
  }

  static double trimError(double error) {
    // Ignore deviations of roughly 10^-5
    return (unsigned long)(error * (1UL<<16)) / (double)(1UL<<16);
  }

private:

  const Derived& derived() const {
    return static_cast<const Derived&>(*this);
  }

  std::vector<Gene> gt{};
  size_t origin = (size_t)(~0);
  unsigned long gen = (unsigned long)(~0);

}; // class CandidateBase


namespace internal {

template<class CandidateBase>
class FullPrinter {

public:

  FullPrinter(const CandidateBase& ref_): ref(ref_) { }

  friend std::ostream& operator<< (std::ostream& os, const FullPrinter& p) {
    return p.ref.print_full(os);
  }

private:

  const CandidateBase& ref;

}; // class FullPrinter<Derived>

} // namespace internal

} // namespace QGA
