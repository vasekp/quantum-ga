// For colourful printing of fitness and generation
template<class Candidate>
class BriefPrinter {

public:

  BriefPrinter(const Candidate& ref_): ref(ref_) { }

  friend std::ostream& operator<< (std::ostream& os, const BriefPrinter& fp) {
    os << Colours::green(fp.ref.fitness());
    if(fp.ref.getGen() != (size_t)(~0))
      os << Colours::blue(" [g", fp.ref.getGen(), "]");
    return os;
  }

private:

  const Candidate& ref;

}; // class BriefPrinter
