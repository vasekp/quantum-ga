namespace QGA {

// Forward, full declaration in GenOpCounter.hpp
template<class> class GenOpCounter;

namespace internal {

  class ConstExprString {

  public:

      constexpr ConstExprString(const char * const str) :
        ntstr(str), len(ntlen(str)) { }

      constexpr operator const char* () const {
        return ntstr;
      }

      constexpr int length() const {
        return len;
      }

  private:

    constexpr int ntlen(const char * const ntstr) {
      return *ntstr == 0 ? 0 : 1 + ntlen(ntstr + 1);
    }

    const char * ntstr;
    int len;

  };

} // namespace internal

template<
  class Candidate,
  class Population = gen::NSGAPopulation<Candidate>>
class CandidateFactory {

  using Gene = typename Candidate::GeneType;

  friend class GenOpCounter<CandidateFactory>;

  // private members declared at bottom

public:

  CandidateFactory(Population& pop_): pop(pop_) { }

  static Candidate genInit() {
    // probability of termination; expLengthIni = expected number of genes
    const double probTerm = 1/Config::expLengthIni;
    std::uniform_real_distribution<> dUni{};
    std::vector<Gene> gtOrig{};
    gtOrig.reserve(Config::expLengthIni);
    do
      gtOrig.push_back(Gene::getRandom());
    while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gtOrig)};
  }

  Candidate getNew() {
    std::uniform_int_distribution<size_t> dUni{0, ops.size() - 1};
    size_t index = dUni(gen::rng);
    return (this->*ops[index].fun)().setOrigin(index);
  }

private:

  const Candidate& get() {
    return pop.NSGASelect(Config::selectBias);
  }

  Candidate mAlterDiscrete() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0)
      return parent;
    auto gtNew = gtOrig;
    std::uniform_real_distribution<> dUni{};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    const double probTerm = 1/Config::expMutationCount;
    do
      gtNew[dPos(gen::rng)].getAnother();
    while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gtNew)};
  }

  Candidate mAlterContinuous() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0)
      return parent;
    auto gtNew = gtOrig;
    std::uniform_real_distribution<> dUni{};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    const double probTerm = 1/Config::expMutationCount;
    do
      gtNew[dPos(gen::rng)].mutate();
    while(dUni(gen::rng) > probTerm);
    return gtNew != gtOrig ? Candidate{std::move(gtNew)} : parent;
  }

  Candidate mAddSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    std::uniform_real_distribution<> dUni{};
    std::uniform_int_distribution<size_t> dPos{0, sz};
    size_t pos = dPos(gen::rng);
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expSliceLength);
    double probTerm = 1/Config::expSliceLength;
    do
      ins.push_back(Gene::getRandom());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz + ins.size());
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos);
    gtNew.insert(gtNew.end(), ins.begin(), ins.end());
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mAddPairs() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    std::uniform_real_distribution<> dUni{};
    std::uniform_int_distribution<size_t> dPos{0, sz};
    size_t pos1 = dPos(gen::rng),
           pos2 = dPos(gen::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expSliceLength);
    double probTerm = 1/Config::expSliceLength;
    do
      ins.push_back(Gene::getRandom());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz + 2*ins.size());
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), ins.begin(), ins.end());
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos1, gtOrig.begin() + pos2);
    for(auto& g : ins)
      g.invert();
    gtNew.insert(gtNew.end(),
        std::make_move_iterator(ins.rbegin()),
        std::make_move_iterator(ins.rend()));
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mMutateAddPair() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    size_t pos = dPos(gen::rng);
    Gene gOrig{gtOrig[pos]};
    gOrig.getAnother();
    Gene gNew{Gene::getRandom()};
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz + 2);
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos);
    gtNew.push_back(gNew);
    gtNew.push_back(std::move(gOrig));
    Gene gNewInv{gNew};
    gNewInv.invert();
    gtNew.push_back(std::move(gNewInv));
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos + 1, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mSwapQubits() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0 || Config::nBit < 2)
      return parent;
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expSliceLength};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::uniform_int_distribution<unsigned> dBit{0, Config::nBit - 2};
    unsigned s1 = dBit(gen::rng),
             s2 = dBit(gen::rng);
    // ensure that the two qubit indices are unequal
    s2 += s2 >= s1;
    std::vector<Gene> gtNew{gtOrig};
    for(size_t pos = pos1; pos < pos2; pos++)
      gtNew[pos].swapQubits(s1, s2);
    return Candidate{std::move(gtNew)};
  }

  Candidate mDeleteSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0)
      return parent;
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expSliceLength};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz - (pos2 - pos1));
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mReplaceSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz == 0)
      return parent;
    std::uniform_real_distribution<> dUni{};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expSliceLength};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expSliceLength);
    double probTerm = 1/Config::expSliceLength;
    do
      ins.push_back(Gene::getRandom());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz - (pos2 - pos1) + ins.size());
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), ins.begin(), ins.end());
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mDeleteUniform() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    std::uniform_real_distribution<> dUni{};
    std::vector<Gene> gtNew{};
    gtNew.reserve(gtOrig.size());
    size_t cnt = 0;
    double prob = double(Config::expMutationCount) / sz;
    for(auto& g : gtOrig)
      if(dUni(gen::rng) >= prob)
        gtNew.push_back(g);
      else
        cnt++;
    return cnt ? Candidate{std::move(gtNew)} : parent;
  }

  Candidate mSplitSwap() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 2};
    std::array<size_t, 4> pos;
    for(auto& p : pos)
      p = dPos(gen::rng);
    std::sort(pos.begin(), pos.end());
    // ensure that pos[1]-pos[0] and pos[3]-pos[2] are nonzero
    pos[1]++, pos[2]++, pos[3] += 2;
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz);
    gtNew.insert(gtNew.end(),
        gtOrig.begin(), gtOrig.begin() + pos[0]);
    gtNew.insert(gtNew.end(),
        gtOrig.begin() + pos[2], gtOrig.begin() + pos[3]);
    gtNew.insert(gtNew.end(),
        gtOrig.begin() + pos[1], gtOrig.begin() + pos[2]);
    gtNew.insert(gtNew.end(),
        gtOrig.begin() + pos[0], gtOrig.begin() + pos[1]);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos[3], gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mReverseSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 2};
    size_t pos1 = dPos(gen::rng),
           pos2 = dPos(gen::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    // ensure that pos2-pos1 is at least 2
    pos2 += 2;
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz);
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    {
      auto prev_end = gtNew.end();
      gtNew.insert(gtNew.end(),
          gtOrig.rbegin() + sz - pos2, gtOrig.rbegin() + sz - pos1);
      auto end = gtNew.end();
      for(auto it = prev_end; it != end; it++)
        it->invert();
    }
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mPermuteSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 2};
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expSliceLength};
    size_t pos1 = dPos(gen::rng),
           len = 2 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gtNew = gtOrig;
    std::shuffle(gtNew.begin() + pos1, gtNew.begin() + pos2, gen::rng);
    return Candidate{std::move(gtNew)};
  }

  Candidate mSwapTwo() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 2};
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expSliceLength};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz - 1 ? sz - 1 : pos1 + len;
    std::vector<Gene> gtNew = gtOrig;
    swap(gtNew[pos1], gtNew[pos2]);
    return Candidate{std::move(gtNew)};
  }
  
  Candidate mMoveGate() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 2};
    size_t pos1 = dPos(gen::rng),
           pos2 = dPos(gen::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    // ensure that pos2-pos1 is at least 1
    pos2 += 1;
    std::bernoulli_distribution dir{};
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz);
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    if(dir(gen::rng)) { // move first to end
      gtNew.insert(gtNew.end(), gtOrig.begin() + pos1 + 1, gtOrig.begin() + pos2);
      gtNew.insert(gtNew.end(), gtOrig.begin() + pos1, gtOrig.begin() + pos1 + 1);
    } else { // move last to beginning
      gtNew.insert(gtNew.end(), gtOrig.begin() + pos2 - 1, gtOrig.begin() + pos2);
      gtNew.insert(gtNew.end(), gtOrig.begin() + pos1, gtOrig.begin() + pos2 - 1);
    }
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mRepeatSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    if(sz < 2)
      return parent;
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    size_t pos1 = dPos(gen::rng),
           pos2 = dPos(gen::rng);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    // ensure that pos2-pos1 is at least 1
    pos2 += 1;
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz + pos2 - pos1);
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos1, gtOrig.begin() + pos2);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos1, gtOrig.begin() + pos2);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate crossoverUniform() {
    auto &parent1 = get(),
         &parent2 = get();
    auto *gt1 = &parent1.genotype(),
         *gt2 = &parent2.genotype();
    // Most of these are here just for clarity and will hopefully be
    // optimized away
    size_t sz1 = gt1->size(),
           sz2 = gt2->size(),
           pos1 = 0, pos2 = 0;
    double pCross1 = std::min(Config::expMutationCount / sz1, 1.0),
           pCross2 = std::min(Config::expMutationCount / sz2, 1.0);
    std::geometric_distribution<size_t> dGeom1{pCross1};
    std::geometric_distribution<size_t> dGeom2{pCross2};
    std::vector<Gene> gtNew{};

    gtNew.reserve(std::max(sz1, sz2));
    for(;;) {
      // Take roughly expLen1 genes from gt1
      size_t upto = pos1 + dGeom1(gen::rng) + 1;
      if(upto >= sz1)
        break; // just use the rest of gt1
      // Skip roughly expLen2 genes from gt2
      pos2 += dGeom2(gen::rng) + 1;
      if(pos2 >= sz2)
        break; // ditto
      gtNew.insert(gtNew.end(), gt1->begin() + pos1, gt1->begin() + upto);
      pos1 = upto;
      // Swap the two
      std::swap(gt1, gt2);
      std::swap(sz1, sz2);
      std::swap(pos1, pos2);
      std::swap(dGeom1, dGeom2);
    }
    // If we get here then either more was requested of gt1 than available or
    // gt2 went empty. In either case, we just take whatever's left and we
    // finish the crossover operation.
    gtNew.insert(gtNew.end(), gt1->begin() + pos1, gt1->end());

    return Candidate{std::move(gtNew)};
  }

  Candidate concat3() {
    auto &parent1 = get(),
         &parent2 = get(),
         &parent3 = get();
    auto &gt1 = parent1.genotype(),
         &gt2 = parent2.genotype(),
         &gt3 = parent3.genotype();
    std::vector<Gene> gtNew{};
    gtNew.reserve(gt1.size() + gt2.size() + gt3.size());
    gtNew.insert(gtNew.end(), gt1.begin(), gt1.end());
    {
      auto start = gtNew.end();
      gtNew.insert(gtNew.end(), gt2.rbegin(), gt2.rend());
      auto end = gtNew.end();
      for(auto it = start; it != end; it++)
        it->invert();
    }
    gtNew.insert(gtNew.end(), gt3.begin(), gt3.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate simplify() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    size_t sz = gtOrig.size();
    if(sz == 0)
      return parent;
    std::vector<Gene> gtNew = gtOrig;
    for(auto& g : gtNew)
      g.simplify();
    return gtNew != gtOrig ? Candidate{std::move(gtNew)} : parent;
  }

  struct GenOp {

    using FunPtr = Candidate (CandidateFactory::*)();
    FunPtr fun;
    internal::ConstExprString name;

  };

  static constexpr std::array<GenOp, 14> ops{{
    { &CandidateFactory::mAlterDiscrete,   "MDiscrete" },
    { &CandidateFactory::mAlterContinuous, "MContns" },
    { &CandidateFactory::mAddSlice,        "AddSlice" },
  //{ &CandidateFactory::mAddPairs,        "AddPairs" },
    { &CandidateFactory::mMutateAddPair,   "MutAddPair" },
    { &CandidateFactory::mSwapQubits,      "SwapQubits" },
    { &CandidateFactory::mDeleteSlice,     "DelShort" },
    { &CandidateFactory::mDeleteUniform,   "DelUnif"  },
    { &CandidateFactory::mReplaceSlice,    "ReplSlice" },
    { &CandidateFactory::mSplitSwap,       "SpltSwp"  },
    { &CandidateFactory::mReverseSlice,    "InvSlice" },
  //{ &CandidateFactory::mPermuteSlice,    "PermSlice" },
  //{ &CandidateFactory::mSwapTwo,         "SwapTwo" },
    { &CandidateFactory::mMoveGate,        "MoveGate" },
    { &CandidateFactory::mRepeatSlice,     "ReptSlice" },
    { &CandidateFactory::crossoverUniform, "C/Over"   },
  //{ &CandidateFactory::concat3,          "Concat3"  },
    { &CandidateFactory::simplify,         "Simplify" }
  }};

  Population& pop;

}; // class CandidateFactory

} // namespace QGA
