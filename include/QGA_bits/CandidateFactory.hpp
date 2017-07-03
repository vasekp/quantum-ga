namespace QGA {

template<
  class Candidate,
  class Population = gen::NSGAPopulation<Candidate>>
class CandidateFactory {

  using Gene = typename Candidate::GeneType;
  // private members declared at bottom

public:

  class Tracker;

  CandidateFactory(Population& pop_, Tracker& trk_): pop(pop_), trk(trk_) { }

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
    auto op = trk.select();
    return (this->*op.first)().setOrigin(op.second);
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
      gtNew[dPos(gen::rng)] = Gene::getRandom();
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
    ins.reserve(2*Config::expMutationCount);
    double probTerm = 1/Config::expMutationCount;
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
    ins.reserve(2*Config::expMutationCount);
    double probTerm = 1/Config::expMutationCount;
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
    gOrig.mutate();
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
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
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
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
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
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expMutationCount);
    double probTerm = 1/Config::expMutationCount;
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
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
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
    std::geometric_distribution<size_t> dGeom{1.0 / Config::expMutationCount};
    size_t pos1 = dPos(gen::rng),
           len = 1 + dGeom(gen::rng),
           pos2 = pos1 + len > sz - 1 ? sz - 1 : pos1 + len;
    std::vector<Gene> gtNew = gtOrig;
    swap(gtNew[pos1], gtNew[pos2]);
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

public:

  class Tracker {

    using CF = CandidateFactory;
    using FunPtr = Candidate (CF::*)();
    // private members declared at bottom

  public:

    struct GenOp {

      FunPtr fun;
      std::string name;
      unsigned long hits;

      GenOp(FunPtr fun_, std::string name_):
        fun(fun_), name(name_), hits(0) { }

    };

    void hit(size_t ix) {
      if(ix >= 0 && ix < count)
        ops[ix].hits++;
    }

    friend std::ostream& operator<< (std::ostream& os, const Tracker& trk) {
      /* Find the longest GenOp name */
      auto max = std::max_element(trk.ops.begin(), trk.ops.end(),
          [](const GenOp& a, const GenOp& b) {
            return a.name.length() < b.name.length();
          });
      auto maxw = max->name.length();

      /* Preserve settings of os */
      auto flags_ = os.flags(std::ios_base::left);

      /* List all op names and probabilities */
      for(auto& op : trk.ops)
        os << std::setw(maxw+3) << op.name + ':' << op.hits << '\n';

      os.flags(flags_);
      return os;
    }

    std::pair<FunPtr, size_t> select() {
      size_t index = dUni(gen::rng);
      return {ops[index].fun, index};
    }

    Tracker(std::vector<GenOp>&& ops_):
    ops(std::move(ops_)), count(ops.size()), dUni(0, count - 1) { }

  private:

    std::vector<GenOp> ops;
    size_t count;
    std::uniform_int_distribution<> dUni;

  }; // class Tracker

  static Tracker getInitTracker() {
    using CF = CandidateFactory;
    std::vector<typename Tracker::GenOp> ops{};
  //ops.push_back({ &CF::mAlterDiscrete,   "MDiscrete" });
    ops.push_back({ &CF::mAlterContinuous, "MutSingle" });
    ops.push_back({ &CF::mAddSlice,        "AddSlice" });
  //ops.push_back({ &CF::mAddPairs,        "AddPairs" });
    ops.push_back({ &CF::mMutateAddPair,   "MutAddPair" });
    ops.push_back({ &CF::mSwapQubits,      "SwapQubits" });
    ops.push_back({ &CF::mDeleteSlice,     "DelShort" });
    ops.push_back({ &CF::mDeleteUniform,   "DelUnif"  });
    ops.push_back({ &CF::mReplaceSlice,    "ReplSlice" });
    ops.push_back({ &CF::mSplitSwap,       "SpltSwp"  });
    ops.push_back({ &CF::mReverseSlice,    "InvSlice" });
  //ops.push_back({ &CF::mPermuteSlice,    "PermSlice" });
    ops.push_back({ &CF::mSwapTwo,         "SwapTwo" });
    ops.push_back({ &CF::mRepeatSlice,     "ReptSlice" });
    ops.push_back({ &CF::crossoverUniform, "C/Over"   });
  //ops.push_back({ &CF::concat3,          "Concat3"  });
    ops.push_back({ &CF::simplify,         "Simplify" });
    return {std::move(ops)};
  }

private:

  Population& pop;
  Tracker& trk;

}; // class CandidateFactory

} // namespace QGA
