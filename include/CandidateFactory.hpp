namespace QGA {


// Defined below CandidateFactory
template<class Candidate, class Population>
class CFSelector;


template<
  class Candidate,
  class Population = gen::NSGAPopulation<Candidate>>
class CandidateFactory {

  using Gene = typename Candidate::GeneType;

public:

  friend class CFSelector<Candidate, Population>;
  using Selector = CFSelector<Candidate, Population>;

private:

  Population& pop;
  Selector& sel;

public:

  CandidateFactory(Population& pop_, Selector& sel_): pop(pop_), sel(sel_) {
    sel.update();
  }

  static Selector getInitSelector() {
    return Selector{};
  }

  static Candidate genInit() {
    // probability of termination; expLengthIni = expected number of genes
    const double probTerm = 1/Config::expLengthIni;
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gtOrig{};
    gtOrig.reserve(Config::expLengthIni);
    do
      gtOrig.push_back(Gene::getNew());
    while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gtOrig)};
  }

  Candidate getNew() {
    auto op = sel.select();
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
    std::uniform_real_distribution<> dUni{0, 1};
    std::uniform_int_distribution<size_t> dPos{0, sz - 1};
    const double probTerm = 1/Config::expMutationCount;
    do
      gtNew[dPos(gen::rng)] = Gene::getNew();
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
    std::uniform_real_distribution<> dUni{0, 1};
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
    std::uniform_real_distribution<> dUni{0, 1};
    size_t pos = gen::rng() % (sz + 1);
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expMutationCount);
    double probTerm = 1/Config::expMutationCount;
    do
      ins.push_back(Gene::getNew());
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
    std::uniform_real_distribution<> dUni{0, 1};
    size_t pos1 = gen::rng() % (sz + 1),
           pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins{};
    ins.reserve(2*Config::expMutationCount);
    double probTerm = 1/Config::expMutationCount;
    do
      ins.push_back(Gene::getNew());
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

  Candidate mDeleteSlice() {
    auto &parent = get();
    auto &gtOrig = parent.genotype();
    auto sz = gtOrig.size();
    std::uniform_real_distribution<> dUni{0, 1};
    if(sz == 0)
      return parent;
    size_t pos1 = gen::rng() % (sz + 1);
    /* Integer with the same distribution in mAddSlice */
    size_t len = 1 + floor(log(dUni(gen::rng))
        / log(1 - 1/Config::expMutationCount));
    size_t pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz - (pos2 - pos1));
    gtNew.insert(gtNew.end(), gtOrig.begin(), gtOrig.begin() + pos1);
    gtNew.insert(gtNew.end(), gtOrig.begin() + pos2, gtOrig.end());
    return Candidate{std::move(gtNew)};
  }

  Candidate mDeleteUniform() {
    auto &parent = get();
    auto gtOrig = parent.genotype();
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gtNew{};
    gtNew.reserve(gtOrig.size());
    size_t cnt = 0;
    for(auto& g : gtOrig)
      if(dUni(gen::rng) >= Config::pChoiceUniform)
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
    std::array<size_t, 4> pos;
    for(auto& p : pos)
      p = gen::rng() % (sz - 1);
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
    size_t pos1 = gen::rng() % (sz - 1),
           pos2 = gen::rng() % (sz - 1);
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

  Candidate crossoverUniform() {
    auto &parent1 = get(),
         &parent2 = get();
    auto &gt1 = parent1.genotype(),
         &gt2 = parent2.genotype();
    size_t sz1 = gt1.size(),
           sz2 = gt2.size(),
           szShorter = std::min(sz1, sz2);
    double pTest1 = (double)szShorter / sz1,
           pTest2 = (double)szShorter / sz2;
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gtNew{};
    gtNew.reserve(sz1 + sz2 - szShorter);
    auto it1 = gt1.begin(),
         it2 = gt2.begin(),
         ie1 = gt1.end(),
         ie2 = gt2.end();
    while(it1 != ie1) {
      gtNew.push_back(*it1++);
      if(dUni(gen::rng) < pTest1) {
        // Whether to make this a possible crossover point.
        // The above check passes every time on the shorter genotype,
        // 1/ratio of the time on the longer one.
        // We might still decide to stay on the same branch, though.
        if(dUni(gen::rng) < Config::pCrossUniform) {
          std::swap(it1, it2);
          std::swap(ie1, ie2);
          std::swap(pTest1, pTest2);
        }
      }
    }
    // now it1 is at its boundary, if there's anything left past it2 we
    // flush it with no more crossovers
    gtNew.insert(gtNew.end(), it2, ie2);
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
    return simplify(get());
  }

public:

  static Candidate simplify(const Candidate& parent) {
    auto &gtOrig = parent.genotype();
    size_t sz = gtOrig.size();
    if(sz == 0)
      return parent;
    std::vector<Gene> gtNew = gtOrig;
    for(auto& g : gtNew)
      g.simplify();
    return gtNew != gtOrig ? Candidate{std::move(gtNew)} : parent;
  }

}; // class CandidateFactory


template<class Candidate, class Population>
class CFSelector {

  using CF = CandidateFactory<Candidate, Population>;
  using FunPtr = Candidate (CF::*)();

  struct GenOp {

    FunPtr fun;
    std::string name;
    double prob;
    unsigned long hits;
    unsigned long thits;

    GenOp(FunPtr fun_, std::string name_):
      fun(fun_), name(name_), prob(1), hits(0), thits(0) { }

  };

  std::vector<GenOp> ops;
  size_t count;

  std::discrete_distribution<size_t> dFun{};

public:

  void hit(size_t ix) {
    if(ix >= 0 && ix < count)
      ops[ix].hits++;
  }

  void update() {
    /* Calculate the probability distribution of GenOps based on prior success
     * rate */
    double denom = 0;
    for(auto& op : ops)
      denom += op.hits / op.prob;
    /* If we're not counting hits the probabilities will stay constant */
    if(denom != 0)
      for(auto& op : ops) {
        op.prob = (1 - Config::heurFactor) * op.prob
          + Config::heurFactor * op.hits / op.prob / denom;
        op.thits += op.hits;
        op.hits = 0;
      }
    std::vector<double> weights(count);
    for(size_t i = 0; i < count; i++)
      weights[i] = ops[i].prob;
    dFun = std::discrete_distribution<size_t>(weights.begin(), weights.end());
  }

  void dump(std::ostream& os) {
    /* Find the longest GenOp name */
    auto max = std::max_element(ops.begin(), ops.end(),
        [](const GenOp& a, const GenOp& b) {
          return a.name.length() < b.name.length();
        });
    auto maxw = max->name.length();
    /* Preserve settings of os */
    auto flags_ = os.flags(std::ios_base::left | std::ios_base::fixed);
    auto precision_ = os.precision(4);
    /* List all op names and probabilities */
    for(auto& op : ops)
      os << std::setw(maxw+3) << op.name + ':'
         << op.prob << "  " << op.thits << '\n';
    os.flags(flags_);
    os.precision(precision_);
  }

  std::pair<FunPtr, size_t> select() {
    size_t index = dFun(gen::rng);
    return {ops[index].fun, index};
  }

  CFSelector() {
    ops.push_back({ &CF::mAlterDiscrete,   "MDiscrete" });
    ops.push_back({ &CF::mAlterContinuous, "MContns" });
    ops.push_back({ &CF::mAddSlice,        "AddSlice" });
    ops.push_back({ &CF::mAddPairs,        "AddPairs" });
    ops.push_back({ &CF::mDeleteSlice,     "DelShort" });
    ops.push_back({ &CF::mDeleteUniform,   "DelUnif"  });
    ops.push_back({ &CF::mSplitSwap,       "SpltSwp"  });
    ops.push_back({ &CF::mReverseSlice,    "InvSlice" });
    ops.push_back({ &CF::crossoverUniform, "C/Over"   });
  //ops.push_back({ &CF::concat3,          "Concat3"  });
    ops.push_back({ &CF::simplify,         "Simplify" });
    count = ops.size();
    double pUniform = 1.0 / count;
    for(auto& op : ops)
      op.prob = pUniform;
    update();
  }

}; // class CFSelector<CandidateFactory>

} // namespace QGA
