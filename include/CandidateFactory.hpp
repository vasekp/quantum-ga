namespace QGA {


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

  static NOINLINE Candidate genInit() {
    // probability of termination; expLengthIni = expected number of genes
    const double probTerm = 1/Config::expLengthIni;
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do
      gt.push_back(Gene::getNew());
    while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gt)};
  }

  NOINLINE Candidate getNew() {
    auto op = sel.select();
    return (this->*op.first)().setOrigin(op.second);
  }

private:

  const Candidate& get() {
    return pop.NSGASelect(Config::selectBias);
  }

  Candidate mAlterDiscrete() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    std::uniform_real_distribution<> dUni{0, 1};
    size_t cnt = 0;
    for(auto& g : gm)
      if(dUni(gen::rng) < Config::pChoiceUniform) {
        g = Gene::getNew();
        cnt++;
      }
    return cnt ? Candidate{std::move(gm)} : p;
  }

  Candidate mAlterContinuous() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    std::uniform_real_distribution<> dUni{0, 1};
    size_t cnt = 0;
    for(auto& g : gm)
      if(dUni(gen::rng) < Config::pChoiceUniform)
        if(g.mutate())
          cnt++;
    return cnt ? Candidate{std::move(gm)} : p;
  }

  Candidate mAddSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    std::uniform_real_distribution<> dUni{0, 1};
    unsigned pos = gen::rng() % (sz + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do
      ins.push_back(Gene::getNew());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(sz + ins.size());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), gt.begin() + pos, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mAddPairs() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    std::uniform_real_distribution<> dUni{0, 1};
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins;
    ins.reserve(2*Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do
      ins.push_back(Gene::getNew());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(sz + 2*ins.size());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), gt.begin() + pos1, gt.begin() + pos2);
    for(auto& g : ins)
      g.invert();
    gm.insert(gm.end(),
        std::make_move_iterator(ins.rbegin()),
        std::make_move_iterator(ins.rend()));
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    std::uniform_real_distribution<> dUni{0, 1};
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1);
    /* Integer with the same distribution in mAddSlice */
    int len = 1 + floor(log(dUni(gen::rng)) / log(1 - 1/Config::expLengthAdd));
    unsigned pos2 = pos1 + len > sz ? sz : pos1 + len;
    std::vector<Gene> gm;
    gm.reserve(sz - (pos2 - pos1));
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mDeleteUniform() {
    auto &p = get();
    auto gt = p.genotype();
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gm;
    gm.reserve(gt.size());
    size_t cnt = 0;
    for(auto& g : gt)
      if(dUni(gen::rng) >= Config::pChoiceUniform)
        gm.push_back(g);
      else
        cnt++;
    return cnt ? Candidate{std::move(gm)} : p;
  }

  Candidate mSplitSwap() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    std::array<unsigned, 4> pos;
    for(int i = 0; i < 4; i++)
      pos[i] = gen::rng() % (sz + 1);
    std::sort(pos.begin(), pos.end());
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos[0]);
    gm.insert(gm.end(), gt.begin() + pos[2], gt.begin() + pos[3]);
    gm.insert(gm.end(), gt.begin() + pos[1], gt.begin() + pos[2]);
    gm.insert(gm.end(), gt.begin() + pos[0], gt.begin() + pos[1]);
    gm.insert(gm.end(), gt.begin() + pos[3], gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate mReverseSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> gm;
    gm.reserve(sz);
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    {
      auto prev_end = gm.end();
      gm.insert(gm.end(), gt.rbegin() + sz - pos2, gt.rbegin() + sz - pos1);
      auto end = gm.end();
      for(auto it = prev_end; it != end; it++)
        it->invert();
    }
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossoverUniform() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype();
    size_t sz1 = gt1.size(),
           sz2 = gt2.size(),
           szS = std::min(sz1, sz2);
    double pTest1 = (double)szS / sz1,
           pTest2 = (double)szS / sz2;
    std::uniform_real_distribution<> dUni{0, 1};
    std::vector<Gene> gm;
    gm.reserve(sz1 + sz2 - szS);
    auto it1 = gt1.begin(),
         it2 = gt2.begin(),
         ie1 = gt1.end(),
         ie2 = gt2.end();
    while(it1 != ie1) {
      gm.push_back(*it1++);
      if(dUni(gen::rng) < pTest1) {
        // happens every time on the shorter genotype,
        // 1/ratio of the time on the longer one
        // We might still decide to stay, though.
        if(dUni(gen::rng) < Config::pCrossUniform) {
          std::swap(it1, it2);
          std::swap(ie1, ie2);
          std::swap(pTest1, pTest2);
        }
      }
    }
    // now it1 is at its boundary, if there's anything left past it2 we
    // flush it with no more crossovers
    gm.insert(gm.end(), it2, ie2);
    return Candidate{std::move(gm)};
  }

  Candidate concat3() {
    auto &p1 = get(),
         &p2 = get(),
         &p3 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype(),
         &gt3 = p3.genotype();
    std::vector<Gene> gm;
    gm.reserve(gt1.size() + gt2.size() + gt3.size());
    gm.insert(gm.end(), gt1.begin(), gt1.end());
    {
      auto start = gm.end();
      gm.insert(gm.end(), gt2.rbegin(), gt2.rend());
      auto end = gm.end();
      for(auto it = start; it != end; it++)
        it->invert();
    }
    gm.insert(gm.end(), gt3.begin(), gt3.end());
    return Candidate{std::move(gm)};
  }

  Candidate simplify() {
    auto &p = get();
    auto &gt = p.genotype();
    size_t sz = gt.size();
    if(!sz)
      return p;
    std::vector<Gene> gm;
    gm.reserve(sz);
    auto it = gt.begin(), end = gt.end();
    gm.push_back(*it);
    auto last = gm.begin();
    size_t cnt = 0;
    for(it++; it != end; it++)
      if(!last->merge(*it)) {
        gm.push_back(*it);
        last++;
      } else {
        // merge == true: success, skip this it
        cnt++;
      }
    for(auto& g : gm)
      if(g.simplify())
        cnt++;
    return cnt ? Candidate{std::move(gm)} : p;
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
    unsigned hits;

    GenOp(FunPtr fun_, std::string name_):
      fun(fun_), name(name_), prob(1), hits(0) { }

  };

  std::vector<GenOp> ops;
  size_t count;

  std::discrete_distribution<> dFun{};

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
      denom += std::sqrt(op.hits);
    /* If we're not counting hits the probabilities will stay constant */
    if(denom != 0)
      for(auto& op : ops)
        op.prob = (1 - Config::heurFactor) * op.prob
          + Config::heurFactor * std::sqrt(op.hits) / denom;
    std::vector<double> weights(count);
    for(size_t i = 0; i < count; i++)
      weights[i] = ops[i].prob;
    dFun = std::discrete_distribution<>(weights.begin(), weights.end());
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
      os << std::setw(maxw+3) << op.name + ':' << op.prob << '\n';
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
    ops.push_back({ &CF::concat3,          "Concat3"  });
    ops.push_back({ &CF::simplify,         "Simplify" });
    count = ops.size();
    double pUniform = 1.0 / count;
    for(auto& op : ops)
      op.prob = pUniform;
    update();
  }

}; // class CFSelector<CandidateFactory>

} // namespace QGA
