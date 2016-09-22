namespace QGA {

template<
  class Candidate,
  class Population = gen::NSGAPopulation<Candidate>>
class CandidateFactory {

  using Gene = typename Candidate::GeneType;
  using GeneFactory = typename Gene::Factory;

  std::uniform_real_distribution<> dUni{0, 1};

  struct GenOp {

    Candidate (CandidateFactory::*func)();
    std::string name;
    double prob = 1;
    unsigned hits = 0;

    GenOp(Candidate (CandidateFactory::*func_)(), std::string name_):
      func(func_), name(name_) { }

  };

public:

  using Counter = std::vector<GenOp>;

private:

  std::discrete_distribution<> dFun{};

  Population& pop;
  Counter& ctr;
  GeneFactory factory{};

public:

  CandidateFactory(Population& _pop, Counter& _ctr): pop(_pop), ctr(_ctr), factory{} {
    /* Calculate the probability distribution of GenOps based on prior success
     * rate */
    double denom = 0;
    for(auto& op : ctr)
      denom += std::sqrt(op.hits);
    /* If we're not counting hits the probabilities will stay constant */
    if(denom != 0)
      for(auto& op : ctr)
        op.prob = (1 - Config::heurFactor) * op.prob
          + Config::heurFactor * std::sqrt(op.hits / denom);
    std::vector<double> weights(ctr.size());
    for(size_t i = 0; i < ctr.size(); i++)
      weights[i] = ctr[i].prob;
    dFun = std::discrete_distribution<>(weights.begin(), weights.end());
  }

  static Counter getInitCounter();

  static NOINLINE Candidate genInit() {
    // probability of termination; expLengthIni = expected number of genes
    const static double probTerm = 1/Config::expLengthIni;
    static thread_local std::uniform_real_distribution<> dUni{0, 1};
    static GeneFactory staticFactory{};

    std::vector<Gene> gt;
    gt.reserve(Config::expLengthIni);
    do
      gt.push_back(staticFactory.getNew());
    while(dUni(gen::rng) > probTerm);
    return Candidate{std::move(gt)};
  }

  NOINLINE Candidate getNew() {
    int index = dFun(gen::rng);
    Candidate c = (this->*ctr[index].func)();
    c.setOrigin(index);
    return c;
  }

  void hit(int ix) {
    if(ix >= 0)
      ctr[ix].hits++;
  }

  void dumpWeights(std::ostream& os) {
    /* Find the longest GenOp name */
    using cmp = typename decltype(ctr)::value_type;
    auto max = std::max_element(ctr.begin(), ctr.end(),
        [](const cmp& v1, const cmp& v2) {
          return v1.name.length() < v2.name.length();
        });
    auto maxw = max->second.length();
    auto _flags = os.flags(std::ios_base::left | std::ios_base::fixed);
    auto _precision = os.precision(4);
    for(auto& op : ctr)
      os << std::setw(maxw+3) << op.name + ':' << op.prob << '\n';
    os.flags(_flags);
    os.precision(_precision);
  }

private:

  const Candidate& get() {
    return pop.NSGASelect(Config::selectBias);
  }

  Candidate mAlterSingle() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    if(!sz)
      return p;
    auto gm = gt;
    unsigned pos = gen::rng() % sz;
    gm[pos] = factory.getNew();
    return Candidate{std::move(gm)};
  }

  Candidate mAddSlice() {
    auto &p = get();
    auto &gt = p.genotype();
    auto sz = gt.size();
    unsigned pos = gen::rng() % (sz + 1);
    std::vector<Gene> ins;
    ins.reserve(Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do
      ins.push_back(factory.getNew());
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
    unsigned pos1 = gen::rng() % (sz + 1),
             pos2 = gen::rng() % (sz + 1);
    if(pos2 < pos1)
      std::swap(pos1, pos2);
    std::vector<Gene> ins;
    ins.reserve(2*Config::expLengthAdd);
    double probTerm = 1/Config::expLengthAdd;
    do
      ins.push_back(factory.getNew());
    while(dUni(gen::rng) > probTerm);
    std::vector<Gene> gm;
    gm.reserve(sz + 2*ins.size());
    gm.insert(gm.end(), gt.begin(), gt.begin() + pos1);
    gm.insert(gm.end(), ins.begin(), ins.end());
    gm.insert(gm.end(), gt.begin() + pos1, gt.begin() + pos2);
    for(auto& g : ins)
      // don't really move, just mark with a && for disposal
      factory.invert(std::move(g));
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
    std::vector<Gene> gm;
    gm.reserve(gt.size());
    for(auto& g : gt)
      if(dUni(gen::rng) >= Config::pDeleteUniform)
        gm.push_back(g);
    return Candidate{std::move(gm)};
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
      auto end = gt.begin() + pos2;
      for(auto it = gt.begin() + pos1; it != end; it++)
        factory.invert(std::move(*it));
    }
    gm.insert(gm.end(), gt.rbegin() + sz - pos2, gt.rbegin() + sz - pos1);
    gm.insert(gm.end(), gt.begin() + pos2, gt.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossover1() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype();
    unsigned pos1 = gen::rng() % (gt1.size() + 1),
             pos2 = gen::rng() % (gt2.size() + 1);
    std::vector<Gene> gm;
    gm.reserve(pos1 + (gt2.size() - pos2));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1);
    gm.insert(gm.end(), gt2.begin() + pos2, gt2.end());
    return Candidate{std::move(gm)};
  }

  Candidate crossover2() {
    auto &p1 = get(),
         &p2 = get();
    auto &gt1 = p1.genotype(),
         &gt2 = p2.genotype();
    unsigned pos1l = gen::rng() % (gt1.size() + 1),
             pos1r = gen::rng() % (gt1.size() + 1),
             pos2l = gen::rng() % (gt2.size() + 1),
             pos2r = gen::rng() % (gt2.size() + 1);
    if(pos1r < pos1l) std::swap(pos1l, pos1r);
    if(pos2r < pos2l) std::swap(pos2l, pos2r);
    std::vector<Gene> gm;
    gm.reserve(gt1.size() - (pos1r - pos1l) + (pos2r - pos2r));
    gm.insert(gm.end(), gt1.begin(), gt1.begin() + pos1l);
    gm.insert(gm.end(), gt2.begin() + pos2l, gt2.begin() + pos2r);
    gm.insert(gm.end(), gt1.begin() + pos1r, gt1.end());
    return Candidate{std::move(gm)};
  }

}; // class CandidateFactory


template<class Candidate, class Population>
typename CandidateFactory<Candidate, Population>::Counter
CandidateFactory<Candidate, Population>::getInitCounter() {
  typename CandidateFactory<Candidate, Population>::Counter ret{};
  ret.emplace_back(&CandidateFactory::mAlterSingle,    "MSingle" );
  ret.emplace_back(&CandidateFactory::mAddSlice,       "AddSlice");
  ret.emplace_back(&CandidateFactory::mAddPairs,       "AddPairs");
  ret.emplace_back(&CandidateFactory::mDeleteSlice,    "DelSlice");
  ret.emplace_back(&CandidateFactory::mDeleteUniform,  "DelUnif" );
  ret.emplace_back(&CandidateFactory::mSplitSwap,      "SpltSwp" );
  ret.emplace_back(&CandidateFactory::mReverseSlice,   "InvSlice");
  ret.emplace_back(&CandidateFactory::crossover1,      "C/Over1" );
  ret.emplace_back(&CandidateFactory::crossover2,      "C/Over2" );
  double pUniform = 1.0 / ret.size();
  for(auto& op : ret)
    op.prob = pUniform;
  return ret;
}

} // namespace QGA
