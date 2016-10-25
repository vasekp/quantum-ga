namespace QGA {

namespace internal {

// Defined below
template<class, class, class...>
class Chooser;

template<class, class...>
class Reader;


/* The Gene type ready for use with a CandidateBase.
 *
 * This holds a shared pointer to GateBase (as typedef'd in GateBase.hpp) and
 * as such can refer to different gates polymorphically.
 *
 * Implements a static getRandom() function which randomly picks from the given
 * Gates, along with several shortcuts to functions of the GateBase.  Other
 * functions like applyTo() or hit() are redirected using a ->. */

template<class Context, class... Gates>
class Gene : GateBase<Context, Gates...>::Pointer {

  using GBase = GateBase<Context, Gates...>;
  using Pointer = typename GBase::Pointer;

public:

  using Counter = typename GBase::Counter;

  template<class Context_>
  using WithContext = Gene<Context_, Gates...>;

  Gene() = default; // Needed in CandidateBase::read()

  Gene(const Pointer& ptr): Pointer(ptr) { }

  Gene(Pointer&& ptr): Pointer(std::move(ptr)) { }

  static Gene getRandom() {
    return {internal::Chooser<
      Pointer, typename Gates::template Template<GBase>...
    >::getNew()};
  }

  const GBase* operator->() const {
    return Pointer::operator->();
  }

  friend std::ostream& operator<<(std::ostream& os, const Gene& g) {
    return os << *g;
  }

  friend std::istream& operator>>(std::istream& is, Gene& g) {
    std::string gene{};
    if(!(is >> gene))
      return is;
    Pointer ptr = internal::Reader<
      typename Gates::template Template<GBase>...
    >::read(gene);
    if(!ptr) {
      is.setstate(std::ios::failbit);
      return is;
    }
    g = {ptr};
    return is;
  }

  /* These helper functions simplify the call pattern of the pointer-passing
   * virtual functions of GateBase. */
  void invert() {
    pointer() = pointer()->invert(pointer());
  }

  void mutate() {
    pointer() = pointer()->mutate(pointer());
  }

  void simplify() {
    pointer() = pointer()->simplify(pointer());
  }

  bool merge(Gene& other) {
    Pointer result = pointer()->merge(pointer(), other);
    if(result) {
      pointer() = result;
      return true;
    } else
      return false;
  }

  /* Two genes are equal iff they point to the same object. This is used in
   * CandidateFactory to check if a genotype has changed at all during a
   * mutation (if not, the parent is returned). */
  bool operator==(const Gene& other) const {
    return pointer() == other.pointer();
  }

private:

  Pointer& pointer() {
    return static_cast<Pointer&>(*this);
  }

  const Pointer& pointer() const {
    return static_cast<const Pointer&>(*this);
  }

}; // class Gene<Context, Gates...>

} // namespace internal


/* The default value for Context (additional data passed to GateBase::applyTo)
 * is void. Other classes can be supplied later using the alias declaration
 * internal::Gene::WithContext. */

template<class... Gates>
using Gene = internal::Gene<void, Gates...>;


namespace internal {

/* Called as Chooser<P, A, B, C, ...>::getNew(i), returns:
 * i = 0: std::make_shared<A>()
 * i = 0: std::make_shared<B>()
 * i = 0: std::make_shared<C>()
 * etc. The classes A, B, C, ... need to be default-constructible and
 * convertible to the class referred to by the pointer P. */

template<class Pointer, class Head, class... Tail>
class Chooser {

  template<class, class, class...>
  friend class Chooser;

  static Pointer getNew(unsigned index) {
    if(index == 0)
      return std::make_shared<Head>();
    else
      return Chooser<Pointer, Tail...>::getNew(index - 1);
  }

public:

  static Pointer getNew() {
    // upper bound inclusive: no need to add 1 for Head
    std::uniform_int_distribution<> dist(0, sizeof...(Tail));
    return getNew(dist(gen::rng));
  }

}; // class Chooser<Pointer, Head, Tail...>

template<class Pointer, class Last>
class Chooser<Pointer, Last> {

public:

#ifdef DEBUG
  static Pointer getNew(unsigned index = 0) {
    if(index == 0)
      return std::make_shared<Last>();
    else
      throw std::logic_error("Index too large in Chooser!");
  }
#else
  static Pointer getNew(unsigned = 0) {
    return std::make_shared<Last>();
  }
#endif

}; // class Chooser<Pointer, Last>


/* Called as Reader<A, B, C, ...>::read(input), this helper tries calling
 * A::read(input), B::read(input), C::read(input)... until the first returns a
 * value convertible to boolean true. */

template<class Head, class... Tail>
class Reader {

public:

  static auto read(std::string input) -> decltype(Head::read(input)) {
    if(auto ret = Head::read(input))
      return ret;
    else
      return Reader<Tail...>::read(input);
  }

}; // class Reader<Head, Tail...>

template<class Last>
class Reader<Last> {

public:

  static auto read(std::string input) -> decltype(Last::read(input)) {
    return Last::read(input);
  }

}; // class Reader<Last>

} // namespace internal

} // namespace QGA
