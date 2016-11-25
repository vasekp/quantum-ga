#include "include/regex.hpp"
#include <regex>

/* See include/regex.hpp */

namespace regex {

  class regex::regex_impl : public std::regex {
  public:
    using std::regex::regex;
  };

  class matches::matches_impl : public std::smatch { };


  regex::regex(std::string expr):
    pImpl(std::make_unique<regex_impl>(expr))
  { }

  regex::~regex() { }

  bool regex::match(std::string searched, matches& ms) {
    return std::regex_match(searched, ms.impl(), *pImpl);
  }

  matches::matches():
    pImpl(std::make_unique<matches_impl>())
  { }

  matches::~matches() { }

  bool matches::matched(size_t index) {
    return impl()[index].matched;
  }

  std::string matches::match(size_t index) {
    return impl()[index].str();
  }

  matches::matches_impl& matches::impl() {
    return *pImpl;
  }

} // namespace regex
