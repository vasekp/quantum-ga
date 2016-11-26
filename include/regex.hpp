#ifndef REGEX_HPP
#define REGEX_HPP

#include <string>
#include <memory>
#include <stddef.h>

/* The purpose of these classes is to allow compiling the std::regex template
 * sources separately and then use the functionality via an opaque library, as
 * recompiling them every time with -O3 is very time-consuming. */

namespace regex {

  class matches {

    class matches_impl;
    std::unique_ptr<matches_impl> pImpl;

  public:

    matches();
    ~matches();
    bool matched(size_t index);
    std::string match(size_t index);

  private:

    matches_impl& impl();

    friend class regex;

  }; // class matches


  class regex {

    class regex_impl;
    std::unique_ptr<regex_impl> pImpl;

  public:

    regex(std::string expr);
    ~regex();
    bool match(std::string searched, matches& ms);

  }; // class regex

} // namespace regex

#endif // !defined REGEX_HPP
