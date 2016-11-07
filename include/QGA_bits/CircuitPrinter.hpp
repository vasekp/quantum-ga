namespace QGA {

namespace internal {

class CircuitPrinter {

  /* This allows names like Φ for quantum gates without distorting
   * the output. */
  class UTF8String {

  public:

    UTF8String(): s(), len(0) { }

    UTF8String(const std::string& s_): s(s_), len(utf8_length(s)) { }

    UTF8String(char c, size_t len_ = 1): s(len_, c), len(len_) { }

    size_t length() const { return len; }

    operator std::string() const { return s; }

    void operator+= (const UTF8String& other) {
      s += other.s;
      len += other.len;
    }

    friend std::ostream& operator<< (std::ostream& os, const UTF8String& utf) {
      return os << utf.s;
    }

    /* Symetrically pad a string to a given length using a specified character */
    UTF8String pad(size_t target, char padding) const& {
      size_t padLeft = (target - len) / 2,
             padRight = target - len - padLeft;
      UTF8String ret{padding, padLeft};
      ret += *this;
      ret += {padding, padRight};
      return ret;
    }

  private:

    size_t utf8_length(const std::string& s) {
      size_t ret = 0;
      for(char c : s)
        if((c & 0xC0) != 0x80)
          ret++;
      return ret;
    }

    std::string s;
    size_t len;

  }; // inner class UTF8String

public:

  CircuitPrinter(unsigned nBit_): nBit(nBit_), lines(2*nBit - 1) {
    // lines[0], lines[2], ..., lines[2*(nBit - 1)]: qubits
    // lines[1], lines[3], ..., lines[2*nBit - 3]: padding space
    for(unsigned i = 0; i < nBit; i++)
      lines[2*i] = '-';
  }

  /* Prints a single gate. */
  void addGate(unsigned line, std::string name) {
    addGates({line}, {name});
  }

  /* Every item in the vector is a pair of:
   *  1. a list of indices,
   *  2. a string to put at those indices.
   * All the gates will be connected by a vertical line. */
  void addGates(std::vector<
                  std::pair<std::vector<unsigned>, std::string>
                > pairs) {
    std::vector<unsigned> qubits{};
    std::vector<UTF8String> names{};
    for(auto& pair : pairs)
      for(unsigned qubit : pair.first) {
        qubits.push_back(qubit);
        names.push_back(pair.second);
      }
    addGates(qubits, names);
  }

  /* Prints a gate spanning all lines of output. */
  void addBarrierGate(std::string name) {
    addBroadGate(0, nBit - 1, {name});
  }

  template<class Gate>
  void print(const Gate& g) {
    g->printOn(*this);
  }

  friend std::ostream& operator<< (std::ostream& os, CircuitPrinter& p) {
    p.alignAll();
    for(auto& line : p.lines)
      os << line << '\n';
    return os;
  }

  /* Normally operator<< would accept a const CircuitPrinter&, but we're
   * modifying it using alignAll(). The other common situation is printing a
   * temporary variable, which would not be accepted by the above function.
   * This trivial wrapper gives the temporary a name (p) and thus converts the
   * rvalue reference to an lvalue reference. */
  friend std::ostream& operator<< (std::ostream& os, CircuitPrinter&& p) {
    return os << p;
  }

private:

  void addGates(std::vector<unsigned> qubits, std::vector<UTF8String> names) {
    /* Find the longest gate name to be printed */
    unsigned maxLen = 0;
    for(auto& name : names)
      if(name.length() > maxLen)
        maxLen = name.length();

    /* Find the total span of the gates being printed */
    auto minmax = std::minmax_element(qubits.begin(), qubits.end());
    align(*minmax.first, *minmax.second);

    /* Print the gate names or connecting | on even (qubit) lines */
    for(unsigned i = *minmax.first; i <= *minmax.second; i++) {
      auto it = std::find(qubits.begin(), qubits.end(), i);
      unsigned pos = it - qubits.begin();
      lines[2*i] += (it != qubits.end()
          ? names[pos]      // found
          : UTF8String('|') // not found
        ).pad(maxLen, '-');
      lines[2*i] += '-';
    }

    /* Print the connecting | on odd (padding) lines */
    for(unsigned i = 2*(*minmax.first) + 1; i < 2*(*minmax.second); i += 2) {
      lines[i] += UTF8String('|').pad(maxLen, ' ');
      lines[i] += ' ';
    }
  }

  void addBroadGate(unsigned lineFrom, unsigned lineTo, UTF8String name) {
    alignAll();
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++)
      lines[i] += '[';
    lines[lineFrom + lineTo] += name; // arithmetic average of 2*lineX
    alignAll(' ');
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++)
      lines[i] += ']';
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i += 2)
      lines[i] += '-';
  }

  void align(unsigned lineFrom, unsigned lineTo, char c = '-') {
    /* Find the longest line wihnin the span so far */
    unsigned max = 0;
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++) {
      if(lines[i].length() > max)
        max = lines[i].length();
    }

    /* Pad the lines with a dash (even: qubit) or a space (odd: padding) */
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++)
      if(lines[i].length() < max)
        lines[i] += {i&1 ? ' ' : c, max - lines[i].length()};
  }

  void alignAll(char c = '-') {
    align(0, nBit - 1, c);
  }

  unsigned nBit;
  std::vector<UTF8String> lines;

}; // class CircuitPrinter

} // namespace internal

} // namespace QGA