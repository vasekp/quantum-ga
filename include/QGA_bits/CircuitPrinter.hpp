namespace QGA {

namespace internal {

class CircuitPrinter {

public:

  CircuitPrinter(unsigned nBit_): nBit(nBit_), lines(2*(nBit - 1)) {
    // lines[0], lines[2], ..., lines[2*(nBit - 1)]: qubits
    // lines[1], lines[3], ..., lines[2*nBit - 3]: padding space
    for(unsigned i = 0; i < nBit; i++)
      lines[2*i] = "-";
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
    std::vector<std::string> names{};
    for(auto& pair : pairs)
      for(unsigned qubit : pair.first) {
        qubits.push_back(qubit);
        names.push_back(pair.second);
      }
    addGates(qubits, names);
  }

  /* Prints a gate spanning all lines of output. */
  void addBarrierGate(std::string name) {
    alignAll();
    for(auto& line : lines)
      line += '[';
    lines[nBit - 1] += name;
    alignAll(' ');
    for(auto& line : lines)
      line += ']';
    for(unsigned i = 0; i < nBit; i++)
      lines[2*i] += '-';
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

  void addGates(std::vector<unsigned> qubits, std::vector<std::string> names) {
    /* Find the longest gate name to be printed */
    unsigned maxLen = 0;
    for(auto& name : names)
      if(length_utf8(name) > maxLen)
        maxLen = length_utf8(name);

    /* Find the total span of the gates being printed */
    auto minmax = std::minmax_element(qubits.begin(), qubits.end());
    align(*minmax.first, *minmax.second);

    /* Print the gate names or connecting | on even (qubit) lines */
    for(unsigned i = *minmax.first; i <= *minmax.second; i++) {
      auto it = std::find(qubits.begin(), qubits.end(), i);
      unsigned pos = it - qubits.begin();
      if(it != qubits.end())
        lines[2*i] = lines[2*i] + pad(names[pos], maxLen, '-') + "-";
      else
        lines[2*i] = lines[2*i] + pad("|", maxLen, '-') + "-";
    }

    /* Print the connecting | on odd (padding) lines */
    for(unsigned i = 2*(*minmax.first) + 1; i < 2*(*minmax.second); i += 2)
      lines[i] = lines[i] + pad("|", maxLen, ' ') + ' ';
  }

  void align(unsigned lineFrom, unsigned lineTo) {
    /* Find the longest line wihnin the span so far */
    unsigned max = 0;
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++) {
      if(length_utf8(lines[i]) > max)
        max = length_utf8(lines[i]);
    }

    /* Pad the lines with a dash (even: qubit) or a space (odd: padding) */
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++)
      if(length_utf8(lines[i]) < max)
        lines[i] += std::string(max - length_utf8(lines[i]), i&1 ? ' ' : '-');
  }

  void alignAll() {
    align(0, nBit - 1);
  }

  /* Symetrically pad a string to a given length using a specified character */
  std::string pad(std::string s, unsigned len, char padding) {
    unsigned padLeft = (len - length_utf8(s)) / 2,
             padRight = len - length_utf8(s) - padLeft;
    return std::string(padLeft, padding) + s + std::string(padRight, padding);
  }

  /* This allows names like Φ for quantum gates without distorting
   * the output. */
  size_t length_utf8(std::string s) {
    size_t len = 0;
    for(char c : s)
      if((c & 0xC0) != 0x80)
        len++;
    return len;
  }

  unsigned nBit;
  std::vector<std::string> lines;

}; // class CircuitPrinter

} // namespace internal

} // namespace QGA
