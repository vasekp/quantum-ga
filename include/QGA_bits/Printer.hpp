class Printer {

public:

  Printer(unsigned nBit_): nBit(nBit_), lines(2*nBit - 1) {
    for(unsigned i = 0; i < nBit; i++)
      lines[2*i] = "-";
  }

  void addGate(unsigned line, std::string name) {
    addGates({line}, {name});
  }

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

  friend std::ostream& operator<< (std::ostream& os, Printer& p) {
    p.alignAll();
    for(auto& line : p.lines)
      os << line << '\n';
    return os;
  }

private:

  void addGates(std::vector<unsigned> qubits, std::vector<std::string> names) {
    unsigned maxLen = 0;
    for(auto& name : names)
      if(length_utf8(name) > maxLen)
        maxLen = length_utf8(name);
    auto minmax = std::minmax_element(qubits.begin(), qubits.end());
    align(*minmax.first, *minmax.second);
    for(unsigned i = *minmax.first; i <= *minmax.second; i++) {
      auto it = std::find(qubits.begin(), qubits.end(), i);
      unsigned pos = it - qubits.begin();
      if(it != qubits.end())
        lines[2*i] = lines[2*i] + pad(names[pos], maxLen, '-') + "-";
      else
        lines[2*i] = lines[2*i] + pad("|", maxLen, '-') + "-";
    }
    for(unsigned i = 2*(*minmax.first) + 1; i < 2*(*minmax.second); i += 2)
      lines[i] = lines[i] + pad("|", maxLen, ' ') + ' ';
  }

  void align(unsigned lineFrom, unsigned lineTo, char c = '-') {
    unsigned max = 0;
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++) {
      if(length_utf8(lines[i]) > max)
        max = length_utf8(lines[i]);
    }
    for(unsigned i = 2*lineFrom; i <= 2*lineTo; i++)
      if(length_utf8(lines[i]) < max)
        lines[i] += std::string(max - length_utf8(lines[i]), i&1 ? ' ' : c);
  }

  void alignAll(char c = '-') {
    align(0, nBit - 1, c);
  }

  std::string pad(std::string s, unsigned len, char padding) {
    unsigned padLeft = (len - length_utf8(s)) / 2,
             padRight = len - length_utf8(s) - padLeft;
    return std::string(padLeft, padding) + s + std::string(padRight, padding);
  }

  size_t length_utf8(std::string s) {
    size_t len = 0;
    for(char c : s)
      if((c & 0xC0) != 0x80)
        len++;
    return len;
  }

  unsigned nBit;
  std::vector<std::string> lines;

};
