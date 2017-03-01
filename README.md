# Introduction

The purpose of this project is to provide a platform for automated discovery of [quantum circuits](https://en.wikipedia.org/wiki/Quantum_circuit) by means of a [genetic evolutionary algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm). Defining a simulated population to consist of initially randomly designed quantum circuits and an environment making those individual circuits providing a result closer to a designated goal more likely to survive and procreate, the selection pressure along with the random alterations (mutations, cross-overs) gradually eliminates unsuitable circuits and retains, improves and combines elements of those circuits which perform better. Ideally, this randomly driven process produces (over the course of many generations) solutions very close to perfect.

Currently two real-world applications are demonstrated in this repository:

1. Quantum search circuit discovery, aiming to evolve towards [Grover's algorithm](https://arxiv.org/abs/quant-ph/9605043),

2. Quantum factoring circuit discovery, aiming to produce the core of [Shor's algorithm](https://arxiv.org/abs/quant-ph/9508027).

However, the purpose of this document is to explain how these two goals (with known results) only illustrate the functionality of a more general platform, and guide through the process of implementing evolutions with other, custom goals.

## Read more about:

* [Installation](https://github.com/vasekp/quantum-ga/blob/master/manual/Installation.md)
* [Running](https://github.com/vasekp/quantum-ga/blob/master/manual/Running.md)
* [Extending](https://github.com/vasekp/quantum-ga/blob/master/manual/Extending.md)
* [Customizing](https://github.com/vasekp/quantum-ga/blob/master/manual/Customizing.md)
