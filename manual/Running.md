# Running

After successful installation, run `./search` or `./fourier` without further parameters. The genetic evolution will commence immediately, starting with a random population. The evaluation of one generation usually takes several hundredths of seconds so the output flows quickly and only summative information is displayed for an overview of the progress and is colour-coded for clarity. On a more concrete example:

![aa](http://i.imgur.com/PV3Bj5q.png)

* **A** denotes the current generation number,
* **B** shows the population size (after culling redundant copies of equal candidates),
* **C** shows the properties of the best-so-far candidate (nondominated and with minimal error),
* **D** shows the size of the nondominated front (the internal population, or archive),
* **E** displays the newest addition to the front,
* **F** is a text-based visualisation of the candidate summarized in **C**.

Due to the nature of the problem and to the nature of the search, there is no quick and guaranteed way of finding the perfect circuit for a given task. A great part of the configuration space needs to be explored before exploiting the discovered features to approach an optimal solution. For this reason a strategy of taking many small steps in a large number of generations has been chosen over a small number of generations employing very elaborate genetic operations. Thus the search may run for a few thousand generations if a perfect solution is required (typically between 1000 and 2000 for the two benchmark problems on 3 qubits, taking 1 to 2 minutes of run time on a modern 4-core processor).

There is no hard-coded termination condition. Surpassing a given error bound and checking for stalled evolution have been considered and rejected. Instead, the user is given a liberty of interrupting and examining the evolution at any point, and to a very limited extent, direct intervention into the evolution is allowed as well. Also see below for termination of the program.

## Interruptions

To pause or stop the evolution and examine the results, hit `Ctrl+C`. The following menu appears:

```
Computation stopped. Choose action:
a: abort,
c: continue,
d: diagnose / list current results,
e: evaluate a candidate in full,
f: filter the front on fitness,
i: inject a candidate,
l: list 20 random candidates,
p: pretty-print a candidate as a circuit,
r: restart,
t: format a candidate as a LuaLaTeX Q-circuit,
q: quit after this generation.
```

User input is then expected in the form of a single lower-case letter followed by `Enter`. If, at this point, you want to exit the program, use **q** (natural stop) or **a** (forced exit) or `Ctrl+C` for a second time (ditto).

The main go-to option to examine the nondominated front is **d**. This lists all the nondominated candidate circuits sorted by their error from highest to lowest. (A high-error candidate can still be a member of the front when it's nondominated in other fitness aspects, e.g., number of gates. An empty circuit usually appears at the top of the list.) Some summative information follows, like the total number of evaluated candidates and time taken.

The list provided by **d** illustrates the format of textual encoding of candidate solutions that is accepted where a candidate is input explicitly by the user (**e**, **i**, **p**, or **t**). The details depend on the problem but the syntax follows a general pattern:

```
SWAP23 Y1(0.6156π) Y2(-0.0021π) P14(0.2028π) Y4(-0.3333π) SWAP14
```

Listing the gates from left to right, i.e., in the order they are applied on the initial state, first comes the name of each, followed by the 1-based qubit indices it acts on. This can be a single digit optionally followed by control qubit indices enclosed in brackets `[`...`]`, or, if all the affected qubits are treated equally, like in a control-Z gate or a swap gate, just a set of the qubit indices. (This assumes no more than 9 qubit lines will be needed, otherwise the current system needs a redesign.) Finally, if a gate has one or more continuous angle parameters, these appear in round parentheses as a multiple of π. Although this is always formatted with a fixed precision and the symbol for π, neither is required when parsing user input. However, no extra space is permitted.

A final choice that deserves some attention in this manual is **f**. Like **d**, this allows to examine the current front, but with the added benefit of filtering on an upper bound of some fitness aspects. For example, in the search problem, one may be interested only in those solutions which use 2 or less oracle calls and don't surpass an error of `0.5`. One could surely list the whole front and ignore solutions which don't qualify but the **f** choice simplifies this. Given that the latter value comes first and the former third within the fitness vector, a specification of such filter would be

```
0.5 _ 2
```

(Parts of fitness which come after the last one we're interested in can be left out.)
