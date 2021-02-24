# KPartiteKClique
Iterate over all k-cliques of a k-partite graph

# Building

This should just build fine using `g++` and probably also with other compilers.

About 25 percent performance gain can achieved using built-in popcount.
To use it, you can compile with `-mpopcnt -mbmi` or `-march=native`.
The first will lead to illegal instructions, if not available.
The second is safe, if the compiler supports it.
In particular the GCC suite supports it.

# Usage

## bitCLQ

We provide an implementation of `bitCLQ` that was first described in

- Grunert, Tore & Irnich, Stefan & Zimmermann, Hans-Jürgen & Schneider, Markus & Wulfhorst, Burkhard. (2001). Cliques in k-partite Graphs and their Application in Textile Engineering

and the follow up

- Mirghorbani, M. & Krokhmal, P.. (2013). On finding k-cliques in k-partite graphs. Optimization Letters. 7. 10.1007/s11590-012-0536-y

## KPartiteKClique -- kpkc

We provide a different pivot strategy with the class `KPartiteKClique`.

The pivot selection of this algorithm was first described in

Konc,J. and Janezic,D. An improved branch and bound algorithm for the maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

https://gitlab.com/janezkonc/mcqd

Our implementation exploits that we already have a k-partite graph and
are looking for k-cliques only.

# Algorithm

The algorithm `KPartiteKClique` is described [here](ALGORITHM.md)

# Python wrapper

A python wrapper is avaible [here](https://github.com/kliem/PyKPartiteKClique).
