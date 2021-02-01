# KPartiteKClique
Iterate over all k-cliques of a k-partite graph

# Building

This should just build fine using `g++` and probably also with other compilers.

About 25 percent performance gain can achieved using build in popcount.

To use it, you can compile with `-mpopcnt` or `-march=native`.
The first will lead to illegal instructions, if not available.
The second is safe, if the compiler supports it.
In particular the GCC suite supports it.

# Algorithm

The algorithm is described [here](ALGORITHM.md)

# See also

The pivot selection of this algorithm was first described in

Konc,J. and Janezic,D. An improved branch and bound algorithm for the maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

https://gitlab.com/janezkonc/mcqd
