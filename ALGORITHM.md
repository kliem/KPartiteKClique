# General Approach

To iterate over all k-cliques of a k-partite graph, there is a very
natural approach that is also described in

- Grunert, Tore & Irnich, Stefan & Zimmermann, Hans-Jürgen & Schneider, Markus & Wulfhorst, Burkhard. (2001). Cliques in k-partite Graphs and their Application in Textile Engineering

and the follow up

- Mirghorbani, M. & Krokhmal, P.. (2013). On finding k-cliques in k-partite graphs. Optimization Letters. 7. 10.1007/s11590-012-0536-y

This approach goes as follows:

1. Select a vertex.
2. Iterate over all (k-1)-cliques of the induced subgraph of the
   neighbors.
3. Remove the vertex and repeat.

The crucial ingredient to this algorithm is the pivot strategy when
selecting the vertex in 1.
Both references above agree in this strategy and the more recent
publication adds the use of bitsets:

*Select an arbitrary vertex in the part with minimal size*

However, this approach seems to perform only well on graphs with
somewhat equally distributed edges.
In this case the degrees of the vertices do not differ significantly.

Instead we use the pivot selection strategy described in

- Konc,J. and Janezic,D. An improved branch and bound algorithm for the maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

*Select a vertex with minimal degree*

This approach might seem counterintuitive, as we start with the vertex
that is least likely to be contained in a k-clique.
However, this allows to reduce the recursion depth significantly.
While selecting a vertex in the smallest part completly fails in some
instances (did not finish in 24 hours), the implementation mcqd of the last reference finishes in 12
seconds.

mcqd completely ignores the given k-partition
and just computes the clique number by coloring the graph (but likely
not optimal).

We specialize their pivot strategy for a given k-partition and for
iterating over the k-cliques instead of finding the size of a maximum
clique.

# Implementation Details

We also represent incidences by bitsets as done by Mirghorbani et. al.
This allows architecture dependent optimization.
