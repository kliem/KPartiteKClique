# General Approach

## FindClique

To iterate over all k-cliques of a k-partite graph, there is a very
natural approach that is also described in

- Grünert, Tore & Irnich, Stefan & Zimmermann, Hans-Jürgen & Schneider, Markus & Wulfhorst, Burkhard. (2001). Cliques in k-partite Graphs and their Application in Textile Engineering

The following follow up proposes to use bitsets as a data structure and
improves upon the implementation:

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

We have implemented this algorithm `FindClique` as free software, as we could not find an available implementation.

While the later paper is an important improvement, we use the name
`FindClique` as proposed in the earlier paper. It is the same algorithm
and bitsets are a common data structure for graphs.


## KPartiteKClique

`FindClique` performs especially well on graphs with
somewhat equally distributed edges.

However, we have encountered graphs that appear infeasible with this approach.
Those graphs arise from point configurations.
They have far from equally distributed edges.

We use the pivot selection strategy described in

- Konc,J. and Janežič,D. An improved branch and bound algorithm for the maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.

*Select a vertex with minimal degree*

This approach might seem counterintuitive, as we start with the vertex
that is least likely to be contained in a k-clique.
However, this allows to reduce the recursion depth significantly.

`mcqd` completely ignores the given k-partition
and just computes the clique number by coloring the graph (but likely
not optimal).

We specialize their pivot strategy for a given k-partition and for
iterating over the k-cliques instead of finding the size of a maximum
clique.

# Implementation Details

We also represent incidences by bitsets as done by Mirghorbani et. al.
