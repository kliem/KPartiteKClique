# The used classes

KPartiteKClique
- public:
  - `void KPartiteKClique(bool** incidences, int n_vertices, int* parts, int k)`
    Initializes vertices.
    Initializes k-graphs (the first full, the others empty).
  - `int next()`
  - destructor
  - `int* k_clique`
- private:
  - `int* parts`
  - `int k`
  - `current_depth`
  - `int n_vertices`
  - `Vertex* vertices`
  - `KPartiteGraph* recursive graphs`
  - `int traceback`
    traces back until it finds a valid graph,
    if not possible return 0

KPartiteGraph
- public:
  - `void KPartiteGraph(int n_vertices, bool fill, int* parts, int k)`
  - `operator assign`
    `copy active_vertices`
    `copy vertices`
    `check for trivial cases and sort vertices`
    `resizes the vector`
    `removes old selection`
    `checks for validity`
  - `int next()`
- private:
  - `int* parts`
  - int k
  - `int selected_k`
  - `int* k_clique`
  - `Bitset active_vertices`
  - `vector<Vertex> vertices`

Bitset
- public:
  - `void Bitset(int n_vertices)`
  - `void Bitset(*bool set_bits, int n_vertices)`
  - `intersection`
  - `intersection_assign`
  - `count(int start, int stop)`
  - `intersection_count(int start, int stop, Bitset other)`


# This is how it works:

## From the users point of view

- Initialize `KPartiteKClique`
- Call `next` until exhausted.
  Returns `1` on sucess and `0` otherwise.
- `k_clique` points to the indices of the `k_clique`.
