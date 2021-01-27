#if __POPCNT__
    #include <immintrin.h>
#endif

#include "k_partite_k_clique.h"

#include <algorithm>
#include <cstdlib>
#include <cassert>

const int ALIGNMENT = 8;

// Bitsets

void Bitset::intersection_assign(Bitset& l, Bitset& r){
    // Assumes all of same length.
    for (int i=0; i<limbs; i++)
        data[i] = l[i] & r[i];
}

inline int popcount(uint64_t i){
#if (__POPCNT__) && (INTPTR_MAX == INT64_MAX)
    return _mm_popcnt_u64(i);
#else
    i = i - ((i >> 1) & 0x5555555555555555ULL);
    i = (i & 0x3333333333333333ULL) + ((i >> 2) & 0x3333333333333333ULL);
    return ( ((i + (i >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
#endif
}

inline int Bitset::intersection_count(Bitset& r, int start, int stop){
    /*
    Count the number of set bits in ``this & r``
    in ``range(start, stop)``.
    */
    int counter = 0;
    // The easy part, count any complete ``uint64_t``.
    for (int i=start/64 + 1; i< stop/64; i++)
        counter += popcount(data[i] & r[i]);

    uint64_t start_limb = data[start/64] & r[start/64];
    if (start % 64)
        start_limb &= ~(((uint64_t) -1) >> (64 - (start % 64)));

    uint64_t end_limb;
    if (stop/64 < limbs){
        end_limb = data[stop/64] & r[stop/64];
        if (stop % 64)
            end_limb &= (((uint64_t) -1) >> (64 - (stop % 64)));
    }

    if (start/64 == stop/64){
        // The start limb is the end limb.
        counter += popcount(start_limb & end_limb);
    } else {
        if (stop/64 < limbs){
            counter += popcount(start_limb) + popcount(end_limb);
        } else {
            // There is no end limb.
            counter += popcount(start_limb);
        }
    }
    return counter;
}

bool Bitset::is_empty(int start, int stop){
    /*
    Currently not used.
    */
    for (int i=start/64 + 1; i< stop/64; i++){
        if (data[i]) return false;
    }

    uint64_t start_limb = data[start/64];
    if (start % 64)
        start_limb &= ~(((uint64_t) -1) >> (64 - (start % 64)));

    uint64_t end_limb;
    if (stop/64 < limbs){
        end_limb = data[stop/64];
        if (stop % 64)
            end_limb &= (((uint64_t) -1) >> (64 - (stop % 64)));
    }

    if (start/64 == stop/64)
        return 0 == start_limb & end_limb;

    if (stop/64 == 0)
        return 0 == start_limb;

    return 0 == start_limb | end_limb;
}

void Bitset::set(int index){
    data[index/64] |= ((uint64_t) 1) << (index % 64);
}

void Bitset::unset(int index){
    data[index/64] &= ~(((uint64_t) 1) << (index % 64));
}

bool Bitset::has(int index){
    return data[index/64] & (((uint64_t) 1) << (index % 64));
}

Bitset::Bitset(int n_vertices, bool fill){
    /*
    Initalize bitset.

    Fill if ``fill``.
    */
    shallow = false;
    limbs = (n_vertices-1)/64 + 1;
    allocate(n_vertices);
    if (!fill)
        return;
    for(int i=0; i<n_vertices/64; i++){
        data[i] = -1;
    }
    if (n_vertices % 64)
        data[n_vertices/64] = (((uint64_t) -1) >> (64 - (n_vertices % 64)));
}

void Bitset::allocate(int n_vertices){
    limbs = ((n_vertices-1)/(ALIGNMENT*8) + 1)*(ALIGNMENT/8);
#if MEM_DBG
    cout << "limbs " << limbs << " n_vertices " << n_vertices << endl;
#endif
    data = (uint64_t*) aligned_alloc(ALIGNMENT, limbs*8);
}

Bitset::~Bitset(){
    if (!shallow){
#if MEM_DBG
        cout << "freeing a bitset" << (size_t) data << endl;
#endif
        free(data);
    }
}

Bitset::Bitset(const bool* set_bits, int n_vertices){
    shallow = false;
    allocate(n_vertices);

    for(int i=0; i < (n_vertices-1)/64 + 1; i++)
        data[i] = 0;

    for(int i=0; i<n_vertices; i++){
        if (set_bits[i])
            set(i);
    }
}

Bitset::Bitset(const Bitset& obj){
#if MEM_DBG
    cout << "create a shallow bitset" << endl;
#endif
    data = obj.data;
    limbs = obj.limbs;
    shallow = true;
}

// Vertex

inline void KPartiteKClique::Vertex::set_weight(){
    int counter = 0;
    int tmp;
    if (!problem->current_graph().active_vertices->has(index)){
        weight = 0;
        return;
    }
    for (int i=0; i<get_k(); i++){
        tmp = intersection_count(get_active_vertices(), i);
        counter += tmp;
        if (!tmp){
            weight = 0;
            return;
        }
    }
    weight = counter;
}

inline void KPartiteKClique::KPartiteGraph::pop_last_vertex(){
    Vertex& v = vertices.back();
    part_sizes[v.part] -= 1;
    active_vertices[0].unset(v.index);
    vertices.pop_back();
}

inline KPartiteKClique::Vertex* KPartiteKClique::KPartiteGraph::last_vertex(){
    if (!vertices.size())
        return NULL;
    Vertex& v = vertices.back();

    // Remove all vertices,
    // that are no longer
    // an option.
    while (!v.weight){
#if DBG
        cout << "actual remove of vertex " << v.index << endl;
#endif
        pop_last_vertex();
        if (!vertices.size())
            return NULL;
        v = vertices.back();
    }
#if DBG
    cout << "last vertex is " << v.index << " with weight " << v.weight << endl;
#endif
    return &v;
}

// KPartiteGraph

bool KPartiteKClique::KPartiteGraph::is_valid(){
    for (int i=0; i<get_k(); i++){
        if (part_sizes[i] == 0)
            return false;
    }
    return true;
}

void KPartiteKClique::KPartiteGraph::init(KPartiteKClique* problem, bool fill){
    active_vertices = new Bitset(problem[0].n_vertices, fill);
    vertices = vector<Vertex>();
    part_sizes = new int[problem[0].k];
    int counter = 0;
    for (int i=0; i<problem[0].k; i++){
        part_sizes[i] = problem[0].parts[i+1];
        counter += part_sizes[i];
    }
    this->problem = problem;
}

KPartiteKClique::KPartiteGraph::KPartiteGraph(){
    active_vertices = NULL;
    part_sizes = NULL;
}

KPartiteKClique::KPartiteGraph::~KPartiteGraph(){
#if MEM_DBG
    cout << "deleting KPartiteGraph " << problem->k << (size_t) active_vertices << endl;
#endif
    delete active_vertices;
    delete[] part_sizes;
}

bool KPartiteKClique::KPartiteGraph::select(KPartiteKClique::KPartiteGraph& next){
    Vertex* v = last_vertex();
    if (!v)
        return false;

    // Copy the current sizes.
    for (int i=0; i<get_k(); i++)
        next.part_sizes[i] = part_sizes[i];

    // pelect v.
    problem[0]._k_clique[v->part] = v->index;
    v->intersection(next.active_vertices[0], active_vertices[0]);
#if DBG
    cout << "select the vertex " << v->index << endl;
#endif

    // v may no longer be selected.
    // In current not, because we have removed it.
    // In next not, because it is selected already.
    pop_last_vertex();
    next.vertices.assign(vertices.begin(), vertices.end());

    // Note that above, the part size of ``_k_clique[v.part]``
    // current is one smaller than next.
    // This is intentional, as in next the vertex was selected, not
    // removed.

    // Raise the current
    // depth, such that the
    // weights get set
    // accordingly.
    problem[0].current_depth += 1;

    next.set_weights();
    sort(next.vertices.begin(), next.vertices.end());
#if DBG
    for (Vertex& v1: next.vertices){
        cout << v1.index << " " << v1.weight << endl;
    }
#endif

    return true;
}

// KPartiteKClique

bool KPartiteKClique::traceback(){
    while (current_depth >= 1){
        current_depth -= 1;
        if (current_graph().is_valid())
            return true;
    }
    return false;
}

bool KPartiteKClique::next(){
    /*
    Set the next clique.
    Return whether there is a next clique.
    */
    while (true){
        if (current_depth < k-1){
#if DBG
            cout << current_depth << " depth" << endl;
#endif
            if (!current_graph().select(next_graph())){
#if DBG
                cout << "need to trace back" << endl;
#endif
                if (!traceback())
                    // Out of options.
                    return false;
            }
        } else {
#if DBG
            cout << "to the end" << endl;
#endif
            Vertex* vpt = current_graph().last_vertex();
            if (!vpt){
#if DBG
                cout << "found no vertex" << endl;
#endif
                if (!traceback()){
#if DBG
                    cout << "out of options" << endl;
#endif
                    // Out of options.
                    return false;
                }
            } else {
#if DBG
                cout << "found something" << endl;
                cout << "it is vertex " << vpt[0].index << endl;
#endif
                _k_clique[vpt[0].part] = vpt[0].index;
                current_graph().pop_last_vertex();
#if DBG
                cout << "will return true" << endl;
#endif
                return true;
            }
        }
    }
}

KPartiteKClique::KPartiteKClique(){
    _k_clique = NULL;
    parts = NULL;
    all_vertices = NULL;
    recursive_graphs = NULL;
}

KPartiteKClique::KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k){
    assert(k>0);

    current_depth = 0;

    _k_clique = new int[k];
    parts = new int[k+1];
    for (int i=0; i<k; i++){
        parts[i] = first_per_part[i];
    }
    parts[k] = n_vertices;
    this->n_vertices = n_vertices;
    this->k = k;

    recursive_graphs = new KPartiteGraph[k];
    for (int i=0; i<k; i++)
        recursive_graphs[i].init(this, i==0);

    all_vertices = new Vertex[n_vertices];
    int current_part = 0;
    for (int i=0; i<n_vertices; i++){
        while ((current_part < k-1) && (i >= parts[current_part + 1]))
            current_part += 1;
        all_vertices[i].init(this, incidences[i], n_vertices, current_part, i);
    }

    recursive_graphs[0].vertices.assign(all_vertices, all_vertices + n_vertices);
    recursive_graphs[0].set_weights();
}


KPartiteKClique::~KPartiteKClique(){
#if MEM_DBG
    cout << "hello" << endl;
#endif
    delete[] _k_clique;
    delete[] parts;
    delete[] all_vertices;
    delete[] recursive_graphs;
}
