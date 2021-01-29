#if __POPCNT__
    #include <immintrin.h>
#endif

#include "k_partite_k_clique.h"

#include <algorithm>
#include <cstdlib>
#include <cassert>

const int ALIGNMENT = 8;

// Bitset helpers.

inline int popcount(uint64_t i){
#if (__POPCNT__) && (INTPTR_MAX == INT64_MAX)
    return _mm_popcnt_u64(i);
#else
    i = i - ((i >> 1) & 0x5555555555555555ULL);
    i = (i & 0x3333333333333333ULL) + ((i >> 2) & 0x3333333333333333ULL);
    return ( ((i + (i >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
#endif
}

inline uint64_t lower_n_bits(int n){
    return ((uint64_t) -1) >> (64 - n);
}

inline uint64_t one_set_bit(int n){
    return ((uint64_t) 1) << (n % 64);
}


// Bitsets

Bitset::Bitset(int n_vertices, int k, int* first_per_part, bool fill){
    /*
    Initalize bitset.

    Fill if ``fill``.
    */
    allocate(n_vertices, k, first_per_part);
    if (!fill)
        return;

    // Fill.
    for(int i=0; i<limbs; i++){
        data[i] = -1;
    }
    // Remove trailing bits.
    for(int i=0; i<k; i++){
        int n_trailing_bits = (first_per_part[i+1] - first_per_part[i]) % 64;
        if (n_trailing_bits)
            data[first_limb_per_part[i+1] - 1] = lower_n_bits(n_trailing_bits);
    }
}

Bitset::Bitset(const bool* set_bits, int n_vertices, int k, int* first_per_part){
    /*
    Initialize bitset with the given bits.
    */
    allocate(n_vertices, k, first_per_part);

    for(int i=0; i < limbs; i++)
        data[i] = 0;

    int part = -1;
    for(int i=0; i<n_vertices; i++){
        if (i == first_per_part[part + 1])
            part += 1;
        if (set_bits[i])
            set(i, part);
    }
}

Bitset::Bitset(const Bitset& obj){
#if MEM_DBG
    cout << "create a shallow bitset" << endl;
#endif
    data = obj.data;
    limbs = obj.limbs;
    first_limb_per_part = obj.first_limb_per_part;
    offset = obj.offset;
    shallow = true;
}

Bitset::~Bitset(){
    if (!shallow){
#if MEM_DBG
        cout << "freeing a bitset" << (size_t) data << endl;
#endif
        delete[] data;
        delete[] first_limb_per_part;
        delete[] offset;
    }
}

void Bitset::intersection_assign(Bitset& l, Bitset& r){
    // Assumes all of same length.
    for (int i=0; i<limbs; i++)
        data[i] = l[i] & r[i];
}

inline int Bitset::intersection_count(Bitset& r, int part){
    /*
    Count the number of set bits in ``this & r``
    in ``range(start, stop)``.
    */
    int counter = 0;
    int i;
    for (i = first_limb_per_part[part]; i < first_limb_per_part[part + 1]; i++)
        counter += popcount(data[i] & r[i]);
    return counter;
}

inline bool Bitset::intersection_is_empty(Bitset& r, int part){
    /*
    Return whether the part is empty when intersected with ``r``.
    */
    int i;
    for (i = first_limb_per_part[part]; i < first_limb_per_part[part + 1]; i++){
        if (data[i] & r[i])
            return false;
    }
    return true;
}

bool Bitset::is_empty(int part){
    /*
    Currently not used.
    */
    int counter = 0;
    int i;
    for (i = first_limb_per_part[part]; i < first_limb_per_part[part + 1]; i++){
        if (data[i])
            return false;
    }
    return true;
}

void Bitset::set(int index, int part){
    index += offset[part];
    data[index/64] |= one_set_bit(index % 64);
}

void Bitset::unset(int index, int part){
    index += offset[part];
    data[index/64] &= ~one_set_bit(index % 64);
}

bool Bitset::has(int index, int part){
    index += offset[part];
    return data[index/64] & one_set_bit(index % 64);
}

void Bitset::allocate(int n_vertices, int k, int* first_per_part){
    shallow = false;
    offset = new int[k+1];
    first_limb_per_part = new int[k+1];
    first_limb_per_part[0] = 0;
    offset[0] = 0;
    for (int i = 0; i < k; i++){
        int partsize = first_per_part[i+1] - first_per_part[i];
        int n_limbs = (partsize - 1)/64 + 1;
        first_limb_per_part[i+1] = first_limb_per_part[i] + n_limbs;
        offset[i+1] = 64*n_limbs - partsize + offset[i];
    }
    limbs = first_limb_per_part[k];
#if MEM_DBG
    cout << "limbs " << limbs << " n_vertices " << n_vertices << endl;
#endif
    data = new uint64_t[limbs];
}


// Vertex

inline bool KPartiteKClique::Vertex::set_weight(){
    // The weight is the number of vertices that are still available when
    // selecting this vertex.
    // However, when selecting the vertex no longer allows a k-clique,
    // the weight is always set to 0.
    //
    // Return ``true`` if and only if this vertex is newly removed.
    int counter = 0;
    int tmp;
    Bitset& active_vertices = get_active_vertices();
    if (!active_vertices.has(index, part)){
        weight = 0;
        return false;
    }
    if (problem->current_depth > 5){
        weight = 1;
        return false;
    }
    for (int i=0; i<get_k(); i++){
        tmp = intersection_count(active_vertices, i);
        counter += tmp;
        if (!tmp){
            // This vertex would not allow for a k-clique anymore.
            weight = 0;
            active_vertices.unset(index, part);
            return true;
        }
    }
    weight = counter;
    return false;
}

inline bool KPartiteKClique::Vertex::is_valid(){
    // The weight is the number of vertices that are still available when
    // selecting this vertex.
    // However, when selecting the vertex no longer allows a k-clique,
    // the weight is always set to 0.
    //
    // Return ``true`` if and only if this vertex is newly removed.
    return true;
    if (problem->current_depth <= 10)
        return true;
    int counter = 0;
    int tmp;
    Bitset& active_vertices = get_active_vertices();
    if (!active_vertices.has(index, part)){
        return false;
    }
    /*
    this->set_weight();
    return weight != 0;
    */
    for (int i=0; i<get_k(); i++){
        if (intersection_is_empty(active_vertices, i))
            return false;
    }
    return true;
}


// KPartiteGraph

inline void KPartiteKClique::KPartiteGraph::pop_last_vertex(){
    Vertex& v = vertices.back();
    part_sizes[v.part] -= 1;
    active_vertices->unset(v.index, v.part);
    vertices.pop_back();
}

inline KPartiteKClique::Vertex* KPartiteKClique::KPartiteGraph::last_vertex(){
    /*
    Get the last vertex, which is (possibly) a valid choice.

    Pop all vertices that are no longer valid choices.
    */
    if (!vertices.size())
        return NULL;
    Vertex& v = vertices.back();
    if (v.weight && !v.is_valid())
        v.weight = 0;

    // Remove all vertices,
    // that are no longer
    // a valid choice.
    while (!v.weight){
#if DBG
        cout << "actual remove of vertex " << v.index << endl;
#endif
        pop_last_vertex();
        if (!vertices.size())
            return NULL;
        v = vertices.back();
        if (v.weight && !v.is_valid())
            v.weight = 0;
    }
#if DBG
    cout << "last vertex is " << v.index << " with weight " << v.weight << endl;
#endif
    return &v;
}

bool KPartiteKClique::KPartiteGraph::is_valid(){
    for (int i=0; i<get_k(); i++){
        if (part_sizes[i] == 0)
            return false;
    }
    return true;
}

void KPartiteKClique::KPartiteGraph::init(KPartiteKClique* problem, bool fill){
    active_vertices = new Bitset(problem->n_vertices, problem->k, problem->parts, fill);
    part_sizes = new int[problem->k];
    for (int i=0; i < problem->k; i++){
        part_sizes[i] = problem->parts[i+1] - problem->parts[i];
    }
    this->problem = problem;
}

KPartiteKClique::KPartiteGraph::KPartiteGraph(){
    active_vertices = NULL;
    part_sizes = NULL;
    vertices = vector<Vertex>();
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

    // select v.
    problem->_k_clique[v->part] = v->index;
    v->intersection(*next.active_vertices, *active_vertices);
#if DBG
    cout << "select the vertex " << v->index << endl;
#endif

    int part = v->part;

    // v may no longer be selected.
    // In current not, because we have removed it.
    // In next not, because it is selected already.
    pop_last_vertex();
    next.vertices.assign(vertices.begin(), vertices.end());

    // Note that above, the part size of ``_k_clique[v.part]``
    // current is one smaller than next.
    // This is intentional, as in next the vertex was selected, not
    // removed.

    // Apply the new knowledge.
    if (part_sizes[part] == 1){
        this->set_weights();
        sort(vertices.begin(), vertices.end());
    }
    if (part_sizes[part] == 0)
        vertices.resize(0);

    // Raise the current
    // depth, such that the
    // weights get set
    // accordingly.
    problem->current_depth += 1;

    next.set_weights();
    if (problem->current_depth < 5 && next.set_weights()){
        next.set_weights();
    }
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
                cout << "it is vertex " << vpt->index << endl;
#endif
                _k_clique[vpt->part] = vpt->index;
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

    recursive_graphs->vertices.assign(all_vertices, all_vertices + n_vertices);
    recursive_graphs->set_weights();
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
