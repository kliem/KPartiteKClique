#if __POPCNT__ || __BMI__
    #include <immintrin.h>
#endif

#include "kpkc.h"

#include <algorithm>
#include <csignal>
#include <iostream>
#include <stdlib.h>
#include <stdexcept>

// Handle keyboard
// interrupts

volatile sig_atomic_t kpkc_interrupted = 0;

void interrupt_signal_handler(int signal) {
    kpkc_interrupted = 1;
}

struct sigaction prev_action;
struct sigaction prev_action2;

struct sigaction sigIntHandler;

#define RESTORE_SIGNALS sigaction(SIGINT, &prev_action, NULL); \
                        sigaction(SIGALRM, &prev_action2, NULL);

#define REGISTER_SIGNALS \
    sigIntHandler.sa_handler = interrupt_signal_handler; \
    sigemptyset(&sigIntHandler.sa_mask); \
    sigIntHandler.sa_flags = 0; \
    sigaction(SIGINT, &sigIntHandler, &prev_action); \
    sigaction(SIGALRM, &sigIntHandler, &prev_action2);

#define CHECK_FOR_INTERRUPT             \
    if (kpkc_interrupted) {                              \
        kpkc_interrupted = 0; \
        RESTORE_SIGNALS \
        throw runtime_error("computation with kpkc was interrupted"); \
}

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

inline int first_in_limb(uint64_t i){
    // Return the position of the first bit.
    //
    // Assumes that ``i`` is nonzero.
#if (__BMI__) && (INTPTR_MAX == INT64_MAX)
    return _tzcnt_u64(i);
#else
    int output = 63;
    (i & 0x00000000FFFFFFFF) ? output -= 32 : (i >>= 32);
    (i & 0x000000000000FFFF) ? output -= 16 : (i >>= 16);
    (i & 0x00000000000000FF) ? output -=  8 : (i >>=  8);
    (i & 0x000000000000000F) ? output -=  4 : (i >>=  4);
    (i & 0x0000000000000003) ? output -=  2 : (i >>=  2);
    if (i & 0x0000000000000001) output -=  1;
    return output;
#endif
}

inline uint64_t one_set_bit(int n){
    return ((uint64_t) 1) << (n % 64);
}

// Bitsets

Bitset::Bitset(int n_vertices, bool fill){
    /*
    Initalize bitset.

    Fill if ``fill``.
    */
    allocate(n_vertices);
    if (!fill)
        return;

    // Fill.
    for(int i=0; i<n_vertices/64; i++){
        data[i] = -1;
    }
    // Remove trailing bits.
    if (n_vertices % 64)
        data[n_vertices/64] = lower_n_bits(n_vertices % 64);
}

Bitset::Bitset(const bool* set_bits, int n_vertices){
    /*
    Initialize bitset with the given bits.
    */
    allocate(n_vertices);

    for(int i=0; i < (n_vertices-1)/64 + 1; i++)
        data[i] = 0;

    for(int i=0; i<n_vertices; i++){
        if (set_bits[i])
            set(i);
    }
}

Bitset::~Bitset(){
    delete[] data;
}

inline void Bitset::intersection_assign(Bitset& l, Bitset& r){
    // Assumes all of same length.
    for (int i=0; i<limbs; i++)
        data[i] = l[i] & r[i];
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
        // Remove the lower bits.
        start_limb &= ~lower_n_bits(start % 64);

    if (stop % 64 == 0)
        return counter + popcount(start_limb);

    uint64_t end_limb = 0;
    if (stop/64 < limbs){
        end_limb = data[stop/64] & r[stop/64];
        // Remove the upper bits.
        end_limb &= lower_n_bits(stop % 64);
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

inline int Bitset::count(int start, int stop){
    /*
    Count the number of set bits in ``this``
    in ``range(start, stop)``.
    */
    int counter = 0;
    // The easy part, count any complete ``uint64_t``.
    for (int i=start/64 + 1; i< stop/64; i++)
        counter += popcount(data[i]);

    uint64_t start_limb = data[start/64];
    if (start % 64)
        // Remove the lower bits.
        start_limb &= ~lower_n_bits(start % 64);

    if (stop % 64 == 0)
        return counter + popcount(start_limb);

    uint64_t end_limb = 0;
    if (stop/64 < limbs){
        end_limb = data[stop/64];
        // Remove the upper bits.
        end_limb &= lower_n_bits(stop % 64);
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

inline int Bitset::first(int start){
    /*
    Return the first bit in ``this``
    in ``range(start, stop)``.

    It assumes that the range is valid and that there is at least on non-zero bit.
    */
    uint64_t start_limb = data[start/64];
    if (start % 64)
        // Remove the lower bits.
        start_limb &= ~lower_n_bits(start % 64);

    int counter = (start/64)*64;
    if (start_limb)
        return counter + first_in_limb(start_limb);
    else
        counter += 64;

    // The easy part, count any complete ``uint64_t``.
    for (int i=start/64 + 1; i< limbs; i++){
        if (data[i])
            return counter + first_in_limb(data[i]);
        else
            counter += 64;
    }

    return limbs * 64;
}

void Bitset::set(int index){
    data[index/64] |= one_set_bit(index % 64);
}

void Bitset::unset(int index){
    data[index/64] &= ~one_set_bit(index % 64);
}

bool Bitset::has(int index){
    return data[index/64] & one_set_bit(index % 64);
}

void Bitset::allocate(int n_vertices){
    limbs = ((n_vertices-1)/64+ 1);
    data = new uint64_t[limbs];
}

// KPartiteKClique_base

KPartiteKClique_base::KPartiteKClique_base(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k){
    if (k <= 0) throw invalid_argument("k must be at least 1");

    current_depth = 0;

    _k_clique = new int[k];
    parts = new int[k+1];
    for (int i=0; i<k; i++){
        parts[i] = first_per_part[i];
    }
    parts[k] = n_vertices;

    for (int i=0; i<k; i++){
        if (parts[i+1] - parts[i] == 0)
            throw invalid_argument("parts may not be empty");
    }

    this->n_vertices = n_vertices;
    this->k = k;

    all_vertices = new Vertex_template[n_vertices];
    int current_part = 0;
    for (int i=0; i<n_vertices; i++){
        while ((current_part < k-1) && (i >= parts[current_part + 1]))
            current_part += 1;
        Vertex_template tmp(this, incidences[i], n_vertices, current_part, i);
        swap(tmp, all_vertices[i]);
    }
}

KPartiteKClique_base::~KPartiteKClique_base(){
    delete[] _k_clique;
    delete[] parts;
    delete[] all_vertices;
}

KPartiteKClique_base::KPartiteGraph* KPartiteKClique_base::current_graph(){
    throw runtime_error("a derived class must implement this");
}

KPartiteKClique_base::KPartiteGraph* KPartiteKClique_base::next_graph(){
    throw runtime_error("a derived class must implement this");
}

bool KPartiteKClique_base::backtrack(){
    /*
    Go the the last valid graph.

    If none exists, return false.
    */
    while (current_depth >= 1){
        current_depth -= 1;
        if (current_graph()->permits_another_choice())
            return true;
    }
    return false;
}

bool KPartiteKClique_base::next(){
    /*
    Set the next clique.
    Return whether there is a next clique.
    */
    REGISTER_SIGNALS

    while (true){
        if (current_depth < k-1){

            // Note that the interrupt can also be abused to pause.
            CHECK_FOR_INTERRUPT

            if (!current_graph()->select(next_graph())){
                if (!backtrack()){
                    // Out of options.
                    RESTORE_SIGNALS
                    return false;
                }
            }
        } else {
            if (!current_graph()->select()){
                if (!backtrack()){
                    // Out of options.
                    RESTORE_SIGNALS
                    return false;
                }
            } else {
                RESTORE_SIGNALS
                return true;
            }
        }
    }
}

// KPartiteKClique

KPartiteKClique::KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k, const int prec_depth)
        : KPartiteKClique_base::KPartiteKClique_base(incidences, n_vertices, first_per_part, k){
    this->prec_depth = prec_depth;

    recursive_graphs = new KPartiteGraph[k];
    for (int i=0; i<k; i++){
        KPartiteGraph tmp(this, i==0);
        swap(tmp, recursive_graphs[i]);
    }

    // Assign vertices to the first graph.
    recursive_graphs->vertices.assign(all_vertices, all_vertices + n_vertices);

    // Set weights twice, if there is new knowledge already.
    if (recursive_graphs->set_weights())
        recursive_graphs->set_weights();

    sort(recursive_graphs->vertices.begin(), recursive_graphs->vertices.end());
}

KPartiteKClique::~KPartiteKClique(){
    delete[] recursive_graphs;
}

KPartiteKClique_base::KPartiteGraph* KPartiteKClique::current_graph(){
    return (KPartiteKClique_base::KPartiteGraph*) &(recursive_graphs[current_depth]);
}

KPartiteKClique_base::KPartiteGraph* KPartiteKClique::next_graph(){
    return (KPartiteKClique_base::KPartiteGraph*) &(recursive_graphs[current_depth + 1]);
}

// FindClique

FindClique::FindClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k)
        : KPartiteKClique_base::KPartiteKClique_base(incidences, n_vertices, first_per_part, k){

    recursive_graphs = new KPartiteGraph[k];
    for (int i=0; i<k; i++){
        KPartiteGraph tmp(this, i==0);
        swap(tmp, recursive_graphs[i]);
    }

    // Take care of trivial parts.
    // ``set_part_sizes`` assumes that parts of size 1 were selected
    // already, which would not be the case, if the part is 1 to start
    // with.
    for (int i=0; i<k; i++){
        if (parts[i+1] - parts[i] == 1)
            current_graph()->part_sizes[i] = 2;

    }

    recursive_graphs->set_part_sizes();
}

FindClique::~FindClique(){
    delete[] recursive_graphs;
}

// Vertex_template

KPartiteKClique_base::Vertex_template::Vertex_template(KPartiteKClique_base* problem, const bool* incidences, int n_vertices, int part, int index){
    bitset = new Bitset(incidences, n_vertices);
    this->part = part;
    this->index = index;
    this->problem = problem;

    // Set each vertex adjacent to itself.
    // This is important, so that after selecting a vertex
    // the corresponding part will have one ``active_vertex``.
    bitset->set(index);
}

inline KPartiteKClique_base::Vertex_template::~Vertex_template(){
    delete bitset;
}

void KPartiteKClique_base::swap(KPartiteKClique_base::Vertex_template& a, KPartiteKClique_base::Vertex_template& b){
    std::swap(a.bitset, b.bitset);
    std::swap(a.part, b.part);
    std::swap(a.index, b.index);
    std::swap(a.problem, b.problem);
}

// KPartiteGraph in KPartiteKClique_base

KPartiteKClique_base::KPartiteGraph::KPartiteGraph(KPartiteKClique_base* problem, bool fill){
    active_vertices = new Bitset(problem->n_vertices, fill);
    part_sizes = new int[problem->k+1];
    for (int i=0; i < problem->k; i++){
        part_sizes[i] = problem->parts[i+1] - problem->parts[i];
    }
    this->problem = problem;
}

KPartiteKClique_base::KPartiteGraph::KPartiteGraph(){
    active_vertices = NULL;
    part_sizes = NULL;
}

KPartiteKClique_base::KPartiteGraph::~KPartiteGraph(){
    delete active_vertices;
    delete[] part_sizes;
}

bool KPartiteKClique_base::KPartiteGraph::permits_another_choice(){
    throw runtime_error("a derived class must implement this");
}

bool KPartiteKClique_base::KPartiteGraph::select(KPartiteKClique_base::KPartiteGraph* next){
    throw runtime_error("a derived class must implement this");
}

bool KPartiteKClique_base::KPartiteGraph::select(){
    throw runtime_error("a derived class must implement this");
}

// Vertex

inline KPartiteKClique::Vertex::Vertex(const KPartiteKClique_base::Vertex_template& obj){
    bitset = obj.bitset;
    part = obj.part;
    weight = -1;
    index = obj.index;
    problem = (KPartiteKClique*) obj.problem;

    if (1 != bitset->count(get_parts()[part], get_parts()[part + 1]))
        throw invalid_argument("the graph is not k-partite");
}

inline KPartiteKClique::Vertex::Vertex(const Vertex& obj){
    bitset = obj.bitset;
    weight = obj.weight;
    part = obj.part;
    index = obj.index;
    problem = obj.problem;
}

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
    if (!active_vertices.has(index)){
        weight = 0;
        return false;
    }
    if (problem->current_depth > problem->prec_depth){
        weight = 1;
        return false;
    }
    for (int i=0; i<get_k(); i++){
        tmp = intersection_count(active_vertices, i);
        counter += tmp;
        if (!tmp){
            // This vertex would not allow for a k-clique anymore.
            weight = 0;
            active_vertices.unset(index);
            return true;
        }
    }
    weight = counter;
    return false;
}

// KPartiteGraph in KPartiteKClique

void KPartiteKClique_base::swap(KPartiteKClique_base::KPartiteGraph& a, KPartiteKClique_base::KPartiteGraph& b){
    std::swap(a.active_vertices, b.active_vertices);
    std::swap(a.part_sizes, b.part_sizes);
    std::swap(a.problem, b.problem);
}

void KPartiteKClique::swap(KPartiteKClique::KPartiteGraph& a, KPartiteKClique::KPartiteGraph& b){
    KPartiteKClique_base::swap(*((KPartiteKClique_base::KPartiteGraph *) &a), *((KPartiteKClique_base::KPartiteGraph *) &b));
    std::swap(a.problem, b.problem);
    std::swap(a.vertices, b.vertices);
}

void FindClique::swap(FindClique::KPartiteGraph& a, FindClique::KPartiteGraph& b){
    KPartiteKClique_base::swap(*((KPartiteKClique_base::KPartiteGraph *) &a), *((KPartiteKClique_base::KPartiteGraph *) &b));
    std::swap(a.problem, b.problem);
    std::swap(a.selected_part, b.selected_part);
}

KPartiteKClique::KPartiteGraph::KPartiteGraph() : KPartiteKClique_base::KPartiteGraph::KPartiteGraph(){
    vertices = vector<Vertex>();
}

KPartiteKClique::KPartiteGraph::KPartiteGraph(KPartiteKClique* problem, bool fill) : KPartiteKClique_base::KPartiteGraph::KPartiteGraph((KPartiteKClique_base*) problem, fill){
    vertices = vector<Vertex>();
    this->problem = problem;
}

inline void KPartiteKClique::KPartiteGraph::pop_last_vertex(){
    Vertex& v = vertices.back();
    part_sizes[v.part] -= 1;
    active_vertices->unset(v.index);
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

    // Remove all vertices,
    // that are no longer
    // a valid choice.
    while (!v.weight){
        pop_last_vertex();
        if (!vertices.size())
            return NULL;
        v = vertices.back();
    }
    return &v;
}

inline bool KPartiteKClique::KPartiteGraph::permits_another_choice(){
    for (int i=0; i<get_k(); i++){
        if (part_sizes[i] == 0)
            return false;
    }
    return true;
}

bool KPartiteKClique::KPartiteGraph::select(KPartiteKClique_base::KPartiteGraph* next2){
    /*
    Select the last (valid) vertex of the current graph set up the next graph
    to be all vertices connected to that last vertex.

    Return false, if there are no vertices left.
    */
    KPartiteGraph* next = (KPartiteGraph*) next2;
    Vertex* v = last_vertex();
    if (!v)
        return false;

    // Copy the current sizes.
    for (int i=0; i<get_k(); i++)
        next->part_sizes[i] = part_sizes[i];

    // Select v.
    problem->_k_clique[v->part] = v->index;
    intersection(*next->active_vertices, *v, *active_vertices);

    int part = v->part;

    // v may no longer be selected.
    // In current not, because we have removed it.
    // In next not, because it is selected already.
    pop_last_vertex();
    next->vertices.assign(vertices.begin(), vertices.end());

    // Note the part size of ``_k_clique[v.part]``:
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

    next->set_weights();

    // When setting weights, we also discover that some vertices are not
    // longer possible.
    // This information can be used to call set weights again and get
    // preciser results.
    // This appears to pay off if we are not too deep in the recursion
    // tree.
    if (problem->current_depth < problem->prec_depth && next->set_weights())
        next->set_weights();

    sort(next->vertices.begin(), next->vertices.end());

    return true;
}

bool KPartiteKClique::KPartiteGraph::select(){
    /*
    Select the last (valid) vertex of the current graph and return
    whether it exists.

    It is assumed that there is no next graph.
    */
    assert(problem->current_depth == problem->k -1);
    Vertex* v = last_vertex();
    if (!v)
        return false;

    // Select v.
    problem->_k_clique[v->part] = v->index;
    int part = v->part;

    // v may no longer be selected.
    // In current not, because we have removed it.
    pop_last_vertex();

    if (part_sizes[part] == 0)
        vertices.resize(0);

    return true;
}

// KPartiteGraph in FindClique

FindClique::KPartiteGraph::KPartiteGraph(FindClique* problem, bool fill) : KPartiteKClique_base::KPartiteGraph::KPartiteGraph((KPartiteKClique_base*) problem, fill){
    this->problem = problem;
}

inline bool FindClique::KPartiteGraph::set_part_sizes(){
    /*
    Set the sizes of the parts.

    Return false, if any part has size 0, otherwise true
    */
    int i;
    int min_so_far = problem->n_vertices;
    selected_part = -1;
    for(i=0; i < problem->k; i++){
        if (part_sizes[i] != 1){
            int j = count(i);
            part_sizes[i] = j;
            if (j == 0){
                // this part is empty; need to backtrack
                selected_part = -2;
                return false;
            }
            if (j == 1){
                // this part has a unique choice
                selected_part = i;
                return true;
            } else if (j < min_so_far){
                min_so_far = j;
                selected_part = i;
            }
        }
    }
    return true;
}

bool FindClique::KPartiteGraph::select(KPartiteKClique_base::KPartiteGraph* next2){
    /*
    Select the first vertex in the smallest part.
    Return false, if there are no vertices left.
    */
    KPartiteGraph* next = (KPartiteGraph*) next2;
    assert(selected_part != -1); // Should not be called, if we found a clique already.
    if (!part_sizes[selected_part])
        return false;

    // Copy the current sizes.
    for (int i=0; i<get_k(); i++)
        next->part_sizes[i] = part_sizes[i];

    next->part_sizes[selected_part] = 1;

    // Select v.
    int v = first(selected_part);
    if (v == -1) return false;

    problem->_k_clique[selected_part] = v;
    intersection(*next->active_vertices, problem->all_vertices[v], *active_vertices);

    // v may no longer be selected on the next call.
    pop_vertex(selected_part, v);

    // Raise the current
    // depth, such that the
    // parts get set
    // accordingly.
    problem->current_depth += 1;

    return next->set_part_sizes();
}

bool FindClique::KPartiteGraph::select(){
    /*
    Select the first vertex in the smallest part.
    Return false, if there are no vertices left.

    It is assumed that there is no next graph.
    */
    assert(selected_part != -1); // Should not be called, if we found a clique already.
    assert(problem->current_depth == problem->k -1);
    if (!part_sizes[selected_part])
        return false;

    // Select v.
    int v = first(selected_part);
    if (v == -1) return false;

    problem->_k_clique[selected_part] = v;

    // v may no longer be selected on the next call.
    pop_vertex(selected_part, v);

    return true;
}

KPartiteKClique_base::KPartiteGraph* FindClique::current_graph(){
    return (KPartiteKClique_base::KPartiteGraph*) &(recursive_graphs[current_depth]);
}

KPartiteKClique_base::KPartiteGraph* FindClique::next_graph(){
    return (KPartiteKClique_base::KPartiteGraph*) &(recursive_graphs[current_depth + 1]);
}

inline bool FindClique::KPartiteGraph::permits_another_choice(){
    return ((selected_part >= 0) && (part_sizes[selected_part]));
}

