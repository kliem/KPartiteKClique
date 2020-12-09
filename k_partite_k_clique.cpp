#include "k_partite_k_clique.h"

#include <algorithm>
#include <cstdlib>
#include <cassert>

const int ALIGNMENT = 8;

void Bitset::intersection_assign(Bitset& l, Bitset& r){
    // Assumes all of same length.
    for (int i=0; i<limbs; i++)
        data[i] = l[i] & r[i];
}

int popcount(uint64_t i){
    i = i - ((i >> 1) & 0x5555555555555555ULL);
    i = (i & 0x3333333333333333ULL) + ((i >> 2) & 0x3333333333333333ULL);
    return ( ((i + (i >> 4)) & 0x0f0f0f0f0f0f0f0fULL) * 0x0101010101010101ULL ) >> 56;
}

int Bitset::intersection_count(Bitset& r, int start, int stop){
    int counter = 0;
    for (int i=start/64 + 1; i< stop/64; i++)
        counter += popcount(data[i] & r[i]);

    uint64_t start_limb = data[start/64] & r[start/64];
    if (start % 64)
        start_limb &= ((uint64_t) -1) >> (64 - (start % 64));

    uint64_t end_limb;
    if (stop/64 < limbs){
        end_limb = data[stop/64] & r[stop/64];
        if (stop % 64)
            end_limb &= ~(((uint64_t) -1) >> (64 - (stop % 64)));
    }

    if (start/64 == stop/64){
        counter += popcount(start_limb & end_limb);
    } else {
        if (stop/64 < limbs){
            counter += popcount(start_limb) + popcount(end_limb);
        } else {
            counter += popcount(start_limb);
        }
    }
    return counter;
}

void Bitset::set(int index){
    data[index/64] |= ((uint64_t) 1) << (index % 64);
}

void Bitset::unset(int index){
    data[index/64] &= ~(((uint64_t) 1) << (index % 64));
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
    cout << "limbs " << limbs << " n_vertices " << n_vertices << endl;
    data = (uint64_t*) aligned_alloc(ALIGNMENT, limbs*8);
}

Bitset::~Bitset(){
    if (!shallow){
        cout << "freeing a bitset" << (size_t) data << endl;
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
    data = obj.data;
    limbs = obj.limbs;
    shallow = true;
}

void KPartiteKClique::KPartiteGraph::Vertex::set_weight(){
    int counter = 0;
    int tmp;
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

void KPartiteKClique::KPartiteGraph::pop_last_vertex(){
    Vertex& v = vertices.back();
    part_sizes[v.part] -= 1;
    active_vertices[0].unset(v.index);
    vertices.pop_back();
}

KPartiteKClique::KPartiteGraph::Vertex* KPartiteKClique::KPartiteGraph::last_vertex(){
    if (!vertices.size())
        return NULL;
    Vertex& v = vertices.back();

    // Remove all vertices,
    // that are no longer
    // an option.
    while (!v.weight){
        pop_last_vertex();
        if (!vertices.size())
            return NULL;
        v = vertices.back();
    }
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
    active_vertices = new Bitset(problem[0].n_vertices, fill);
    vertices = vector<Vertex>();
    part_sizes = new int[problem[0].k];
    this->problem = problem;
}

KPartiteKClique::KPartiteGraph::KPartiteGraph(){
    active_vertices = NULL;
    part_sizes = NULL;
}

KPartiteKClique::KPartiteGraph::~KPartiteGraph(){
    cout << "hello1 " << problem->k << (size_t) active_vertices << endl;
    delete active_vertices;
    delete[] part_sizes;
}


bool KPartiteKClique::KPartiteGraph::select(KPartiteKClique::KPartiteGraph& next){
    KPartiteGraph::Vertex* v = last_vertex();
    if (!v)
        return false;

    // Copy the current sizes.
    for (int i=0; i<get_k(); i++)
        next.part_sizes[i] = part_sizes[i];

    // Select v.
    problem[0]._k_clique[v->part] = v->index;
    v->intersection(next.active_vertices[0], active_vertices[0]);
    cout << "select the vertex " << v->index << endl;;

    // v may no longer be selected.
    // In current not, because we have removed it.
    // In next not, because it is selected already.
    pop_last_vertex();
    next.vertices = vertices;

    // Note that above, the part size of ``_k_clique[v.part]``
    // current is one smaller than next.
    // This is intentional, as in next the vertex was selected, not
    // removed.

    next.set_weights();
    sort(next.vertices.begin(), next.vertices.end());

    problem[0].current_depth += 1;

    return true;
}

bool KPartiteKClique::traceback(){
    while (current_depth >= 0){
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
            cout << current_depth << " depth" << endl;
            if (!current_graph().select(next_graph())){
                cout << "need to trace back" << endl;
                if (!traceback())
                    // Out of options.
                    return false;
            }
        } else {
            cout << "to the end" << endl;
            KPartiteGraph::Vertex* vpt = current_graph().last_vertex();
            if (!vpt){
                cout << "found no vertex" << endl;
                if (!traceback())
                    cout << "out of options" << endl;
                    // Out of options.
                    return false;
            } else {
                cout << "found something" << endl;
                _k_clique[vpt[0].part] = vpt[0].index;
                current_graph().pop_last_vertex();
                cout << "will return true" << endl;
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

    _k_clique = new int[k];
    parts = new int[k];
    for (int i=0; i<k; i++){
        parts[i] = first_per_part[i];
    }
    this->n_vertices = n_vertices;
    this->k = k;

    recursive_graphs = new KPartiteGraph[k];
    for (int i=0; i<k; i++)
        recursive_graphs[i].init(this, i==0);

    all_vertices = new KPartiteGraph::Vertex[n_vertices];
    int current_part = 0;
    for (int i=0; i<n_vertices; i++){
        while ((current_part < k-1) && (i >= parts[current_part + 1]))
            current_part += 1;
        all_vertices[i].init(recursive_graphs, incidences[i], n_vertices, current_part, i);
    }

    recursive_graphs[0].vertices.assign(all_vertices, all_vertices + n_vertices);
    recursive_graphs[0].set_weights();
}


KPartiteKClique::~KPartiteKClique(){
    cout << "hello" << endl;
    delete[] _k_clique;
    delete[] parts;
    delete[] all_vertices;
    delete[] recursive_graphs;
}
