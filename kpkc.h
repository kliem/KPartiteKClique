#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <algorithm>

// namespace k-partite-k-clicque
namespace kpkc
{
inline uint64_t one_set_bit(int n){
    return ((uint64_t) 1) << (n % 64);
}

class Bitset;

class Bitset {
    public:
        Bitset() = default;
        Bitset(int n_vertices, bool fill=false);
        Bitset(const bool* set_bits, int n_vertices);
        Bitset(const Bitset&) = delete;
        Bitset(Bitset&&) = default;
        Bitset& operator=(const Bitset&) = delete;
        Bitset& operator=(Bitset&&) = default;

        void unset(int index){
            data[index/64] &= ~one_set_bit(index % 64);
        }

        bool has(int index){
            return data[index/64] & one_set_bit(index % 64);
        }

        void set(int index);
        int intersection_count(Bitset& r, int start, int stop);
        int count(int start, int stop);
        int first(int start);
        void intersection_assign(Bitset& l, Bitset& r);
    protected:
        std::vector<uint64_t> data;
        int limbs{0};
        void allocate(int n_vertices);
        uint64_t& operator[](int i){ return data[i]; }
};

enum Algorithm { kpkc, FindClique };

template<Algorithm A>
class KPartiteKClique {
    class Vertex_template;
    class Vertex;
    class KPartiteGraph;

    public:
        KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k, const int prec_depth=5);
        virtual ~KPartiteKClique() = default;
        const int* k_clique(){ return _k_clique.data(); }
        bool next();

        KPartiteKClique(const KPartiteKClique&) = delete;
        KPartiteKClique& operator=(const KPartiteKClique&) = delete;
    private:
        std::vector<int> _k_clique;
        std::vector<int> parts;
        int k;
        int current_depth;
        int n_vertices;
        int prec_depth;
        std::vector<Vertex_template> all_vertices;
        std::vector<KPartiteGraph> recursive_graphs;
        KPartiteGraph* current_graph(){
            return &(recursive_graphs[current_depth]);
        }
        KPartiteGraph* next_graph(){
            return &(recursive_graphs[current_depth + 1]);
        }
        bool backtrack();
        void finish_init();
};

template class KPartiteKClique<kpkc>;
template class KPartiteKClique<FindClique>;

template<Algorithm A>
class KPartiteKClique<A>::Vertex_template {
    // Takes care of the memory allocation for vertices.

    friend KPartiteKClique<A>;
    friend KPartiteKClique<A>::Vertex;

    friend void intersection(Bitset& c, Vertex_template& l, Bitset& r){
        c.intersection_assign(l.bitset, r);
    }

    public:
        int index;  // The index in the original graph.
        int part;  // The part in the orginal graph.

        Vertex_template(const Vertex_template&) = delete;
        Vertex_template& operator=(const Vertex_template&) = delete;
        Vertex_template(Vertex_template&&) = default;
        Vertex_template& operator=(Vertex_template&&) = default;
        Vertex_template(KPartiteKClique<A>* problem, const bool* incidences, int n_vertices, int part, int index);

    private:
        Bitset bitset;
        KPartiteKClique<A>* problem;
};

template<Algorithm A>
class KPartiteKClique<A>::KPartiteGraph {
    public:
        Bitset active_vertices;
        std::vector<int> part_sizes;
        std::vector<Vertex> vertices;
        KPartiteGraph() = default;
        KPartiteGraph(KPartiteKClique<A>* problem, bool fill);

        bool permits_another_choice();
        bool select(KPartiteGraph* next);
        bool select();

        // kpkc
        Vertex* last_vertex();
        void pop_last_vertex();
        bool set_weights();

        // FindClique
        int selected_part;
        bool set_part_sizes();
        int count(int start, int stop){
            return active_vertices.count(start, stop);
        }
        int count(int part){
            return count(get_parts()[part], get_parts()[part+1]);
        }
        int first(int part);
        void pop_vertex(int part, int vertex);

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
        KPartiteGraph(KPartiteGraph&&) = default;
        KPartiteGraph& operator=(KPartiteGraph&&) = default;
    private:
        const int* get_parts() { assert(problem); return problem->parts.data(); }
        const int get_k() { assert(problem); return problem->k; }
        KPartiteKClique<A>::KPartiteGraph* current_graph(){ return problem->current_graph(); }
        KPartiteKClique<A>::KPartiteGraph* next_graph(){ return problem->next_graph(); }
        KPartiteKClique<A>* problem{nullptr};
};

template<>
inline KPartiteKClique<FindClique>::Vertex* KPartiteKClique<FindClique>::KPartiteGraph::last_vertex() { return NULL; }

template<>
inline void KPartiteKClique<FindClique>::KPartiteGraph::pop_last_vertex() {}

template<>
inline bool KPartiteKClique<FindClique>::KPartiteGraph::set_weights() { return false; }

template<>
inline bool KPartiteKClique<kpkc>::KPartiteGraph::set_part_sizes() { return false;}

template<>
inline int KPartiteKClique<kpkc>::KPartiteGraph::count(int, int) { return 0;}

template<>
inline int KPartiteKClique<kpkc>::KPartiteGraph::count(int) { return 0;}

template<>
inline int KPartiteKClique<kpkc>::KPartiteGraph::first(int) { return 0;}

template<>
inline void KPartiteKClique<kpkc>::KPartiteGraph::pop_vertex(int, int) {}

template<>
class KPartiteKClique<FindClique>:: Vertex {};

template<>
class KPartiteKClique<kpkc>::Vertex {
    // Vertex is a shallow copy of Vertex_template.

    friend bool operator<(const Vertex& l, const Vertex& r){
        // The lower the weight, the higher the obstruction when
        // selecting this vertex.
        // We want to select vertices with high obstruction first
        // and put the last (to be popped first).
        return l.weight > r.weight;
    }
    friend void intersection(Bitset& c, Vertex& l, Bitset& r){
        c.intersection_assign(*(l.bitset), r);
    }

    public:
        int index;  // The index in the original graph.
        int part;  // The part in the orginal graph.
        int weight;  // The higher, the higher the likelihood of a k-clique with this vertex.

        Vertex() = default;

        Vertex(Vertex_template& obj){
            bitset = &(obj.bitset);
            part = obj.part;
            weight = -1;
            index = obj.index;
            problem = obj.problem;

            if (1 != bitset->count(get_parts()[part], get_parts()[part + 1]))
                throw std::invalid_argument("the graph is not k-partite");
        }

        Vertex(const Vertex& obj){
            bitset = obj.bitset;
            weight = obj.weight;
            part = obj.part;
            index = obj.index;
            problem = obj.problem;
        }

        bool set_weight(){
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
        int intersection_count(Bitset& r, int start, int stop){
            return bitset->intersection_count(r, start, stop);
        }

        int intersection_count(Bitset& r, int part){
            return intersection_count(r, get_parts()[part], get_parts()[part+1]);
        }

    private:
        Bitset* bitset;
        KPartiteKClique* problem;
        const int* get_parts() { return problem->parts.data(); }
        const int get_k() { return problem->k; }
        Bitset& get_active_vertices() { return problem->current_graph()->active_vertices; }
};

// KPartiteGraph functions that are inlined

template<>
inline bool KPartiteKClique<kpkc>::KPartiteGraph::set_weights(){
    bool new_knowledge = false;
    for(Vertex& v: vertices)
        new_knowledge |= v.set_weight();
    return new_knowledge;
}

template<>
inline void KPartiteKClique<kpkc>::KPartiteGraph::pop_last_vertex(){
    Vertex& v = vertices.back();
    part_sizes[v.part] -= 1;
    active_vertices.unset(v.index);
    vertices.pop_back();
}

template<>
inline KPartiteKClique<kpkc>::Vertex* KPartiteKClique<kpkc>::KPartiteGraph::last_vertex(){
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

template<>
inline int KPartiteKClique<FindClique>::KPartiteGraph::count(int start, int stop){
    return active_vertices.count(start, stop);
}

template<>
inline int KPartiteKClique<FindClique>::KPartiteGraph::count(int part){
    return count(get_parts()[part], get_parts()[part+1]);
}

template<>
inline int KPartiteKClique<FindClique>::KPartiteGraph::first(int part){
    int the_first = active_vertices.first(get_parts()[part]);
    return (the_first < get_parts()[part+1]) ? the_first : -1;
}

template<>
inline void KPartiteKClique<FindClique>::KPartiteGraph::pop_vertex(int part, int vertex){
    active_vertices.unset(vertex);
    part_sizes[part] -= 1;
}

template<>
inline bool KPartiteKClique<FindClique>::KPartiteGraph::set_part_sizes(){
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

template<>
inline bool KPartiteKClique<kpkc>::KPartiteGraph::permits_another_choice(){
    for (int i=0; i<get_k(); i++){
        if (part_sizes[i] == 0)
            return false;
    }
    return true;
}

template<>
inline bool KPartiteKClique<FindClique>::KPartiteGraph::permits_another_choice(){
    return ((selected_part >= 0) && (part_sizes[selected_part]));
}

template<>
inline bool KPartiteKClique<kpkc>::KPartiteGraph::select(KPartiteGraph* next){
    /*
    Select the last (valid) vertex of the current graph set up the next graph
    to be all vertices connected to that last vertex.

    Return false, if there are no vertices left.
    */
    Vertex* v = last_vertex();
    if (!v)
        return false;

    // Copy the current sizes.
    for (int i=0; i<get_k(); i++)
        next->part_sizes[i] = part_sizes[i];

    // Select v.
    problem->_k_clique[v->part] = v->index;
    intersection(next->active_vertices, *v, active_vertices);

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

template<>
bool KPartiteKClique<kpkc>::KPartiteGraph::select(){
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

template<Algorithm A>
bool KPartiteKClique<A>::backtrack(){
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

}

#endif
