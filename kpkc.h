#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <cassert>
#include <stdexcept>

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

class KPartiteKClique_base {
    protected:
        class Vertex_template;
        class KPartiteGraph;
    public:
        KPartiteKClique_base(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k);
        virtual ~KPartiteKClique_base() = default;
        const int* k_clique(){ return _k_clique.data(); }
        bool next();

        KPartiteKClique_base(const KPartiteKClique_base&) = delete;
        KPartiteKClique_base& operator=(const KPartiteKClique_base&) = delete;
    protected:
        std::vector<int> _k_clique;
        std::vector<int> parts;
        int k;
        int current_depth;
        int n_vertices;
        std::vector<Vertex_template> all_vertices;
        virtual KPartiteGraph* current_graph() = 0;
        virtual KPartiteGraph* next_graph() = 0;
        bool backtrack();
};

class KPartiteKClique : public KPartiteKClique_base {
    friend KPartiteKClique_base::Vertex_template;
    class Vertex;
    class KPartiteGraph;

    public:
        KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k, const int prec_depth=5);
        KPartiteKClique(const KPartiteKClique&) = delete;
        KPartiteKClique& operator=(const KPartiteKClique&) = delete;
    private:
        using KPartiteKClique_base::_k_clique;
        using KPartiteKClique_base::all_vertices;
        int prec_depth;
        std::vector<KPartiteGraph> recursive_graphs;
        KPartiteKClique_base::KPartiteGraph* current_graph() override;
        KPartiteKClique_base::KPartiteGraph* next_graph() override;
};

class FindClique : public KPartiteKClique_base {
    class KPartiteGraph;

    public:
        FindClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k);

        FindClique(const FindClique&) = delete;
        FindClique& operator=(const FindClique&) = delete;
    private:
        std::vector<KPartiteGraph> recursive_graphs;
        KPartiteKClique_base::KPartiteGraph* current_graph() override;
        KPartiteKClique_base::KPartiteGraph* next_graph() override;
};

class KPartiteKClique_base::Vertex_template {
    // Takes care of the memory allocation for vertices.

    friend KPartiteKClique;
    friend KPartiteKClique::Vertex;
    friend FindClique;

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
        Vertex_template(KPartiteKClique_base* problem, const bool* incidences, int n_vertices, int part, int index);

    protected:
        Bitset bitset;
        KPartiteKClique_base* problem;
};

class KPartiteKClique_base::KPartiteGraph {
    public:
        Bitset active_vertices;
        std::vector<int> part_sizes;
        KPartiteGraph() = default;
        KPartiteGraph(KPartiteKClique_base* problem, bool fill);

        virtual bool permits_another_choice() = 0;
        virtual bool select(KPartiteGraph* next) = 0;
        virtual bool select() = 0;

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
        KPartiteGraph(KPartiteGraph&&) = default;
        KPartiteGraph& operator=(KPartiteGraph&&) = default;
    protected:
        const int* get_parts() { assert(problem); return problem->parts.data(); }
        const int get_k() { assert(problem); return problem->k; }
        KPartiteKClique_base::KPartiteGraph* current_graph(){ return problem->current_graph(); }
        KPartiteKClique_base::KPartiteGraph* next_graph(){ return problem->next_graph(); }
        KPartiteKClique_base* problem{nullptr};
};

class KPartiteKClique::Vertex {
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

        Vertex(KPartiteKClique_base::Vertex_template& obj){
            bitset = &(obj.bitset);
            part = obj.part;
            weight = -1;
            index = obj.index;
            problem = (KPartiteKClique*) obj.problem;

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

class KPartiteKClique::KPartiteGraph : KPartiteKClique_base::KPartiteGraph {
    public:
        std::vector<Vertex> vertices;
        KPartiteGraph();
        KPartiteGraph(KPartiteKClique* problem, bool fill);

        bool permits_another_choice() override;
        bool select(KPartiteKClique_base::KPartiteGraph* next2) override;
        bool select() override;

        Vertex* last_vertex();
        void pop_last_vertex(){
            Vertex& v = vertices.back();
            part_sizes[v.part] -= 1;
            active_vertices.unset(v.index);
            vertices.pop_back();
        }
        bool set_weights(){
            bool new_knowledge = false;
            for(Vertex& v: vertices)
                new_knowledge |= v.set_weight();
            return new_knowledge;
        }

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
        KPartiteGraph(KPartiteGraph&&) = default;
        KPartiteGraph& operator=(KPartiteGraph&&) = default;
    protected:
        KPartiteKClique* problem;
};

class FindClique::KPartiteGraph : KPartiteKClique_base::KPartiteGraph {
    public:
        using KPartiteKClique_base::KPartiteGraph::active_vertices;
        int selected_part;
        KPartiteGraph() : KPartiteKClique_base::KPartiteGraph() {}
        KPartiteGraph(FindClique* problem, bool fill);

        bool permits_another_choice() override;
        bool select(KPartiteKClique_base::KPartiteGraph* next2) override;
        bool select() override;

        bool set_part_sizes();
        int count(int start, int stop){
            return active_vertices.count(start, stop);
        }
        int count(int part){
            return count(get_parts()[part], get_parts()[part+1]);
        }
        int first(int part){
            int the_first = active_vertices.first(get_parts()[part]);
            return (the_first < get_parts()[part+1]) ? the_first : -1;
        }
        void pop_vertex(int part, int vertex){
            active_vertices.unset(vertex);
            part_sizes[part] -= 1;
        }

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
        KPartiteGraph(KPartiteGraph&&) = default;
        KPartiteGraph& operator=(KPartiteGraph&&) = default;
    protected:
        FindClique* problem;
};

inline bool KPartiteKClique_base::backtrack(){
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
