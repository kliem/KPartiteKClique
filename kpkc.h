#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <cassert>

// namespace k-partite-k-clicque
namespace kpkc
{
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
        void unset(int index);
        bool has(int index);
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

        Vertex() {}
        Vertex(KPartiteKClique_base::Vertex_template& obj);
        Vertex(const Vertex& obj);
        bool set_weight();
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
        void pop_last_vertex();
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
}

#endif
