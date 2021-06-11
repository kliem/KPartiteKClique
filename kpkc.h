#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <cassert>
using namespace std;

class Bitset;

class Bitset {
    public:
        Bitset(int n_vertices, bool fill=false);
        Bitset(const bool* set_bits, int n_vertices);
        Bitset(const Bitset&) = delete;
        Bitset& operator=(const Bitset&) = delete;
        ~Bitset();
        void unset(int index);
        bool has(int index);
        void set(int index);
        int intersection_count(Bitset& r, int start, int stop);
        int count(int start, int stop);
        int first(int start);
        void intersection_assign(Bitset& l, Bitset& r);
    protected:
        uint64_t* data;
        int limbs;
        void allocate(int n_vertices);
        inline uint64_t& operator[](int i){ return data[i]; }
};

class KPartiteKClique_base {
    protected:
        class Vertex_template;
        class KPartiteGraph;
    public:
        KPartiteKClique_base(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k);
        virtual ~KPartiteKClique_base();
        const int* k_clique(){ return _k_clique; }
        bool next();

        KPartiteKClique_base(const KPartiteKClique_base&) = delete;
        KPartiteKClique_base& operator=(const KPartiteKClique_base&) = delete;
    protected:
        int* _k_clique;
        int* parts;
        int k;
        int current_depth;
        int n_vertices;
        Vertex_template* all_vertices;
        virtual KPartiteGraph* current_graph();
        virtual KPartiteGraph* next_graph();
        bool backtrack();
        void swap(Vertex_template& a, Vertex_template& b);
        void swap(KPartiteGraph& a, KPartiteGraph& b);
};

class KPartiteKClique : public KPartiteKClique_base {
    friend KPartiteKClique_base::Vertex_template;
    class Vertex;
    class KPartiteGraph;

    public:
        KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k, const int prec_depth=5);
        ~KPartiteKClique();

        KPartiteKClique(const KPartiteKClique&) = delete;
        KPartiteKClique& operator=(const KPartiteKClique&) = delete;
    private:
        using KPartiteKClique_base::_k_clique;
        using KPartiteKClique_base::all_vertices;
        int prec_depth;
        KPartiteGraph* recursive_graphs;
        KPartiteKClique_base::KPartiteGraph* current_graph();
        KPartiteKClique_base::KPartiteGraph* next_graph();
        void swap(KPartiteGraph& a, KPartiteGraph& b);
};

class FindClique : public KPartiteKClique_base {
    class KPartiteGraph;

    public:
        FindClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k);
        ~FindClique();

        FindClique(const FindClique&) = delete;
        FindClique& operator=(const FindClique&) = delete;
    private:
        KPartiteGraph* recursive_graphs;
        KPartiteKClique_base::KPartiteGraph* current_graph();
        KPartiteKClique_base::KPartiteGraph* next_graph();
        void swap(KPartiteGraph& a, KPartiteGraph& b);
};

class KPartiteKClique_base::Vertex_template {
    // Takes care of the memory allocation for vertices.

    friend KPartiteKClique;
    friend KPartiteKClique::Vertex;
    friend FindClique;
    friend void KPartiteKClique_base::swap(Vertex_template& a, Vertex_template& b);

    inline friend void intersection(Bitset& c, Vertex_template& l, Bitset& r){
        c.intersection_assign(*(l.bitset), r);
    }

    public:
        int index;  // The index in the original graph.
        int part;  // The part in the orginal graph.

        Vertex_template() { bitset = NULL; }
        Vertex_template(const Vertex_template&) = delete;
        Vertex_template(KPartiteKClique_base* problem, const bool* incidences, int n_vertices, int part, int index);
        ~Vertex_template();

    protected:
        Bitset* bitset;
        KPartiteKClique_base* problem;
};

class KPartiteKClique_base::KPartiteGraph {
    friend void KPartiteKClique_base::swap(KPartiteGraph& a, KPartiteGraph& b);
    public:
        Bitset* active_vertices;
        int* part_sizes;
        KPartiteGraph();
        KPartiteGraph(KPartiteKClique_base* problem, bool fill);
        virtual ~KPartiteGraph();

        virtual bool permits_another_choice();
        virtual bool select(KPartiteGraph* next);
        virtual bool select();

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
    protected:
        inline const int* get_parts() { assert(problem); return problem->parts; }
        inline const int get_k() { assert(problem); return problem->k; }
        KPartiteKClique_base::KPartiteGraph* current_graph(){ return problem->current_graph(); }
        KPartiteKClique_base::KPartiteGraph* next_graph(){ return problem->next_graph(); }
        KPartiteKClique_base* problem;
};

class KPartiteKClique::Vertex {
    // Vertex is a shallow copy of Vertex_template.

    inline friend bool operator<(const Vertex& l, const Vertex& r){
        // The lower the weight, the higher the obstruction when
        // selecting this vertex.
        // We want to select vertices with high obstruction first
        // and put the last (to be popped first).
        return l.weight > r.weight;
    }
    inline friend void intersection(Bitset& c, Vertex& l, Bitset& r){
        c.intersection_assign(*(l.bitset), r);
    }

    public:
        int index;  // The index in the original graph.
        int part;  // The part in the orginal graph.
        int weight;  // The higher, the higher the likelihood of a k-clique with this vertex.

        Vertex() {}
        Vertex(const KPartiteKClique_base::Vertex_template& obj);
        Vertex(const Vertex& obj);
        bool set_weight();
        inline int intersection_count(Bitset& r, int start, int stop){
            return bitset->intersection_count(r, start, stop);
        }
        inline int intersection_count(Bitset& r, int part){
            return intersection_count(r, get_parts()[part], get_parts()[part+1]);
        }

    private:
        Bitset* bitset;
        KPartiteKClique* problem;
        inline const int* get_parts() { return problem->parts; }
        inline const int get_k() { return problem->k; }
        inline Bitset& get_active_vertices() { return *(problem->current_graph()->active_vertices); }
};

class KPartiteKClique::KPartiteGraph : KPartiteKClique_base::KPartiteGraph {
    friend void KPartiteKClique::swap(KPartiteGraph& a, KPartiteGraph& b);
    public:
        vector<Vertex> vertices;
        KPartiteGraph();
        KPartiteGraph(KPartiteKClique* problem, bool fill);
        ~KPartiteGraph() {}

        bool permits_another_choice();
        bool select(KPartiteKClique_base::KPartiteGraph* next2);
        bool select();

        Vertex* last_vertex();
        void pop_last_vertex();
        inline bool set_weights(){
            bool new_knowledge = false;
            for(Vertex& v: vertices)
                new_knowledge |= v.set_weight();
            return new_knowledge;
        }

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
    protected:
        KPartiteKClique* problem;
};

class FindClique::KPartiteGraph : KPartiteKClique_base::KPartiteGraph {
    friend void FindClique::swap(KPartiteGraph& a, KPartiteGraph& b);
    public:
        using KPartiteKClique_base::KPartiteGraph::active_vertices;
        int selected_part;
        KPartiteGraph() : KPartiteKClique_base::KPartiteGraph() {}
        KPartiteGraph(FindClique* problem, bool fill);
        ~KPartiteGraph() {}

        bool permits_another_choice();
        bool select(KPartiteKClique_base::KPartiteGraph* next2);
        bool select();

        bool set_part_sizes();
        inline int count(int start, int stop){
            return active_vertices->count(start, stop);
        }
        inline int count(int part){
            return count(get_parts()[part], get_parts()[part+1]);
        }
        inline int first(int part){
            int the_first = active_vertices->first(get_parts()[part]);
            return (the_first < get_parts()[part+1]) ? the_first : -1;
        }
        inline void pop_vertex(int part, int vertex){
            active_vertices->unset(vertex);
            part_sizes[part] -= 1;
        }

        KPartiteGraph(const KPartiteGraph&) = delete;
        KPartiteGraph& operator=(const KPartiteGraph&) = delete;
    protected:
        FindClique* problem;
};

#endif
