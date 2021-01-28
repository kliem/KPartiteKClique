#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;
#define DBG 0
#define MEM_DBG 0

class Bitset;

class Bitset {

    public:
        Bitset(int n_vertices, int k, int *first_per_part, bool fill=false);
        Bitset(const bool* set_bits, int n_vertices, int k, int *first_per_part);
        Bitset(const Bitset& obj);
        ~Bitset();
        int first(int start=0);  // curently neither used nor defined
        void unset(int index, int part);
        bool has(int index, int part);
        void set(int index, int part);
        int intersection_count(Bitset& r, int part);
        bool intersection_is_empty(Bitset& r, int part);
        bool is_empty(int part);  // currently not used
        void intersection_assign(Bitset& l, Bitset& r);
    private:
        uint64_t* data;
        int limbs;
        int* first_limb_per_part;
        int* offset;
        void allocate(int n_vertices, int k, int* first_per_part);
        uint64_t& operator[](int i){ return data[i]; }
        bool shallow;
};

class KPartiteKClique {
        class Vertex {
            friend bool operator<(const Vertex& l, const Vertex& r){
                // The lower the weight, the higher the obstruction when
                // selecting this vertex.
                // We want to select vertices with high obstruction first
                // and put the last (to be popped first).
                return l.weight > r.weight;
            }
            public:
                int index;  // The index in the original graph.
                int part;  // The part in the orginal graph.
                Vertex(){
                    is_shallow = true;
                }
                void init(KPartiteKClique* problem, const bool* incidences, int n_vertices, int part, int index){
                    bitset = new Bitset(incidences, n_vertices, problem->k, problem->parts);
                    is_shallow = false;
                    weight = -1;
                    this->part = part;
                    this->index = index;
                    this->problem = problem;

                    // Set each vertex adjacent to itself.
                    // This is important, so that after selecting a vertex
                    // the corresponding part will have one ``active_vertex``.
                    bitset->set(index, part);
                };
                Vertex(const Vertex& obj){
                    // Make a shallow copy.
                    bitset = obj.bitset;
                    is_shallow = true;
                    weight = obj.weight;
                    part = obj.part;
                    index = obj.index;
                    problem = obj.problem;
                }
                ~Vertex(){
                    if (!is_shallow){
#if MEM_DBG
                        cout << "deleteing a vertex" << (size_t) bitset << endl;
#endif
                        delete bitset;
                    }
                }
                bool set_weight();
                bool is_valid();
                int weight;
                void intersection(Bitset& c, Bitset& r){
                    // c = this & r.
                    c.intersection_assign(*bitset, r);
                }
                inline int intersection_count(Bitset& r, int part){
                    return bitset->intersection_count(r, part);
                }
                inline bool intersection_is_empty(Bitset& r, int part){
                    return bitset->intersection_is_empty(r, part);
                }

            private:
                bool is_shallow;
                Bitset* bitset;
                KPartiteKClique* problem;
                const int* get_parts() { return problem->parts; }
                const int get_k() { return problem->k; }
                Bitset& get_active_vertices() { return *(problem->current_graph()).active_vertices; }
        };

    class KPartiteGraph {
        public:
            vector<Vertex> vertices;
            Bitset* active_vertices;

            KPartiteGraph();
            void init(KPartiteKClique* problem, bool fill);
            ~KPartiteGraph();
            Vertex* last_vertex();
            void pop_last_vertex();
            bool is_valid();
            inline bool set_weights(){
                bool new_knowledge = false;
                for(Vertex& v: vertices){
#if DBG
                    cout << "set weight of " << v->index << endl;
#endif
                    new_knowledge |= v.set_weight();
#if DBG
                    cout << "weight is " << v->weight << endl;
#endif
                }
                return new_knowledge;
            }
            bool select(KPartiteGraph& next);
            inline int count_active_vertices(int part){
                return active_vertices->intersection_count(*active_vertices, part);
            }
        private:
            const int* get_parts() { assert(problem); return problem->parts; }
            const int get_k() { assert(problem); return problem->k; }
            KPartiteGraph& current_graph(){ return problem->current_graph(); }
            KPartiteGraph& next_graph(){ return problem->next_graph(); }
            int* part_sizes;
            KPartiteKClique* problem;
    };

    public:
        const int* k_clique(){ return _k_clique; }
        KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k);
        KPartiteKClique();
        ~KPartiteKClique();
        bool next();
    private:
        int* _k_clique;
        int* parts;
        int k;
        int current_depth;
        int n_vertices;
        Vertex* all_vertices;
        KPartiteGraph* recursive_graphs;
        KPartiteGraph& current_graph(){ return recursive_graphs[current_depth]; }
        KPartiteGraph& next_graph(){ return recursive_graphs[current_depth + 1]; }
        bool traceback();
};

#endif
