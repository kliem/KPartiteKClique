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
        Bitset(const Bitset& obj){
            // Not defined.
            assert(0);
        }
        ~Bitset();
        void unset(int index);
        bool has(int index);
        void set(int index);
        int intersection_count(Bitset& r, int start, int stop);
        int count(int start, int stop);
        int first(int start, int stop);
        void intersection_assign(Bitset& l, Bitset& r);
    private:
        uint64_t* data;
        int limbs;
        void allocate(int n_vertices);
        inline uint64_t& operator[](int i){ return data[i]; }
};

class KPartiteKClique {
    protected:
        class Vertex {
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

                Vertex();
                Vertex(const Vertex& obj);
                void init(KPartiteKClique* problem, const bool* incidences, int n_vertices, int part, int index);
                ~Vertex();
                bool set_weight();
                inline int intersection_count(Bitset& r, int start, int stop){
                    return bitset->intersection_count(r, start, stop);
                }
                inline int intersection_count(Bitset& r, int part){
                    return intersection_count(r, get_parts()[part], get_parts()[part+1]);
                }

            private:
                bool is_shallow;
                Bitset* bitset;
                KPartiteKClique* problem;
                inline const int* get_parts() { return problem->parts; }
                inline const int get_k() { return problem->k; }
                inline Bitset& get_active_vertices() { return *(problem->current_graph()).active_vertices; }
        };

        class KPartiteGraph {
            public:
                vector<Vertex> vertices;
                Bitset* active_vertices;
                int* part_sizes;
                int selected_part; // Needed for bitCLQ

                KPartiteGraph();
                KPartiteGraph(const KPartiteGraph& obj){
                    // Not defined.
                    assert(0);
                }
                void init(KPartiteKClique* problem, bool fill);
                ~KPartiteGraph();
                Vertex* last_vertex();
                void pop_last_vertex();
                bool is_valid();
                inline bool set_weights(){
                    bool new_knowledge = false;
                    for(Vertex& v: vertices)
                        new_knowledge |= v.set_weight();
                    return new_knowledge;
                }
                bool select(KPartiteGraph& next);
                bool select_bitCLQ(KPartiteGraph& next);

                // Used by bitCLQ.
                bool set_part_sizes();
                inline int count(int start, int stop){
                    return active_vertices->count(start, stop);
                }
                inline int count(int part){
                    return count(get_parts()[part], get_parts()[part+1]);
                }
                inline int first(int start, int stop){
                    return active_vertices->first(start, stop);
                }
                inline int first(int part){
                    return first(get_parts()[part], get_parts()[part+1]);
                }
                inline void pop_vertex(int part, int vertex){
                    active_vertices->unset(vertex);
                    part_sizes[part] -= 1;
                }
            private:
                inline const int* get_parts() { assert(problem); return problem->parts; }
                inline const int get_k() { assert(problem); return problem->k; }
                inline KPartiteGraph& current_graph(){ return problem->current_graph(); }
                inline KPartiteGraph& next_graph(){ return problem->next_graph(); }
                KPartiteKClique* problem;
        };

    public:
        const int* k_clique(){ return _k_clique; }
        KPartiteKClique(const bool* const* incidences, const int n_vertices, const int* first_per_part, const int k, const int prec_depth=5);
        KPartiteKClique(const KPartiteKClique& obj){
            // Not defined.
            assert(0);
        }
        KPartiteKClique();
        ~KPartiteKClique();
        bool next();
    protected:
        int* _k_clique;
        int* parts;
        int k;
        int current_depth;
        int n_vertices;
        int prec_depth;
        Vertex* all_vertices;
        KPartiteGraph* recursive_graphs;
        KPartiteGraph& current_graph(){ return recursive_graphs[current_depth]; }
        KPartiteGraph& next_graph(){ return recursive_graphs[current_depth + 1]; }
        bool traceback();
};

class bitCLQ : public KPartiteKClique {
    public:
        bool next();
};

#endif
