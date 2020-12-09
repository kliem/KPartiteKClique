#ifndef KPartiteKClique_header
#define KPartiteKClique_header

#include <cstdint>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;

class Bitset;

class Bitset {

    public:
        Bitset(int n_vertices){
            Bitset(n_vertices, false);
        }
        Bitset(int n_vertices, bool fill);
        Bitset(const bool* set_bits, int n_vertices);
        ~Bitset();
        Bitset(const Bitset& obj);
        int first(int start=0);
        void unset(int index);
        void set(int index);
        int intersection_count(Bitset& r, int start, int stop);
        void intersection_assign(Bitset& l, Bitset& r);
    private:
        uint64_t* data;
        int limbs;
        void allocate(int n_vertices);
        uint64_t& operator[](int i){ return data[i]; }
        bool shallow;
};

class KPartiteKClique {
    class KPartiteGraph {
        public:
            class Vertex {
                friend bool operator<(const Vertex& l, const Vertex& r){
                    return l.weight > r.weight;
                }
                public:
                    int index;
                    int part;
                    Vertex(){
                        is_shallow = true;
                    }
                    void init(KPartiteGraph* graph, const bool* incidences, int n_vertices, int part, int index){
                        bitset = new Bitset(incidences, n_vertices);
                        is_shallow = false;
                        weight = -1;
                        this->part = part;
                        this->index = index;
                        this->graph = graph;

                        // Set each vertex adjacent to itself.
                        // This is important, so that after selecting a vertex
                        // the corresponding part will have one ``active_vertex``.
                        bitset->set(index);
                    };
                    Vertex(const Vertex& obj){
                        bitset = obj.bitset;
                        is_shallow = true;
                        weight = -1;
                        part = obj.part;
                        index = obj.index;
                        graph = obj.graph;
                    }
                    ~Vertex(){
                        if (!is_shallow){
                            cout << "deleteing a vertex" << (size_t) bitset << endl;
                            delete bitset;
                        }
                    }
                    void set_weight();
                    int weight;
                    void intersection(Bitset& c, Bitset& r){
                        // c = this & r.
                        c.intersection_assign(bitset[0], r);
                    }
                    int intersection_count(Bitset& r, int start, int stop){
                        return bitset[0].intersection_count(r, start, stop);
                    }
                    int intersection_count(Bitset& r, int part){
                        return intersection_count(r, get_parts()[part], get_parts()[part+1]);
                    }

                private:
                    const int* get_parts() { return graph[0].get_parts(); }
                    const int get_k() { return graph[0].get_k(); }
                    Bitset& get_active_vertices() { return graph[0].active_vertices[0]; }
                    bool is_shallow;
                    Bitset* bitset;
                    KPartiteGraph* graph;
            };
            vector<Vertex> vertices;
            KPartiteGraph();
            void init(KPartiteKClique* problem, bool fill);
            ~KPartiteGraph();
            Vertex* last_vertex();
            void pop_last_vertex();
            bool is_valid();
            void set_weights(){
                for(auto v : vertices)
                    v.set_weight();
            }
            bool select(KPartiteGraph& next);
        private:
            const int* get_parts() { assert(problem); return problem->parts; }
            const int get_k() { assert(problem); return problem->k; }
            Bitset* active_vertices;
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
        KPartiteGraph::Vertex* all_vertices;
        KPartiteGraph* recursive_graphs;
        KPartiteGraph& current_graph(){ return recursive_graphs[current_depth]; }
        KPartiteGraph& next_graph(){ return recursive_graphs[current_depth + 1]; }
        bool traceback();
};

#endif
