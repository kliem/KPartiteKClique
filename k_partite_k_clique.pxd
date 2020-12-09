# distutils: depends = k_partite_k_clique.cpp k_partite_k_clique.h
# distutils: extra_compile_args=-O3 -march=native -std=c++11
# distutils: language = c++

from libcpp cimport bool

cdef extern from "k_partite_k_clique.h":
    cdef cppclass KPartiteKClique:
        KPartiteKClique(bool **, int n_vertices, int* first_per_part, int k)
        KPartiteKClique()
        bool next()
        const int* k_clique()
