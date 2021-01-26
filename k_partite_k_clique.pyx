# distutils: depends = k_partite_k_clique.cpp k_partite_k_clique.h
# distutils: sources = k_partite_k_clique.cpp
# distutils: extra_compile_args=-O3 -march=native -std=c++11
# distutils: language = c++

from libcpp cimport bool
from cysignals.signals      cimport sig_on, sig_off

cdef extern from "k_partite_k_clique.h":
    cdef cppclass KPartiteKClique:
        KPartiteKClique(bool **, int n_vertices, int* first_per_part, int k)
        KPartiteKClique()
        bool next()
        const int* k_clique()

from sage.ext.memory_allocator cimport MemoryAllocator
from libcpp cimport bool

def KPartiteKClique_iter(G, parts, benchmark=False):
    """
    Iterates over all k-cliques
    """
    cdef int n = G.order()
    cdef int i, j

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef bool ** incidences = <bool **> mem.allocarray(n, sizeof(bool *))
    for i in range(n):
        incidences[i] = <bool *> mem.calloc(n, sizeof(bool))

    assert isinstance(parts, list)  # We will allow more flexibility later.

    cdef int k = len(parts)
    cdef int* first_per_part = <int*> mem.allocarray(k, sizeof(int))

    cdef int counter = 0
    for i in range(k):
        first_per_part[i] = counter
        counter += len(parts[i])

    def id_to_part(index):
        for i in range(k):
            if first_per_part[i] > index:
                break
        else:
            i += 1
        i -= 1
        return i

    def id_to_vertex(index):
        i = id_to_part(index)
        return parts[i][index - first_per_part[i]]

    cdef dict vertex_to_id = {id_to_vertex(i):i for i in range(n)}

    cdef int ui, vi

    for u in G:
        ui = vertex_to_id[u]
        for v in G.neighbors(u):
            vi = vertex_to_id[v]
            if id_to_part(ui) == id_to_part(vi):
                raise ValueError("not a k-partite graph")
            incidences[ui][vi] = True
            incidences[vi][ui] = True

    if benchmark:
        yield []

    cdef KPartiteKClique * K = new KPartiteKClique(incidences, n, first_per_part, k)

    try:
        sig_on()
        foo = K.next()
        sig_off()
        while foo:
            yield [id_to_vertex(K.k_clique()[i]) for i in range(k)]
            sig_on()
            foo = K.next()
            sig_off()
    finally:
        del K
