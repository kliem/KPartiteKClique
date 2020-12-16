from sage.all import *

def get_random_k_partite_graph(n_parts, min_part_size, max_part_size, dens1, dens2):
    r"""
    Randomly generated k-partite graphs according to Grünert et al.
    """
    B = list(range(n_parts))
    # Construct nodes
    offset = 0
    sizes = [randint(min_part_size, max_part_size) for _ in range(n_parts)]
    parts = []
    for b in B:
        parts.append(list(range(offset, offset+sizes[b])))
        offset += sizes[b]

    vertex_weight_generators = [dens1 + random()*(dens2 - dens1) for _ in range(offset)]

    G = Graph()
    G.add_vertices(range(sum(sizes)))

    for b in range(n_parts):
        for i in parts[b]:
            for b1 in range(b+1, n_parts):
                for i1 in parts[b1]:
                    if random() < (vertex_weight_generators[i] + vertex_weight_generators[i1])/2:
                        G.add_edge([i,i1])

    return G, parts
