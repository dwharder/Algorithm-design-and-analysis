# Graph representations: memory and access patterns

This repository provides a collection of graph classes that illustrate how
edges may be stored in memory and how they may be accessed in C++.

The focus is not on graph algorithms, but rather on:

- how different graph representations organize edge data in memory,
- how edges are accessed and modified, and
- the trade-offs between these representations.

All graphs:

- assume a fixed number of vertices N,
- use zero-based indexing, and
- store nonnegative floating-point weights, where
  - w_{i,j} > 0 indicates the presence of an edge,
  - w_{i,j} = 0 indicates that no edge exists.

Edge access is uniformly provided through:

    g.edge(i, j)

which:

- returns the weight of the edge (v_i, v_j),
- returns 0 if no edge exists, and
- in the non-const case, returns a proxy that validates assignments.

---

# Overview of representations

The following classes are included:

- Directed_graph_matrix
- Undirected_graph_calculated
- Undirected_graph_row_pointers
- Directed_graph_list
- Undirected_graph_list
- Directed_graph_hashed
- Undirected_graph_hashed

Each class demonstrates a different way of storing edges.

---

# Dense representations

## Directed_graph_matrix

A dense directed graph is stored using a contiguous N x N array.

- Each edge (i, j) is stored at index i*N + j
- The diagonal entries are unused (loops are not allowed)
- Memory usage:
  Θ(|V|^2)

Access:

- Direct indexing provides constant-time access
- No search is required

This representation is simple and efficient for dense graphs, but wastes memory
for sparse graphs.

---

## Undirected_graph_calculated

An undirected graph is stored using only the lower-triangular portion
of the adjacency matrix.

- Each unordered pair {i, j} is mapped to a unique index
- The index is computed at each access
- Memory usage:
  Θ(|V|^2 / 2)

Access:

- Requires a small amount of arithmetic to compute the index
- Still constant time

This reduces memory compared to the full matrix but introduces a small
computational overhead.

---

## Undirected_graph_row_pointers

This representation improves on the calculated version by precomputing row offsets.

- The lower-triangular data is stored contiguously
- An additional array of pointers stores the start of each row
- Memory usage:
  Θ(|V|^2 / 2)

Access:

- No index calculation is required
- Direct pointer indexing is used

This trades a small amount of additional memory for slightly faster access.

---

# Sparse representations

## Directed_graph_list

A sparse directed graph is stored as an array of adjacency lists.

- Each vertex stores a vector of pairs (j, w_{i,j})
- Memory usage:
  Θ(|V| + |E|)

Access:

- Finding an edge requires a linear search over outgoing edges
- Efficient when the number of edges per vertex is small

This representation is well suited for sparse graphs where iteration over
neighbors is common.

---

## Undirected_graph_list

An undirected graph is stored using adjacency lists with a canonical ordering.

- Each edge {i, j} is stored once
- Indices are normalized to ensure consistent storage
- Memory usage:
  Θ(|V| + |E|)

Access:

- Linear search within a small list
- Efficient when vertex degrees are small

This avoids duplication of edges while maintaining compact storage.

---

## Directed_graph_hashed

A sparse directed graph is stored using a hash table.

- Keys are ordered pairs (i, j)
- Values are edge weights
- Memory usage:
  Θ(|E|)

Access:

- Expected constant-time lookup, insertion, and removal

This representation is useful when fast edge lookup is required, but it is less
efficient for iterating over all neighbors of a vertex.

---

## Undirected_graph_hashed

An undirected graph is stored using a hash table with canonicalized keys.

- Each edge {i, j} is stored once in a normalized form
- Memory usage:
  Θ(|E|)

Access:

- Expected constant-time lookup
- Indices are normalized before access

This provides fast edge queries while avoiding duplication.

---

# Proxy-based access

All classes use a proxy object for non-const access:

    g.edge(i, j) = value;

The proxy ensures that:

- assigned weights are finite and nonnegative,
- assigning 0 removes the edge,
- invalid assignments raise an exception.

This allows all graph representations to share a consistent interface while
enforcing correctness.

---

# Summary

These classes demonstrate the trade-offs between:

- memory usage
- access time
- iteration efficiency

| Representation | Memory                           | Lookup          | Iteration      |
|----------------|----------------------------------|-----------------|----------------|
| Matrix         | Θ(&#124;V&#124;^2)               | Θ(1)            | inefficient    |
| List           | Θ(&#124;V&#124; + &#124;E&#124;) | Θ(deg(i))       | efficient      |
| Hashed         | Θ(&#124;E&#124;)                 | Θ(1) (expected) | less efficient |

No single representation is optimal in all cases; the choice depends on the
structure of the graph and the operations being performed.
