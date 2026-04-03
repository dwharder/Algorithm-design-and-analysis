// Directed_graph_matrix
//     Directed_graph_matrix implements a dense directed graph without loops using
//     a contiguous N × N adjacency matrix. Each edge (v_i, v_j) is stored explicitly,
//     and edge(i, j) provides access to its weight. The class enforces finite,
//     nonnegative weights via a proxy object, with a weight of 0 indicating that no
//     edge exists. This representation offers constant-time access to any edge but
//     requires Theta(|V|^2) memory, making it most suitable for dense graphs.
//
// Undirected_graph_calculated
//     Undirected_graph_calculated implements a dense undirected graph without loops
//     using a compact contiguous array that stores only the lower-triangular portion
//     of the adjacency matrix. Indices are computed on demand to map each unordered
//     pair {i, j} to a unique storage location. Edge access is provided through
//     edge(i, j), and assignments are validated via a proxy to ensure finite,
//     nonnegative weights, with 0 indicating no edge. This representation reduces
//     storage by roughly half compared to a full matrix while retaining constant-time access.
//
// Undirected_graph_row_pointers
//     Undirected_graph_row_pointers implements a dense undirected graph without loops
//     using a compact lower-triangular storage scheme augmented with precomputed row
//     pointers. These pointers allow direct indexing into the contiguous array,
//     reducing arithmetic overhead compared to computed indexing. The class provides
//     edge(i, j) access with proxy-based assignment enforcing finite, nonnegative
//     weights, where 0 indicates no edge. This design trades a small amount of extra
//     memory for slightly faster access compared to the calculated-index variant.
//
// Directed_graph_list
//     Directed_graph_list implements a sparse directed graph using adjacency lists,
//     where each vertex stores a vector of outgoing edges (j, w_{i,j}). Edge access
//     is provided through edge(i, j), which performs a linear search over the outgoing
//     edges of vertex i. Assignments are handled via a proxy that enforces finite,
//     nonnegative weights, with 0 indicating that an edge should be removed. This
//     representation requires Theta(|V| + |E|) memory and is well suited for sparse graphs
//     where iteration over neighbors is common.
//
// Undirected_graph_list
//     Undirected_graph_list implements a sparse undirected graph using adjacency
//     lists, storing each edge {i, j} exactly once in a canonical order. Access via
//     edge(i, j) normalizes indices to ensure consistent lookup and modification.
//     Assignments are validated through a proxy, enforcing finite, nonnegative weights,
//     with 0 indicating removal of the edge. This representation requires
//     Theta(|V| + |E|) memory and is efficient for sparse graphs where traversal of
//     neighboring vertices is the dominant operation.
//
//
// Directed_graph_hashed
//     Directed_graph_hashed implements a sparse directed graph using a hash table
//     that maps ordered vertex pairs (i, j) to edge weights. This allows expected
//     constant-time lookup, insertion, and removal of edges. Access is provided via
//     edge(i, j), with proxy-based assignment enforcing finite, nonnegative weights,
//     where 0 removes the edge. This representation is suitable when fast edge
//     existence queries are required, though it may be less efficient for iterating
//     over all neighbors of a vertex.
//
//
// Undirected_graph_hashed
//     Undirected_graph_hashed implements a sparse undirected graph using a hash table
//     that maps unordered vertex pairs {i, j} (stored in canonical form) to edge
//     weights. Access through edge(i, j) ensures consistent normalization of indices,
//     and assignments are validated via a proxy enforcing finite, nonnegative weights,
//     with 0 removing the edge. This representation provides expected constant-time
//     edge lookup and modification, making it suitable for applications requiring
//     fast access to individual edges in sparse graphs.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

  ////////////////////////////////////////////////
 // A directed graph using an adjacency matrix //
////////////////////////////////////////////////
template <typename T, std::size_t N>
class Directed_graph_matrix {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Directed_graph_matrix;

        entry( Directed_graph_matrix &graph, size_type i, size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Directed_graph_matrix *p_graph_;
        size_type i_;
        size_type j_;
    };

    [[nodiscard]] T edge( size_type i, size_type j ) const noexcept {
        assert( (i < N) && (j < N) && (i != j) );
        return edges_[i*N + j];
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) noexcept {
        assert( (i < N) && (j < N) && (i != j) );
        return entry{ *this, i, j };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative."
            };
        }

        edges_[i*N + j] = value;
    }

    std::array<T, N*N> edges_{};
};

   ///////////////////////////////////////////////////
  // An undirected graph using an adjacency matrix //
 // where entries are calculated.                 //
///////////////////////////////////////////////////

template <typename T, std::size_t N>
class Undirected_graph_calculated {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Undirected_graph_calculated;

        entry( Undirected_graph_calculated &graph,
               size_type i,
               size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Undirected_graph_calculated *p_graph_;
        size_type i_;
        size_type j_;
    };

    [[nodiscard]] T edge( size_type i, size_type j ) const noexcept {
        assert( (i < N) && (j < N) && (i != j) );

        size_type const index{
            (i > j) ? (i*(i - 1))/2 + j
                    : (j*(j - 1))/2 + i
        };

        return edges_[index];
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) noexcept {
        assert( (i < N) && (j < N) && (i != j) );
        return (j < i) ? entry{ *this, i, j } : entry{ *this, j, i };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );
        assert( i > j );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative."
            };
        }

        size_type const index{
            (i > j) ? (i*(i - 1))/2 + j
                    : (j*(j - 1))/2 + i
        };

        edges_[index] = value;
    }

    std::array<T, (N*(N - 1))/2> edges_{};
};

   /////////////////////////////////////////////////////////
  // An undirected graph using an adjacency matrix where //
 //  row ponters are determined in the constructor.     //
/////////////////////////////////////////////////////////

template <typename T, std::size_t N>
class Undirected_graph_row_pointers {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Undirected_graph_row_pointers;

        entry( Undirected_graph_row_pointers &graph,
               size_type i,
               size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Undirected_graph_row_pointers *p_graph_;
        size_type i_;
        size_type j_;
    };

    Undirected_graph_row_pointers() {
        for ( size_type i{ 0 }; i < N; ++i ) {
            edge_rows_[i] = edges_.data() + (i*(i - 1))/2;
        }
    }

    [[nodiscard]] T edge( size_type i, size_type j ) const noexcept {
        assert( (i < N) && (j < N) && (i != j) );

        return (i > j)
            ? edge_rows_[i][j]
            : edge_rows_[j][i];
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) noexcept {
        assert( (i < N) && (j < N) && (i != j) );
        return (j < i) ? entry{ *this, i, j } : entry{ *this, j, i };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );
        assert( i > j );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative."
            };
        }

        if ( i > j ) {
            edge_rows_[i][j] = value;
        } else {
            edge_rows_[j][i] = value;
        }
    }

    std::array<T*, N> edge_rows_{};
    std::array<T, (N*(N - 1))/2> edges_{};
};

  ///////////////////////////////////////////////
 // A directed graph using an adjacency list. //
///////////////////////////////////////////////

template <typename T, std::size_t N>
class Directed_graph_list {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;
    using edge_type = std::pair<size_type, T>;

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Directed_graph_list;

        entry( Directed_graph_list &graph, size_type i, size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Directed_graph_list *p_graph_;
        size_type i_;
        size_type j_;
    };

    Directed_graph_list() = default;

    // The call edge( i, j ) returns the weight of the edge (v_i, v_j).
    [[nodiscard]] T edge( size_type i, size_type j ) const {
        assert( (i < N) && (j < N) );
        auto const &row{ edges_[i] };

        auto itr{ std::find_if(
            row.begin(),
            row.end(),
            [j]( edge_type const &e ) {
                return e.first == j;
            }
        ) };

        return (itr == row.end()) ? T{ 0 } : itr->second;
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) {
        assert( (i < N) && (j < N) );
        return entry{ *this, i, j };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative; a weight of 0 removes the edge."
            };
        }

        auto &row{ edges_[i] };

        auto itr{ std::find_if(
            row.begin(),
            row.end(),
            [j]( edge_type const &e ) {
                return e.first == j;
            }
        ) };

        if ( value == T{ 0 } ) {
            if ( itr != row.end() ) {
                row.erase( itr );
            }

            return;
        }

        if ( itr == row.end() ) {
            row.emplace_back( j, value );
        } else {
            itr->second = value;
        }
    }

    std::array<std::vector<edge_type>, N> edges_{};
};

  //////////////////////////////////////////////////
 // An undirected graph using an adjacency list. //
//////////////////////////////////////////////////

template <typename T, std::size_t N>
class Undirected_graph_list {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;
    using edge_type = std::pair<size_type, T>;

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Undirected_graph_list;

        entry( Undirected_graph_list &graph, size_type i, size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Undirected_graph_list *p_graph_;
        size_type i_;
        size_type j_;
    };

    Undirected_graph_list() = default;

    // The call edge( i, j ) returns the weight of the edge (v_i, v_j).
    [[nodiscard]] T edge( size_type i, size_type j ) const {
        assert( (i < N) && (j < N) );

        if ( i < j ) {
            std::swap( i, j );
        }

        auto const &row{ edges_[i] };

        auto itr{ std::find_if(
            row.begin(),
            row.end(),
            [j]( edge_type const &e ) {
                return e.first == j;
            }
        ) };

        return (itr == row.end()) ? T{ 0 } : itr->second;
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) {
        assert( (i < N) && (j < N) );
        return (j < i) ? entry{ *this, i, j } : entry{ *this, j, i };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );
        assert( i > j );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative; a weight of 0 removes the edge."
            };
        }

        auto &row{ edges_[i] };

        auto itr{ std::find_if(
            row.begin(),
            row.end(),
            [j]( edge_type const &e ) {
                return e.first == j;
            }
        ) };

        if ( value == T{ 0 } ) {
            if ( itr != row.end() ) {
                row.erase( itr );
            }

            return;
        }

        if ( itr == row.end() ) {
            row.emplace_back( j, value );
        } else {
            itr->second = value;
        }
    }

    std::array<std::vector<edge_type>, N> edges_{};
};

   /////////////////////////////////////////
  // A directed graph using a hash table //
 //  for the adjacency matrix.          //
/////////////////////////////////////////

template <typename T, std::size_t N>
class Directed_graph_hashed {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;
    using key_type = std::pair<size_type, size_type>;

    struct pair_hash {
        [[nodiscard]] std::size_t operator()(
            key_type const &p
        ) const noexcept {
            // The constant 0x9e3779b9 is derived from the golden ratio phi
            // and is commonly used in hash-combine functions due to its
            // good bit-mixing properties.
            return std::hash<size_type>{}( p.first )
                 ^ (std::hash<size_type>{}( p.second )
                 + 0x9e3779b9u
                 + (p.first << 6)
                 + (p.first >> 2));
        }
    };

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Directed_graph_hashed;

        entry( Directed_graph_hashed &graph, size_type i, size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Directed_graph_hashed *p_graph_;
        size_type i_;
        size_type j_;
    };

    Directed_graph_hashed() = default;

    [[nodiscard]] T edge( size_type i, size_type j ) const {
        assert( (i < N) && (j < N) );

        auto const itr{ edges_.find( key_type{ i, j } ) };
        return (itr == edges_.end()) ? T{ 0 } : itr->second;
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) {
        assert( (i < N) && (j < N) );
        return entry{ *this, i, j };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative; a weight of 0 removes the edge."
            };
        }

        key_type const key{ i, j };

        if ( value == T{ 0 } ) {
            edges_.erase( key );
            return;
        }

        edges_[key] = value;
    }

    std::unordered_map<key_type, T, pair_hash> edges_{};
};

   ////////////////////////////////////////////
  // An undirected graph using a hash table //
 //  for the adjacency matrix.             //
////////////////////////////////////////////

template <typename T, std::size_t N>
class Undirected_graph_hashed {
    static_assert( std::floating_point<T> );

public:
    using size_type = std::size_t;
    using key_type = std::pair<size_type, size_type>;

    struct pair_hash {
        [[nodiscard]] std::size_t operator()(
            key_type const &p
        ) const noexcept {
            // The constant 0x9e3779b9 is derived from the golden ratio phi
            // and is commonly used in hash-combine functions due to its
            // good bit-mixing properties.
            return std::hash<size_type>{}( p.first )
                 ^ (std::hash<size_type>{}( p.second )
                 + 0x9e3779b9u
                 + (p.first << 6)
                 + (p.first >> 2));
        }
    };

    class entry {
    public:
        [[nodiscard]] operator T() const {
            return std::as_const( *p_graph_ ).edge( i_, j_ );
        }

        entry &operator=( T const &value ) {
            p_graph_->set_edge( i_, j_, value );
            return *this;
        }

        entry &operator=( entry const &other ) {
            return *this = static_cast<T>( other );
        }

    private:
        friend class Undirected_graph_hashed;

        entry( Undirected_graph_hashed &graph, size_type i, size_type j ):
            p_graph_{ &graph },
            i_{ i },
            j_{ j } {
        }

        Undirected_graph_hashed *p_graph_;
        size_type i_;
        size_type j_;
    };

    Undirected_graph_hashed() = default;

    [[nodiscard]] T edge( size_type i, size_type j ) const {
        assert( (i < N) && (j < N) );

        if ( i < j ) {
            std::swap( i, j );
        }

        auto const itr{ edges_.find( key_type{ i, j } ) };
        return (itr == edges_.end()) ? T{ 0 } : itr->second;
    }

    [[nodiscard]] entry edge( size_type i, size_type j ) {
        assert( (i < N) && (j < N) );
        return (j < i) ? entry{ *this, i, j } : entry{ *this, j, i };
    }

private:
    void set_edge( size_type i, size_type j, T value ) {
        assert( (i < N) && (j < N) && (i != j) );
        assert( i > j );

        if ( !std::isfinite( value ) || (value < T{ 0 }) ) {
            throw std::logic_error{
                "Edge weights must be finite and nonnegative; a weight of 0 removes the edge."
            };
        }

        key_type const key{ i, j };

        if ( value == T{ 0 } ) {
            edges_.erase( key );
            return;
        }

        edges_[key] = value;
    }

    std::unordered_map<key_type, T, pair_hash> edges_{};
};
