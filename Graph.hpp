#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#include <thrust/iterator/transform_iterator.h>
#include <functional>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  /** Type of the nodes uid. */
  typedef unsigned uid_type;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;
  /** Type of the nodes' value. */
  typedef V node_value_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;
  /** Type of the edge' value. */
  typedef E edge_value_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    n_nodes = 0;
    n_edges = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   *
   * Size : 8 Bytes
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return (graph_->node_position(uid_));
    }

    /** Return a modifiable node position. */
    Point& position() {
      return (graph_->node_position(uid_));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return (graph_->nodes_.at(uid_).idx_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if(graph_ == n.graph_)
        return uid_ < n.uid_;
      else
        return graph_ < n.graph_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value.
     * Complexity: O(1).
     */
    node_value_type& value() {
      return (graph_->node_value(uid_));
    }

    /** Return this node's value.
     * Complexity: O(1).
     */
    const node_value_type& value() const {
      return (graph_->node_value(uid_));
    }

    /** Return this node's degree. 
     * Complexity: O(1).
     */
    size_type degree() const {
      return (graph_->adj_).at(uid_).size();
    }

    /** Return an iterator to the first incident of this node.
     * @post If result_iterator != edge_end(),
     *       *result_iterator.node1() == @a *this
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(*this, 0);
    }

    /** Returns an iterator corresponding to the past-the-end 
     * incident edge of the node.
     * Should not be dereferenced.
     * If the node does not have any neighbors, result_iterator == edge_begin().
     *
     * Complexity: O(1).
     */
    incident_iterator edge_end() const {
      return IncidentIterator(*this, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // Pointer back to the graph's nodes
    Graph* graph_;
    // This node unique identifier
    uid_type uid_;
    /** Private constructor for valid Node objects. */
    Node(const Graph* graph, uid_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_nodes;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The value of the node (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    Node n = Node(this, next_uid());

    // add node's position, value and idx
    nodes_.push_back(NodeInfo(position, value, n_nodes));
    idx2uid_.push_back(n.uid_);

    n_nodes ++;

    // initialize vector of adjacent nodes
    adj_.push_back(std::vector<EdgeInfo>());

    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(n.graph_ == this && n.index() < n_nodes)
      return (n.uid_ == idx2uid_[n.index()]);
    else
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    uid_type uid = idx2uid_.at(i);
    return Node(this, uid);
  }

  /** Remove the node @a n and all incident edges.
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph: n.index() < num_nodes
   * @post old num_nodes() == new num_nodes() + 1
   * @post old n_edges() == new num_edges() + old n.degree()
   * @post all elements with index i != @a n.index() and i != num_nodes()-1
   *       are not affected. In particular their index stay unchanged
   * @return The index of the next node if it exists (old n.index())
   *         or new num_nodes().
   *
   * @note invalidates index old n.index() and the last index (num_nodes()-1)
   *
   * Complexity: O(n.degree()), worst case: O(num_nodes())
   */
  size_type remove_node(const Node& n) {
    return remove_node(n.uid_);
  }

  /** Remove the node with position @a node_it and all incident edges.
   * @param[in] node_it Iterator pointing to the node to remove
   * @pre node_it != node_end()
   * @post old num_nodes() == new num_nodes() + 1
   * @post old n_edges() == new num_edges() + old (*node_it).degree()
   * @post all iterators in the range [node_begin(), @a node_it) are not affected
   * @post all elements with position i > @a node_it
   *       and i < num_nodes()-1 are not affected
   * @return an iterator to the next node, or node_end().
   *
   * @note invalidates previously defined iterators in the range
   *       [@a node_it, node_end())
   *
   * Complexity: O((*node_it).degree()), worst case: O(num_nodes()).
   */
  node_iterator remove_node(node_iterator node_it) {
    uid_type uid = idx2uid_.at(node_it.i_);
    size_type idx = remove_node(uid);
    return(NodeIterator(this, idx));
  }



  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   *
   * Size : 12 Bytes
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, nid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, nid2_);
    }

    /** Return the distance between the two nodes */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return the value of this edge.
     * Complexity: O(max node1().degree(), node2().degree()),
     *             worst case: O(num_nodes()).
     */
    edge_value_type& value() {
      return graph_->edge_value(nid1_, nid2_);
    }

    /** Return the value of this edge.
     * Complexity: O(max node1().degree(), node2().degree()),
     *             worst case: O(num_nodes()).
     */
    const edge_value_type& value() const {
      return graph_->edge_value(nid1_, nid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_ != e.graph_)
        return false;

      return ((nid1_ == e.nid1_ && nid2_ == e.nid2_)
              || (nid1_ == e.nid2_ && nid2_ == e.nid1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      uid_type m1 = std::min(nid1_, nid2_);
      uid_type m2 = std::min(e.nid1_, e.nid2_);
      // first compare graph adresses
      // then min ids, and finally max ids if necessary
      if(graph_ != e.graph_)
        return graph_ < e.graph_;
      else if(m1 != m2)
        return m1 < m2;
      else
        return std::max(nid1_, nid2_) < std::max(e.nid1_, e.nid2_);
    } 

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // A pointer back to the graph
    Graph* graph_;
    // IDs of the nodes connected by this edge
    uid_type nid1_, nid2_;

    /** Private constructor through Nodes, for valid Edge objects. */
    Edge(const Graph* graph, Node n1, Node n2)
      : graph_(const_cast<Graph*>(graph)), nid1_(n1.uid_), nid2_(n2.uid_) {
    }

    /** Private constructor through Node Ids, for valid Edge objects. */
    Edge(const Graph* graph, uid_type id1, uid_type id2)
      : graph_(const_cast<Graph*>(graph)), nid1_(id1), nid2_(id2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return n_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
   */
  Edge edge(size_type i) const {
    std::cout << "WARNING: Graph::edge() is deprecated, ";
    std::cout << "should better use edge iterators instead" << std::endl;

    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(min(a.degree(), b.degree()); worst case: O(num_nodes())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // small optimization : search on node with smallest degree
    if(a.degree() < b.degree())
      return (find_adj_node(b.uid_, a.uid_) != a.degree());

    else
      return (find_adj_node(a.uid_, b.uid_) != b.degree());
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @param[in] a,b The two nodes the edge connects
   * @param[in] value The value of the edge (optional)
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   * @post if old has_edge(a,b) == false, then result_edge.value() == @a value
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(min(a.degree(), b.degree()); worst case: O(num_nodes())
   */
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& value = edge_value_type()) {
    Edge e(this, a, b);
    if(has_edge(a,b))
      return e;

    else{
      // add each node in the adjacency list of the other one
      // twin_idx is the degree before insertion
      adj_.at(a.uid_).push_back(EdgeInfo(b.uid_, b.degree(), value));
      adj_.at(b.uid_).push_back(EdgeInfo(a.uid_, a.degree()-1, value));

      ++n_edges;

      return e;
    }
  }

  /** Remove the edge connecting nodes @a n1 and @a n2.
   * @param[in] n1,n2 The nodes whose connecting edge should be removed
   * @pre @a n1 and @a n2 are two valid node of the graph
   * @post if result == 1, old num_edges() == new num_edges() + 1
   * @return The number of edges removed (0 or 1).
   *
   * @note Invalidates all iterators on edges in the range
   *       [@a e_it, edge_end())
   *
   * Complexity: O(n1.degree() + n2.degree()),
   *             worst case: O(num_nodes()).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    return remove_edge(n1.uid_, n2.uid_);
  }


  /** Remove edge @a e.
   * @param[in] e The edge to remove
   * @pre The nodes of @a e are valid nodes of the graph
   * @post if result == 1, old num_edges() == new num_edges() + 1
   * @return The number of edges removed (0 or 1).
   *
   * @note Invalidates all iterators on edges in the range
   *       [@a e_it, edge_end())
   *
   * Complexity: O(e.node1().degree() + e.node2().degree()),
   *             worst case: O(num_nodes()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.nid1_, e.nid2_);
  }

  /** Remove edge at position @a e_it.
   * @param[in] e_it Iterator pointing to the edge to remove
   * @pre edge_begin() <= e_it < edge_end()
   * @post old num_edges() == new num_edges() + 1
   * @post all elements in the range [edge_begin(), @a e_it) are not affected
   * @post all elements with position i > @a e_it are moved to position i-1
   * @return An iterator to the next edge, or edge_end().
   *
   * @note Invalidates all iterators on edges in the range
   *       [@a e_it, edge_end())
   *
   * Complexity: O(*e_it.node1().degree() + *e_it.node2().degree()) in average,
   *             worst case: O(num_nodes() + num_edge()).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    // e_it is invalid after remove_edge
    e_it.fix();
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    n_nodes = 0;
    n_edges = 0;
    nodes_.clear();
    idx2uid_.clear();
    adj_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Derefence the NodeIterator to return the Node it corresponds to.
     * @pre @a this != node_end().
     *
     * Complexity: O(1).
     */
    Node operator*() const {
      return graph_->node(i_);
    }

    /** Increment the iterator.
     * @pre @a this != node_end().
     * @post If result_iterator != node_end(),
     *       *old_iterator.index() + 1 == *result_iterator.index()
     *
     * Complexity: O(1).
     */
    NodeIterator& operator++() {
      ++ i_;
      return *this;
    }

    /** Test whether this iterator and @a ni are equal.
     * Two iterators that are not node_end() are equal if and only if
     * the nodes they point to are equal : 
     *   it1 == it2 <=> *it1 == *it2
     *
     * Complexity: O(1).
     */
    bool operator==(const NodeIterator& ni) const {
      if(graph_ == ni.graph_)
        return (i_ == ni.i_);
      else
        return false;
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    /* A pointer to the graph whose nodes it iterates through. */
    Graph* graph_;
    /* Index of the node this iterator represents. */
    size_type i_;

    /* A valid constructor for the Graph class. */
    NodeIterator(const Graph* graph, size_type i)
      : graph_(const_cast<Graph*>(graph)), i_(i) {
    }

  };
  
  struct node_of_uid : public std::unary_function<uid_type,Node>
  {
       node_of_uid(Graph* g) : g_(g) {
	   }
       Node operator()(uid_type uid) { return Node(g_, uid); }
     private:
       Graph* g_ = NULL; 
  };
  
  struct NodeIterator2 : thrust::transform_iterator<node_of_uid,std::vector<uid_type>, Node> {
      // Import super class's constructors
      using NodeIterator2::transform_iterator::transform_iterator;
	  
	  NodeIterator2(const graph_type* g, uid_type uid)
	    :  NodeIterator::transform_iterator(uid, node_of_uid(g))
	    {}  
  };
  
  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns an iterator to the first node of the graph.
   * @post result_iterator->index() == 0.
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns an iterator corresponding to the past-the-end node of the graph.
   * Should not be dereferenced.
   * If the graph does not have any node, result_iterator == node_begin().
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Dereference the iterator to return the edge it corresponds to.
     * @pre @a this != edge_end().
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      uid_type id2 = (graph_->adj_).at(i1_).at(i2_).adj_uid_;
      return Edge(graph_, i1_, id2);
    }

    /** Increment the iterator.
     * Increment such that result_iterator can be dereferenced if it
     * is not edge_end().
     * Iterating from edge_begin() to edge_end() will visit each edge
     * exactly once.
     * @post If result_iterator != edge_end(),
     *       *result_iterator.node1() < *result_iterator.node2().
     *
     * Complexity: O(1) if graph is not too sparse, 
     *             worst case: O(num_nodes() + num_edges()).
     *             In any case, iterating through all the edges will not
     *             take more than O(num_nodes() + num_edges()).
     */
    EdgeIterator& operator++() {
      ++ i2_;
      fix();
      return *this;
    }

    /** Test whether this iterator equals @a ei.
     *
     * Complexity: O(1).
     */
    bool operator==(const EdgeIterator& ei) const {
      return ((i1_ == ei.i1_) && (i2_ == ei.i2_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /* A pointer to the graph whose edges it iterates through. */
    Graph* graph_;
    /* First index in adj_ */
    uid_type i1_;
    /* Second index in adj_ */
    size_type i2_;

    /* Valid constructor for the Graph class.
     * @post If result_iterator != edge_end(),
     *       *result_iterator.node1() < *result_iterator.node2().
     *
     * Complexity: O(1) if graph is not too sparse, 
     *             worst case: O(num_nodes() + num_edges()).
     */
    EdgeIterator(const Graph* graph, uid_type i1, size_type i2)
      : graph_(const_cast<Graph*>(graph)), i1_(i1), i2_(i2) {
      fix();
    }

    /** Helper function to ensure (graph_->adj_)[i1_][i2_] is defined.
     * and that i1_ < (graph_->adj_)[i1_][i2_].
     * @post If *this != edge_end(),
     *  (graph_->adj_)[i1_][i2_] is valid
     *  and i1_ < (graph_->adj_)[i1_][i2_] 
     *      ie. *this.node1().id < *this.node2().id
     */
    void fix() {
      if(i1_ >= graph_->next_uid()){
        // then it should be edge_end()
        i1_ = graph_->next_uid();
        i2_ = 0;
      }
      else if(i2_ >= (graph_->adj_).at(i1_).size()){
        ++ i1_;
        i2_ = 0;
        fix();
      }
      else if(i1_ >= (graph_->adj_).at(i1_).at(i2_).adj_uid_){
        ++ i2_;
        fix();
      }
      // else everything is fine
    }
  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns an iterator to the first edge of the graph.
   * @post If result_iterator != edge_end(),
   *       result_iterator->node1() < result_iterator->node2().
   *
   * Complexity: O(1) if graph is not too sparse,
   *             worst case: O(num_nodes()).
   */
  edge_iterator edge_begin() const {
    EdgeIterator ei(this, 0, 0);
    return ei;
  }

  /** Returns an iterator corresponding to the past-the-end edge of the graph.
   * Should not be dereferenced.
   * If the graph does not have any edge, result_iterator == edge_begin().
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, next_uid(), 0);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereference the iterator to return an incident edge.
     * @pre this != edge_end().
     * @post result_edge.node1() == node that spawns the incident_iterator.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      uid_type id2 = (graph_->adj_).at(nid_).at(i_).adj_uid_;
      return Edge(graph_, nid_, id2);
    }

    /** Increment the iterator.
     * @pre this != edge_end().
     *
     * Complexity: O(1).
     */
    IncidentIterator& operator++() {
      ++ i_;
      return *this;
    }

    /** Test whether this is equal to @a ii.
     * Two IncidentIterator are equal if they iterates through the incident edges
     * of the same node, and if they point to the same incident edge.
     *
     * Complexity: O(1).
     */
    bool operator==(const IncidentIterator& ii) const {
      return((graph_ == ii.graph_) &&
             (nid_ == ii.nid_) &&
             (i_ == ii.i_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    /* A pointer back to the graph. */
    Graph* graph_;
    /* The uid of the node whose edges this iterates. */
    uid_type nid_;
    /* Index of adjacent node. */
    size_type i_;

    /* A valid constructor for the Graph class. */
    IncidentIterator(const Node& n, size_type i)
      : graph_(const_cast<Graph*>(n.graph_)), nid_(n.uid_), i_(i) {
    }
  };

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // A structure to store a node's attribute
  struct NodeInfo {
    Point pos_;
    node_value_type val_;
    size_type idx_;

    // constructor
    NodeInfo(Point pos, node_value_type val, size_type idx)
      : pos_(pos), val_(val), idx_(idx) {}
  };

  // A structure to store an edge's attribute
  struct EdgeInfo {
    // uid of the second node of the edge
    uid_type adj_uid_;
    // idx of the first node in the second node's adjacency list
    size_type twin_idx_;
    // the edge's value
    edge_value_type val_;

    // constructor
    EdgeInfo(uid_type uid, size_type twin_idx, edge_value_type val)
      : adj_uid_(uid), twin_idx_(twin_idx), val_(val) {}
  };

  // The container for the nodes' info, the index is uid
  std::vector<NodeInfo> nodes_;

  // A map from a node's idx to its uid
  std::vector<uid_type> idx2uid_;
 
  // The container for edges, stored as an adjacency list
  // First vector has length next_uid() and is ordered by the nodes' id
  // Second vector is not sorted
  std::vector<std::vector<EdgeInfo>> adj_;
  /** Representation invariant:
   * for all valid edges (n1,n2),
   * twin_idx1 and twin_idx2 are such that
   *   adj_.at(n1.uid).at(twin_idx2).ajd_uid_ == n2.uid
   *   adj_.at(n1.uid).at(twin_idx2).twin_idx_ == twin_idx1
   *   adj_.at(n2.uid).at(twin_idx1).ajd_uid_ == n1.uid
   *   adj_.at(n2.uid).at(twin_idx1).twin_idx_ == twin_idx2
   * adj_.at(n1.uid).at(twin_idx2).val_ is valid iff n1.uid < n2.uid
   */

  // The number of nodes in the graph
  // n_nodes == idx2uid_.size()
  size_type n_nodes;

  // The number of undirected edges in the graph
  size_type n_edges;


  /** Return the next uid to attribute to a new node.
   * @post result == nodes_.size() == adj_.size()
   * Complexity: O(1).
   */
  uid_type next_uid() const {
    return nodes_.size();
  }

  /** Return the position of the node with uid @a id.
   * @pre 0 <= @a id < next_uid()
   *
   * Complexity: O(1).
   */  
  Point& node_position(uid_type id) {
    return nodes_.at(id).pos_;
  }

  /** Return the value of the node with uid @a id.
   * @pre 0 <= @a id < next_uid()
   *
   * Complexity: O(1).
   */
  node_value_type& node_value(uid_type id) {
      return nodes_.at(id).val_;
  }

  /** Return the value of the edge connecting nodes @a uid1 and @a uid2.
   * @pre (Node(uid1), Node(uid2)) is a valid edge of the graph.
   *
   * Complexity: O(max degree(uid1),degree(uid2)), worst case: O(num_nodes()).
   */
  edge_value_type& edge_value(uid_type uid1, uid_type uid2) {
    // edge value is relevent only for uid1 < uid2 in adj_
    size_type idx = find_adj_node(std::max(uid1, uid2), std::min(uid1, uid2));
    return adj_.at(std::min(uid1, uid2)).at(idx).val_;
  }

  /** Find the index of an adjacent node @a uid in adj_.at(@a from_uid).
   * @pre @uid is in adj_.at(@a from_uid),
   *      ie. (Node(uid), Node(from_uid)) is a valid edge of the graph
   *
   * Complexity: O(degree(from_uid)), worst case: O(num_nodes())
   */
  size_type find_adj_node(uid_type uid, uid_type from_uid) const {
    size_type idx = 0;
    size_type d = adj_.at(from_uid).size();
    while(idx < d && adj_.at(from_uid).at(idx).adj_uid_ != uid)
      ++idx;
    return idx;
  }

  /** Remove the node with id @a uid and all incident edges.
   * @pre This uid represents a valid node:
   *       - uid < next_uid()
   *       - Node(uid).index() < num_nodes
   * @post old num_nodes() == new num_nodes() + 1
   * @post old n_edges() == new num_edges() + old degree(uid)
   * @post all elements with index i != Node(@a uid).index() and i != num_nodes()-1
   *       are not affected. In particular their index stay unchanged
   * @return The index of the next node if it exists (old Node(@a uid).index())
   *         or new num_nodes().
   *
   * @note invalidates indexes
   *        - Node(@a uid).index()
   *        - num_nodes() (last index)
   *
   * Complexity: O(degree()), worst case: O(num_nodes()).
   */
  size_type remove_node(uid_type uid) {
    // change idx and 'remove' from nodes in O(1)
    size_type idx = nodes_.at(uid).idx_;
    if(idx < n_nodes-1) {
      // set idx to uid that will become invalid
      nodes_.at(uid).idx_ = n_nodes-1;
      // switch with node that has the last idx
      uid_type uid2 = idx2uid_.back();
      nodes_.at(uid2).idx_ = idx;
      idx2uid_.at(idx) = uid2;
    }
    // remove last element
    idx2uid_.pop_back();

    // remove all incident edges
    for(EdgeInfo ei : adj_.at(uid)) {
      remove_adj_node(ei.twin_idx_, ei.adj_uid_);
      --n_edges;
    }

    // reset adjencecy list of node uid
    adj_.at(uid) = std::vector<EdgeInfo>();

    --n_nodes;

    return idx;
  }

  /** Remove node at position @a idx in the adjacency list of
   * node with uid @a from_uid.
   * @pre from_uid < adj_.size() == next_uid()
   * @pre idx < adj_.at(from_uid).size() == Node(from_uid).degree()
   * @post new adj_.at(from_uid).size() == old adj_.at(from_uid).size() -1
   * Potentially invalidates previously defined indeces in adj_.at(from_uid)
   * and 'twin indeces'.
   *
   * Complexity: O(1)
   */
  void remove_adj_node(size_type idx, uid_type from_uid) {
    if(idx < adj_.at(from_uid).size() - 1) {
      // need to swap nodes to remove in O(1)
      adj_.at(from_uid).at(idx) = adj_.at(from_uid).back();

      // need to update twin_idx of the new node at position idx
      uid_type uid2 = adj_.at(from_uid).at(idx).adj_uid_;
      size_type idx2 = adj_.at(from_uid).at(idx).twin_idx_;
      adj_.at(uid2).at(idx2).twin_idx_ = idx;
    }
    // remove last element in O(1)
    adj_.at(from_uid).pop_back();
  }


  /** Remove the edge connecting nodes with uid @a uid1 and @a uid2.
   * @pre @a uid1 and @a uid2 correpond to valid nodes of the graph
   * @post if result == 1, old n_edges == new n_edges + 1
   * @return number of edges removed (0 or 1)
   * 
   * @note Invalidates all iterators on edges in the range
   *       [@a e_it, edge_end())
   *
   * Complexity: O(degree(uid1) + degree(uid2)),
   *             worst case: O(num_nodes()).
   */
  size_type remove_edge(uid_type uid1, uid_type uid2) {
    size_type idx1 = find_adj_node(uid1, uid2);
    size_type idx2 = find_adj_node(uid2, uid1);

    if(idx1 == adj_.at(uid2).size() || idx2 == adj_.at(uid1).size())
      return 0;

    remove_adj_node(idx1, uid2);
    remove_adj_node(idx2, uid1);

    -- n_edges;
    
    return 1;
  }

};

#endif // CME212_GRAPH_HPP
