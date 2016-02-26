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
  
  /** Type of values stored in node */
  typedef V node_value_type;
  
  /** Type of values stored in edge */
  typedef E edge_value_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

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
  
  typedef unsigned uid_type;
  
  typedef unsigned index_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
   : g_nodes(), g_edges(), g_evalues(), g_values(), g_indices(), g_idx2uid(){
	   g_num_nodes = 0;
	   g_num_edges = 0;
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
     * 
     * Invalid node graph pointer initialized using nullptr
     * Invalid node index not initialized to value
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
	  uid_type uid = n_graph->g_idx2uid[n_index];
      return n_graph->g_nodes[uid];
    }
    
    Point& position() {
		uid_type uid = n_graph->g_idx2uid[n_index];
		return n_graph->g_nodes[uid];
	}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return n_index;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n_index == n.n_index && n_graph == n.n_graph){
		  return true;
      }
      return false;
    }
    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * 
     * Global ordering is first by graph pointer ordering, 
     * then by index if nodes of same graph
     */
    bool operator<(const Node& n) const {
      if ((n_graph <n.n_graph) or (n_graph == n.n_graph and n_index < n.n_index)){
		  return true;
	  }
	  return false;
    }
    
    /** Return the value stored in this node. */
    node_value_type& value() {
		//std::cout <<"INDEX " << n_index <<std::endl;
		uid_type uid = n_graph->g_idx2uid[n_index];
		//std::cout<<"Returning Value "<<uid<<std::endl;
		return const_cast<graph_type*>(n_graph)->g_values[uid];
	}
	
	/** Return the value stored in this node */
	const node_value_type& value() const {
		uid_type uid = n_graph->g_idx2uid[n_index];
		return n_graph->g_values[uid];
	}
	
	/** Return the degree of this node (number of edges) */
	size_type degree() const {
		uid_type uid = n_graph->g_idx2uid[n_index];
		return n_graph->connectivity(uid);
	}
	
	/** Return a begin iterator that iterates over edges of node */
	incident_iterator edge_begin() const {
		return IncidentIterator(n_graph,*this,0);
	}
	
	/** Return an end iterator that iterates for edges of node */
	incident_iterator edge_end() const {
		return IncidentIterator(n_graph,*this,degree());
	}

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	
    graph_type* n_graph = nullptr;
    index_type n_index;
    
    /** Construct a valid Node object
     * 
     * @return node object with given graph and index
     * 
     */
    Node(const graph_type* node_graph, index_type node_index) {
		n_index = node_index;
		n_graph = const_cast<graph_type*>(node_graph);	
	}
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return g_idx2uid.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    g_nodes.push_back(position);
    g_values.push_back(v);
    ++g_num_nodes;
    std::vector<size_type> empty_size_vector;
    g_edges.push_back(empty_size_vector);
    g_evalues.resize(g_num_nodes);
    for (size_type i = 0; i < g_num_nodes; ++i) {
		g_evalues[i].resize(g_num_nodes);
	}
	g_idx2uid.push_back(g_nodes.size()-1);
	//std::cout << "$$$$$$$$$$$$$$$$$$$$$$" <<std::endl;
	//for (int i = 0; i < g_idx2uid.size(); ++i) {
	//std::cout << "$$ " << g_idx2uid[i] <<std::endl;
    //}
    //std::cout<<"yo" << node(0).value()<<std::endl;
    //std::cout<<g_idx2uid[0]<<std::endl;
    //std::cout<<g_values[1]<<std::endl;
    return Node(this, g_idx2uid.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.n_graph == this){
		return true;
	}
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(index_type i) const {
    return Node(this, i);
  }
  
  size_type remove_node(const Node& n) {
	  for (size_type k = 0; k < g_edges[n.index()].size();) {
		  const node_type n2 = node(g_edges[n.index()][k]);
		  if (has_edge(n,n2)) {
			  remove_edge(n,n2);
		  } else {
			  ++k;
		  }
	  }
	  g_idx2uid.erase(g_idx2uid.begin() + n.index());
	  --g_num_nodes;
	  //remove edges
	  return 1;
  }
  
  node_iterator remove_node(node_iterator n_it) {
       g_idx2uid.erase(g_idx2uid.begin() + (*n_it).index());
	  --g_num_nodes;
	  //remove edges
	  return n_it;
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(e_graph,e_node1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(e_graph,e_node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	    if (e_graph == e.e_graph and ((e_node1 == e.e_node1 and 
	    e_node2 == e.e_node2) or (e_node1 == e.e_node2 and e_node2 == e.e_node1))){
			return true;
	    }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if ((e_graph <e.e_graph) or (e_graph == e.e_graph and node1().index() < e.node1().index())) {
	        return true;
	    }
      return false;
    }
    
    double length() const {
		return norm(e_graph->node(e_node1).position()-e_graph->node(e_node2).position());
	}
	
	edge_value_type& value() {
		index_type i_min = std::min(e_node1,e_node2);
		index_type i_max = std::max(e_node1,e_node2);
		return const_cast<graph_type*>(e_graph)->g_evalues[i_min][i_max];
	}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* e_graph = nullptr;
    index_type e_node1, e_node2;
    
    /** Construct a valid Edge object
     * 
     * @return edge object with given graph and nodes
     * 
     */
    Edge(const graph_type* edge_graph, index_type edge_node1, index_type edge_node2) {
		e_graph = const_cast<graph_type*>(edge_graph);
		e_node1 = edge_node1;
		e_node2 = edge_node2;
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return g_num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(index_type i) const {
    auto ei = edge_begin();
    for (index_type j = 0; j < i; ++j) {
		++ei;
	}
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
	  if (g_edges[a.index()].size() != 0) { 
	      for (index_type i = 0; i < g_edges[a.index()].size(); ++i) {
		      if (g_edges[a.index()][i] == b.index()) {
		          return true;
		      }
	      }
      }
      return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& v = edge_value_type()) {
    edge_type new_edge = Edge(this, a.index(), b.index());
    if (!has_edge(a, b)) {
	    g_edges[a.index()].push_back(b.index());
	    g_edges[b.index()].push_back(a.index());
	    g_evalues[std::min(a.index(),b.index())][std::max(a.index(),b.index())] = v;
	    ++g_num_edges;
	}
    return new_edge;
  }
  
  size_type remove_edge(const Node& n1, const Node& n2) {
	  if (has_edge(n1,n2)) {
		  index_type i1 = n1.index();
		  index_type i2 = n2.index();
		  for (index_type k = 0; k < g_edges[i1].size(); ++k) {
			  if (g_edges[i1][k] == i2) {
				  g_edges[i1].erase(g_edges[i1].begin() + k);
			  }
		  }
		  for (index_type k = 0; k < g_edges[i2].size(); ++k) {
			  if (g_edges[i2][k] == i1) {
				  g_edges[i2].erase(g_edges[i2].begin() + k);
			  }
		  }
		  --g_num_edges;
      }
	  return 0;
  }
  
  size_type remove_edge(const Edge& e) {
	  if (has_edge(e.node1(),e.node2())) {
		  index_type i1 = e.node1().index();
		  index_type i2 = e.node2().index();
		  for (index_type k = 0; k < g_edges[i1].size(); ++k) {
			  if (g_edges[i1][k] == i2) {
				  g_edges[i1].erase(g_edges[i1].begin() + k);
			  }
		  }
		  
		  for (index_type k = 0; k < g_edges[i2].size(); ++k) {
			  if (g_edges[i2][k] == i1) {
				  g_edges[i2].erase(g_edges[i2].begin() + k);
			  }
		  }
		  --g_num_edges;
      }
	  return 0;
  }
  
  edge_iterator remove_edge(edge_iterator e_it) {
	  edge_type e = (*e_it);
	  if (has_edge(e.node1(),e.node2())) {
		  index_type i1 = e.node1().index();
		  index_type i2 = e.node2().index();
		  for (index_type k = 0; k < g_edges[i1].size(); ++k) {
			  if (g_edges[i1][k] == i2) {
				  g_edges[i1].erase(g_edges[i1].begin() + k);
			  }
		  }
		  
		  for (index_type k = 0; k < g_edges[i2].size(); ++k) {
			  if (g_edges[i2][k] == i1) {
				  g_edges[i2].erase(g_edges[i2].begin() + k);
			  }
		  }
		  --g_num_edges;
      }
	  return 0;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    g_nodes.clear();
    g_edges.clear();
    g_evalues.clear();
    g_values.clear();
    g_indices.clear();
    g_idx2uid.clear();
    g_num_nodes = 0;
    g_num_edges = 0;
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
    
    /** Return node corresponding to this iterator */
    Node operator*() const {
		return Node(ni_graph, ni_index);
	}
	
	/** Increment iterator */
	node_iterator& operator++() {
		++ni_index;
		return *this;
	}
	
	/** Test equality of iterators */
	bool operator==(const node_iterator& ni) const {
		if (ni_index == ni.ni_index and ni_graph == ni.ni_graph) {
			return true;
		}
		return false;
	}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    index_type ni_index;
    graph_type* ni_graph;
    
    NodeIterator(const graph_type* current_graph, index_type index) {
		ni_graph = const_cast<graph_type*>(current_graph);
		ni_index = index;
	}
    
  };
  
  /** Return begin node iterator of current graph */
  node_iterator node_begin() const {
	  return NodeIterator(this, 0);
  }
  
  /** Return end node iterator of current graph */
  node_iterator node_end() const {
	  return NodeIterator(this, g_num_nodes);
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
    
    /** Return edge of current edge iterator */
    Edge operator*() const {
	    return Edge(ei_graph, ei_node1_i, ei_graph->index_node2(*this));
	}
	
	/** Increment current edge iterator */
	EdgeIterator& operator ++() {
		do {
		    ++ei_node2_p;
		    if (ei_node2_p == ei_graph->connectivity(*this)) {
			    ei_node2_p = 0;
			    do {
				    ++ei_node1_i;
				    if (ei_node1_i == ei_graph->graph_size()) {
						return *this;
					}
			    } while (ei_graph->connectivity(*this) == 0);
		    }
		} while (ei_node1_i < ei_graph->index_node2(*this));
		return *this;
	}
	
	/** Test equality of two edge iterators */
	bool operator==(const EdgeIterator& eit) const {
		if (ei_graph == eit.ei_graph and ei_node1_i == eit.ei_node1_i and ei_node2_p == eit.ei_node2_p) {
			return true;
		}
		return false;
	}

   private:
    friend class Graph;
    graph_type* ei_graph;
    index_type ei_node1_i;
    index_type ei_node2_p;
    
    EdgeIterator(const graph_type* current_graph, index_type node1_i, index_type node2_p) {
		ei_graph = const_cast<graph_type*>(current_graph);
		ei_node1_i = node1_i;
		ei_node2_p = node2_p;
	}
  };
  
  /** Return begin edge iterator of current graph */
  edge_iterator edge_begin() const {
	  index_type i = 0;
	  while (g_edges[i].size() == 0) {
		  ++i;
	  }
	  return edge_iterator(this,i,0);
  }
  
  /** Return end edge iterator of current graph */
  edge_iterator edge_end() const {
	  return edge_iterator(this,g_num_nodes,0);
  }
  
  /** Return index of node two of current edge */
  index_type index_node2(const edge_iterator& eit) const {
	  return g_edges[eit.ei_node1_i][eit.ei_node2_p];
  }
  
  /** Return index of node two of current edge */
  index_type index_node2(const incident_iterator& iit) const {
	  return g_edges[iit.iit_node][iit.iit_edge_p];
  }
  
  /** Return connectivity of node 1 of current edge */
  size_type connectivity(const edge_iterator& eit) const {
	  return g_edges[eit.ei_node1_i].size();
  }
  
  /** Return connectivity of node n */
  size_type connectivity(const index_type n) const {
	  return g_edges[n].size();
  }
  
  /** Return total number of nodes of graph */
  size_type graph_size() {
	  return g_num_nodes;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    
    /** Return edge of current iterator */
    Edge operator*() const {
		return Edge(iit_graph,iit_node,iit_graph->index_node2(*this));
	}
	
	/** Increment current iterator */
	incident_iterator& operator++() {
		++iit_edge_p;
		return *this;
	}
	
	/** Test equality of two iterators */
	bool operator==(const incident_iterator& iit) const {
		if (iit_node == iit.iit_node and iit_edge_p == iit.iit_edge_p and iit_graph == iit.iit_graph) {
			return true;
		}
		return false;
	}

   private:
    friend class Graph;
    graph_type* iit_graph;
    index_type iit_node;
    index_type iit_edge_p;
    
    IncidentIterator(const graph_type* current_graph, const Node& n, size_type edge_p) {
		iit_graph = const_cast<graph_type*>(current_graph);
		iit_node = n.index();
		iit_edge_p = edge_p;
	}
  };

 private:
  
  std::vector<Point> g_nodes;
  std::vector<std::vector<index_type>> g_edges;
  std::vector<std::vector<edge_value_type>> g_evalues;
  std::vector<node_value_type> g_values;
  std::vector<size_type> g_indices;
  size_type g_num_nodes;
  size_type g_num_edges;
  std::vector<uid_type> g_idx2uid;

};

#endif // CME212_GRAPH_HPP
