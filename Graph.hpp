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
template <typename V>
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
  
  /** Type of user-specified value. */
  typedef V node_value_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
   : g_nodes(), g_values(), g_edges() {
	   g_node_begin = node_iterator(Node(this,0));
	   g_node_end = node_iterator(Node(this,0));
	   
	   g_edge_begin = edge_iterator(this,0,0);
	   g_edge_end = edge_iterator(this,0,0);  
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
  class Node : private totally_ordered<Node>{
   public:
   
    /** Define node value type. */
    node_value_type& value() {
        return n_graph->g_values[n_index];
    }

    const node_value_type& value() const {
        return n_graph->g_values[n_index];
    }

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
      return n_graph->g_nodes[n_index];
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

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	
    const graph_type* n_graph = nullptr;
    size_type n_index;
    
    /** Construct a valid Node object
     * 
     * @return node object with given graph and index
     * 
     */
    Node(const graph_type* node_graph, size_type node_index) {
		n_index = node_index;
		n_graph = node_graph;
	}
    
  };
  
  /** Define node iterator class
   * 
   * 
   */
   
   class node_iterator : private totally_ordered<node_iterator> {
	   public:
	     
	     node_iterator() {
		 }
		 
	     /** Return node that iterator is representing */
	     Node operator*() const {
			 return ni_node;
		 }
		 
		 /** Increment node iterator to next node in graph */
		 node_iterator& operator++() {
			 ++ni_node.n_index;
			 return *this;
		 }
		 
		 /** Compare equality of two node iterators */
		 bool operator==(const node_iterator& nit) const {
			 return ni_node == nit.ni_node;
		 }
	   
	   private:
	     friend class Graph;
	     node_type ni_node;
	     
	     node_iterator(node_type no) {
			 ni_node = no;
		 }
	   
   };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return g_nodes.size();
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
  Node add_node(const Point& position, const node_value_type& n = node_value_type()) {
    g_nodes.push_back(position);
    g_values.push_back(n);
    ++g_node_end;
    std::cout<<g_nodes[1] <<std::endl;
    return Node(this, g_nodes.size()-1);
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
  Node node(size_type i) const {
    return Node(this, i);
  }
  
  /** Return node iterator to first element in sequence
   * 
   * 
   * 
   */
  node_iterator node_begin() const {
	  return g_node_begin;
  }
   
  /** Return node iterator to one past last element in sequence
   * 
   * 
   * 
   */
   node_iterator node_end() const {
	   return g_node_end;
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

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const graph_type* e_graph = nullptr;
    size_type e_node1, e_node2;
    
    /** Construct a valid Edge object
     * 
     * @return edge object with given graph and nodes
     * 
     */
    Edge(const graph_type* edge_graph, size_type edge_node1, size_type edge_node2) {
		e_graph = edge_graph;
		e_node1 = edge_node1;
		e_node2 = edge_node2;
	}
	
  };
  
  /** Defenition of edge iterator class
   * 
   * 
   * 
   */
    
  class edge_iterator : private totally_ordered<edge_iterator> {
	  public:
	  
	    edge_iterator() {
		}
		
	    /** Return edge referenced by this edge iterator */
	    Edge operator*() const {
			return Edge(ei_graph, node1_index, ei_graph->index_node2(*this));
		}
		
		/** Increment edge iterator to next edge in sequence */
		edge_iterator& operator++() {
			while(ei_graph->index_node2(*this) < node1_index) {
			    if (node2_pos == ei_graph->connectivity((*this).node1_index) - 1 ) {
				    node1_index += 1;
				    node2_pos = 0;
			    }
			    else {
				    node2_pos += 1;
			    }
		    }
		    return *this;
		}
		
		/** Compare equality of two edge iterators */
		bool operator==(const edge_iterator& eit) const {
			if (node1_index == eit.node1_index and node2_pos == eit.node2_pos and ei_graph == eit.ei_graph) {
				return true;
			}
			return false;
		}
	    
	  
	  private:
	    friend class Graph;
	    size_type node1_index;
	    size_type node2_pos;
	    const graph_type* ei_graph;
	    
	    edge_iterator(const graph_type* edge_graph, size_type n1_i, size_type n2_pos) {
			ei_graph = edge_graph;
			node1_index = n1_i;
			node2_pos = n2_pos;
		}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
	  size_type edge_count = 0;
      for (auto ei = edge_begin(); ei != edge_end(); ++ei) {
		  edge_count += 1;
	  }
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
	  auto ei = edge_begin();
      for (size_type j = 0; j < i; j++) {
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
    if (std::find(g_edges[a.index()].begin(), g_edges[a.index()].end(), b.index()) != g_edges[a.index()].end()) {
		return true;
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
  Edge add_edge(const Node& a, const Node& b) {
    edge_type new_edge = Edge(this, a.index(), b.index());
    if (!has_edge(a, b)) {
	    g_edges[a.index()].push_back(b.index());
	    g_edges[b.index()].push_back(a.index());
	}
    return new_edge;
  }
  
  /** Return edge iterator of first edge in sequence
   * 
   * 
   * 
   */
  edge_iterator edge_begin() const {
	  return g_edge_begin;
  }
  
  /** Return edge iterator of on past last edge in sequence
   * 
   * 
   * 
   */
   edge_iterator edge_end() const {
	   return g_edge_end;
   }
   
   /** Return index of the second node refered to by edge iterator */
   size_type index_node2(const edge_iterator& eit) const {
	   return g_edges[eit.node1_index][eit.node2_pos];
   }
   
   /** Return number of edges connected to current node */
   size_type connectivity(const size_type& no_index) const {
	   return g_edges[no_index].size();
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    g_nodes.clear();
    g_edges.clear();
    g_values.clear();
  }

 private:
  std::vector<Point> g_nodes;
  std::vector<node_value_type> g_values;
  std::vector<std::vector<size_type>> g_edges;
  node_iterator g_node_begin;
  node_iterator g_node_end;
  edge_iterator g_edge_begin;
  edge_iterator g_edge_end;

};

#endif // CME212_GRAPH_HPP
