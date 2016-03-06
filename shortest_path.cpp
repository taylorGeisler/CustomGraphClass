/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

#include <limits>
#include <queue>


/** Comparator that compares the distance from a given point p.
 * 
 * Not employed in current implementation
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };
   
   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    (void) node1; (void) node2;    // Quiet compiler warning
    return false;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int,int>& g, const Point& point) {

  // Calculate root node
  int init_val = -1; // Initialize values to large number
  float closest_dist = std::numeric_limits<float>::max();
  Graph<int>::size_type root;
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
	  float dist = norm((*ni).position()-point);
	  Graph<int>::size_type node_index = (*ni).index();
	  g.node(node_index).value() = init_val;
	  if (dist < closest_dist) {
		  root = (*ni).index();
		  closest_dist = dist;
	  }
  }
  
  std::cout<< root <<std::endl;

  // Breadth-First Search
  std::vector<Graph<int>::size_type> parent (g.size(),g.size());
  
  std::queue<Graph<int>::size_type> Q;
  g.node(root).value() = 0;
  Q.push(root);
  
  while (!(Q.empty())) {
	  Graph<int>::size_type current = Q.front();
	  Q.pop();
	  
	  for (auto eit = g.node(current).edge_begin(); !(eit == g.node(current).edge_end()); ++eit ) {
		  if ((*eit).node2().value() == -1) {
			  (*eit).node2().value() = (*eit).node1().value() + 1;
			  Q.push((*eit).node2().index());
		  }
	  }
  }
  
  // Calculate longest path
  int longest_path = 0;
  for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
	  Graph<int>::node_type ni_node = *ni;
	  if (ni_node.value() > longest_path) {
		  longest_path = ni_node.value();
	  }
  }
  return longest_path;
}



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  viewer.launch();

  // Use shortest_path_lengths to set the node values to the path lengths
  Point root_point = Point(-1,0,1);
  int longest_path = shortest_path_lengths(graph, root_point);
  std::cout << "longest path  ";
  std::cout << longest_path << std::endl;
  // Construct a Color functor and view with the SDLViewer
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), [longest_path](Graph<int>::node_type n){ float c = (float)n.value()/(float)(longest_path+1);
	  if (c < 0.0) {c = 0.0;};return CME212::Color::make_heat(-(c-1));}, node_map);
	  
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  return 0;
}
