/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"



// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel = Point(0,0,0);       //< Node velocity
  double mass;     //< Node mass
};

struct EdgeData {
	double L_init;
	double K;
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData,EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

struct con_wall {
	  void operator() (GraphType& g, double t) {
		  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			  if (dot((*it).position(),Point(0,0,1)) < zmax) {
				  (*it).position()[2] = zmax;
				  (*it).value().vel[2] = 0;
			  }
		  }
		  (void)t;
	  }
	  double zmax = -0.75;
};

struct con_sphere {
	  void operator() (GraphType& g, double t) {
		  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			  double a = norm((*it).position()-center);
			  if (a < radius) {
				  Point dir = (*it).position() - center;
				  dir = dir/norm(dir);
				  (*it).position() = center + radius*dir;
				  (*it).value().vel = (*it).value().vel - dot((*it).value().vel,dir)*dir;
			  }
		  }
		  (void)t;
	  }
	  Point center = Point(0.5,0.5,-0.5);
	  double radius = 0.15;
};

struct con_sphere_rem {
	  void operator() (GraphType& g, double t) {
		  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
			  double a = norm((*it).position()-center);
			  if (a < radius) {
			  }
		  }
		  (void)t;
	  }
	  Point center = Point(0.5,0.5,-0.5);
	  double radius = 0.15;
};

template <typename C1, typename C2>
struct CombinedConstraint {
	C1 c1;
	C2 c2;
	void operator()(GraphType& g, double t) {
		c1(g,t);
		c2(g,t);
	}
};

template <typename C1, typename C2>
CombinedConstraint<C1,C2> make_combined_constraint(C1 constraint1, C2 constraint2) {
	CombinedConstraint<C1,C2> constraint;
	constraint.c1 = constraint1;
	constraint.c2 = constraint2;
	return constraint;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)) {	
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
	}
  }
  //con_sphere hi;
  //hi(g,t);
  CombinedConstraint<con_wall,con_sphere> ZZZ = make_combined_constraint(con_wall(), con_sphere());
  ZZZ(g, t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}



template <typename F1, typename F2>
struct CombinedForce {
	F1 f1;
	F2 f2;
	template <typename NODE>
	Point operator()(NODE n, double t) {
		return f1(n,t) + f2(n,t);
	}
};

template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1 force1, F2 force2) {
	CombinedForce<F1,F2> force;
	force.f1 = force1;
	force.f2 = force2;
	return force;
}

template <typename F1, typename F2, typename F3>
CombinedForce<F1,CombinedForce<F2,F3>> make_combined_force(F1 force1, F2 force2, F3 force3) {
	CombinedForce<F2,F3> force23;
	CombinedForce<F1,CombinedForce<F2,F3>> force;
	force23.f1 = force2;
	force23.f2 = force3;
	force.f1 = force1;
	force.f2 = force23;
	return force;
}

struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t.
   *
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
		//Gravity Force
		double m = n.value().mass;
		return Point(0,0,-grav)*m;
		(void)t;
	}
};

struct MassSpringForce {
  /** Return the mass spring force applying to @a n at time @a t.
   *
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
	    Point F = Point(0,0,0);
		//Spring Forces
		//Iterate over each edge
		for (auto eit = n.edge_begin(); !(eit == n.edge_end()); ++eit ) {
			Point F_unit_dir = ((*eit).node1().position() - (*eit).node2().position())/(*eit).length();
			double F_mag = (*eit).value().K*((*eit).value().L_init-(*eit).length());
			F += F_unit_dir*F_mag;
			(void)t;
		}
		return F;
	}
};

struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   *
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
		//Damping force
		return n.value().vel*-c_;
		(void)t;
	}
	
    double c_ = 1/100;
};


int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));
  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);

    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  (*it).value().mass = 1.0/double(graph.size());
  }
  
  for (auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit) {
	  (*eit).value().K = 100;
	  (*eit).value().L_init = (*eit).length();
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    //damp_force.c = 1/graph.size();
    //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce()));
    symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce()));
    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CME212::sleep(0.001);
  }
  return 0;
}
