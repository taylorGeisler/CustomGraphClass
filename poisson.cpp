/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>


// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!

struct GraphValue {
	bool b_;
	double v_;
};

typedef Graph<GraphValue,char> GraphType;

struct GraphSymmetricMatrix {
	GraphSymmetricMatrix(const GraphType& g) : s_(g.size()), g_(&g) {
	}
	template <typename VectorIn , typename VectorOut , typename Assign > 
    void mult(const VectorIn& v, VectorOut& w, Assign) const {
		assert(size(v) == size(w));
		assert(size(v) == s_);
		
		// Iterate over all nodes
	    for (auto it = g_->node_begin(); it != g_->node_end(); ++it) { 
			double store = 0;
			for (auto it1 = g_->node_begin(); it1 != g_->node_end(); ++it1) {
				if (it == it1 and (*it).value().b_) {
					store += v[(*it1).index()];
				}
				else if (it != it1 and ((*it).value().b_ or (*it1).value().b_)) {
				}
				else {
					if (it == it1) {
						store -= v[(*it1).index()]*double((*it).degree());
					}
					else if (g_->has_edge(*it,*it1)) {
						store += v[(*it1).index()];
					}
					else {
					}
				}
			}
			Assign::apply(w[(*it).index()],store);
		}
    }
    
    /** Matvec forwards to MTL’s lazy mat_cvec_multiplier operator */ 
    template <typename Vector > 
    mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector> 
    operator*(const Vector& v) const {
        return {*this, v};
	}
    
    std::size_t s_;
	
	private:
	const GraphType* g_;
};

inline std::size_t size(const GraphSymmetricMatrix& A) {
	return A.s_*A.s_;
}

inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
	return A.s_;
}

inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
	return A.s_;
}

namespace mtl {
namespace ashape {

/** Define IdentityMatrix to be a non-scalar type. */
template<>
struct ashape_aux<GraphSymmetricMatrix> {
	typedef nonscal type;
};
}    //end namespace ashape

/** IdentityMatric implements the Collection concept
 * with value_type and size_type */
 template<>
 struct Collection <GraphSymmetricMatrix> {
	 typedef double value_type;
	 typedef unsigned size_type;
};
}  // end namespace mtl

/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  // HW3: YOUR CODE HERE
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
	  if (bb.contains((*it).position())) {
		  g.remove_node(it);
	  }
  }
  return;
}

double f_fun(const Point pnt) {
	  return 5*cos(norm_1(pnt));
}

double g_fun(const Point pnt) {
	if (norm_inf(pnt) == 1.0) {
		return 0;
	}
	else if (norm_inf(pnt - Point(0.6, 0.6, 0)) < 0.2 or norm_inf(pnt - Point(-0.6, -0.6, 0)) < 0.2 
	  or norm_inf(pnt - Point(-0.6, 0.6, 0)) < 0.2 or norm_inf(pnt - Point(0.6, -0.6, 0)) < 0.2) {
        return -0.2;
	}
	else if (Box3D(Point(-0.6,-0.2,-1), Point(0.6,0.2,1)).contains(pnt)) {
	  return 1;
    }
    return 0;
}

/// Class for iteration control that cyclically prints residual
  template <class Real, class OStream = std::ostream>
  class visual_iteration : public itl::basic_iteration<Real> 
  {
      typedef itl::basic_iteration<Real> super;
      typedef visual_iteration self;

      void print_resid()
      {
	  if (!this->my_quite && this->i % cycle == 0)
	      if (multi_print || this->i != last_print) { // Avoid multiple print-outs in same iteration
		  out << "iteration " << this->i << ": resid " << this->resid() 
		      // << " / " << this->norm_r0 << " = " << this->resid() / this->norm_r0 << " (rel. error)"
		      << std::endl;
		  last_print= this->i;
		  
		  //for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	          //(*it).value().v_ = x[(*it).index()];
          //}
  
		  //CME212::SDLViewer viewer;
		  //auto node_map = viewer.empty_node_map(graph);
		  
		  //viewer.clear();
		  //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
		  //viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
		
		  //viewer.center_view();
	      }
      }

    public:
  
      template <class Vector>
      visual_iteration(const Vector& r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100,
		       OStream& out = std::cout)
	: super(r0, max_iter_, tol_, atol_), cycle(cycle_), last_print(-1), multi_print(false), out(out)
      {}

      visual_iteration(Real r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100,
		       OStream& out = std::cout)
	: super(r0, max_iter_, tol_, atol_), cycle(cycle_), last_print(-1), multi_print(false), out(out)
      {}
      

      bool finished() { return super::finished(); }

      template <typename T>
      bool finished(const T& r) 
      {
	  bool ret= super::finished(r);
	  print_resid();
	  return ret;
      }

      inline self& operator++() { ++this->i; return *this; }
      
      inline self& operator+=(int n) { this->i+= n; return *this; }

      operator int() const { return error_code(); }

      /// Whether the residual is printed multiple times in iteration
      bool is_multi_print() const { return multi_print; }

      /// Set whether the residual is printed multiple times in iteration
      void set_multi_print(bool m) { multi_print= m; }

      int error_code() const 
      {
	  if (!this->my_suppress)
	      out << "finished! error code = " << this->error << '\n'
		  << this->iterations() << " iterations\n"
		  << this->resid() << " is actual final residual. \n"
		  << this->relresid() << " is actual relative tolerance achieved. \n"
		  << "Relative tol: " << this->rtol_ << "  Absolute tol: " << this->atol_ << '\n'
		  << "Convergence:  " << pow(this->relresid(), 1.0 / double(this->iterations())) << std::endl;
	  return this->error;
      }
    protected:
      int        cycle, last_print;
      bool       multi_print;
      OStream&   out;
  };


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }
  
  Point a;

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());

  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.
  
  // Mark Edges
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  if ((*it).degree() != 4) {
		  (*it).value().b_ = true;
	  } else if (Box3D(Point(-0.8,-0.8,-1), Point(-0.4,-0.4,1)).contains((*it).position()) or
	      Box3D(Point(0.4,-0.8,-1), Point(0.8,-0.4,1)).contains((*it).position()) or
	      Box3D(Point(-0.8,0.4,-1), Point(-0.4,0.8,1)).contains((*it).position()) or
	      Box3D(Point(0.4,0.4,-1), Point(0.8,0.8,1)).contains((*it).position()) or
	      Box3D(Point(-0.6,-0.2,-1), Point(0.6,0.2,1)).contains((*it).position())) {
			  (*it).value().b_ = true;
	  } else {
		  (*it).value().b_ = false;
	  }
  }
  
  // Make A
  GraphSymmetricMatrix A(graph);
  
  // Make x
  mtl::dense_vector<double> x(graph.size());
  for (auto it = x.begin(); it != x.end(); ++it)
      *it = 0;
  
  // Make b
  mtl::dense_vector<double> b(graph.size());
  for (auto it = b.begin(); it != b.end(); ++it)
      *it = 0;
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  if ((*it).value().b_) {
		 b[(*it).index()] = g_fun((*it).position());
	  } else {
		  b[(*it).index()] = h*h*f_fun((*it).position());
		  for (auto iit = (*it).edge_begin(); iit != (*it).edge_end(); ++iit) {
			  b[(*it).index()] -= g_fun((*iit).node2().position());
		  }
	  }
  }

  // Preconditioner
  itl::pc::identity<GraphSymmetricMatrix> P(A);
  
  // Iterator
  // Termination criterion: r < 1e-6 * b or N iterations
  //itl::cyclic_iteration<double> iter(b, 500, 1.e-10, 0, 50);
  visual_iteration<double> iter(b, 500, 1.e-10, 0, 50);
  
  // Solve for x
  itl::cg(A, x, b, P, iter);
  
   for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	  (*it).value().v_ = x[(*it).index()];
  }
  
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  return 0;
}
