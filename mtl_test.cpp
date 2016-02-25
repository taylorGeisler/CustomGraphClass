/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <math.h>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL

struct IdentityMatrix {
    IdentityMatrix(std::size_t s) : s_(s) {}
    /** Helper function to perform multiplication. Allows for delayed
     * evaluation of results.
     * Assign::apply(a, b) resolves to an assignment operation such as * a += b, a -= b, or a = b.
     * @pre @a size(v) == size(w) */
    template <typename VectorIn , typename VectorOut , typename Assign > 
    void mult(const VectorIn& v, VectorOut& w, Assign) const {
		assert(size(v) == size(w));
		assert(size(v) == s_);
		
	    for (std::size_t i = 0; i < s_; i++) {
			Assign::apply(w[i],v[i]);
		}
    }
	
    /** Matvec forwards to MTL’s lazy mat_cvec_multiplier operator */ 
    template <typename Vector > 
    mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector> 
    operator*(const Vector& v) const {
        return {*this, v};
	}
	
	std::size_t s_;
	
	private:
};

inline std::size_t size(const IdentityMatrix& A) {
	return A.s_*A.s_;
}

inline std::size_t num_rows(const IdentityMatrix& A) {
	return A.s_;
}

inline std::size_t num_cols(const IdentityMatrix& A) {
	return A.s_;
}

/** Traits that MTL uses to determine properties of our Identity Matrix. */
namespace mtl {
namespace ashape {

/** Define IdentityMatrix to be a non-scalar type. */
template<>
struct ashape_aux<IdentityMatrix> {
	typedef nonscal type;
};
}    //end namespace ashape

/** IdentityMatric implements the Collection concept
 * with value_type and size_type */
 template<>
 struct Collection <IdentityMatrix> {
	 typedef double value_type;
	 typedef unsigned size_type;
};
}  // end namespace mtl

int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver
  
  const std::size_t size = 10;
  
  IdentityMatrix A(size);

  mtl::dense_vector<double> x(size), b(size);
  std::cout << A.s_<<std::endl;
  iota(b);
  
  // Termination criterion: r < 1e-6 * b or N iterations
  itl::noisy_iteration<double> iter(b, 500, 1.e-6);
  
  // Create an ILU(0) preconditioner
  itl::pc::identity<IdentityMatrix> P(A);
  
  itl::cg(A, x, b, P, iter);
  
  std::cout << "x is " << x << std::endl;
 
  return 0;
}
