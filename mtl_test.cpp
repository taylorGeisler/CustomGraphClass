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
	/** Conpute the product of a vector with this identity matrix
	 */
    template <typename Vector>
    Vector operator*(const Vector& x) const {
		return x;
	}
	
	int s_;
	
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

  return 0;
}
