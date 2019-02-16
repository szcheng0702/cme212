/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL


class IdentityMatrix{
public:
	IdentityMatrix(int s1) :s(s1) {}

	/** Helper function to perform multiplication. Allows for delayed
	 *  evaluation of results.
	 *  Assign:: apply(a,b) resolves to an assignment operator such as
	    a+=b, a-=b, or a=b
	 *  @pre @a size(v)==size(w) */
	template<typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v, VectorOut& w, Assign) const {
		assert(int(size(v))==s);
		assert(size(v)==size(w));//The VectorIn is exactly VectorOut

		for(int i=0;i<s;i++){
			Assign::apply(w[i],v[i]);
			
		}

	}
	
	/** Matvec forwards to MTL's lazy mat_cvec_multiplier operator */
	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector> operator * (const Vector& v) const {
		return { *this, v };
	}

	int s;
};

/** The number of elements in the matrix */
inline std::size_t size(const IdentityMatrix& A) {
	return A.s*A.s;
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const IdentityMatrix& A) {
	return A.s;
}

/** The number of columns in the matrix */
inline std::size_t num_cols(const IdentityMatrix &A) {
	return A.s;
}


/** Traits that MTL uses to determine properties of our IdentityMatrix. */
namespace mtl {
	namespace ashape {
		/** Define IdentityMatrix to be a non-scalar type. */
		template<>
		struct ashape_aux<IdentityMatrix> {
			typedef nonscal type;
		};
	}// end namespace ashape

	 /** IdentityMatrix implements the Collection concept
	 *  with value_type and size_type */
	template<>
	struct Collection<IdentityMatrix> {
		typedef double value_type;
		typedef unsigned size_type;
	};
}//end namespace mtl

int main()
{
  // HW3: YOUR CODE HERE
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

	using namespace mtl;
	using namespace itl;

	const int N = 10;

	IdentityMatrix I(10);

	pc::identity<IdentityMatrix> P(I);

	//set b such that x==1 is solution; start with x==0
	dense_vector<double> x(N, 1.0), b(N);
	b = I * x; x = 0;

	//Termination criterion: r<1e-6* b or N iterations
	noisy_iteration<double> iter(b, 500, 1.e-6);

	cg(I, x, b, P, iter);


  return 0;
}
