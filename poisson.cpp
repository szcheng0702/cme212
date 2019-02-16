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
 * Launches an SFML_Viewer to visualize the solution.
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"


// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
using GraphType = Graph<char,char>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;

/** @class GraphSymmetricMatrix
* @brief A template for a matrix to solve a Poisson Problem for a graph.
*
* Users can solve the Poisson Problem, given the graph and 
* the boundary condition.
*/
class GraphSymmetricMatrix {
public:
	GraphSymmetricMatrix(GraphType g1,std::vector<bool> boundary) :s(g1.num_nodes()),g(g1),btag(boundary){}

	/** Helper function to perform multiplication. Allows for delayed
	*  evaluation of results.
	*  Assign:: apply(a,b) resolves to an assignment operator such as
	a+=b, a-=b, or a=b
	*  @pre @a size(v)==size(w) 
	*  @pre @a size(v)==@g.num_nodes()
	*  @pre @a btag.size()=size(v) */
	template<typename VectorIn, typename VectorOut, typename Assign>
	void mult(const VectorIn& v, VectorOut& w, Assign) const {
		
		assert(size(v)== g.num_nodes());
		assert(size(v) == size(w));
        

		for (unsigned int i = 0; i<s; i++) {
			if(btag[i]){
				Assign::apply(w[i], v[i]);
				//when i is on the boundary
			}
			else {
				//if i is not on the boundary, Aii=-deg(ni)
				auto update = -(double)g.node(i).degree()*v[i];
				auto n=g.node(i);
				for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
					auto n2id = (*it).node2().index();
					if (!btag[n2id])//if not in the boundary
						update += v[n2id];//Aij=1

				}
				Assign::apply(w[i], update);
			}

		}

	}

	/** Matvec forwards to MTL's lazy mat_cvec_multiplier operator */
	template <typename Vector>
	mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector> operator * (const Vector& v) const {
		return { *this, v };
	}

	unsigned int s;

private:
	GraphType g;
	std::vector<bool> btag;//vector of booleans to tag boundary
};

/** The number of elements in the matrix */
inline std::size_t size(const GraphSymmetricMatrix& A) {
	return A.s*A.s;
}

/** The number of rows in the matrix */
inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
	return A.s;
}

/** The number of columns in the matrix */
inline std::size_t num_cols(const GraphSymmetricMatrix &A) {
	return A.s;
}


/** Traits that MTL uses to determine properties of our GraphSymmetricMatrix. */
namespace mtl {
	namespace ashape {
		/** Define GraphSymmetricMatrix to be a non-scalar type. */
		template<>
struct ashape_aux<GraphSymmetricMatrix> {
	typedef nonscal type;
};
	}// end namespace ashape

	 /** GraphSymmetricMatrix implements the Collection concept
	 *  with value_type and size_type */
	template<>
	struct Collection<GraphSymmetricMatrix> {
		typedef double value_type;
		typedef unsigned size_type;
	};
}//end namespace mtl


/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
* @param[in,out] g  The Graph to remove nodes from
* @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
* @post For all i, 0 <= i < @a g.num_nodes(),
*        not bb.contains(g.node(i).position())
*/
void remove_box(GraphType& g, const Box3D& bb) {
	// HW3: YOUR CODE HERE
	//(void) g; (void) bb;   //< Quiet compiler
	for (unsigned i = 0; i < g.num_nodes(); i++) {
		if (bb.contains(g.node(i).position())) {
			g.remove_node(g.node(i));
		}
	}
	return;
}


struct NodeColor {
	/** Constructor for the node color, given a @a mtl::dense_vector<double> u_ 
	 *  Inside the constructor, we first find the element inside @a u_
	 *  with maximum absolute value*/
	NodeColor(mtl::dense_vector<double> u) {
		u_ = u;
		double maxelem = std::abs(u_[0]);
		for (unsigned i = 1; i < size(u_); ++i) {
			if (std::abs(u[i]) > maxelem)
				maxelem = std::abs(u[i]);
		}
		maximum = maxelem;
	}

	/** Return the color applying to Node @a n 
	 *  Each node is assocaited with a fraction, which is the value 
	 *  obtained by the absolute value of u_[n.index()] divided by the 
	    maximum absolute value*/
	CME212::Color operator()(NodeType n) {
		unsigned id=n.index();
		double fraction;
		CME212::Color result;
		if (maximum!=0) {
			fraction = std::abs(u_[id])/ maximum;
			result = CME212::Color::make_heat(fraction);
	    }
	    else{
	    	/* If the denominator is 0 
			* we introduce a small perturbation and do the similar handling as above */
	    	fraction = (std::abs(u_[id])+0.01)/ (maximum+0.01);
			result = CME212::Color::make_heat(fraction);
	    }
		return result;
	}

private:
	 mtl::dense_vector<double> u_;
	 double maximum;
};

struct NodePosition {
	NodePosition(mtl::dense_vector<double> u): u_(u) {

	}

	/** Return the new position to Node @a n
	*  x, y coordinates of the new position is set to its current value
	*  while the z coordinate of the position that we are going to plot is
	*  @a u_[n.index()] */
	Point operator()(NodeType n) {
		unsigned id=n.index();
		Point p=n.position();
		return Point(p.x,p.y,u_[id]);
	}
private:
 mtl::dense_vector<double> u_;

};

namespace itl {
	/** @class itl::visual_iteration
	* @brief Class which ouputs the residual periodically, and plot the current solution
	* at every iteration
	*
	* visual_iteration objects are used in linear solvers.
	*/
	template<class GraphType, class Vector, class Real, class Ostream=std::ostream>
	class visual_iteration :public cyclic_iteration<Real,Ostream>
	{

		typedef cyclic_iteration<Real,Ostream> super;
	
		/* Launch the solver*/
		CME212::SFML_Viewer viewer;

		/*initialize node_map */
		std::map<typename GraphType::node_type,unsigned>node_map=viewer.empty_node_map(g);

		/* Plot the current solution inside the SFML_viewer solver*/
		void plot_sol() {

			//Clear
			viewer.clear();
			node_map.clear();

			//Draw Nodes 
			viewer.add_nodes(g.node_begin(), g.node_end(), NodeColor(*u), NodePosition(*u), node_map);


			//Draw Edges
			viewer.add_edges(g.edge_begin(), g.edge_end(), node_map);


			// Center the view and enter the event loop for interactivity 
			viewer.center_view();
			//viewer.event_loop();

			
		}

	public:
		/* Constructor for visual_iteration*/
		visual_iteration(const GraphType& g_, Vector* u_, const Vector& r0, int max_iter_, Real tol_, Real atol_=Real(0),
			             int cycle_=100, Ostream& out=std::cout)
			:super(r0,max_iter_,tol_,atol_,cycle_,out),g(g_), u(u_){}
		
		/* Return a boolean to indicate whether the solver stops or not
		 * based on the same convergence criterion as cyclic_iteration */
		bool finished() {
			return super::finished();
		}


		/* Return a boolean to indicate whether the solver stops or not
		* When the boolean is false, keep plot the current solution in SFML_viewer*/
		template <typename T>
		bool finished(const T& r) {
			bool ret=super::finished(r);
			plot_sol();
			if(ret){
				viewer.event_loop();
			}
			return ret;
		}



	protected:
		GraphType g;
		Vector* u; //pass u by pointer so that u can keep updated at each iteration
		
	};



}//end namespace itl


int main(int argc, char** argv)
{
	// Check arguments
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
		exit(1);
	}

	// Define an empty Graph
	GraphType graph;

	{
		// Create a nodes_file from the first input argument
		std::ifstream nodes_file(argv[1]);
		// Interpret each line of the nodes_file as a 3D Point and add to the Graph
		std::vector<NodeType> node_vec;
		Point p;
		while (CME212::getline_parsed(nodes_file, p))
			node_vec.push_back(graph.add_node(2 * p - Point(1, 1, 0)));

		// Create a tets_file from the second input argument
		std::ifstream tets_file(argv[2]);
		// Interpret each line of the tets_file as four ints which refer to nodes
		std::array<int, 4> t;
		while (CME212::getline_parsed(tets_file, t)) {
			graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
			graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
			graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
			graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
		}
	}

	// Get the edge length, should be the same for each edge
	auto it = graph.edge_begin();
	assert(it != graph.edge_end());
	double h = norm((*it).node1().position() - (*it).node2().position());

	// Make holes in our Graph
	remove_box(graph, Box3D(Point(-0.8 + h, -0.8 + h, -1), Point(-0.4 - h, -0.4 - h, 1)));
	remove_box(graph, Box3D(Point(0.4 + h, -0.8 + h, -1), Point(0.8 - h, -0.4 - h, 1)));
	remove_box(graph, Box3D(Point(-0.8 + h, 0.4 + h, -1), Point(-0.4 - h, 0.8 - h, 1)));
	remove_box(graph, Box3D(Point(0.4 + h, 0.4 + h, -1), Point(0.8 - h, 0.8 - h, 1)));
	remove_box(graph, Box3D(Point(-0.6 + h, -0.2 + h, -1), Point(0.6 - h, 0.2 - h, 1)));

	// HW3: YOUR CODE HERE
	// Define b using the graph, f, and g.
	// Construct the GraphSymmetricMatrix A using the graph
	// Solve Au = b using MTL.

	//vector to tag whether the node is a boundary node or not
	std::vector<bool>boundary(graph.num_nodes());

	//Define the boundary bounding box
	Box3D Box=Box3D(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));

	//Define b
	mtl::dense_vector<double> b(graph.num_nodes(), 0.0);

	//Fill in boundary points first
	for (unsigned int i = 0; i < graph.num_nodes(); ++i){
		auto n = graph.node(i);
		Point d1 = n.position() - Point(0.6, 0.6, 0);
		Point d2 = n.position() - Point(-0.6, -0.6, 0);
		Point d3 = n.position() - Point(0.6, -0.6, 0);
		Point d4 = n.position() - Point(-0.6, 0.6, 0);
		if (norm_inf(n.position()) == 1) {
			b[i] = 0;
			boundary[i] = 1;
		}
		else if ((norm_inf(d1) < 0.2) || (norm_inf(d2) < 0.2) || (norm_inf(d3) < 0.2) || (norm_inf(d4) < 0.2)) {
			b[i] = -0.2;
			boundary[i] = 1;
		}
		else if (Box.contains(n.position())) {
			b[i] = 1;
			boundary[i] = 1;
		}
		else
			boundary[i] = 0;
	}

	//Fill in points that are not in the boundary
	for (unsigned int i = 0; i < graph.num_nodes(); ++i) {
		if (!boundary[i]) {
			auto n = graph.node(i);
			double f = 5 * cos(norm_1(n.position()));
			double g_sum = 0;
			for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
				unsigned j = (*it).node2().index();
				if (boundary[j])
					g_sum += b[j];
			}
			b[i] = h * h*f - g_sum;

		}
	}


	GraphSymmetricMatrix A(graph,boundary);

	itl::pc::identity<GraphSymmetricMatrix> P(A);

	//start with x==0
	mtl::dense_vector<double> u(graph.num_nodes(), 0.0);

	auto uptr=&u;

	//Termination criterion: r<1e-6* b or N iterations
	itl::visual_iteration<GraphType,mtl::dense_vector<double>,double> 
		iter(graph,uptr,b, 1000, 1.e-10,0.0,50);

	cg(A, u, b, P, iter);

    return 0;
}
