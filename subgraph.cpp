/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */
#include <fstream>
#include <iterator>
//#include<queue>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
// Define our types
using GraphType = Graph<int>;
using NodeType = typename GraphType::node_type;
using NodeIter = typename GraphType::node_iterator;
using size_type = unsigned int;

template <typename Pred, typename It>
class filter_iterator : private equality_comparable<filter_iterator<Pred, It>>
{
public:
	// Get all of the iterator traits and make them our own
	using value_type = typename std::iterator_traits<It>::value_type;
	using pointer = typename std::iterator_traits<It>::pointer;
	using reference = typename std::iterator_traits<It>::reference;
	using difference_type = typename std::iterator_traits<It>::difference_type;
	using iterator_category = typename std::input_iterator_tag;
	using self_type = filter_iterator<Pred, It>;

	// Constructor
	filter_iterator(const Pred& p, const It& first, const It& last)
		: p_(p), it_(first), end_(last) {
		// HW1 #4: YOUR CODE HERE
		while ((!p_(*it_)) && (it_ != end_)) {
			++it_;
		}			

	}

	// HW1 #4: YOUR CODE HERE
	// Supply definitions AND SPECIFICATIONS for:
	value_type operator*() const {
		return *it_;
	}

	filter_iterator& operator++(){
		++it_;
		while ((!p_(*it_)) && (it_ != end_)) {
			++it_;
		}
		return *this;
    }

	bool operator==(const self_type& x) const {
		if ((it_ !=x.it_)||(end_!=x.end_)) {
			return false;
		}
		else {
			if (it_ != end_) {
				if ((p_(*it_)) && (p_(*(x.it_))))
					return true;
				else
					return false;
			}
			else
				return true;
		}

	}

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators. This deduces the type of
 * the predicate function and the iterator so the user doesn't have to write it.
 * This also allows the use of lambda functions as predicates.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

// HW1 #4: YOUR CODE HERE
// Specify and write an interesting predicate on the nodes.
// Explain what your predicate is intended to do and test it.
// If you'd like you may create new nodes and tets files.

/** Predicate to test whether the node is a isolated node (no edges) or not
 *  based on the nodes degree. 
 *  If the degree of a node (number of incident edges) is 0, then it is an isolated edge
 *
 *  Given a NODE @a node, the overloading operator ()returns a boolean to indicate whether 
 *  the node is an isolated node or not. If the node is isolated, the boolean returns false 
 *  and it is removed from the subgraph.
 *  It has template typename NODE, which takes any nodetype input argument
 *  **/
struct IsolatedPredicate {
	template <typename NODE>
	bool operator()(const NODE& node) const {
		if (node.degree() == 0)
			return false;
		else
			return true;
	}
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) const {
    return n.position().x < 0;
  }
};


/** combining subgraph and color **/
/*
//Euclidean Distance Comparator
struct EDComp {
	EDComp(Point P) :P_(P) {}
	bool operator()(const NodeType& n1, const NodeType& n2) {
		double dist_1 = norm(n1.position() - P_);
		double dist_2 = norm(n2.position() - P_);
		return dist_1 < dist_2;
	}
private:
	Point P_;
};

NodeIter nearest_node(const GraphType& g, const Point& point)
{
	// HW1 #3: YOUR CODE HERE

	if (g.num_nodes() == 0)
		return g.node_end();
	else {
		NodeIter nstnodeit = std::min_element(g.node_begin(), g.node_end(), EDComp(point));
		return nstnodeit;
	}

}
*/
/** Update a graph with the shortest path lengths from a root node.
* @param[in,out] g     Input graph
* @param[in,out] root  Root node to start the search.
* @return The maximum path length found.
*
* @post root.value() == 0
* @post Graph has modified node values indicating the minimum path length
*           to the root.
* @post Graph nodes that are unreachable from the root have value() == -1.
*
* This sets all nodes' value() to the length of the shortest path to
* the root node. The root's value() is 0. Nodes unreachable from
* the root have value() -1.
*/
/*
bool PathComp(const NodeType& n1, const NodeType& n2) {
	return n1.value() < n2.value();
}

int shortest_path_lengths(GraphType& g, NodeType& root)
{
	// HW1 #3: YOUR CODE HERE

	for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
		(*ni).value() = -1;//Assign the distance of all unvisited nodes to be -1
	}

	std::queue<size_type> Q; //queue of the index of next visited node

	Q.push(root.index());
	root.value() = 0;

	while (!Q.empty()) {
		size_type currid = Q.front();
		Q.pop();
		NodeType curr = g.node(currid);

		for (auto icr = curr.edge_begin(); icr != curr.edge_end(); ++icr) {
			if ((*icr).node2().value() == -1) {// if the node is unvisited
				(*icr).node2().value() = curr.value() + 1;
				Q.push((*icr).node2().index());
			}
		}
	}

	NodeIter furthest = std::max_element(g.node_begin(), g.node_end(), PathComp);
	int lengths = (*furthest).value();
	return lengths;
}
*/
int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }


  // Construct a Graph
  GraphType graph;
  std::vector<NodeType> nodes;

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


  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  // HW1 #4: YOUR CODE HERE
  // Use the filter_iterator to plot an induced subgraph.
  auto node_map = viewer.empty_node_map(graph);
  //IsolatedPredicate predic;
  SlicePredicate predic;
  auto start = make_filtered(graph.node_begin(), graph.node_end(), predic);
  auto end = make_filtered(graph.node_end(), graph.node_end(), predic);

  //add nodes
 viewer.add_nodes(start, end, node_map);
  
  //add edges

  for (auto ni = start; ni != end; ++ni) {
	  NodeType n = *ni;
	  viewer.add_edges(n.edge_begin(), n.edge_end(), node_map);
  }

/*
  struct colorfn {
	  colorfn(int shortest) :shortest_(shortest) {
	  }
	  CME212::Color operator()(NodeType n) {
		  float fraction = 1 - ((float)n.value() / (float)shortest_);
		  CME212::Color result = CME212::Color::make_heat(fraction);

		  return result;
	  }
  private:
	  int shortest_;
  };

  NodeType root = *(nearest_node(graph, Point(-1, 0, 1)));
  int shortest_lgth = shortest_path_lengths(graph, root);
  viewer.add_nodes(start, end, colorfn(shortest_lgth),node_map);

  for (auto ni = start; ni != end; ++ni) {
	  NodeType n = *ni;
	  viewer.add_edges(n.edge_begin(), n.edge_end(), node_map);
  }

*/  
  // Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
