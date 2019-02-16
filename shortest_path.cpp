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
#include <queue>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;
using InciIter = typename GraphType::incident_iterator;
using size_type = unsigned int;

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */

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

	if(g.num_nodes()==0)
		return g.node_end();
	else {
		NodeIter nstnodeit = std::min_element(g.node_begin(), g.node_end(), EDComp(point));
		return nstnodeit;
	}
 
}

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
				(*icr).node2().value() = curr.value() +1;
				Q.push((*icr).node2().index());
			}
		}
	}

	NodeIter furthest = std::max_element(g.node_begin(), g.node_end(), PathComp);
	int lengths = (*furthest).value();
	return lengths;
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
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

  // Launch the SFML_Viewer
  CME212::SFML_Viewer viewer;

  // HW1 #3: YOUR CODE HERE
  // Use nearest_node and shortest_path_lengths to set the node values
  // Construct a Color functor and view with the SFML_Viewer

  struct colorfn {
	  colorfn(int shortest) :shortest_(shortest) {
	  }
	  CME212::Color operator()(NodeType n) {
		  float fraction = 1.0-((float)n.value()/(float)shortest_);
		  CME212::Color result = CME212::Color::make_heat(fraction);
		  
		  return result;
	  }
  private:
	  int shortest_;
  };
  NodeType root = *(nearest_node(graph, Point(-1, 0, 1)));
  int shortest_lgth = shortest_path_lengths(graph, root);

  //Draw Nodes
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), colorfn(shortest_lgth),node_map);

  //Draw Edges
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  
// Center the view and enter the event loop for interactivity
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
