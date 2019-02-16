#include <vector>
#include <fstream>
#include <queue>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int>;
using NodeType = typename GraphType::node_type;
using NodeIter = typename GraphType::node_iterator;
using InciIter = typename GraphType::incident_iterator;
using size_type = unsigned int;

NodeIter nearest_node(const GraphType& g, const Point& point);

