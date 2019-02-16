#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V,typename E>
class Graph {
 private:
	 
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
	 unsigned size_;
	 unsigned num_edges_;

	 //internal struct
	 struct inciedges;

	 struct proxynodes;

	 struct proxyedges;

	 //STL containers
	 std::vector <proxynodes> Nvec;
	 std::vector <proxyedges> Evec;
     std::vector <unsigned int> i2u;//map of id to the unique id of node, which maps to the exact location in the Nvec
	 std::vector<unsigned int> ei2eu;//similar as above, but is a map between edge id to the unique edge id



 public:
//
  // PUBLIC TYPE DEFINITIONS
  //

  
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Type of node values.*/
  typedef V node_value_type;
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
	  size_ = 0;
	  num_edges_ = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node:private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
		return graph_->Nvec[uindex_].P;
		assert(false);
    }

	/** Return the reference of the position so that we can modify it **/
	Point& position() {
		return graph_->Nvec[uindex_].P;
	}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
		return graph_->Nvec[uindex_].i;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	node_value_type& value() {
		return graph_->Nvec[uindex_].val;
	}

	const node_value_type& value() const {
		return graph_->Nvec[uindex_].val;
	}


	size_type degree() const {
		return graph_->Nvec[uindex_].adjN.size();
	}
    
	incident_iterator edge_begin() const {
		return IncidentIterator(graph_, this, 0);
	}

	incident_iterator edge_end() const {
		return IncidentIterator(graph_,this,this->degree());
	}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning

	  if ((graph_==n.graph_)&& (uindex_ == n.uindex_))
		  return true;
	  else
          return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
    
	  if (uindex_ < n.uindex_)
		  return true;
	  else
		  return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	graph_type* graph_;
	size_type uindex_;

	//private constructor
	Node(const graph_type* graph, size_type id)
		:graph_(const_cast<graph_type*>(graph)), uindex_(id) {

	}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size_;
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nodeval =node_value_type()) {
    // HW0: YOUR CODE HERE
	  std::vector<inciedges> adjN{};
	  Nvec.emplace_back(size_,position,nodeval,adjN);
	  i2u.push_back(Nvec.size()-1);
	  ++size_;
	  //Now the now node index is just size_-1
  
    return Node(this,Nvec.size()-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
	  Node N1 = node(n.index());
	  if (N1==n)
		  return true;
	  else
		  return false;
   // (void) n;            // Quiet compiler warning
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
	if (i < size()) {
		return Node(this, i2u[i]);
	}
    else
		return Node();        // Invalid node
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
		return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
		return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
		if ((node1() == e.node1()) && (node2() == e.node2()))
			return true;
		else if ((node2() == e.node1()) && (node1() == e.node2()))
			return true;
		else
			return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
		if (euid_ < e.euid_)
			return true;
		else if (node1() < e.node1())
			return true;
		else if (graph_!=e.graph_)
			return true;
		else
			return false;
    }

	/* Length of the edge, set it to be the Euclidean distance of the 1st Edge*/
	double length() const{
		Point p1 = node1().position();
		Point p2 = node2().position();
		return norm(p1 - p2);
	}
	
	edge_value_type& value() {
		return graph_->Evec[euid_].eval;
	}

	const edge_value_type& value() const {
		return graph_->Evec[euid_].eval;
	}
	
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	graph_type* graph_;
	size_type euid_;
	size_type uid1_;
	size_type uid2_;


	//Private Constructor
	Edge(const graph_type* graph, size_type euid, size_type uid1, size_type uid2)
		:graph_(const_cast<graph_type*>(graph)), euid_(euid),uid1_(uid1), uid2_(uid2) {

	}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
	  return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
	  if (i < num_edges_) {
		  size_type uniqueid = ei2eu[i];
		  return Edge(this, uniqueid, Evec[uniqueid].id1, Evec[uniqueid].id2);
	  }
	  else//if i exceeds the num of edges
		  return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	size_type Aid = a.uindex_;
	for (unsigned i = 0; i < Nvec[Aid].adjN.size(); ++i) {
		if (Nvec[Aid].adjN[i].id2 == b.uindex_)
			return true;
	}
	return false;

  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edgeval = edge_value_type()) {
    // HW0: YOUR CODE HERE
    //(void) a, (void) b;   // Quiet compiler warning
	  if (!(a==b)) {
		  if (!has_edge(a, b)) {
			  Evec.emplace_back(num_edges_,a.uindex_, b.uindex_, edgeval);
			  ei2eu.push_back(Evec.size()-1);//Add the map from ei to the unique eid
			  Nvec[a.uindex_].adjN.emplace_back(Evec.size()-1, b.uindex_);//Add b as the adjacent node of a
			  Nvec[b.uindex_].adjN.emplace_back(Evec.size()-1, a.uindex_);//Add a as the adjacent node of b
			  num_edges_++;
			  return Edge(this, Evec.size()-1,a.uindex_, b.uindex_);
		  }
		  else
			  return Edge(this, Evec.size()-1, a.uindex_, b.uindex_);
	  }
	else
	 	return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
	  size_ = 0;
	  num_edges_ = 0;
	  Nvec.clear();
	  Evec.clear();
	  i2u.clear();
	  ei2eu.clear();

  }

  //
  // Node Iterator
  //

  struct get_node{
    get_node(graph_type* g):graph_(const_cast<graph_type*>(g)){}
    Node operator()(size_type i) const{
      return graph_->node(i);
    }
  private:
    const graph_type* graph_;
  };

  typedef thrust::counting_iterator<size_type> CountingIter;

  /** @class Graph::NodeIterator
  * @brief Iterator class for nodes. A forward iterator. 
  * This class is inherited from the transform_iterator constructed by counting_iterator
  * and get_node functor defined above*/
  struct NodeIterator:thrust::transform_iterator<get_node,CountingIter, Node>{
    using super_t=thrust::transform_iterator<get_node,CountingIter,Node>;
    NodeIterator(){}//Default (invalid) constructor
    NodeIterator(const super_t& ti): super_t{ti}{}//KLUDGE Conversion Constructor

  private:
    friend class Graph;
    NodeIterator(size_type i, const graph_type* g)
       :super_t{thrust::make_transform_iterator(thrust::make_counting_iterator(i),
        get_node(const_cast<graph_type*>(g)))}{}
  };

  node_iterator node_begin() const{
    return NodeIterator(0,this);
  }

  node_iterator node_end() const{
    return NodeIterator(this->num_nodes(),this);
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	Edge operator*() const {
		size_type currid = node_->uindex_;
		return Edge(graph_, graph_->Nvec[currid].adjN[iid_].eid, currid,graph_->Nvec[currid].adjN[iid_].id2);
	}
    
	IncidentIterator& operator++() {
		iid_++;
		return *this;
	}

	bool operator==(const IncidentIterator& x) const {
		if ((graph_ == x.graph_) && (node_ == x.node_) && (iid_ == x.iid_))
			return true;
		else
			return false;
	}

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
	graph_type* graph_;
	node_type* node_;
	size_type iid_;

	//private constructor
	IncidentIterator(const graph_type* graph,const node_type* node,size_type iid) :
		graph_(const_cast<graph_type*>(graph)),node_(const_cast<node_type*>(node)),iid_(iid) {}

  };

  //
  // Edge Iterator
  //

  
  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	
	Edge operator*() const {
		return graph_->edge(eid_);
	}
    
	EdgeIterator& operator++() {
		eid_++;
		return *this;
	}

	bool operator==(const EdgeIterator& x) const {
		if ((graph_ == x.graph_) && (eid_ == x.eid_))
			return true;
		else 
			return false;
	}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	graph_type* graph_;
	size_type eid_;

	//private constructor
	EdgeIterator(const graph_type* graph, size_type eid) :
		graph_(const_cast<graph_type*>(graph)), eid_(eid) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
	  return EdgeIterator(this, 0);
  }
  
  edge_iterator edge_end() const {
	  return EdgeIterator(this,this->num_edges());

  };


  /** Remove an edge from the graph, and return 1 if one of the edges is removed, 0 if no edge is removed
  * @pre @a n1 @a n2 are valid nodes of this graph
  * @return 1 if @pre has_edge(n1,n2), 0 otherwise
  * @post has_edge(n1,n2) == false
  * @post If the edge connecting n1 and n2 is removed, new num_edges() == old num_edges()-1.
  *       Else,                                        new num_edges() == old num_edges()
  *
  * Can invalidate vector ei2eu indexes and vector adjN indexes for this @a n1 @a n2 
  * -- in other words,
  * old ei2eu indexes may not be equal new indexes and the same for adjN.
  * Must not invalidate the indexes for vector Evec.
  *
  * Complexity: No more than O(num_nodes()+num_edges())
  */

  size_type remove_edge(const Node& n1, const Node& n2) {
	  if (has_edge(n1, n2)) {
		  for (size_type i = 0; i < n1.degree(); ++i) {
			  if (Nvec[n1.uindex_].adjN[i].id2 == n2.uindex_) {

				  size_type uniqeid = Nvec[n1.uindex_].adjN[i].eid;
				  size_type realeid = Evec[uniqeid].ei;

				  //change the adjcent nodes
				  //For the adjacency of n1
				  Nvec[n1.uindex_].adjN.erase(Nvec[n1.uindex_].adjN.begin() + i);
				  
				  //For the adjacency of n2
				  for (size_type j = 0; j < n2.degree(); ++j) {
					  if (Nvec[n2.uindex_].adjN[j].id2 == n1.uindex_) {
						  Nvec[n2.uindex_].adjN.erase(Nvec[n2.uindex_].adjN.begin() + j);
					  }
				  }

				  ei2eu.erase(ei2eu.begin() + realeid);//remove the mapping

				  //change the realid
				  for (size_type k = realeid; k < ei2eu.size();++k) {
					  Evec[ei2eu[k]].ei = k;

				  }
				  num_edges_ -= 1;
				  return 1;
			  }
		  }
	  }

	return 0;
  }

  /** Remove an edge from the graph, and return 1 if one of the edges is removed, 0 if no edge is removed
  * @pre @a e is a valid edge of this graph
  * @post has_edge(e.node1(),e.node2()) == false
  * @post If @a e is removed, new num_edges() == old num_edges()-1.
  *       Else,               new num_edges() == old num_edges()
  *
  * Can invalidate vector ei2eu indexes and vector adjN indexes for this @a e.node1() and @a e.node2()
  * -- in other words,
  * old ei2eu indexes may not be equal new indexes and the same for adjN of @a e.node1) and @a e.node2().
  * Must not invalidate the indexes for vector Evec.
  *
  * Complexity: No more than O(num_nodes()+num_edges())
  */
  size_type remove_edge(const Edge& e) {
	  Node N1 = e.node1();
	  Node N2 = e.node2();
	  size_type num = remove_edge(N1, N2);
	  return num;
  }

  /** Remove an edge from the graph, and return the new iterator followed by the erased edge
  * @pre @a e_it is a edge iterator
  * @post e_it is not a valid edge iterator
  * @post If @a (*e_it) is removed, new num_edges() == old num_edges()-1.
  *       Else,                     new num_edges() == old num_edges()
  *
  * Can invalidate vector ei2eu indexes and vector adjN indexes for this @a (*e_it).node1() and @a (*e_it).node2()
  * -- in other words,
  * old ei2eu indexes may not be equal new indexes and the same for adjN of @a (*e_it).node1) and @a (*e_it).node2().
  * Must not invalidate the indexes for vector Evec.
  *
  * Complexity: No more than O(num_nodes()+num_edges())
  */

  edge_iterator remove_edge(edge_iterator e_it) {
	  auto e = *(e_it);
	  remove_edge(e);
	  edge_iterator nextit = e_it;
	  return nextit;
  }

  /** Remove an node from the graph, and return 1 if one of the nodes is removed, 0 if no node is removed
  * @pre @a n is a node
  * @return 1 if n is a valid node of this graph, 0 if n is not a valid node of this graph
  * @post has_node(n) == false
  * @post If @a n is removed, new num_nodes() == old num_nodes()-1.
  *       Else,               new num_nodes() == old num_nodess()
  *
  * Can invalidate vector i2u indexes and vector adjN for this @a n -- in other words, 
  * old i2u indexes may not be equal new indexes and adjN will become empty. 
  * Must not invalidate the indexes for vector Nvec.
  *
  * Complexity: No more than O(num_nodes())
  */

  size_type remove_node(const Node& n) {
	  if (has_node(n)) {
		  //Remove all edges incident to it

		  while (n.degree() > 0) {
			  auto e = *(n.edge_begin());
			  remove_edge(e);
		  }

		  i2u.erase(i2u.begin() + n.index());

		  for (size_type k = n.index(); k < i2u.size(); ++k) {
			  Nvec[i2u[k]].i = k;

		  }
		  size_ -= 1;
		  return 1;
	  }
	  else
		  return 0;
  }

  /** Remove an node from the graph, and return the new location followed by the erased element
  * @pre @a n_it is a node iterator
  * @return the new iterator followed by the erase element
  * @post n_it is not a valid node iterator
  * @post If @a n is removed, new num_nodes() == old num_nodes()-1.
  *       Else,               new num_nodes() == old num_nodes()
  *
  * Can invalidate vector i2u indexes and vector adjN for this @a n -- in other words,
  * old i2u indexes may not be equal new indexes and adjN will become empty.
  * Must not invalidate the indexes for vector Nvec.
  *
  * Complexity: No more than O(num_nodes())
  */

  node_iterator remove_node(node_iterator n_it) {
	  auto n = *(n_it);
	  remove_node(n);
	  node_iterator nextit = n_it;
	  return nextit;
	  
  }


 
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

	 struct inciedges {
		 size_type eid;
		 size_type id2;

		 inciedges(size_type eid_1, size_type id2_1) :eid(eid_1), id2(id2_1) {}
	 };

	 struct proxynodes {
		 size_type i;
		 Point P;
		 node_value_type val;
		 std::vector<inciedges> adjN;

		 proxynodes(size_type reali,Point P1, node_value_type val1,std::vector<inciedges> adjN1) 
			 :i(reali),P(P1), val(val1),adjN(adjN1) {}

	 };

	 struct proxyedges {
		 size_type ei;
		 size_type id1;
		 size_type id2;
		 edge_value_type eval;

		 proxyedges(size_type realei,size_type id_1, size_type id_2,edge_value_type eval1) 
			 :ei(realei),id1(id_1), id2(id_2), eval(eval1) {}
	 };


};

#endif // CME212_GRAPH_HPP
