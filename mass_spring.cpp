/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <thread>
#include <thrust/for_each.h>
#include <thrust/system/omp/execution_policy.h>

#include "CME212/BoundingBox.hpp"
#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "SpaceSearcher.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Spring constant
  double length;     //< Initial edge length
  EdgeData() : K(100.0), length(1.0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using NodeIterator=typename GraphType::node_iterator;
using Edge = typename GraphType::edge_type;


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
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(@a t)
 *           where @a t is the current time.
 *           @a constraint must apply the given constraint on
 *           all nodes at time @a t.
 */

/* Functor to update the position for a given single node
 * given update step size dt. 
 * This is later used in symp_euler_step function */
struct updatepos{
  updatepos(double dt1):dt(dt1){}
  
  __host__ __device__
  Point operator()(Node n){
    if(n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0))
        n.position()+=n.value().vel*dt;
    return n.position();
  }
private:
  double dt;
};

/* Functor to update the velocity for a given single node
 * given force, time, and the update step size
 * This is later used in symp_euler_step function */
template<typename F>
struct updatevel{
  updatevel(F f1,double t1, double dt1)
    : force(f1),t(t1),dt(dt1){}

  __host__ __device__
  Point operator()(Node n){
    if(n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0))
        n.value().vel+=force(n,t)*(dt/n.value().mass);
    return n.value().vel;
  }
private:
  F force;
  double t;
  double dt;
};

/* symp-euler_step function using parallelism */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  thrust::for_each(thrust::omp::par,g.node_begin(),g.node_end(),updatepos(dt));
  constraint(g,t);
  thrust::for_each(thrust::omp::par,g.node_begin(),g.node_end(),updatevel<F>(force,t,dt));
  

  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    //Set force at end-points to 0
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    {
      return Point(0, 0, 0);
    }
    //Compute gravitational force
    Point f_grav = n.value().mass * Point(0, 0, -grav);
    //Compute spring force
    Point f_spr = Point(0, 0, 0);
    for(auto incid_itr = n.edge_begin(); incid_itr != n.edge_end(); ++incid_itr)
    {
      auto edge_curr = *incid_itr;
      Point edge_vec = Point(0, 0, 0);
      if(edge_curr.node1() == n)
        edge_vec = n.position() - edge_curr.node2().position();
      else
        edge_vec = n.position() - edge_curr.node1().position();
      f_spr -= edge_curr.value().K*edge_vec*(norm(edge_vec) - edge_curr.value().length)/norm(edge_vec);
    }
    (void) t; // silence warnings
    return f_grav + f_spr;
  }
};


/** Force function object for spring force. */
struct SpringForce {
  /** Return the spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point f_spr = Point(0, 0, 0);
    for(auto incid_itr = n.edge_begin(); incid_itr != n.edge_end(); ++incid_itr)
    {
      auto edge_curr = *incid_itr;
      Point edge_vec = Point(0, 0, 0);
      if(edge_curr.node1() == n)
        edge_vec = n.position() - edge_curr.node2().position();
      else
        edge_vec = n.position() - edge_curr.node1().position();
      f_spr -= edge_curr.value().K*edge_vec*(norm(edge_vec) - edge_curr.value().length)/norm(edge_vec);
    }
    (void) t; // silence warnings
    return f_spr;
  }
};

/** Force function object for gravitational force. */
struct GravitationalForce {
  /** Return the gravitational force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Force function object for damping force. */
struct DampingForce {
  double damp_coeff;
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    return -n.value().vel*damp_coeff;
  }
};


/** Add two forces together.
 * @tparam force_type_1 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 * @tparam force_type_2 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 */
template<typename force_type_1, typename force_type_2>
struct AddForce {
  force_type_1 force_1;
  force_type_2 force_2;
  /** Return the sum of the two forces @a n at time @a t.
   * @param[in]     n      node in the graph
   * @param[in]     t      the current time (useful for time-dependent forces)
   * @return the sum of the two forces @a n at time @a t.
   *
   * @tparam NODE is a class representing the graph nodes.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    return force_1(n, t) + force_2(n, t);
  }
};

/** Combine two forces together.
 * @tparam force_type_1 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 * @tparam force_type_2 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 */
template<typename force_type_1, typename force_type_2>
 /** Return the sum of the two forces.
  * @param[in]     force_1      force of type force_type_1
  * @param[in]     force_2      force of type force_type_2
  * @return a new force containing the sum of the two forces.
  *
  */
AddForce<force_type_1, force_type_2> make_combined_force(force_type_1 force_1, force_type_2 force_2){
  return {force_1, force_2};
};

/** Combine three forces together.
 * @tparam force_type_1 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 * @tparam force_type_2 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 * @tparam force_type_3 is a function object called as @a force(n, @a t)
 *                      where n is a node of the graph and @a t is the current time.
 *                      @a force must return a Point representing the force vector on
 *                      Node n at time @a t.
 */
template<typename force_type_1, typename force_type_2, typename force_type_3>
  /** Return the sum of the two forces.
  * @param[in]     force_1      force of type force_type_1
  * @param[in]     force_2      force of type force_type_2
  * @param[in]     force_3      force of type force_type_3
  * @return a new force containing the sum of the three forces.
  *
  */
AddForce<AddForce<force_type_1, force_type_2>, force_type_3> make_combined_force(force_type_1 force_1, force_type_2 force_2, force_type_3 force_3){
  return make_combined_force(make_combined_force(force_1, force_2), force_3);
};


/** Apply two non-contradicting constraints together.
 * @tparam constr_type_1 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 * @tparam constr_type_2 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 */
template<typename constr_type_1, typename constr_type_2>
struct ApplyConstraints {
  constr_type_1 constr_1;
  constr_type_2 constr_2;
  /** Apply the two constraints at time @a t.
   * @param[in]     g   the graph on which to apply constraints
   * @param[in]     t   the current time (useful for time-dependent forces)
   */
  void operator()(GraphType& g, double t) {
    constr_1(g, t);
    constr_2(g, t);
  }
};

/** Combine two constraints together.
 * @tparam constr_type_1 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 * @tparam constr_type_2 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 */
template<typename constr_type_1, typename constr_type_2>
 /** Return the sum of the two forces.
  * @param[in]     constr_1   constraint of type constr_type_1
  * @param[in]     constr_2   constraint of type constr_type_2
  * @return a new constraint imposing both the constraints.
  *
  */
ApplyConstraints<constr_type_1, constr_type_2> make_combined_constraints(constr_type_1 constr_1, constr_type_2 constr_2){
  return {constr_1, constr_2};
};

/** Combine two constraints together.
 * @tparam constr_type_1 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 * @tparam constr_type_2 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 * @tparam constr_type_3 is a function object called as @a constraint(@a t)
 *                       where @a t is the current time.
 *                       @a constraint must apply the given constraint on
 *                       all nodes at time @a t.
 */
template<typename constr_type_1, typename constr_type_2, typename constr_type_3>
 /** Return the sum of the two forces.
  * @param[in]     constr_1   constraint of type constr_type_1
  * @param[in]     constr_2   constraint of type constr_type_2
  * @param[in]     constr_3   constraint of type constr_type_3
  * @return a new constraint imposing all three constraints.
  *
  */
ApplyConstraints<ApplyConstraints<constr_type_1, constr_type_2>, constr_type_3> make_combined_constraints(constr_type_1 constr_1, constr_type_2 constr_2, constr_type_3 constr_3){
  return make_combined_constraints(make_combined_constraints(constr_1, constr_2), constr_3);
};

/** Constraint function object for plane constraint. */
struct ConstraintPlane
{
  /** Applying plane constraint to @a g at time @a t. */
  void operator()(GraphType& g, double t){
    (void) t;
    for(auto node_itr = g.node_begin(); node_itr != g.node_end(); ++node_itr)
    {
      auto node_curr = *node_itr;
      if(node_curr.position().z < -0.75)
      {
        node_curr.value().vel.z = 0;
        node_curr.position().z = -0.75;
      }
    }
  }
};

/** Constraint function object for sphere constraint. */
struct ConstraintSphere
{
  /** Applying sphere constraint to @a g at time @a t. */
  void operator()(GraphType& g, double t){
    (void) t;
    Point sph_center = Point(0.5, 0.5, -0.5);
    double sph_rad = 0.15;
    for(auto node_itr = g.node_begin(); node_itr != g.node_end(); ++node_itr)
    {
      auto node_curr = *node_itr;
      if(norm(node_curr.position() - sph_center) != 0 && norm(node_curr.position() - sph_center) < sph_rad)
      {
        Point Ri = (node_curr.position() - sph_center)/norm(node_curr.position() - sph_center);
        node_curr.value().vel -= dot(node_curr.value().vel, Ri)*Ri;
        node_curr.position() = sph_center + (node_curr.position() - sph_center)*sph_rad/norm(node_curr.position() - sph_center);
      }
      else if(norm(node_curr.position() - sph_center) == 0)
      {
        node_curr.value().vel.x = 0;
        node_curr.position() = sph_center + Point(0.15, 0, 0);
      }
    }
  }
};

/** Constraint function object for sphere constraint with tearing. */
struct ConstraintSphereTear
{
  /** Applying sphere constraint to @a g at time @a t. */
  void operator()(GraphType& g, double t){
    (void) t;
    Point sph_center = Point(0.5, 0.5, -0.5);
    double sph_rad = 0.15;
    for(auto node_itr = g.node_begin(); node_itr != g.node_end(); ++node_itr)
    {
      auto node_curr = *node_itr;
      if(norm(node_curr.position() - sph_center) < sph_rad)
        g.remove_node(node_curr);
    }
  }
};

/* Functor to remove a certain component of the velocity 
*  when certain constraints are satisfied*/
struct rmvel {
	//Constructor
	rmvel(Point c, Node n2_1, double r2) :center(c), n2(n2_1), radius2(r2) {}

	Point operator()(Node n) {
		Point r = center - n2.position();
		double l2 = normSq(r);
		if (n != n2 && l2 < radius2) {
			n.value().vel -= (dot(r, n.value().vel) / 12)*r;
		}
		return n.value().vel;
	}
private:
	Point center;
	Node n2;
	double radius2;
};

/* Functor to find the radius2, given a previous found radius2
*  and the next edge. This is later passed in std::accumulate function 
*  to find the radius2 after loop through all incident edges*/
struct get_r2 {
	get_r2(Point c):center(c) {}


	double operator()(double previous, Edge e2) {
		double radius2 = std::min(previous, normSq(e2.node2().position()-center));
		return radius2;
	}
private:
	Point center;
};

/* Functor to map a node to its position. This will be useful when 
*  define SpaceSearcher and constructing a point iterator by NodeIterator.*/
struct get_position {
	Point operator()(Node n) {
		return n.position();
	}
};

/* Function to return a graph_spacesearcher by constructing a large enough boundingbox
*  and set t2p to be get_position() functor */
SpaceSearcher<Node> graph_spacesearcher(const GraphType& g){
  thrust::transform_iterator<get_position,NodeIterator,Point>iter_begin
        (g.node_begin(),get_position());
      thrust::transform_iterator<get_position,NodeIterator,Point>iter_end
        (g.node_end(),get_position());
  Box3D largebox(iter_begin,iter_end);

  //Relax it by adding two large points to this box
  largebox |=Point(100,100,100);
  largebox |=Point(-100,-100,-100);

  SpaceSearcher<Node> ss (largebox,g.node_begin(), g.node_end(), get_position());
 
  return ss;
}

/* Functor to check selfcollision constraints and remove some components of velocities
*  given a single node*/
struct do_selfcollision {
  //Constructor
    do_selfcollision(GraphType* g, SpaceSearcher<Node> ss1)
     :graph_(const_cast<GraphType*>(g)),ss(ss1){
      
    }

  __host__ __device__
	void operator()(Node n) {
		const Point& center = n.position();
		double r2 = std::numeric_limits<double>::max();
		double rad2=std::accumulate(n.edge_begin(), n.edge_end(), r2, get_r2(center));
		rad2 *= 0.9;
		double r = sqrt(rad2);

    /* Since the boundingbox is not a sphere, we want to leave more spaces for relax 
      and check the velocity removal condition in rmvel*/
		Point lower = center - Point(r, r, r);
		Point upper = center + Point(r, r, r);
		Box3D bb(lower, upper);

		thrust::for_each(ss.begin(bb), ss.end(bb), rmvel(center,n,rad2));
	}
private:
	GraphType* graph_;
  SpaceSearcher<Node> ss;
};

/* Apply do_selfcollision to all nodes in the graph*/
struct SelfCollisionConstraint {
	void operator()(GraphType& g, double)const {
    
		thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), do_selfcollision(&g,graph_spacesearcher(g)));

	}
};


int main(int argc, char** argv)
{
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

  //Set initial values for nodes and edges
  for(auto node_itr = graph.node_begin(); node_itr != graph.node_end(); ++node_itr)
  {
    auto node = *node_itr;
    node.value().mass = (1.0/graph.num_nodes())/graph.num_nodes();
    node.value().vel = Point(0, 0, 0);
  }
  for(auto edge_itr = graph.edge_begin(); edge_itr != graph.edge_end(); ++edge_itr)
  {
    auto edge = *edge_itr;
    edge.value().K = 100.0/graph.num_nodes();
    edge.value().length = edge.length();
  }

  // Add forces
  //double damp_coeff = 1.0/graph.num_nodes();
  //DampingForce f_damp{damp_coeff};
  SpringForce f_spr;
  GravitationalForce f_grav;
  auto f_added = make_combined_force(f_spr, f_grav);

  // Combine constraints
  ConstraintPlane plane_constr;
  //ConstraintSphere sphere_constr;
  SelfCollisionConstraint sc_constr; 
  //auto combined_constraints = make_combined_constraints(plane_constr, sphere_constr);
  auto combined_constraints = make_combined_constraints(plane_constr, sc_constr);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // We want viewer interaction and the simulation at the same time
  // Viewer is thread-safe, so launch the simulation in a child thread
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 1.0/graph.num_nodes();
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, f_added, combined_constraints);

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}