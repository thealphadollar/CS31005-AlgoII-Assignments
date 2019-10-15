#include<iostream>
#include<memory>
#ifndef VoronoiDiagram_hpp
#define VoronoiDiagram_hpp


#include <iostream>
#include <vector>
#include <limits>
#include <cmath>

#define POINT_EPSILON 1.0e-6

class Point2D {

	struct Point2D_XY_Compare {
		bool operator()(const Point2D &p1, const Point2D &p2) {
			return (p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y));
		}
	};

public:
    double x, y;
    int id;
    const static double Inf;
	static Point2D_XY_Compare xy_compare;
    
    Point2D(double x = 0.0, double y = 0.0);
//    Point2D(double x = 0.0, double y = 0.0, int id=0);
    Point2D(const Point2D &point);
    
    friend double dotProduct(const Point2D &p1, const Point2D &p2);
    friend double crossProduct(const Point2D &p1, const Point2D &p2);
    
    friend Point2D operator+(const Point2D &p1, const Point2D &p2);
    friend Point2D operator-(const Point2D &p1, const Point2D &p2);
    friend Point2D operator/(const Point2D &p1, const Point2D &p2);
    friend Point2D operator*(const Point2D &p, double value);
    friend Point2D operator*(double value, const Point2D &p);
    friend Point2D operator/(const Point2D &p, double value);
	friend Point2D operator-(const Point2D &p);

    friend std::ostream &operator<<(std::ostream &stream, const Point2D &p);
    friend std::vector<Point2D> &operator<<(std::vector<Point2D> &v, const Point2D &p);
    
    Point2D &operator-=(const Point2D &p);
    Point2D &operator+=(const Point2D &p);
    Point2D &operator*=(double value);
    Point2D &operator/=(double value);
    
    Point2D normalized();
    void normalize();
    double norm();
    double norm2();
    
    Point2D getRotated90CW();
    Point2D getRotated90CCW();

	static bool isLeftTurn(const Point2D &p1, const Point2D &p2, const Point2D &p3);
	static bool isRightTurn(const Point2D &p1, const Point2D &p2, const Point2D &p3);
    
    double operator[](int i);
    
    void setX(double x);
    void setY(double y);
    void setId(int id);
    
    bool isVertical();
    bool isHorizontal();
    bool isValid();
    
};

double dotProduct(const Point2D &p1, const Point2D &p2);
double crossProduct(const Point2D &p1, const Point2D &p2);

bool equal(const Point2D &p1, const Point2D &p2, double EPSILON = POINT_EPSILON);
bool equal(double v1, double v2, double EPSILON = POINT_EPSILON);

// PARABOLA


/**
 
 Calculate number of intersection points between two parabolas with foci `f1` and `f2` and with given `directrix`
 
 */
int intersectionPointsNum(const Point2D &f1, const Point2D &f2, double directrix);


/**

 Find intersection points of two parabolas with foci `f1` and `f2` and with given `directrix`
 
 */
std::vector<Point2D> findIntersectionPoints(const Point2D &f1, const Point2D &f2, double directrix);







#include <time.h>
#include <iostream>
#include<memory>
#include <limits>
#include <iomanip>
#include <vector>
#include <random>
#include <cassert>




// DCEL

namespace DCEL {

    class Vertex;
    class HalfEdge;
    
    typedef std::shared_ptr<HalfEdge> HalfEdgePtr;
    typedef std::shared_ptr<Vertex> VertexPtr;


    class Vertex {
    public:
        
        Point2D point;
        HalfEdgePtr edge; // The edge points towards this vertex [-->o]
        
        Vertex(const Point2D &pos, HalfEdgePtr incident_edge = nullptr);
        
        inline double x() { return point.x; }
        inline double y() { return point.y; }

    };

    
    class HalfEdge {
    public:
        
        int l_index, r_index;
        
        VertexPtr vertex;
        HalfEdgePtr twin;
        HalfEdgePtr next;
        HalfEdgePtr prev;
        
        HalfEdge(int _l_index, int _r_index, VertexPtr _vertex = nullptr);
        
        inline VertexPtr vertex0() { return vertex; }
        inline VertexPtr vertex1() { return twin->vertex; }
        inline bool is_finite() {
            return vertex != nullptr && twin->vertex != nullptr;
        }
        
        // Iterators around vertex
        HalfEdgePtr vertexNextCCW();
        HalfEdgePtr vertexNextCW();
    };

    
    std::pair<HalfEdgePtr, HalfEdgePtr> make_twins(int left_index, int right_index);
    
    
    std::pair<HalfEdgePtr, HalfEdgePtr> make_twins(const std::pair<int,int> &indices);
    
    
    void connect_halfedges(HalfEdgePtr p1, HalfEdgePtr p2);
    
}




class Event;


namespace beachline {
    
    using namespace DCEL;
    
    class BLNode;
    typedef std::shared_ptr<BLNode> BLNodePtr;

    class BLNode {
    public:
        
        // Height of the tree
        int height;
        
        // Pointer to a position of a sweepline
        double *sweepline;
        
        // Pointer to a vector of input points
        const std::vector<Point2D> *points;
        
        // Indices of the points
        std::pair<int, int> indices;
        
        // Pointers to left, right children and parent node
        BLNodePtr left, right, parent;
        
        // Pointer to a circle event for a leaf node or halfedge for an internal node
        std::shared_ptr<Event> circle_event;
        std::shared_ptr<HalfEdge> edge;
        
        // Constructor
        BLNode(const std::pair<int,int>& _indices,
               double* _sweepline = nullptr,
               const std::vector<Point2D>* _points = nullptr,
               BLNodePtr _left = nullptr,
               BLNodePtr _right = nullptr,
               BLNodePtr _parent = nullptr,
               int _height = 1);
        
        // Pointers to a next and previous arc-nodes
        BLNodePtr next, prev;
        
        // Leaf is defined as <p_i,p_i>
        inline bool is_leaf() {
            return indices.first == indices.second;
        }
        
        // Returns id of the node (only for leafs)
        inline int get_id() {
            return indices.first;
        }
        
        // Check if indices are equal
        inline bool has_indices(int a, int b) {
            return indices.first == a && indices.second == b;
        }
        
        // Check if indices are equal
        inline bool has_indices(const std::pair<int,int> &p) {
            return indices.first == p.first && indices.second == p.second;
        }
        
        // Return x-coordinate of:
        //  - in case of leaf node - corresponding focus of parabola;
        //  - in case of internal node - breakpoint;
        double value();
        
    };
    
    
    /**
     Connect as a list
     */
    void connect(BLNodePtr prev, BLNodePtr next);


    /**
     Check if the node is a root node
     */
    bool is_root(BLNodePtr node);


    /**
     Get height of the node
     */
    int get_height(BLNodePtr node);


    /**
     Update height of the node
     */
    void update_height(BLNodePtr node);


    /**
     Get balance of the node (difference between the height of left and right subtrees)
     */
    int get_balance(BLNodePtr node);


    /**
     Performs rotation of a tree around `node` such that it goes to the left subtree
     */
    BLNodePtr rotate_left(BLNodePtr node);


    /**
     Performs rotation of a tree around `node` such that it goes to the right subtree
     */
    BLNodePtr rotate_right(BLNodePtr node);


    /**
     Find a leaf in a tree such that x is under the parabolic arc,
     which corresponds to this leaf.
     */
    BLNodePtr find(BLNodePtr root, double x);
    

    /**
     Replace a leaf `node` with a new subtree, which has root `new_node`.
     The function rebalances the tree and returns the pointer to a new root node.
     */
    BLNodePtr replace(BLNodePtr node, BLNodePtr new_node);
    
    
    /**
     Remove a disappearing arc related to a circle event.
     The function rebalances the tree and returns the pointer to a new root node.
     */
    BLNodePtr remove(BLNodePtr leaf);;
    
    
    /**
     Returns breakpoints for a given arc
     */
    std::pair<BLNodePtr, BLNodePtr> breakpoints(BLNodePtr leaf);
    
    
    BLNodePtr make_subtree(int index, int index_behind, double *sweepline,
                           const std::vector<Point2D> *points,
                           std::vector<HalfEdgePtr> &edges);
    
    
    BLNodePtr make_simple_subtree(int index, int index_behind, double *sweepline,
                                  const std::vector<Point2D> *points,
                                  std::vector<HalfEdgePtr> &edges);
    
    
    bool _validate(BLNodePtr node);
    
    
    bool _check_balance(BLNodePtr node);


    /**
     Print tree
     */
    void print_tree(BLNodePtr root, int width = 7);
    
}


namespace bl = beachline;


void build_voronoi(const std::vector<Point2D> &points,
                   std::vector<bl::HalfEdgePtr> &halfedges,
                   std::vector<bl::VertexPtr> &vertices,
                   std::vector<bl::HalfEdgePtr> &faces);

//std::vector<bl::HalfEdgePtr> init
//

#endif /* VoronoiDiagram_hpp */
