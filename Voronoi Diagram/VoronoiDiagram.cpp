
#include "VoronoiDiagram.hpp"
#include <cmath>
#include <queue>

#define BREAKPOINTS_EPSILON 1.0e-5
#define _DEBUG_
#define CIRCLE_CENTER_EPSILON 1.0e-4


#if defined(_WIN64) || defined(_WIN32)      
    #define isnan(x) _isnan(x)
#endif

using namespace std;

/* =================================

Function definitions for 2D Points

================================= */


const double Point2D::Inf = std::numeric_limits<double>::infinity();

Point2D::Point2D_XY_Compare Point2D::xy_compare = Point2D::Point2D_XY_Compare();

Point2D::Point2D(double _x, double _y) : x(_x), y(_y) {
}
//Point2D::Point2D(double _x,double _y,int _id) :  x(_x), y(_y),id(_id) {  }
Point2D::Point2D(const Point2D &point) : x(point.x), y(point.y) ,id(point.id) {
}

double dotProduct(const Point2D &p1, const Point2D &p2) {
    return p1.x * p2.x + p1.y * p2.y;
}

double crossProduct(const Point2D &p1, const Point2D &p2) {
    return p1.x * p2.y - p1.y * p2.x;
}

Point2D operator+(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x + p2.x, p1.y + p2.y);
}

Point2D operator-(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x - p2.x, p1.y - p2.y);
}

Point2D operator/(const Point2D &p1, const Point2D &p2) {
    return Point2D(p1.x / p2.x, p1.y / p2.y);
}

Point2D operator*(const Point2D &p, double value) {
    return Point2D(p.x * value, p.y * value);
}

Point2D operator*(double value, const Point2D &p) {
    return Point2D(p.x * value, p.y * value);
}

Point2D operator/(const Point2D &p, double value) {
    return Point2D(p.x / value, p.y / value);
}

Point2D operator-(const Point2D &p) {
	return Point2D(-p.x, -p.y);
}

std::ostream &operator<<(std::ostream &stream, const Point2D &p) {
    stream << "(" << p.x << "," << p.y << ")";
    return stream;
}

std::vector<Point2D> &operator<<(std::vector<Point2D> &v, const Point2D &p) {
    v.push_back(p);
    return v;
}

Point2D &Point2D::operator-=(const Point2D &p) {
    x -= p.x;
    y -= p.y;
    return *this;
}

Point2D &Point2D::operator+=(const Point2D &p) {
    x += p.x;
    y += p.y;
    return *this;
}

Point2D &Point2D::operator*=(double value) {
    x *= value;
    y *= value;
    return *this;
}

Point2D &Point2D::operator/=(double value) {
    x /= value;
    y /= value;
    return *this;
}

double Point2D::operator[](int i) {
    if (i==0) return x;
    else return y;
}

void Point2D::setX(double x) {
    this->x = x;
}

void Point2D::setY(double y) {
    this->y = y;
}

void Point2D::setId(int id){
    this->id = id;

}    

bool Point2D::isVertical() {
    return (y == Inf && !isnan(x) && x != Inf);
}

bool Point2D::isHorizontal() {
    return (x == Inf && !isnan(y) && y != Inf);
}

bool Point2D::isValid() {
    if (x == Inf && y == Inf)
        return false;
    return (!isnan(x) && !isnan(y));
}

Point2D Point2D::normalized() {
    return (*this) / this->norm();
}

void Point2D::normalize() {
    double n = norm();
    x /= n;
    y /= n;
}

double Point2D::norm() {
    return sqrt(x * x + y * y);
}

double Point2D::norm2() {
    return x *x + y * y;
}

Point2D Point2D::getRotated90CW() {
    return Point2D(y, -x);
}

Point2D Point2D::getRotated90CCW() {
    return Point2D(-y, x);
}

bool Point2D::isLeftTurn(const Point2D &p1, const Point2D &p2, 
						 const Point2D &p3) {
	return (crossProduct(p2 - p1, p3 - p2) > 0.0);
}

bool Point2D::isRightTurn(const Point2D &p1, const Point2D &p2, 
						  const Point2D &p3) {
	return (crossProduct(p2 - p1, p3 - p2) < 0.0);
}

bool equal(const Point2D &p1, const Point2D &p2, double EPSILON) {
    return (fabs(p1.x - p2.x) < EPSILON && fabs(p1.y - p2.y) < EPSILON);
}

bool equal(double v1, double v2, double EPSILON) {
    return fabs(v1 - v2) < EPSILON;
}


/* =================================

Function definitions for Circle

================================= */


bool findCircleCenter(const Point2D &p1, const Point2D &p2, const Point2D &p3, Point2D &center) {
    
    Point2D u1 = (p1 - p2).normalized(), u2 = (p3 - p2).normalized();
    
    double cross = crossProduct(u1, u2);
    
    // check if vectors are collinear
    if (fabs(cross) < CIRCLE_CENTER_EPSILON) {
        return false;
    }
    
    // get cental points
    Point2D pc1 = 0.5 * (p1 + p2), pc2 = 0.5 * (p2 + p3);
    
    // get free components
    double b1 = dotProduct(u1, pc1), b2 = dotProduct(u2, pc2);
    
    // calculate the center of a circle
    center.x = (b1 * u2.y - b2 * u1.y) / cross;
    center.y = (u1.x * b2 - u2.x * b1) / cross;
    
    return true;
}


/* =================================

Function definitions for DCEL

================================= */

namespace DCEL {
    
   
    Vertex::Vertex(const Point2D &pos, HalfEdgePtr incident_edge) : point(pos), edge(incident_edge)  { }
    HalfEdge::HalfEdge(int _l_index, int _r_index, VertexPtr _vertex) :
        l_index(_l_index), r_index(_r_index), vertex(_vertex) {}
    HalfEdgePtr HalfEdge::vertexNextCCW() {
        return twin->prev;
    }
    HalfEdgePtr HalfEdge::vertexNextCW() {
        return next->twin;
    }
    std::pair<HalfEdgePtr, HalfEdgePtr> make_twins(int left_index, int right_index) {
        
        HalfEdgePtr h = std::make_shared<HalfEdge>(left_index, right_index);
        HalfEdgePtr h_twin = std::make_shared<HalfEdge>(right_index, left_index);
        
        h->twin = h_twin;
        h_twin->twin = h;
        
        return std::make_pair(h, h_twin);
    }
    std::pair<HalfEdgePtr, HalfEdgePtr> make_twins(const std::pair<int,int> &indices) {
        
        return make_twins(indices.first, indices.second);
    }
    void connect_halfedges(HalfEdgePtr p1, HalfEdgePtr p2) {
        p1->next = p2;
        p2->prev = p1;
    }
}



/* =================================

Function definitions for Parabola

================================= */

int intersectionPointsNum(const Point2D &f1, const Point2D &f2, double directrix) {
    if (fabs(f1.x - f2.x) < POINT_EPSILON && fabs(f1.y - f2.y) < POINT_EPSILON) {
        return -1;
    }
    if (fabs(f1.y - f2.y) < POINT_EPSILON)
        return 1;
    return 2;
}

std::vector<Point2D> findIntersectionPoints(const Point2D &f1, const Point2D &f2, double d) {
    std::vector<Point2D> result;
    if (fabs(f1.x - f2.x) < POINT_EPSILON) {
        double y = 0.5 * (f1.y + f2.y), D = sqrt(d * d - d * (f1.y + f2.y) + f1.y * f2.y);
        result.push_back(Point2D(f1.x - D, y));
        result.push_back(Point2D(f1.x + D, y));
    } else if (fabs(f1.y - f2.y) < POINT_EPSILON) { 
        double x = 0.5 * (f1.x + f2.x);
        result.push_back(Point2D(x, 0.5 * ((x - f1.x) * (x - f1.x) + f1.y * f1.y  - d * d) / (f1.y - d)));
    } else {
        
        double D = 2. * sqrt(pow(f1.x - f2.x, 2) * (d - f1.y) * (d - f2.y) * (pow(f1.x - f2.x, 2) + pow(f1.y - f2.y, 2)));
        double T = -2. * d * pow(f1.x - f2.x, 2) + (f1.y + f2.y) * (pow(f2.x - f1.x, 2) + pow(f2.y - f1.y, 2));
        double Q = 2. * pow(f1.y - f2.y, 2);
        
        double y1 = (T - D) / Q, y2 = (T + D) / Q;
        double x1 = 0.5 * (f1.x * f1.x - f2.x * f2.x + (2 * y1 - f2.y - f1.y) * (f2.y - f1.y)) / (f1.x - f2.x);
        double x2 = 0.5 * (f1.x * f1.x - f2.x * f2.x + (2 * y2 - f2.y - f1.y) * (f2.y - f1.y)) / (f1.x - f2.x);
        
        if (x1 > x2) {
            std::swap(x1, x2);
            std::swap(y1, y2);
        }
        result.push_back(Point2D(x1, y1));
        result.push_back(Point2D(x2, y2));
    }
    return result;
}


/* =================================

Function definitions for Beachline

================================= */

namespace bl = beachline;
namespace beachline {

    BLNode::BLNode(const std::pair<int,int>& _indices,
                   double* _sweepline ,
                   const std::vector<Point2D>* _points,
                   BLNodePtr _left,
                   BLNodePtr _right,
                   BLNodePtr _parent,
                   int _height) : indices(_indices), left(_left), right(_right),
                                  parent(_parent), height(_height),
                                  sweepline(_sweepline), points(_points),
                                  next(nullptr), prev(nullptr) {}

    double BLNode::value() {
        if (points == nullptr)
            return std::numeric_limits<double>::infinity();
        if (is_leaf()) {
            return (*points)[indices.first].x;
        } else {
            Point2D p1 = (*points)[indices.first], p2 = (*points)[indices.second];
            
            std::vector<Point2D> ips = findIntersectionPoints(p1, p2, *sweepline);
            if (ips.size() == 2) {
                if (p1.y < p2.y) {
                    return ips[0].x;
                } else {
                    return ips[1].x;
                }
            } else {
                return ips[0].x;
            }
        }
    }

    /**
     Connect as a list
     */
    void connect(BLNodePtr prev, BLNodePtr next) {
        prev->next = next;
        next->prev = prev;
    }
    /**
     Check if the node is a root node
     */
    bool is_root(BLNodePtr node) {
        return node->parent == nullptr;
    }
    /**
     Get height of the node
     */
    int get_height(BLNodePtr node) {
        if (node == nullptr) return 0;
        return node->height;
    }
    /**
     Update height of the node
     */
    void update_height(BLNodePtr node) {
        if (node == nullptr)
            return;
        node->height = std::max(get_height(node->left), get_height(node->right)) + 1;
    }
    /**
     Get balance of the node (difference between the height of left and right subtrees)
     */
    int get_balance(BLNodePtr node) {
        return get_height(node->left) - get_height(node->right);
    }

    BLNodePtr rotate_left(BLNodePtr node) {
        
        if (node == nullptr)
            return nullptr;
        
        if (node->right == nullptr)
            return node;
        
        // get right node, which becomes a new root node
        BLNodePtr rnode = node->right;
        
        // establish connections with a root node if threre is one
        if (!is_root(node)) {
            if (node->parent->left == node) {
                node->parent->left = rnode;
            } else {
                node->parent->right = rnode;
            }
        }
        rnode->parent = node->parent;
        
        // connect right subtree of the left child as a left subtree of `node`
        node->right = rnode->left;
        if (rnode->left != nullptr) {
            rnode->left->parent = node;
        }
        
        // connect `node` as a right child of it's child
        rnode->left = node;
        node->parent = rnode;
        
        // update height attribute
        update_height(node);
        update_height(rnode);
        update_height(rnode->parent);
        
        return rnode;
    }
    /**
     Performs rotation of a tree around `node` such that it goes to the right subtree
     */
    BLNodePtr rotate_right(BLNodePtr node) {
        
        if (node == nullptr)
            return nullptr;
        
        if (node->left == nullptr)
            return node;
        
        // left node becomes root node of subtree
        BLNodePtr lnode = node->left;
        
        // establish connections with a root node if threre is one
        if (!is_root(node)) {
            if (node->parent->left == node) {
                node->parent->left = lnode;
            } else {
                node->parent->right = lnode;
            }
        }
        lnode->parent = node->parent;
        
        // connect right subtree of the left child as a left subtree of `node`
        node->left = lnode->right;
        if (lnode->right != nullptr) {
            lnode->right->parent = node;
        }
        
        // connect `node` as a right child of it's child
        lnode->right = node;
        node->parent = lnode;
        
        // update height attribute
        update_height(node);
        update_height(lnode);
        update_height(lnode->parent);
        
        return lnode;
    }

    BLNodePtr find(BLNodePtr root, double x) {
        if (root == nullptr) {
            return nullptr;
        }
        BLNodePtr node = root;
        while (!node->is_leaf()) {
            if (node->value() < x) {
                node = node->right;
            } else {
                node = node->left;
            }
        }
        return node;
    }

    BLNodePtr replace(BLNodePtr node, BLNodePtr new_node) {
        
        if (node == nullptr) {
            return new_node;
        }
        
        // Find x-coordinate
        double x = new_node->value();
        
        // Get a parent node
        BLNodePtr parent_node = node->parent;
        

        new_node->parent = parent_node;
        if (parent_node != nullptr) {
            if (parent_node->value() < x) {
                parent_node->right = new_node;
            } else {
                parent_node->left = new_node;
            }
        }
        
        // Rebalance the tree
        node = new_node;
        while (parent_node != nullptr) {
            update_height(parent_node);
            int balance = get_balance(parent_node);
            if (balance > 1) { // left subtree is higher than right subtree by more than 1
                if (parent_node->left != nullptr && !parent_node->left->is_leaf() && get_balance(parent_node->left) < 0) { // @TODO ensure that
                    parent_node->left = rotate_left(parent_node->left);
                }
                parent_node = rotate_right(parent_node);
            } else if (balance < -1) { // right subtree is lower than left subtree by more than 1
                if (parent_node->right != nullptr && !parent_node->right->is_leaf() && get_balance(parent_node->right) > 0) {
                    parent_node->right = rotate_right(parent_node->right);
                }
                parent_node = rotate_left(parent_node);
            }
            
            //_validate(parent_node);
            
            node = parent_node;
            parent_node = parent_node->parent;
        }
        
        //_check_balance(node);
        
        return node;
    }


    BLNodePtr remove(BLNodePtr leaf) {
        
        if (leaf == nullptr)
            return nullptr;
        
        BLNodePtr parent = leaf->parent, grandparent = parent->parent;
        std::pair<int,int> bp1(leaf->prev->get_id(), leaf->get_id());
        std::pair<int,int> bp2(leaf->get_id(), leaf->next->get_id());
        std::pair<int,int> other_bp;
        
        assert(leaf->next != nullptr);
        assert(leaf->prev != nullptr);
        assert(parent != nullptr);
        assert(grandparent != nullptr);
        
        assert(parent->has_indices(bp1) || parent->has_indices(bp2));
        
        if (parent->has_indices(bp1)) {
            other_bp = bp2;
        } else if (parent->has_indices(bp2)) {
            other_bp = bp1;
        }
        
        BLNodePtr other_subtree;
        if (parent->left == leaf)
            other_subtree = parent->right;
        else
            other_subtree = parent->left;
        
        other_subtree->parent = grandparent;
        if (grandparent->left == parent) {
            grandparent->left = other_subtree;
        } else {
            grandparent->right = other_subtree;
        }
        
        BLNodePtr new_root = grandparent;
        // Go up and rebalance the whole tree
        while (grandparent != nullptr) {
            if (grandparent->has_indices(other_bp))
                grandparent->indices = std::make_pair(leaf->prev->get_id(), leaf->next->get_id());
            // update height of a node
            update_height(grandparent);
            // calculate balance factor of a node
            int balance = get_balance(grandparent);
            if (balance > 1) { // left subtree is higher than right subtree by more than 1
                if (grandparent->left != nullptr && !grandparent->left->is_leaf() && get_balance(grandparent->left) < 0) {
                    grandparent->left = rotate_left(grandparent->left);
                }
                grandparent = rotate_right(grandparent);
            } else if (balance < -1) { // right subtree is lower than left subtree by more than 1
                if (grandparent->right != nullptr && !grandparent->right->is_leaf() && get_balance(grandparent->right) > 0) {
                    grandparent->right = rotate_right(grandparent->right);
                }
                grandparent = rotate_left(grandparent);
            }
            
            //_validate(grandparent);
            
            new_root = grandparent;
            grandparent = grandparent->parent;
        }
        
        // Connect previous with next leaf
        connect(leaf->prev, leaf->next);
        
        //_check_balance(new_root);
        
        return new_root;
    }

    std::pair<BLNodePtr, BLNodePtr> breakpoints(BLNodePtr leaf) {
        
        if (leaf == nullptr || leaf->next == nullptr || leaf->prev == nullptr)
            return std::make_pair<BLNodePtr>(nullptr, nullptr);
        
        BLNodePtr parent = leaf->parent, gparent = leaf->parent;
        std::pair<int,int> bp1(leaf->prev->get_id(), leaf->get_id()); // left breakpoint
        std::pair<int,int> bp2(leaf->get_id(), leaf->next->get_id()); // right breakpoint
        std::pair<int,int> other_bp;
        
        bool left_is_missing = true;
        
        if (parent->has_indices(bp1)) {
            other_bp = bp2;
            left_is_missing = false;
        } else if (parent->has_indices(bp2)) {
            other_bp = bp1;
            left_is_missing = true;
        }
        
        // Go up and rebalance the whole tree
        while (gparent != nullptr) {
            if (gparent->has_indices(other_bp)) {
                break;
            }
            gparent = gparent->parent;
        }
        
        if (left_is_missing) {
            return std::make_pair(gparent, parent);
        } else {
            return std::make_pair(parent, gparent);
        }

    }

    BLNodePtr make_subtree(int index, int index_behind, double *sweepline,
                           const std::vector<Point2D> *points,
                           std::vector<HalfEdgePtr> &edges) {
        
        // create nodes corresponding to branching points
        BLNodePtr node1 = std::make_shared<BLNode>(std::make_pair(index_behind, index), sweepline, points);
        BLNodePtr node2 = std::make_shared<BLNode>(std::make_pair(index, index_behind), sweepline, points);
        
        // create leaf nodes
        BLNodePtr leaf1 = std::make_shared<BLNode>(std::make_pair(index_behind, index_behind), sweepline, points);
        BLNodePtr leaf2 = std::make_shared<BLNode>(std::make_pair(index, index), sweepline, points);
        BLNodePtr leaf3 = std::make_shared<BLNode>(std::make_pair(index_behind, index_behind), sweepline, points);
        
        // adjust tree connections
        node1->right = node2;
        node2->parent = node1;
        
        node1->left = leaf1;
        leaf1->parent = node1;
        
        node2->left = leaf2;
        leaf2->parent = node2;
        
        node2->right = leaf3;
        leaf3->parent = node2;
        
        // add halfedges
        std::pair<HalfEdgePtr, HalfEdgePtr> twin_edges = make_twins(index_behind, index);
        node1->edge = twin_edges.first;//second;//first;
        node2->edge = twin_edges.second;//first;//second;
        
        edges.push_back(twin_edges.first);
        edges.push_back(twin_edges.second);
        
        // connect leaf nodes
        connect(leaf1, leaf2);
        connect(leaf2, leaf3);
        
        // reset height of a node
        update_height(node2);
        update_height(node1);
        
        // return the result
        return node1;
    }

    BLNodePtr make_simple_subtree(int index, int index_behind, double *sweepline,
                                  const std::vector<Point2D> *points,
                                  std::vector<HalfEdgePtr> &edges) {
        
        BLNodePtr node, leaf_l, leaf_r;
        
        std::pair<HalfEdgePtr, HalfEdgePtr> twin_edges = make_twins(index_behind, index);
        
        edges.push_back(twin_edges.first);
        edges.push_back(twin_edges.second);
        
        if ((*points)[index].x < (*points)[index_behind].x) {
            // Depends on the point order
            node = std::make_shared<BLNode>(std::make_pair(index, index_behind), sweepline, points);
            leaf_l = std::make_shared<BLNode>(std::make_pair(index, index), sweepline, points);
            leaf_r = std::make_shared<BLNode>(std::make_pair(index_behind, index_behind), sweepline, points);
            node->edge = twin_edges.second;//twin_edges.first;
        } else {
            node = std::make_shared<BLNode>(std::make_pair(index_behind, index), sweepline, points);
            leaf_l = std::make_shared<BLNode>(std::make_pair(index_behind, index_behind), sweepline, points);
            leaf_r = std::make_shared<BLNode>(std::make_pair(index, index), sweepline, points);
            node->edge = twin_edges.first;//twin_edges.second;
        }
        
        node->left = leaf_l;
        node->right = leaf_r;
        
        leaf_l->parent = node;
        leaf_r->parent = node;
        
        connect(leaf_l, leaf_r);
        update_height(node);
        
        return node;
    }

    bool _validate(BLNodePtr node) {
        
        if (node == nullptr)
            return true;
        
        if (node->is_leaf()) {
            if (node->left != nullptr || node->right != nullptr) {
                std::cout << "LEAF NOT A LEAF: " << node->indices.first << ", " << node->indices.second << std::endl;
                return false;
            }
        } else {
            if (node->left == nullptr || node->right == nullptr) {
                std::cout << " BP WITHOUT LEAF: " << node->indices.first << ", " << node->indices.second << std::endl;
                return false;
            }
        }
        return true;
    }

    bool _check_balance(BLNodePtr node) {
        if (node == nullptr) return true;
        if (_check_balance(node->left) && _check_balance(node->right)) {
            if (fabs(get_balance(node)) > 1) {
                
                std::cout << "+unbalanced (" << node->indices.first << ", " << node->indices.second << ")" << std::endl;
                
                return false;
            }
        }
        return true;
    }


    void print_tree(BLNodePtr root, int width) {
        
        if (root == nullptr)
            return;
        
        std::vector<std::vector<BLNodePtr>> layers(root->height);
        
        layers[0].push_back(root);
        int size = 2;
        for (int i = 1; i < root->height; ++i) {
            layers[i].resize(size);
            for (int j = 0; j < layers[i-1].size(); ++j) {
                if (layers[i-1][j] != nullptr) {
                    layers[i][2*j] = layers[i-1][j]->left;
                    layers[i][2*j+1] = layers[i-1][j]->right;
                }
            }
            size *= 2;
        }
        
        size /= 2;
        for (int i = 0; i < root->height; ++i) {
            for (int j = 0; j < layers[i].size(); ++j) {
                if (layers[i][j] != nullptr)
                    std::cout << std::setw(width * size) << "<" << layers[i][j]->indices.first << ", " << layers[i][j]->indices.second << ">";
                else
                    std::cout << std::setw(width * size) << "      ";
            }
            std::cout << std::endl;
            size /= 2;
        }
    }

}

struct Event {
    enum { SITE = 0, CIRCLE = 1, SKIP = 2, };
    int type;
    Point2D point;
    /*
     Site event attributes:
     */
    int index;
    /*
     Circle event attributes:
     */
    Point2D center;
    bl::BLNodePtr arc;
    Event(int _index = -1, int _type = Event::SKIP, const Point2D &_point = Point2D(0.0, 0.0)) :
    index(_index), type(_type), point(_point), arc(nullptr) {}
    
};


/* =================================

Function definitions for Voronoi Combined with above definitions

================================= */


typedef std::shared_ptr<Event> EventPtr;

struct Point2DComparator {
    bool operator()(const Point2D &p1, const Point2D &p2) {
        return (p1.y == p2.y && p1.x > p2.x) || p1.y > p2.y;
    }
};

struct Point2DComparator2 {
    bool operator()(const Point2D &p1, const Point2D &p2) {
        return (p1.y == p2.y && p1.x < p2.x) || p1.y < p2.y;
    }
};

struct EventPtrComparator {
    Point2DComparator point_cmp;
    bool operator()(const EventPtr &e1, const EventPtr &e2) {
        return point_cmp(e1->point, e2->point);
    }
};

EventPtr checkCircleEvent(bl::BLNodePtr n1, bl::BLNodePtr n2, bl::BLNodePtr n3,const std::vector<Point2D> &points, double sweepline) {
    
    if (n1 == nullptr || n2 == nullptr || n3 == nullptr)
        return nullptr;
    
    Point2D p1 = points[n1->get_id()];
    Point2D p2 = points[n2->get_id()];
    Point2D p3 = points[n3->get_id()];
    Point2D center, bottom;
    
    if (p2.y > p1.y && p2.y > p3.y)
        return nullptr;
    
    if (!findCircleCenter(p1, p2, p3, center))
        return nullptr;
    
    bottom = center;
    bottom.y += (center - p2).norm();
    
    // check circle event
    if (fabs(bottom.y - sweepline) < POINT_EPSILON || sweepline < bottom.y) {
        // create a circle event structure
        EventPtr e = std::make_shared<Event>(-1, Event::CIRCLE, bottom);
        // initialize attributes
        e->center = center;
        e->arc = n2;
        // add reference in the corresponding node
        n2->circle_event = e;
        return e;
    }
    
    return nullptr;
}

void build_voronoi(const std::vector<Point2D> &points,
                   std::vector<bl::HalfEdgePtr> &halfedges,
                   std::vector<bl::VertexPtr> &vertices,
                   std::vector<bl::HalfEdgePtr> &faces) {
    
    // create a priority queue
    std::priority_queue<EventPtr, std::vector<EventPtr>, EventPtrComparator> pq;
    
    // initialize it with all site events
    for (size_t i = 0; i < points.size(); ++i) {
        pq.push(std::make_shared<Event>(static_cast<int>(i), Event::SITE, points[i]));
    }
    
    // initialize vector of halfedges for faces
    faces.resize(points.size(), nullptr);
    
    // create a beachline tree
    bl::BLNodePtr root;
    double sweepline = 0L; // current position of the sweepline
    
    // process events
    while (!pq.empty()) {    
        // extract new event from the queue
        EventPtr e = pq.top(); pq.pop();
        
        // set position of a sweepline
        sweepline = e->point.y;
        
        if (e->type == Event::SITE) { // handle site event
            
            int point_i = e->index;
            if (root == nullptr) { // init empty beachline tree
                root = std::make_shared<bl::BLNode>(std::make_pair(point_i, point_i), &sweepline, &points);
            } else { // if it's not empty
                
                bl::BLNodePtr arc = bl::find(root, e->point.x);
                bl::BLNodePtr subtree, left_leaf, right_leaf;
                
                if (arc->circle_event != nullptr) {
                    EventPtr circle_e = arc->circle_event;
                    circle_e->type = Event::SKIP; // ignore corresponding event
                }
                
                // check number of intersection points
                int isp_num = intersectionPointsNum(points[arc->get_id()], e->point, sweepline);
                
                // different subtrees depending on the number of intersection points
                if (isp_num == 1) {
                    subtree = bl::make_simple_subtree(point_i, arc->get_id(), &sweepline, &points, halfedges);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right;
                } else if (isp_num == 2) {
                    subtree = bl::make_subtree(point_i, arc->get_id(), &sweepline, &points, halfedges);
                    left_leaf = subtree->left;
                    right_leaf = subtree->right->right;
                } else {
                    continue;
                }
                
                if (arc->prev != nullptr)
                    bl::connect(arc->prev, left_leaf);
                
                if (arc->next != nullptr)
                    bl::connect(right_leaf, arc->next);
                
                // Replace old leaf with a subtree and rebalance it
                root = bl::replace(arc, subtree);
                
                // Check circle events
                EventPtr circle_event = checkCircleEvent(left_leaf->prev, left_leaf, left_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(right_leaf->prev, right_leaf, right_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
            }
            
        } else if (e->type == Event::CIRCLE) { // handle circle event
            
            bl::BLNodePtr arc = e->arc, prev_leaf, next_leaf;

            // get breakpoint nodes
            std::pair<bl::BLNodePtr, bl::BLNodePtr> breakpoints = bl::breakpoints(arc);
            
            // recheck if it's a false alarm 1
            if (breakpoints.first == nullptr || breakpoints.second == nullptr) {
                continue;
            }
            
            // recheck if it's a false alarm 2
            double v1 = breakpoints.first->value(), v2 = breakpoints.second->value();
            
            if (fabs(v1 - v2) > BREAKPOINTS_EPSILON) {
                continue;
            }
            
            // create a new vertex and insert into doubly-connected edge list
            bl::VertexPtr vertex = std::make_shared<bl::Vertex>(e->center);
            bl::HalfEdgePtr h_first = breakpoints.first->edge;
            bl::HalfEdgePtr h_second = breakpoints.second->edge;
            
            // store vertex of Voronoi diagram
            vertices.push_back(vertex);
            
            // remove circle event corresponding to next leaf
            if (arc->prev != nullptr && arc->prev->circle_event != nullptr) {
                EventPtr circle_e = arc->prev->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }
            
            // remove circle event corresponding to prev leaf
            if (arc->next != nullptr && arc->next->circle_event != nullptr) {
                EventPtr circle_e = arc->next->circle_event;
                circle_e->type = Event::SKIP; // ignore corresponding event
            }
            
            // store pointers to the next and previous leaves
            prev_leaf = arc->prev;
            next_leaf = arc->next;
            
            // They should not be null
            assert(prev_leaf != nullptr);
            assert(next_leaf != nullptr);
            
            // get node associated with a new edge
            bl::BLNodePtr new_edge_node;
            if (arc->parent == breakpoints.first)
                new_edge_node = breakpoints.second;
            else
                new_edge_node = breakpoints.first;
            
            // remove arc from the beachline
            root = bl::remove(arc);
            
            // make a new pair of halfedges
            std::pair<bl::HalfEdgePtr, bl::HalfEdgePtr> twin_nodes = bl::make_twins(prev_leaf->get_id(), next_leaf->get_id());
            new_edge_node->edge = twin_nodes.first;
            //1/ new_edge_node->edge = twin_nodes.first;

            // connect halfedges
            bl::connect_halfedges(h_second, h_first->twin);
            bl::connect_halfedges(h_first, twin_nodes.first);
            bl::connect_halfedges(twin_nodes.second, h_second->twin);
            
            // halfedges are pointing into a vertex  -----> O <-----
            // not like this <---- O ----->
            // counterclockwise
            h_first->vertex = vertex;
            h_second->vertex = vertex;
            twin_nodes.second->vertex = vertex;
            vertex->edge = h_second;
            
            halfedges.push_back(twin_nodes.first);
            halfedges.push_back(twin_nodes.second);
            
            // check new circle events
            if (prev_leaf != nullptr && next_leaf != nullptr) {
                EventPtr circle_event = checkCircleEvent(prev_leaf->prev, prev_leaf, next_leaf, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
                circle_event = checkCircleEvent(prev_leaf, next_leaf, next_leaf->next, points, sweepline);
                if (circle_event != nullptr) {
                    pq.push(circle_event);
                }
            }
        }
    }
    
    // Fill edges corresponding to faces
    for (size_t i = 0; i < halfedges.size(); ++i) {
        bl::HalfEdgePtr he = halfedges[i];
        if (he->prev == nullptr || faces[he->l_index] == nullptr) {
            faces[he->l_index] = he;
        }
    }
}
