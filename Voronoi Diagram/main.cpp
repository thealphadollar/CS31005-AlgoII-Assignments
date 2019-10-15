#include<bits/stdc++.h>
#include<memory>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include "Point2D.h"
#include "VoronoiDiagram.hpp"
#include "Beachline.hpp"
#define MAX 500
double x_max=0,y_max=0;
double x_min=0,y_min=0;
namespace bl = beachline;

void plot_circle(const Point2D &c, double r, std::ofstream & fout , double shift_x, double shift_y) {

    fout<<"<circle cx=\""<<c.x+shift_x<<"\" "<<"cy=\""<<c.y+shift_y<<"\" "<<"r=\""<<r<<"\" stroke=\"#0066ff\" stroke-width=\"0\" fill=\"#0066ff\" fill-opacity=\"0.4\" />\n";   
}

void plot_point(const Point2D &c, std::ofstream & fout , double shift_x, double shift_y){
        
 fout<<"<circle cx=\""<<c.x+shift_x<<"\" "<<"cy=\""<<c.y+shift_y<<"\" "<<"r=\""<<2<<"\" stroke=\"black\" stroke-width=\"0\" fill=\"#000000\" />\n";
}


double random_number() {
    return (double(rand())*MAX)/ double(RAND_MAX);
}

// random number generator that takes current second as the seed
std::vector<Point2D> randomPoint(int number) {
    srand(static_cast<unsigned int>(time(0)));
    std::vector<Point2D> points;
    for (int i = 0; i < number; ++i) {
        double x = random_number(), y = random_number();
        int id=i+1;
        Point2D p;
        p.x=x;
        p.y=y;
        p.id=id;
        points.push_back(p);
    }
    return points;
}


void initEdgePointsVis(bl::HalfEdgePtr h, std::vector<double> &x, std::vector<double> &y,
                       const std::vector<Point2D> &points) { 
    
    if (h->vertex != nullptr && h->twin->vertex != nullptr) {
        
        x[0] = h->vertex->point.x;
      
        x[1] = h->twin->vertex->point.x;
     
        y[0] = h->vertex->point.y;
        
        y[1] = h->twin->vertex->point.y;
        
      auto  temp1=std::max_element(x.begin(),x.end());
      auto  temp2=std::min_element(x.begin(),x.end());        
        if(x_max< *temp1)
         x_max= (*temp1);
        if(x_min> *temp2)
         x_min= (*temp2);
        temp1=std::max_element(y.begin(),y.end());
        temp2=std::min_element(y.begin(),y.end());        
        if(y_max< *temp1)
         y_max= (*temp1);
        if(y_min> *temp2)
         y_min= (*temp2);         
         
        
    } else if (h->vertex != nullptr) {
        
        x[0] = h->vertex->point.x;
        y[0] = h->vertex->point.y;
        
        Point2D norm = (points[h->l_index] - points[h->r_index]).normalized().getRotated90CCW();
        x[1] = x[0] + norm.x * 1000;
        y[1] = y[0] + norm.y * 1000;
      auto  temp1=std::max_element(x.begin(),x.end());
      auto  temp2=std::min_element(x.begin(),x.end());        
        if(x_max< *temp1)
         x_max= (*temp1);
        if(x_min> *temp2)
         x_min= (*temp2);
        temp1=std::max_element(y.begin(),y.end());
        temp2=std::min_element(y.begin(),y.end());        
        if(y_max< *temp1)
         y_max= (*temp1);
        if(y_min> *temp2)
         y_min= (*temp2);         

        
    } else if (h->twin->vertex != nullptr) {
        
        x[0] = h->twin->vertex->point.x;
        y[0] = h->twin->vertex->point.y;
        
        Point2D norm = (points[h->twin->l_index] - points[h->twin->r_index]).normalized().getRotated90CCW();
        x[1] = x[0] + norm.x * 1000;
        y[1] = y[0] + norm.y * 1000;
      auto  temp1=std::max_element(x.begin(),x.end());
      auto  temp2=std::min_element(x.begin(),x.end());        
        if(x_max< *temp1)
         x_max= (*temp1);
        if(x_min> *temp2)
         x_min= (*temp2);
        temp1=std::max_element(y.begin(),y.end());
        temp2=std::min_element(y.begin(),y.end());        
        if(y_max< *temp1)
         y_max= (*temp1);
        if(y_min> *temp2)
         y_min= (*temp2); 
        
    } else {
        
        Point2D p1 = points[h->l_index], p2 = points[h->r_index];
        
        Point2D norm = (p1 - p2).normalized().getRotated90CCW();
        Point2D c = 0.5 * (p1 + p2);
        
        x[0] = c.x + norm.x *1000 ;
        x[1] = c.x - norm.x *1000 ;
        
        y[0] = c.y + norm.y * 1000;
        y[1] = c.y - norm.y * 1000;
     auto   temp1=std::max_element(x.begin(),x.end());
     auto   temp2=std::min_element(x.begin(),x.end());        
        if(x_max< *temp1)
         x_max= (*temp1);
        if(x_min> *temp2)
         x_min= (*temp2);
        temp1=std::max_element(y.begin(),y.end());
        temp2=std::min_element(y.begin(),y.end());        
        if(y_max< *temp1)
         y_max= (*temp1);
        if(y_min> *temp2)
         y_min= (*temp2); 
    }
}


int main(int argc, const char *argv[]) {
    
    int n;
    // Generate random points
    std::cout<<"Voronoi Visualisation by Shivam Kumar Jha [17CS30033]\n";
    std::cout<<"Number of sites in the diagram:\n";
    std::cin>>n;
    std::vector<Point2D> points = randomPoint(n);
  
   std::vector<bl::HalfEdgePtr> halfedges, faces;
   std::vector<bl::VertexPtr> vertices;

    std::ofstream fout;
    fout.open("t4.svg");

   build_voronoi(points, halfedges, vertices, faces);
    if(n<=1)
     {std::cout<<"No voronoi diagram :-)\n";
       return 0;
     }
     else if(n==2)
     {   for (size_t i = 0; i < halfedges.size(); ++i) {
        bl::HalfEdgePtr h = halfedges[i];
       std::vector<double> x(2, 0.0), y(2, 0.0);
       initEdgePointsVis(h, x, y, points);
    }
    
  double X_MIN=x_min;
  double X_MAX=x_max;
  double Y_MIN=y_min;
  double Y_MAX=y_max;  
  double shift_x=X_MAX-X_MIN;
  double shift_y=Y_MAX-Y_MIN;
   fout<<"<svg xmlns=\"http://www.w3.org/2000/svg\">\n<rect width=\""<<2*shift_x+100<<"\" height=\""<<2*shift_y+100<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />\n";
   for(size_t i = 0; i < points.size(); ++i){
     plot_point(points[i],fout,shift_x,shift_y);
    }
    
      for (size_t i = 0; i < halfedges.size(); ++i) {
        bl::HalfEdgePtr h = halfedges[i];
       std::vector<double> x(2, 0.0), y(2, 0.0);
       initEdgePointsVis(h, x, y, points);
       fout<<"<line x1=\""<<x[0]+shift_x<<"\" y1=\""<<y[0]+shift_y<<"\" x2=\""<<x[1]+shift_x<<"\" y2=\""<<y[1]+shift_y<<"\" style=\"stroke:rgb(255,0,0);stroke-width:1\" />\n";
    } 
      //draw circles manually
      
      double dist2;
      Point2D diff = points[0]-points[1];
      dist2= diff.norm();
      plot_circle(points[0],dist2/2.0,fout,shift_x,shift_y);
      plot_circle(points[1],dist2/2.0,fout,shift_x,shift_y);
      fout<<"</svg>\n"; 
    
     }
    else{
    
   for (size_t i = 0; i < halfedges.size(); ++i) {
        bl::HalfEdgePtr h = halfedges[i];
       std::vector<double> x(2, 0.0), y(2, 0.0);
       initEdgePointsVis(h, x, y, points);
    }
    
  double X_MIN=x_min;
  double X_MAX=x_max;
  double Y_MIN=y_min;
  double Y_MAX=y_max;  
  double shift_x=X_MAX-X_MIN;
  double shift_y=Y_MAX-Y_MIN; 
   fout<<"<svg xmlns=\"http://www.w3.org/2000/svg\">\n<rect width=\""<<2*shift_x+100<<"\" height=\""<<2*shift_y+100<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />\n";
   for(size_t i = 0; i < points.size(); ++i){
     plot_point(points[i],fout,shift_x,shift_y);
    }
    
      for (size_t i = 0; i < halfedges.size(); ++i) {
        bl::HalfEdgePtr h = halfedges[i];
       std::vector<double> x(2, 0.0), y(2, 0.0);
       initEdgePointsVis(h, x, y, points);
       fout<<"<line x1=\""<<x[0]+shift_x<<"\" y1=\""<<y[0]+shift_y<<"\" x2=\""<<x[1]+shift_x<<"\" y2=\""<<y[1]+shift_y<<"\" style=\"stroke:rgb(255,0,0);stroke-width:1\" />\n";
    } 
    
 
    //now plot the circle with least radius for each site
    
    for(size_t i = 0; i < faces.size(); ++i) //each face corresponds to a site 
    { 
      //identify the site 
      Point2D site =points[i];
      bl::HalfEdgePtr h = faces[i];
      bl::HalfEdgePtr h1=h;
      //make it corner edge if possible 
      while(h1->vertex!=nullptr)
      {h1=h1->next;
        if(h1==h)
         {break;}
      } 
      double mindistance=-1;
      double distance;
        do 
        {Point2D othersite = points[h1->r_index];
         Point2D diff = site-othersite;
            distance= diff.norm();
            if(mindistance>distance||mindistance==-1)
             {mindistance = distance;}
            h1=h1->prev;
        }while(h1->twin->vertex!=nullptr&&h1!=h);
        Point2D othersite = points[h1->r_index];
         Point2D diff = site-othersite;
            distance= diff.norm();
            if(mindistance>distance||mindistance==-1)
             {mindistance = distance;}
             
      plot_circle(site,mindistance/2.0,fout,shift_x,shift_y);
    
    
    }
    
std::cout<<"t4.svg has been generated in the current directory!\n";
    
    fout<<"</svg>\n"; 
 
    return 0;
    }
}

