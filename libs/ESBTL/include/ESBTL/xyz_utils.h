// Copyright (c) 2009-2010  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
//This file is part of ESBTL.
//
//ESBTL is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ESBTL is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with ESBTL.  If not, see <http://www.gnu.org/licenses/>.
//
//
//Additional permission under GNU GPL version 3 section 7
//
//If you modify this Library, or any covered work, by linking or
//combining it with CGAL (or a modified version of that library), the
//licensors of this Library grant you additional permission to convey
//the resulting work. Corresponding Source for a non-source form of
//such a combination shall include the source code for the parts of CGAL
//used as well as that of the covered work. 
//
//
//
// Author(s)     :  Sébastien Loriot



#ifndef ESBTL_XYZ_UTILS_H
#define ESBTL_XYZ_UTILS_H

#include <cmath>
#include <math.h>
#include <limits>
#include <ESBTL/constants.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

namespace ESBTL{
  
  /**
   * Basic point type class.
   */
  class Point_3{
  public:
    Point_3(double x,double y,double z):x_(x),y_(y),z_(z){}
    Point_3():x_(0),y_(0),z_(0){}
    
    DECLARE_AND_ACCESS(x,double)
    DECLARE_AND_ACCESS(y,double)
    DECLARE_AND_ACCESS(z,double)
  };

  /** 
   * Computes the square of n. 
   * @tparam NT is the number type of the input.
   * @param n is the number to be squared.
   * @return the square of n.
   */
  template <class NT>
  inline NT square(const NT& n){return n*n;}


  /** 
   * Computes the squared distance between two points. 
   * @tparam Point1 is a point type providing access functions x(),y() and z().
   * @tparam Point2 is a point type providing access functions x(),y() and z().
   * @param p1 is the first point.
   * @param p2 is the second point.
   * @return the squared distance between p1 and p2.
   */  
  template <class Point1,class Point2>
  inline double squared_distance(const Point1& p1,const Point2& p2){
    return square(p1.x()-p2.x())+square(p1.y()-p2.y())+square(p1.z()-p2.z());
  }


  /** 
   * Computes the root mean square deviation (RMSD) between two sets of atoms.
   * @pre The two sets must have the same size, and corresponding atoms must be enumerated in the same order.
   * @tparam Iterator is an iterator over a point type providing access functions x(),y() and z().
   * @param begin1 is the first element of the first range of points.
   * @param end1 is the past-end iterator over elements of the first range of points.
   * @param begin2 is the first element of the second range of points.
   * @param end2 is the past-end iterator over elements of the second range of points.
   * @return The RMSD between two sets of points.
   */  
  template <class Iterator>
  double rms_no_align(Iterator begin1,Iterator end1,Iterator begin2,Iterator end2){
    typedef boost::tuple<Iterator,Iterator>   It_tuple;
    typedef boost::zip_iterator<It_tuple>     Zip_iterator;
    double nb=0;
    double sum_sqd=0;
    for (Zip_iterator it=boost::make_zip_iterator(boost::make_tuple(begin1,begin1));it!=boost::make_zip_iterator(boost::make_tuple(end1,end2));++it){
      sum_sqd+=squared_distance(boost::get<0>(*it),boost::get<1>(*it));
      nb+=1;
    }
    return sqrt(sum_sqd/nb);
  }
  
  /**
   * Computes the bounding box containing a set of points.
   * @tparam Iterator is an iterator over a set of points.
   * @tparam Point_3 is a point type with access method x() ,y() and z().
   * @param begin first iterator of the set.
   * @param end  past end iterator of the set.
   * @return a pair containing the lower corner and the upper corner of the bounding box.
   */
  template <class Point_3,class Iterator>
  std::pair<Point_3,Point_3> 
  bounding_box(Iterator begin,Iterator end){
    double xmin=std::numeric_limits<double>::max(),ymin=std::numeric_limits<double>::max(),zmin=std::numeric_limits<double>::max();
    double xmax=std::numeric_limits<double>::min(),ymax=std::numeric_limits<double>::min(),zmax=std::numeric_limits<double>::min();
    for (Iterator it=begin;it!=end;++it){
      if (it->x() < xmin) xmin=it->x();
      if (it->x() > xmax) xmax=it->x();
      if (it->y() < ymin) ymin=it->y();
      if (it->y() > ymax) ymax=it->y();
      if (it->z() < zmin) zmin=it->z();
      if (it->z() > zmax) zmax=it->z();     
    }
    return std::make_pair(Point_3(xmin,ymin,zmin),Point_3(xmax,ymax,zmax));
  }

  /**
   * Computes a cube containing a set of points.
   * @tparam Iterator is an iterator over a set of points.
   * @tparam Point_3 is a point type with access method x() ,y() and z().
   * @param begin first iterator of the set.
   * @param end  past end iterator of the set.
   * @return a pair containing the lower corner and the edge length of the bounding cube.
   */
  template <class Point_3,class Iterator>
  inline
  std::pair<Point_3,double> 
  bounding_cube(Iterator begin,Iterator end){
    std::pair<Point_3,Point_3> bbox=bounding_box<Point_3>(begin,end);
    double edge_length=bbox.second.x()-bbox.first.x();
    if (edge_length < bbox.second.y()-bbox.first.y()) edge_length=bbox.second.y()-bbox.first.y();
    if (edge_length < bbox.second.z()-bbox.first.z()) edge_length=bbox.second.z()-bbox.first.z();
    
    return std::make_pair( bbox.first,nextafter(edge_length,std::numeric_limits<double>::max()) );
  }  
    
} //namespace ESBTL


#endif //ESBTL_XYZ_UTILS_H
