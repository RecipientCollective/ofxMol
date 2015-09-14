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



#include <cmath>

#ifndef ESBTL_FCC_LATTICE_H
#define ESBTL_FCC_LATTICE_H

namespace ESBTL {
  /** \cond */
  namespace internal_fcc_lattice{
    template <class Point_3,class FT,class Output_iterator>
    void even_layer(const FT& center_x,const FT& center_y,const FT& layer_z,const FT& unit_move,unsigned size,Output_iterator out){
      const FT x_s=center_x - (size-1)*unit_move;
      const FT y_s=center_y - (size-1)*unit_move;
      
      for (unsigned i=0;i<size;++i)
        for (unsigned j=0;j<size;++j)
          *out++=Point_3(x_s+2*i*unit_move,y_s+2*j*unit_move,layer_z);

      for (unsigned i=0;i<size-1;++i)
        for (unsigned j=0;j<size-1;++j)
          *out++=Point_3(x_s+unit_move+2*i*unit_move,y_s+unit_move+2*j*unit_move,layer_z);
    }


    template <class Point_3,class FT,class Output_iterator>
    void odd_layer(const FT& center_x,const FT& center_y,const FT& layer_z,const FT& unit_move,unsigned size,Output_iterator out){
      const FT x_s=center_x - (size-1)*unit_move;
      const FT y_s=center_y - (size-1)*unit_move;
      
      for (unsigned i=0;i<size;++i)
        for (unsigned j=0;j<size-1;++j)
          *out++=Point_3(x_s+2*i*unit_move,y_s+unit_move+2*j*unit_move,layer_z);

      for (unsigned i=0;i<size-1;++i)
        for (unsigned j=0;j<size;++j)
          *out++=Point_3(x_s+unit_move+2*i*unit_move,y_s+2*j*unit_move,layer_z);
    }

  } //namespace internal_fcc_lattice
  /** \endcond */

/**
 * Computes a cube filled with a face centered cubic lattice of balls.
 *@tparam Point_3 is a point type with access method x() ,y() and z().
 *@tparam Output_iterator is an output iterator of Point_3.
 *@param center is the center of the cube.
 *@param radius is the radius used for the balls to compute the lattice.
 *@param min_edge_length is a lower bound of the edge length of the cube.
 *@param out is an output iterator of Point_3.
 */
  template<class Point_3,class Output_iterator>
  void fcc_lattice(const Point_3& center, double radius, double min_edge_length,Output_iterator out){
    double unit_move=sqrt(2*radius);
    unsigned size = static_cast<unsigned>( ceil( min_edge_length/2./unit_move ) )+1;

    double z_s=center.z() - (size-1)*unit_move;
    for (unsigned i=0;i<size;++i)
      internal_fcc_lattice::even_layer<Point_3>(center.x(),center.y(),z_s+2*i*unit_move,unit_move,size,out);
    for (unsigned i=0;i<size-1;++i)
      internal_fcc_lattice::odd_layer<Point_3>(center.x(),center.y(),z_s+2*i*unit_move+unit_move,unit_move,size,out);    
  }
  
/**
 * Computes a cube filled with a face centered cubic lattice of balls.
 *@tparam Point_3 is a point type with access method x() ,y() and z().
 *@tparam Output_iterator is an output iterator of Point_3.
 *@param cube is a pair containing the lower corner of the cube to be filled and its edge length.
 *@param radius is the radius used for the balls to compute the lattice.
 *@param expand_value is the value of which the edge length of the cube must be increased.
 *@param out is an output iterator of Point_3.
 */
  template<class Point_3,class Output_iterator>
  inline
  void fcc_lattice(const std::pair<Point_3,double>& cube, double radius, double expand_value,Output_iterator out){
    double edge=cube.second;
    Point_3 center(cube.first.x()+edge/2.,cube.first.y()+edge/2.,cube.first.z()+edge/2.);
    return fcc_lattice(center,radius,edge+expand_value,out);
  }
    

} //namespace ESBTL

#endif //ESBTL_FCC_LATTICE_H
