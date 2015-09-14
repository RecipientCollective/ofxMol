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



#ifndef ESBTL_GRID_OF_CUBES_H
#define ESBTL_GRID_OF_CUBES_H


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <list>
#include <map>


namespace ESBTL{

using boost::tuple;
using boost::make_tuple;

/** Traits for using with Grid_of_cubes. Objects are atoms/points, maximum radius is fixed to 3.
  *\tparam Point_ is the underlying point type used.
  *\tparam Iterator_ is an iterator over objects to be stored in the grid.
  */
template<class Point_,class Iterator_>
struct Traits_for_grid{
  typedef double NT;
  typedef Iterator_ Iterator;
  typedef Point_ Point;
  
  /** static function indicating how the grid must be initialized (the bounding box)
    *\param begin first iterator over the range of objects to be inserted into the grid.
    *\param end past-end iterator over the range of objects to be inserted into the grid.
    *\param xlim is the number of cubes in the x-direction. This value will be defined by that function.
    *\param ylim is the number of cubes in the y-direction. This value will be defined by that function.
    *\param zlim is the number of cubes in the z-direction. This value will be defined by that function.
    *\param rmax is an upper bound on the radius of the atoms. This value will be defined by that function.
    *\param lower_corner is a point indicating the lower corner of the grid. This value will be defined by that function.
    */
  static void init_grid(const Iterator& begin,const Iterator& end,int& xlim,int& ylim,int& zlim,double& rmax,Point& lower_corner){
    Point min,max;
    boost::tie(min,max)=ESBTL::bounding_box<Point>(begin,end);
    
    //#warning should not be hardcoded
    rmax = 3.;
    xlim=static_cast<int>( ceil( ( max.x()-min.x() + 2*rmax )/rmax) );
    ylim=static_cast<int>( ceil( ( max.y()-min.y() + 2*rmax )/rmax) );
    zlim=static_cast<int>( ceil( ( max.z()-min.z() + 2*rmax )/rmax) );
    lower_corner=min;
  };
  
  /**
    * static function computing the position of an object inside the grid.
    *\tparam Any_iterator is an iterator type over object to be stored in the grid.
    *\param it is an iterator over an object.
    *\param xlim is the number of cubes in the x-direction.
    *\param ylim is the number of cubes in the y-direction.
    *\param zlim is the number of cubes in the z-direction.
    *\param cubel is the edge length of a cube in the grid.
    *\param lower_corner is a point indicating the lower corner of the grid.
    *\return a tuple giving the cube position of the object in the grid.
    */
  template<class Any_iterator>
  static boost::tuple<int,int,int> 
  locate_cube(Any_iterator it,const int& xlim,const int& ylim,const int& zlim,const double& cubel,const Point& lower_corner){
    int x=static_cast<int>( floor( ( it->x()-lower_corner.x() )/cubel ) );
    int y=static_cast<int>( floor( ( it->y()-lower_corner.y() )/cubel ) );
    int z=static_cast<int>( floor( ( it->z()-lower_corner.z() )/cubel ) );
    
    //handle limit cases: snapping point at the limit of the bounding box and possibly
    //intersecting an element of the cube. This is not the optimal method as we are making more test than needed
    //(we know that some cubes are free from intersection) : maybe allow this in the iterator: just modify the is_outside_grid
    //criteria and allow iterator for *lim and -1
    if (x==xlim) x=xlim-1;
    if (y==ylim) y=ylim-1;
    if (z==zlim) z=zlim-1;
    if (x==-1) x=0;
    if (y==-1) y=0;
    if (z==-1) z=0;
    
    assert(!(x<0) && x < xlim);
    assert(!(y<0) && y < ylim);
    assert(!(z<0) && z < zlim);
    return boost::make_tuple(x,y,z);
  };
  
  /**
    * static function indicating if an object can be in conflict with a grid object.
    *\param it is an iterator over an object.
    *\param xlim is the number of cubes in the x-direction.
    *\param ylim is the number of cubes in the y-direction.
    *\param zlim is the number of cubes in the z-direction.
    *\param cubel is the edge length of a cube in the grid.
    *\param lower_corner is a point indicating the lower corner of the grid.
    *\return a tuple giving the cube position of the object in the grid.
    */  
  template<class Any_iterator>
  static bool 
  is_outside_grid(Any_iterator it,const int& xlim,const int& ylim,const int& zlim,const double& cubel,const Point& lower_corner){
    int x=static_cast<int>( floor( ( it->x()-lower_corner.x() )/cubel ) );
    int y=static_cast<int>( floor( ( it->y()-lower_corner.y() )/cubel ) );
    int z=static_cast<int>( floor( ( it->z()-lower_corner.z() )/cubel ) );
    if (x<-1 || x > xlim) return true;
    if (y<-1 || y > ylim) return true;
    if (z<-1 || z > zlim) return true;
    return false;
  };
};  
  
/** This class offers an intersection detection between atoms using a grid.
  * The bounding box of all the atoms is subdivided using a regular grid of cubes
  * (which edge length is twice the largest radius + epsilon). Testing intersection between two atoms
  * then reduces to test intersection between atoms in the same cube or in a direct neighboring cube.
  * \tparam Traits is a class providing static function indicating how to handle stored objects (see Traits_for_grid for an example.)
  */
template<class Traits>
struct Grid_of_cubes
{
  typedef tuple<int,int,int> Cube_coordinates;
  
  typedef typename Traits::Iterator Object_iterator;
  
  struct Cube_unit{
    typedef typename std::list<Object_iterator>::iterator In_cube_iterator;
    std::list<Object_iterator> objects;
    Cube_unit(Object_iterator V){ objects.push_front(V); }
    void insert(Object_iterator V){ objects.push_front(V); }
    unsigned size() const { return objects.size(); }
    In_cube_iterator begin(){ return objects.begin(); }
    In_cube_iterator end(){ return objects.end(); }
  };

  //Iterate over objects in a cube
  typedef typename Grid_of_cubes<Traits>::Cube_unit::In_cube_iterator In_cube_iterator;
  //use pointer to allow to delete a unit cube
  typedef std::map<Cube_coordinates,Cube_unit*> Cube_container;

  //Data members
  Cube_container cube_container;//contains cubes from 0 to *lim-1 in each direction
  int xlim,ylim,zlim;//nb of cubes in each direction
  typename Traits::NT cube_edge_length;
  typename Traits::Point lower_corner;
  
  
  //functions
  void init(const Object_iterator& it_beg,const Object_iterator& it_end){
    Traits::init_grid(it_beg,it_end,xlim,ylim,zlim,cube_edge_length,lower_corner);
  }
  
  void fill(Object_iterator it_beg,Object_iterator it_end){
    for(Object_iterator it=it_beg;it!=it_end;++it)
      insert_in_cube(it);
  }
  
  Grid_of_cubes(){}
  
  Grid_of_cubes(const Object_iterator& it_beg,const Object_iterator& it_end){
    init(it_beg,it_end);
    fill(it_beg,it_end);
  }

  ~Grid_of_cubes(){
     for(typename Cube_container::iterator it=cube_container.begin();it!=cube_container.end();++it)
       delete (*it).second;
  }
  
  
  unsigned nb_element() const {
    unsigned i=0;
    for(typename Cube_container::const_iterator it=cube_container.begin();it!=cube_container.end();++it)
       i+=it->second->size();
    return i;
  }
  
  Cube_unit* 
  get_cube(const Cube_coordinates& t) {
    typename Cube_container::iterator it=cube_container.find(t);
    return (it==cube_container.end())?(NULL):(it->second);
  }
  
  template <class Any_iterator>
  Cube_coordinates 
  locate_cube(Any_iterator v) const {
    return Traits::locate_cube(v,xlim,ylim,zlim,cube_edge_length,lower_corner);
  }
  
  template <class Any_iterator>
  Cube_unit* 
  get_cube(Any_iterator v) {
    return get_cube( locate_cube(v) );
  }
  
  void insert_in_cube(Object_iterator V){
    Cube_coordinates T=Traits::locate_cube(V,xlim,ylim,zlim,cube_edge_length,lower_corner);
    typename Cube_container::iterator it=cube_container.find(T);
    if (it==cube_container.end())
      cube_container[T]=new Cube_unit(V);
    else
      (*it).second->insert(V);
  };
  
  template <class Any_iterator>
  bool is_outside_grid(Any_iterator v) const {
    return Traits::is_outside_grid(v,xlim,ylim,zlim,cube_edge_length,lower_corner);
  }
  
  bool valid_tuple(const Cube_coordinates& T) const {
    return cube_container.find(T)!=cube_container.end();
  }
  
  class object_iterator;
  
  class iterator{
    friend class Grid_of_cubes<Traits>::object_iterator;    
    friend void Grid_of_cubes<Traits>::erase(const object_iterator& it);      
    protected:
    Cube_coordinates current;
    Grid_of_cubes* grid_ptr;
    Cube_coordinates 
    get_next(){
      assert(current!=make_tuple(-1,-1,-1));//DEBUG
      Cube_coordinates ret=current;
      if (current.get<0>() < grid_ptr->xlim-1)
        ret.get<0>()++;
      else{
        if (current.get<1>() < grid_ptr->ylim-1){
          ret.get<1>()++;
          ret.get<0>()=0;
        }
        else{
          if (current.get<2>() < grid_ptr->zlim-1){
          ret.get<2>()++;
          ret.get<0>()=0;
          ret.get<1>()=0;
          }
          else{
            return Cube_coordinates(-1,-1,-1);
          }
        }
      }
      if (grid_ptr->valid_tuple(ret)) return ret;
      current=ret;
      return get_next();
    }
    public:
    iterator(){};
    iterator(const Cube_coordinates& C):current(C),grid_ptr(NULL){};
    iterator(Grid_of_cubes* Grd):grid_ptr(Grd){current=make_tuple(-1,-1,-1);};
    iterator(Grid_of_cubes* Grd,const Cube_coordinates& C):current(C),grid_ptr(Grd){};
    iterator& operator++(){
      current=get_next();
      assert(grid_ptr->valid_tuple(current) || current==make_tuple(-1,-1,-1));//DEBUG
      return *this;
    }
    Cube_unit& operator*(){
      typename Cube_container::iterator it=grid_ptr->cube_container.find(current);
      assert(current!=make_tuple(-1,-1,-1));//DEBUG
      assert(it!=grid_ptr->cube_container.end());//DEBUG
      return *(*it).second;
    }
    Cube_unit* operator->(){
      typename Cube_container::iterator it=grid_ptr->cube_container.find(current);
      assert(it!=grid_ptr->cube_container.end());//DEBUG
      return (*it).second;
    }
    bool operator==(const iterator& it){
      return this->current==it.current;
    }
    bool operator!=(const iterator& it){
      return this->current!=it.current;
    }
  };
  
  class neighbor_iterator : public iterator
  {
    using iterator::current;
    using iterator::grid_ptr;
    private:
    Cube_coordinates center;
    
    Cube_coordinates get_next(){
      Cube_coordinates ret=current;
      assert(current!=make_tuple(-1,-1,-1));//DEBUG
      if ( (current.template get<0>() < center.get<0>() + 1) && (current.template get<0>() < grid_ptr->xlim - 1) )
        ret.get<0>()++;
      else{
        if ( (current.template get<1>() < center.get<1>() + 1 ) && (current.template get<1>() < grid_ptr->ylim-1) ){
          ret.get<1>()++;
          ret.get<0>()=((center.get<0>() == 0)?(center.get<0>()):(center.get<0>()-1));
        }
        else{
          if ( (current.template get<2>() < center.get<2>()+1) && (current.template get<2>() < grid_ptr->zlim-1) ){
            ret.get<2>()++;
            ret.get<0>()=((center.get<0>() == 0)?(center.get<0>()):(center.get<0>()-1));
            ret.get<1>()=((center.get<1>() == 0)?(center.get<1>()):(center.get<1>()-1));
          }
          else{
            return make_tuple(-1,-1,-1);
          }
        }
      }
      if (ret==center){
        current=ret;
        return get_next();
      }
      if (grid_ptr->valid_tuple(ret)) return ret;
      current=ret;
      return get_next();  
    }
    
    public:
    neighbor_iterator(){}
    neighbor_iterator(const Cube_coordinates& C):iterator(C){}
    neighbor_iterator(const Cube_coordinates& C,Grid_of_cubes* Grd):iterator(Grd),center(C){
      current=grid_ptr->get_first_neighbor(center);
      if (!grid_ptr->valid_tuple(current))
        current=this->get_next();
    } 
    
    Cube_coordinates& get_current(){
      return current;
    }
    
    neighbor_iterator& operator++(){
      current=this->get_next();
      return *this;
    }
  };
  
  iterator begin(){
    iterator ret=iterator(this,make_tuple(0,0,0));
    return (this->valid_tuple(make_tuple(0,0,0)))?(ret):(++ret);
  }
  iterator end(){
    return iterator(this,make_tuple(-1,-1,-1));
  }
  neighbor_iterator nend(){
    return neighbor_iterator(make_tuple(-1,-1,-1));
  }
  
  class object_iterator{
    friend void Grid_of_cubes<Traits>::erase(const object_iterator& it);
    protected:
    iterator current_cube;
    In_cube_iterator current_object;

    public:
    object_iterator(){};
    object_iterator(iterator it,In_cube_iterator itc):current_cube(it),current_object(itc){};
    object_iterator(iterator it):current_cube(it),current_object(NULL){};

    object_iterator& operator++(){
      if (++current_object==(*current_cube).end()){
        ++current_cube;
        if (current_cube!=current_cube.grid_ptr->end())
          current_object=(*current_cube).begin();
      }
      return *this;
    }
    
    //postfix 
    object_iterator operator++(int ){
      object_iterator tmp=*this;
      if (++current_object==(*current_cube).end()){
        ++current_cube;
        if (current_cube!=current_cube.grid_ptr->end())
          current_object=(*current_cube).begin();
      }
      return tmp;
    }    
    
    Object_iterator& operator*(){ return (*current_object);};
    Cube_unit* operator->(){ return *(*current_object);   };
    bool operator==(const object_iterator& it){
      bool ret=(this->current_cube==it.current_cube);
      if (!ret)
        return false;
      return (this->current_cube==this->current_cube.grid_ptr->end())?(true):(this->current_object==it.current_object);
    }
    bool operator!=(const object_iterator& it){
      bool ret=(this->current_cube!=it.current_cube);
      if (ret)
        return true;
      return (this->current_cube==this->current_cube.grid_ptr->end())?(false):(this->current_object!=it.current_object);      
    }
  };
  
  object_iterator last_object(){return object_iterator(this->end());}
  object_iterator first_object(){
    if (this->begin()!=this->end())
      return object_iterator(this->begin(),(*(this->begin())).begin());
    else
      return object_iterator(this->end());
  }
  
  void erase(const object_iterator& it){
		#ifndef NDEBUG
    int k=get_cube(it.current_cube.current)->size();//DEBUG
		#endif
    get_cube(it.current_cube.current)->objects.erase(it.current_object);
    #ifndef NDEBUG
    assert(k==(int) get_cube(it.current_cube.current)->objects.size()+1);//DEBUG
    #endif
    if (get_cube(it.current_cube.current)->objects.size()==0){
      typename Cube_container::iterator Iterator=cube_container.find(it.current_cube.current);
      delete (*Iterator).second;
      cube_container.erase(Iterator);
    }
  }
  
  Cube_coordinates get_first_neighbor(Cube_coordinates& T)
  {
    int a,b,c;
    int pb=0;
    if (T.get<0>()>0)
      a=T.get<0>()-1;
    else{
      a=T.get<0>();
      ++pb;
    }
    if (T.get<1>()>0)
      b=T.get<1>()-1;
    else{
      b=T.get<1>();
      ++pb;
    }
    if (T.get<2>()>0)
      c=T.get<2>()-1;
    else{
      c=T.get<2>();
      ++pb;
    }    
    return (pb==3)?(Cube_coordinates(1,0,0)):(Cube_coordinates(a,b,c));
  }  
  
  neighbor_iterator first_neighbor(const Cube_coordinates& C) {
    return neighbor_iterator(C,this);
  }
  
  neighbor_iterator first_neighbor(iterator& it){
    return neighbor_iterator(it.get_current(),this);
  }
};

} //namespace ESBTL

#endif //ESBTL_GRID_OF_CUBES_H
