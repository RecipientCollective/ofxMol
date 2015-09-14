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



#ifndef SELECTED_ATOM_ITERATOR_H
#define SELECTED_ATOM_ITERATOR_H

#include <ESBTL/iterators.h>

//TODO this seems redundant with boost::filter_iterator
// only advantage is that we can pass functor to iterator

namespace ESBTL{


//iterator over all atoms of a certain type of the model
//is_const specify whether this is a const iterator or not
  
/**
  * Class defining an iterator over atoms of a model, restricted to a certain
  * type defined by a function object.
  * \tparam Model is the type of the model used.
  * \tparam Subset_functor is a function following the concept of \ref atomsel.
  * \tparam is_const indicates the constness of the iterator.
  *
  * If you do not need to pass a state to your functor, you might want to use
  * boost::filter_iterator instead.
  */
template <class Model,class Subset_functor,bool is_const>
class Selected_atom_iterator:
public boost::iterator_facade<
    Selected_atom_iterator<Model,Subset_functor,is_const>,
    typename internal::do_add_const<typename Model::System::Atom,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Model::System::Atom   Atom;
  typedef typename internal::Atoms_iterator_from_model<Model,is_const> Base_iterator;

public:
  //functor may be initialized with info.
  /** Constructor.
    * \param it is the iterator indicating the starting point for the iterator.
    * \param functor is the functor used to select atoms.
    */
  Selected_atom_iterator(Base_iterator it,const Subset_functor& functor):atom_iterator(it),keep(functor)
  {if (!atom_iterator.is_end() && !keep(*atom_iterator)) increment();}
  //default constructor
  Selected_atom_iterator(Base_iterator it):atom_iterator(it)
  {if (!atom_iterator.is_end() && !keep(*atom_iterator)) increment();}

private:
  friend class boost::iterator_core_access;

  void increment(){
    do{
      ++atom_iterator;
    }
    while(!atom_iterator.is_end() && !keep(*atom_iterator));
  }
    
  bool equal(Selected_atom_iterator const& other) const{
    return atom_iterator==other.atom_iterator;
  }
  
  typename internal::do_add_const_to_ref<Atom,is_const>::type dereference() const {
    return *atom_iterator;
  }
    
//data members
  Base_iterator atom_iterator;
  Subset_functor keep;
};  
  

/** Short cut function to construct a Selected_atom_iterator.*/
template <class Model,class Subset_functor,bool is_const>
Selected_atom_iterator<Model,Subset_functor,is_const> 
make_selected_atom_iterator(internal::Atoms_iterator_from_model<Model,is_const> iterator,const Subset_functor& functor){
  return Selected_atom_iterator<Model,Subset_functor,is_const> (iterator,functor);
}

/** Short cut function to construct a Selected_atom_iterator.*/
template <class Subset_functor,class Model,bool is_const>
Selected_atom_iterator<Model,Subset_functor,is_const> 
make_selected_atom_iterator(internal::Atoms_iterator_from_model<Model,is_const> iterator){
  return Selected_atom_iterator<Model,Subset_functor,is_const> (iterator);
}
  
}//ESBTL

#endif //SELECTED_ATOM_ITERATOR_H
