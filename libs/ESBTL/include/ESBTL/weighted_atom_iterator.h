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



#ifndef WEIGHTED_ATOM_ITERATOR_H
#define WEIGHTED_ATOM_ITERATOR_H

#include <ESBTL/iterators.h>

//this may be redundant with boost::transform_iterator. Again the default constructor of the 
//functor is not sufficient for use

namespace ESBTL{
  
//iterator over all atoms of the model and create weighted atoms
/** 
  * Class providing an iterator over atoms of a model that returns a
  * weighted point when dereferencing. This is particularly useful when using CGAL
  * and its algorithms working on weighted points.
  * \tparam Model is the type of model used.
  * \tparam Weighted_atom is the weighted point type returned.
  * \tparam Weight_functor is a function object which operator() associate a radius to an atom (as Weight_of_atoms for example.)
  *
  * An alternative to this class is to use the class boost::transform_iterator, provided your function object does not
  * need only information available at compile time.
  */
template <class Model,class Weighted_atom,class Weight_functor>
class Weighted_atom_iterator:
public boost::iterator_facade<
    Weighted_atom_iterator<Model,Weighted_atom,Weight_functor>,
    Weighted_atom,
    boost::forward_traversal_tag,
    Weighted_atom >
{
private:
  typedef typename Model::System::Atom   Atom;

public:

  /** Constructor.
    * \tparam Functor_parameters is a parameter type for the function object constructor.
    * \param it is the corresponding atom iterator.
    * \param params is a parameter given to the function object constructor.
    */
  template<class Functor_parameters>
  Weighted_atom_iterator(typename Model::Atoms_const_iterator it,Functor_parameters params):atom_iterator(it),weight(params){}

private:
  friend class boost::iterator_core_access;

  void increment(){
    ++atom_iterator;
  }
    
  bool equal(Weighted_atom_iterator const& other) const{
    return atom_iterator==other.atom_iterator;
  }
  
  Weighted_atom dereference() const {
    return Weighted_atom(*atom_iterator,weight(*atom_iterator));
  }
    
//data members
  typename Model::Atoms_const_iterator atom_iterator;
  Weight_functor weight;
};  
  
}//namespace ESBTL

#endif //WEIGHTED_ATOM_ITERATOR_H
