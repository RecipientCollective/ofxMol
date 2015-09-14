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



#ifndef ESBTL_COARSE_GRAIN_H
#define ESBTL_COARSE_GRAIN_H

#include <ESBTL/molecular_system.h>

namespace ESBTL{

/** A coarse-grain atom type 
  * \tparam Atom is the base atom class that the coarse atom represent.
  * \tparam Point is a point type with coordinates const access methods x(), y() and z().
  */
template <class Atom,class Point>
class Coarse_atom: public Point{
  typedef typename Atom::System::Residue Residue;
  //list of atoms that contributes to the coarse atom
  std::list<const Atom*> atoms_;
  unsigned index_;
  const Residue* residue_;
  
public:
  typedef Point Point_3;

  /** Add atom to the set of atom represented.*/
  void add(const Atom& atom){atoms_.push_back(&atom);}
  
  Coarse_atom(const Point_3& p,unsigned i,const Residue& res):Point_3(p),index_(i),residue_(&res){}
  Coarse_atom(unsigned i,const Residue& res):Point_3(),index_(i),residue_(&res){}
    
  Coarse_atom():Point_3(),index_(0),residue_(NULL){}
  Coarse_atom(float x,float y,float z):Point_3(x,y,z),index_(0),residue_(NULL){}    
  
  const Residue& residue() const {return *residue_;}
  
  unsigned index() const {return index_;}
};

/** A coarse-grain residue type 
  * \tparam Residue is a base residue class.
  * \tparam Chain is a chain type.
  * \tparam Coarse_atom_ is the Coarse-grain atom type used.
  */
template <class Residue,class Chain,class Coarse_atom_>
class Coarse_residue:public Residue{
  std::vector<Coarse_atom_> coarse_atoms_container;
  typedef Coarse_residue<Residue,Chain,Coarse_atom_> Self;
public:
  typedef std::vector<Coarse_atom_> Coarse_atom_container;
  typedef Coarse_atom_ Coarse_atom;
  typedef typename Residue::Atom  Atom;

  template<class Line_format>
  Coarse_residue(const Line_format& line_format, const std::string& line,const Chain& ch):
    Residue(line_format,line,ch){}
  Coarse_residue(const std::string& resname,int index,char insc,const Chain& ch):
    Residue(resname,index,insc,ch){}

      
  /**
    * Method to insert the coarse grain atoms into that residue,
    * using a coarse grain atom creator.
    * \tparam A coarse grain atom creator model of the concept of \ref coarse_creator.
    */      
  template <class Coarse_creator>
  int create_coarse_atoms(const Coarse_creator& creator){
    return creator(*this,std::back_inserter(coarse_atoms_container));
  }

  /**
    * Method to insert a coarse grain atoms using a point.
    * \param pt is the coordinates of the pseudo-atom to be added.
    * \param i is the index of the pseudo-atom (when a residue is modeled by several).
    */        
  void add_coarse_atom(const typename Residue::Atom::Point_3& pt,unsigned i=0){
    this->coarse_atoms_container.push_back(Coarse_atom(pt,i,*this));
  }
  
  const Coarse_atom get_coarse_atom(unsigned i) const {
    return coarse_atoms_container[i];
  }
  
  //TODO find another way to to this except a simple copy from the base class
  //iterators
  //---------
  //~ //functions needed by generic iterator in iterators.h 
  //~ static const Atom& dereference (typename Atom_container::const_iterator it) {return it->second;}
  //~ static Atom& dereference (typename Atom_container::iterator it) {return it->second;}
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_residue<Self,true>  Atoms_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_residue<Self,false> Atoms_iterator;
  
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_begin() {return Atoms_iterator(this->atom_container_.begin(),this->atom_container_.end());}
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_end()   {return Atoms_iterator(this->atom_container_.end());}  

  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_begin() const {return Atoms_const_iterator(this->atom_container_.begin(),this->atom_container_.end());}
  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_end()   const {return Atoms_const_iterator(this->atom_container_.end());}
  //-------    
  /** \ingroup grp_iters*/
  typedef typename std::vector<Coarse_atom_>::const_iterator Coarse_atom_const_iterator;
  /** \ingroup grp_iters*/
  typedef typename std::vector<Coarse_atom_>::iterator Coarse_atom_iterator;
  /** \ingroup grp_iters*/
  Coarse_atom_const_iterator coarse_atoms_begin() const {return coarse_atoms_container.begin();}
  /** \ingroup grp_iters*/
  Coarse_atom_const_iterator coarse_atoms_end()   const {return coarse_atoms_container.end();}
  /** \ingroup grp_iters*/
  Coarse_atom_iterator coarse_atoms_begin() {return coarse_atoms_container.begin();}
  /** \ingroup grp_iters*/
  Coarse_atom_iterator coarse_atoms_end()   {return coarse_atoms_container.end();}  
  
  using Residue::dereference;
  static const Coarse_atom& dereference(Coarse_atom_const_iterator it){ return *it;};
  static Coarse_atom& dereference(Coarse_atom_iterator it){ return *it;};
};


/** Inserts a set of coarse grain atoms in a molecule using one residue per such atom.
  * \tparam Input_iterator is an iterator type over point type objects. 
  * This type must be the one used to defined the coarse grain atom type in System.
  * \tparam System is a coarse-grain system.
  * \param first is the first iterator over the range of points.
  * \param last is the past-end iterator over the range of points.
  * \param system is the system in which coarse grain atoms must be inserted.
  * \param modelid is the number of the model of the system in which coarse grain atoms must be inserted.
  * \param chainid is the id of the chain of the model in which coarse grain atoms must be inserted.
  * \param resname is the name of the residues in which coarse grain atoms must be inserted.
  * \param starting_res_index is the residue sequence number of the first residue inserted. 
  * Other coarse grain atoms will use an increasing number starting from this value.
  */
template<class Input_iterator,class System>
void insert_coarse_atoms(Input_iterator first,Input_iterator last,System& system,int modelid=1,char chainid='Z',std::string resname="SOL",int starting_res_index=1){
  typename System::Model& model=system.get_or_create_model(modelid);
  typename System::Chain& chain=model.get_or_create_chain(chainid);
  int res_index=starting_res_index-1;
  char insc=' ';
  for (Input_iterator it=first; it!=last;++it){
    typename System::Residue& residue=chain.get_or_create_residue(resname,++res_index,insc);
    residue.add_coarse_atom(*it,0);
  }
}


/** Object defining iterator types over coarse grain residue of a model
  * \tparam Model is the type of the model coarse grain atoms are read.
  * \ingroup grp_iters
  */
template <class Model>
struct Coarse_atoms_iterators{
  typedef internal::Coarse_atom_from_model<Model,true>  const_iterator;
  typedef internal::Coarse_atom_from_model<Model,false>  iterator;
};

/** Return the first iterator over coarse grain atoms of a model (const version).
  * \ingroup grp_iters*/
template <class Model>
typename Coarse_atoms_iterators<Model>::const_iterator coarse_atoms_begin(const Model& model){
  return typename Coarse_atoms_iterators<Model>::const_iterator(model);
}

/** Return the past-end iterator over coarse grain atoms of a model (const version).
  * \ingroup grp_iters*/
template <class Model>
typename Coarse_atoms_iterators<Model>::const_iterator coarse_atoms_end(const Model& model){
  return typename Coarse_atoms_iterators<Model>::const_iterator(model,true);
}

/** Return the first iterator over coarse grain atoms of a model.
  * \ingroup grp_iters*/
template <class Model>
typename Coarse_atoms_iterators<Model>::iterator  coarse_atoms_begin(Model& model){
  return typename Coarse_atoms_iterators<Model>::iterator(model);
}

/** Return the past-end iterator over coarse grain atoms of a model.
  * \ingroup grp_iters*/
template <class Model>
typename Coarse_atoms_iterators<Model>::iterator coarse_atoms_end(Model& model){
  return typename Coarse_atoms_iterators<Model>::iterator(model,true);
}

/** A default coarse-grain system items gathering wrappers
  * (needed to break the cyclic dependency) 
  * that allow to define a model, chain, residue and atom type
  * from the types Molecular_model, Molecular_chain, Coarse_residue and Molecular_atom respectively.
  */
struct System_items_with_coarse_grain{
  template <class System,class Point_3>
  struct Model_wrapper{
    typedef Molecular_model<System>         Type;
  };
  
  template <class System,class Point_3>
  struct Chain_wrapper{
    typedef Molecular_chain<System>         Type;
  };
  
  template <class System,class Point_3>
  struct Atom_wrapper{
    typedef Molecular_atom<System,Point_3>  Type;
  };  
  
  template <class System,class Point_3>
  class Residue_wrapper{
    typedef Molecular_residue<System>         Residue;
    typedef Molecular_atom<System,Point_3>    Atom;
    typedef Molecular_chain<System>           Chain;
  public:
    typedef Coarse_residue<Residue,Chain,Coarse_atom<Atom,Point_3> > Type;
  };
  
};

}//namespace ESBTL

#endif //ESBTL_COARSE_GRAIN_H
