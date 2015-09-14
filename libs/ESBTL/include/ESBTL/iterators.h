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



#ifndef ESBTL_ITERATORS_H
#define ESBTL_ITERATORS_H

#include <boost/iterator/iterator_facade.hpp>

namespace ESBTL{
/** \cond */  
namespace internal{
  
template <class T,bool b>
struct do_add_const_to_ref;

template <class T>
struct do_add_const_to_ref<T,true>{typedef const T& type;};

template <class T>
struct do_add_const_to_ref<T,false>{typedef T& type;};

template <class T,bool b>
struct do_add_const;

template <class T>
struct do_add_const<T,true>{typedef const T type;};

template <class T>
struct do_add_const<T,false>{typedef T type;};


template <class T,bool b>
struct do_add_const_to_ptr;

template <class T>
struct do_add_const_to_ptr<T,true>{typedef const T* type;};

template <class T>
struct do_add_const_to_ptr<T,false>{typedef T* type;};


template <class T,bool b>
struct do_get_const_iterator;

template <class T>
struct do_get_const_iterator<T,true>{typedef typename T::const_iterator type;};

template <class T>
struct do_get_const_iterator<T,false>{typedef typename T::iterator type;};
  

//In the following we define
// Chains_iterator_from_model 
// Residues_iterator_from_model
// Atoms_iterator_from_model
// Residues_iterator_from_chain
// Atoms_iterator_from_chain
// Atoms_iterator_from_residue
// Model_iterator_from_system
// Coarse_atom_from_model


//iterator over all chains of a model
//-----------------------------------
template <class Model,bool is_const>
class Chains_iterator_from_model:
public boost::iterator_facade<
    Chains_iterator_from_model<Model,is_const>,
    typename do_add_const<typename Model::Chain,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Model::System::Chain   Chain;
public:
  Chains_iterator_from_model(typename do_get_const_iterator<typename Model::Chain_container,is_const>::type it):
    chain_it(it){}
  
private:
  friend class boost::iterator_core_access;

  void increment(){++chain_it;}
  bool equal(Chains_iterator_from_model const& other) const{ return other.chain_it==chain_it; }
  typename do_add_const_to_ref<Chain,is_const>::type dereference() const {return Model::dereference(chain_it);}
  
//data member
  typename do_get_const_iterator<typename Model::Chain_container,is_const>::type chain_it;
};

//iterator over all residues of a chain
template <class Chain,bool is_const>
class Residues_iterator_from_chain:
public boost::iterator_facade<
    Residues_iterator_from_chain<Chain,is_const>,
    typename do_add_const<typename Chain::Residue,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Chain::Residue   Residue;
public:
  Residues_iterator_from_chain(typename do_get_const_iterator<typename Chain::Residue_container,is_const>::type it):
    res_it(it){}
  Residues_iterator_from_chain(){}
private:
  friend class boost::iterator_core_access;

  void increment(){++res_it;}
  bool equal(Residues_iterator_from_chain const& other) const{ return other.res_it==res_it; }
  typename do_add_const_to_ref<Residue,is_const>::type dereference() const {return Chain::dereference(res_it);}
  
//data member
  typename do_get_const_iterator<typename Chain::Residue_container,is_const>::type res_it;
};

//iterator over all residues of a model
template <class Model,bool is_const>
class Residues_iterator_from_model:
public boost::iterator_facade<
    Residues_iterator_from_model<Model,is_const>,
    typename do_add_const<typename Model::Chain::Residue,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Model::Chain::Residue   Residue;
  Chains_iterator_from_model<Model,is_const>                    chain_it;
  Chains_iterator_from_model<Model,is_const>                    chain_end;
  Residues_iterator_from_chain<typename Model::Chain,is_const>  res_it;
  Residues_iterator_from_chain<typename Model::Chain,is_const>  res_end;
public:
  Residues_iterator_from_model(typename do_add_const_to_ref<Model,is_const>::type model,bool is_end=false):
    chain_it(is_end?model.chains_end():model.chains_begin()),chain_end(model.chains_end())
  {
    if (chain_it!=chain_end){
      res_it=chain_it->residues_begin();
      res_end=chain_it->residues_end();
    }
  }
private:
  friend class boost::iterator_core_access;

  void increment(){
    ++res_it;
    if (res_it==res_end){
      do{
        ++chain_it;
        if (chain_it!=chain_end){
          res_it=chain_it->residues_begin();
          res_end=chain_it->residues_end();
        }
        else{break;}
      }while (res_it==res_end);
    }
  }
  bool equal(Residues_iterator_from_model const& other) const{
    return chain_it==other.chain_it && ( chain_it==chain_end || other.res_it==res_it ); 
  }
  typename do_add_const_to_ref<Residue,is_const>::type dereference() const {return *res_it;}
};

//iterator over all atoms of a residue
template <class Residue,bool is_const>
class Atoms_iterator_from_residue:
public boost::iterator_facade<
    Atoms_iterator_from_residue<Residue,is_const>,
    typename do_add_const<typename Residue::Atom,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Residue::Atom   Atom;
public:
  //for end
  Atoms_iterator_from_residue(typename do_get_const_iterator<typename Residue::Atom_container,is_const>::type it):
    atm_it(it),atm_end(it){}
  //for begin
  Atoms_iterator_from_residue(
    typename do_get_const_iterator<typename Residue::Atom_container,is_const>::type it,
    typename do_get_const_iterator<typename Residue::Atom_container,is_const>::type end
  ):atm_it(it),atm_end(end)
  {}

  Atoms_iterator_from_residue(){}
private:
  friend class boost::iterator_core_access;

  void increment(){
    ++atm_it; 
  }
  
  bool equal(Atoms_iterator_from_residue const& other) const{ return other.atm_it==atm_it; }
  typename do_add_const_to_ref<Atom,is_const>::type dereference() const {return Residue::dereference(atm_it);}
  
//data member
  typename do_get_const_iterator<typename Residue::Atom_container,is_const>::type atm_it;
  typename do_get_const_iterator<typename Residue::Atom_container,is_const>::type atm_end;
};

//iterator over all atoms of a chain
template <class Chain,bool is_const>
class Atoms_iterator_from_chain:
public boost::iterator_facade<
    Atoms_iterator_from_chain<Chain,is_const>,
    typename do_add_const<typename Chain::Residue::Atom,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Chain::Residue::Atom                       Atom;

  Residues_iterator_from_chain<Chain,is_const>                   res_it;
  Residues_iterator_from_chain<Chain,is_const>                   res_end;
  Atoms_iterator_from_residue<typename Chain::Residue,is_const>  atm_it;
  Atoms_iterator_from_residue<typename Chain::Residue,is_const>  atm_end;
public:
  Atoms_iterator_from_chain(){}
  Atoms_iterator_from_chain(typename do_add_const_to_ref<Chain,is_const>::type chain,bool is_end=false):
    res_it(is_end?chain.residues_end():chain.residues_begin()),res_end(chain.residues_end())
  {
    if (res_it!=res_end){
      atm_it=res_it->atoms_begin();
      atm_end=res_it->atoms_end();
      make_valid();
    }
  }
private:
  friend class boost::iterator_core_access;

  void increment(){
    ++atm_it;
    make_valid();
  }
  
  void make_valid(){
    while(atm_it==atm_end){
      ++res_it;
      if (res_it==res_end) break;
      atm_it=res_it->atoms_begin();
      atm_end=res_it->atoms_end();
    }    
  }
  
  bool equal(Atoms_iterator_from_chain const& other) const{
    return res_it==other.res_it && ( res_it==res_end || other.atm_it==atm_it ); 
  }
  typename do_add_const_to_ref<Atom,is_const>::type dereference() const {return *atm_it;}
};




//iterator over all atoms of a model
template <class Model,bool is_const>
class Atoms_iterator_from_model:
public boost::iterator_facade<
    Atoms_iterator_from_model<Model,is_const>,
    typename do_add_const<typename Model::Chain::Residue::Atom,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Model::Chain::Residue::Atom                       Atom;

  Chains_iterator_from_model<Model,is_const>                   chain_it;
  Chains_iterator_from_model<Model,is_const>                   chain_end;
  Atoms_iterator_from_chain<typename Model::Chain,is_const>  atm_it;
  Atoms_iterator_from_chain<typename Model::Chain,is_const>  atm_end;
public:
  Atoms_iterator_from_model(typename do_add_const_to_ref<Model,is_const>::type model,bool is_end=false):
    chain_it(is_end?model.chains_end():model.chains_begin()),chain_end(model.chains_end())
  {
    if (chain_it!=chain_end){
      atm_it=chain_it->atoms_begin();
      atm_end=chain_it->atoms_end();
      make_valid();
    }
  }
  
  bool is_end(){
    return chain_it==chain_end;
  }
  
private:
  friend class boost::iterator_core_access;

  void increment(){
    ++atm_it;
    make_valid();
  }
  
  void make_valid(){
    while (atm_it==atm_end){
      ++chain_it;
      if (chain_it==chain_end) break;
      atm_it=chain_it->atoms_begin();
      atm_end=chain_it->atoms_end();
    }    
  }
  
  bool equal(Atoms_iterator_from_model const& other) const{
    return chain_it==other.chain_it && ( chain_it==chain_end || other.atm_it==atm_it ); 
  }
  typename do_add_const_to_ref<Atom,is_const>::type dereference() const {return *atm_it;}
};

//iterator over all model of a system
//-----------------------------------
template <class System,bool is_const>
class Model_iterator_from_system:
public boost::iterator_facade<
    Model_iterator_from_system<System,is_const>,
    typename do_add_const<typename System::Model,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename System::Model Model;
public:
  Model_iterator_from_system(typename do_get_const_iterator<typename System::Model_container,is_const>::type it):
    model_it(it){}
  
private:
  friend class boost::iterator_core_access;

  void increment(){++model_it;}
  bool equal(Model_iterator_from_system const& other) const{ return other.model_it==model_it; }
  typename do_add_const_to_ref<Model,is_const>::type dereference() const {return System::dereference(model_it);}
  
//data member
  typename do_get_const_iterator<typename System::Model_container,is_const>::type model_it;
};


//iterator over all coarse_atom of a model
//-----------------------------------
template <class Model,bool is_const>
class Coarse_atom_from_model:
public boost::iterator_facade<
    Coarse_atom_from_model<Model,is_const>,
    typename do_add_const<typename Model::Chain::Residue::Coarse_atom,is_const>::type,
    boost::forward_traversal_tag>
{
private:
  typedef typename Model::Chain::Residue::Coarse_atom Coarse_atom;
public:
  Coarse_atom_from_model(typename do_add_const_to_ref<Model,is_const>::type model,bool is_end=false):
    res_it(is_end?model.residues_end():model.residues_begin()),res_end(model.residues_end())
    {
      if (res_it!=res_end){
        atm_it=res_it->coarse_atoms_begin();
        atm_end=res_it->coarse_atoms_end();
        make_valid();
      }
    }
  
private:
  friend class boost::iterator_core_access;

  void increment(){
    ++atm_it;
    make_valid();
  }
  
  void make_valid(){
    while(atm_it==atm_end){
      ++res_it;
      if (res_it==res_end) break;
      atm_it=res_it->coarse_atoms_begin();
      atm_end=res_it->coarse_atoms_end();
    }    
  }
  bool equal(Coarse_atom_from_model const& other) const{ return other.res_it==res_it && (res_it==res_end || other.atm_it==atm_it); }
  typename do_add_const_to_ref<Coarse_atom,is_const>::type dereference() const {return Model::Chain::Residue::dereference(atm_it);}
  
//data member
  Residues_iterator_from_model<Model,is_const> res_it;
  Residues_iterator_from_model<Model,is_const> res_end;
 typename do_get_const_iterator <typename Model::Chain::Residue::Coarse_atom_container,is_const>::type  atm_it;
 typename do_get_const_iterator <typename Model::Chain::Residue::Coarse_atom_container,is_const>::type atm_end;
};

  
}//namespace internal
/** \endcond */
}//namespace ESBTL

#endif //ESBTL_ITERATORS_H
