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



#ifndef ESBTL_COARSE_CREATORS_H
#define ESBTL_COARSE_CREATORS_H

#include <limits>

namespace ESBTL{

/**
  * \defgroup coarse_creator Coarse grain creator
  * This module gathers all objects provided to create coarse grain atoms from a residue.
  * Such an object must be a class of the following type:
  * \code
  *  template <class Residue>
  *  struct Coarse_creator_concept{
  *    //a string indicating how the coarse grained model is defined.
  *    static const std::string info(){ ... }
  *
  *    // Creates coarse grain atoms from a residue.
  *    // Output_iterator is an output iterator of coarse grain atoms.
  *    // the integer returned is the number of coarse grain atoms created.
  *    template <class Output_iterator>
  *    int operator()(const Residue& res,Output_iterator out) const{ ... }
  *  };
  * \endcode
  */  


/** 
  * Creator of a coarse grain model creating up to two pseudo-atoms: one as the barycenter of backbone atoms 
  * (if relevant) and one as the barycenter of non-barycenter atoms. Backbone atoms
  * are identified using the global function ESBTL::is_backbone.
  * \ingroup coarse_creator
  * \tparam Residue is the coarse grain residue type used.
  */
template <class Residue>
class Coarse_creator_two_barycenters{
  typedef typename Residue::Coarse_atom         Coarse_atom;
  typedef typename Residue::Atom                Atom;
  typedef double                                FT;
  typedef typename Residue::Atom::Point_3       Point_3;
  
  static const std::string  info_;

  void update_coordinates(FT* X,const Atom& atom) const{
    X[0]+=atom.x();
    X[1]+=atom.y();
    X[2]+=atom.z();
    X[3]+=1;
  }
  
public:

  static const std::string info(){return info_;}

  //returns the number of atoms created
  template <class Output_iterator>
  int operator()(const Residue& res,Output_iterator out) const{
    Coarse_atom bb(0,res);
    Coarse_atom sc(1,res);
    FT Xbb[4]={0,0,0,0};
    FT Xsc[4]={0,0,0,0};
    unsigned tot_res=0;
    for(typename Residue::Atoms_const_iterator
        at_it=res.atoms_begin();
        at_it!=res.atoms_end();
        ++at_it)
    {
      ++tot_res;
      if ( is_backbone(*at_it) ){
        bb.add(*at_it);
        update_coordinates(Xbb,*at_it);
      }
      else{
        sc.add(*at_it);
        update_coordinates(Xsc,*at_it);
      }
    }
    if (tot_res==0)
      return 0;
    unsigned ret=0;
    if (Xbb[3]!=0){ //to handle case of a ligand for example.
      static_cast<Point_3&>(bb)=Point_3(Xbb[0]/Xbb[3],Xbb[1]/Xbb[3],Xbb[2]/Xbb[3]);
      *out++=bb;
      ++ret;
    }
    if (Xsc[3]!=0){
      static_cast<Point_3&>(sc)=Point_3(Xsc[0]/Xsc[3],Xsc[1]/Xsc[3],Xsc[2]/Xsc[3]);
      *out++=sc;
      return ret+1;
    }
    else{assert(res.residue_name()=="GLY");}
    assert(ret==1);
    return 1;
  }
};

template<class Residue>
const std::string Coarse_creator_two_barycenters<Residue>::info_="Coarse model using one atom as barycenter of backbone and one atom as barycenter of side chain atoms";


/** 
  * Creator of a coarse grain model using one pseudo-atom centered at
  * the atom closest to the barycenter of the side-chain and the C-alpha of a residue.
  * These atoms are identified using the global function ESBTL::is_side_chain_or_CA.
  * \ingroup coarse_creator
  * \tparam Residue is the coarse grain residue type used.
  * \tparam FT_ is the number type of the coordinates.
  */
template <class Residue,class FT_>
class Coarse_creator_closest_to_barycenter{
  typedef typename Residue::Coarse_atom         Coarse_atom;
  typedef typename Residue::Atom                Atom;
  typedef FT_                                   FT;
  
  static const std::string  info_;

  void update_coordinates(FT* X,const Atom& atom) const{
    X[0]+=atom.x();
    X[1]+=atom.y();
    X[2]+=atom.z();
    X[3]+=1;
  }
  
  static inline FT square(const FT& n){return n*n;}
  
  template <class Point_3>
  static inline FT squared_distance(const Point_3& p1,const Point_3& p2){
    return square(p1.x()-p2.x())+square(p1.y()-p2.y())+square(p1.z()-p2.z());
  }
  
public:

  static const std::string info(){return info_;}

  //returns the number of atoms created
  template <class Output_iterator>
  int operator()(const Residue& res,Output_iterator out) const {
    Coarse_atom c_atm(0,res);
    FT X_c_atm[4]={0,0,0,0};
    std::list<const Atom*> candidate_list;
    for(typename Residue::Atoms_const_iterator
        at_it=res.atoms_begin();
        at_it!=res.atoms_end();
        ++at_it)
    {
      if ( is_side_chain_or_CA(*at_it) ){
        c_atm.add(*at_it);
        update_coordinates(X_c_atm,*at_it);
        candidate_list.push_back( &(*at_it) );
      }
    }
    assert(X_c_atm[3]!=0);
    Atom barycenter=Atom(X_c_atm[0]/X_c_atm[3],X_c_atm[1]/X_c_atm[3],X_c_atm[2]/X_c_atm[3]);
    FT min_dist=std::numeric_limits<FT>::max();
    
    typename std::list<const Atom*>::iterator min_it=candidate_list.end();
    for (typename std::list<const Atom*>::iterator it=candidate_list.begin();it!=candidate_list.end();++it){
      FT curr_dist=squared_distance(barycenter,*(*it));
      if (curr_dist < min_dist){
        min_it=it;
        min_dist=curr_dist;
      }
    }
    assert(min_it!=candidate_list.end());
    static_cast<typename Coarse_atom::Point_3&>(c_atm)=static_cast<const typename Atom::Point_3&>(*(*min_it));

        
    *out++=c_atm;
    return 1;
  }
};

template<class Residue,class FT_>
const std::string 
Coarse_creator_closest_to_barycenter<Residue,FT_>::info_=
"Coarse model using one pseudo-atom centered at the atom closest to the barycenter of the side-chain and the C-alpha of each residue";





} //namespace ESBTL

#endif //ESBTL_COARSE_CREATORS_H
