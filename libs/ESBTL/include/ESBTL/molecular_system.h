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



#ifndef ESBTL_MOLECULAR_SYSTEM_H
#define ESBTL_MOLECULAR_SYSTEM_H

#include <map>
#include <list>
#include <vector>
#include <cassert>
#include <string>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <ESBTL/iterators.h>
#include <ESBTL/constants.h>

/** \defgroup grp_iters Iterators 
  * Iterators and functions offering iteration possibilities are gathered on this page.
  */

namespace ESBTL{

/** 
  * A class representing a molecular system.
  * \tparam Items is a class gathering wrapper defining model, chain , residue and atom types (see Default_system_items or System_items_with_coarse_grain for example).
  * \tparam Point is a point type with coordinates const access methods x(), y() and z().
  */
template <class Items,class Point>
class Molecular_system{
  typedef Molecular_system<Items,Point>                                                 Self;
public:
  typedef typename Items::template   Model_wrapper<Molecular_system,Point>::Type        Model;
  typedef typename Items::template   Chain_wrapper<Molecular_system,Point>::Type        Chain;
  typedef typename Items::template    Atom_wrapper<Molecular_system,Point>::Type         Atom;
  typedef typename Items::template Residue_wrapper<Molecular_system,Point>::Type      Residue;
  typedef std::map<int,Model>                                                   Model_container;   

  //iterators
  //---------
  static const Model& dereference (typename Model_container::const_iterator it) {return it->second;}
  static Model& dereference (typename Model_container::iterator it) {return it->second;}

  /** \ingroup grp_iters*/
  typedef internal::Model_iterator_from_system<Self,false>                        Models_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Model_iterator_from_system<Self,true>                   Models_const_iterator;

  /** \ingroup grp_iters*/
  Models_iterator models_begin() {return Models_iterator(model_container_.begin()); }
  /** \ingroup grp_iters*/
  Models_iterator models_end()   {return Models_iterator(model_container_.end()); }
  /** \ingroup grp_iters*/
  Models_const_iterator models_begin() const {return Models_const_iterator(model_container_.begin());}
  /** \ingroup grp_iters*/
  Models_const_iterator models_end()   const {return Models_const_iterator(model_container_.end());  }
  //---------
  
  Molecular_system(int index,std::string name="no_name"):name_(name),index_(index),alternate_location_(' '){}

  template<class Line_format>
  void interpret_line(const Line_format& line_format,const std::string& line,int current_model){
    Model& model=get_or_create_model(current_model);
    Chain& chain=model.get_or_create_chain(line_format,line);
    Residue& residue=chain.get_or_create_residue(line_format,line);
    residue.add_atom(line_format,line);
  }

  bool has_no_model() const{
    return model_container_.empty();
  }

  bool has_model(int i) const{
    return model_container_.find(i)!=model_container_.end();
  }  
  
  size_t number_of_models() const {
    return model_container_.size();
  }
  
  void set_altloc(char c) {alternate_location_=c;}
  
  //member variables
  //DECLARE_AND_ACCESS(model_container,Model_container)
  DECLARE_AND_ACCESS(name,std::string)
  DECLARE_AND_ACCESS(index,int)
  DECLARE_AND_ACCESS(alternate_location,char)
  
  //create the model if it does not exist
  Model& get_or_create_model(int i){
    typename Model_container::iterator itm=model_container_.find(i);
    if (itm==model_container_.end()){
      itm=model_container_.insert( std::make_pair(i,Model(i,*this)) ).first;
    }
    return itm->second;
  }
  
  Model& get_model(int i){
    typename Model_container::iterator itm=model_container_.find(i);
    if (itm==model_container_.end()){std::cerr <<"Cannot find model " << i << std::endl; exit(EXIT_FAILURE);}
    return itm->second;    
  }

  //Model_container& model_container(){return model_container_;}
private:
  Model_container model_container_;
};

/** A class representing a model.
  * \tparam System_ is a system (like ESBTL::Molecular_system for example).
  */
template<class System_>
class Molecular_model{
 
public:
  typedef System_                          System;
  typedef typename System::Chain           Chain;
  typedef std::map<char,Chain>             Chain_container;
private:
  typedef Molecular_model<System_>         Self;
  const System& system_;
  Chain_container chain_container_;
public:
  const System& system() const {return system_;}
  
  Molecular_model(int nbm,const System& sys):system_(sys),model_number_(nbm){}

  template <class Line_format>
  Chain& get_or_create_chain(const Line_format& line_format,const std::string& line){
    char ch=line_format.get_chain_identifier(line);
    typename Chain_container::iterator itch=chain_container_.find(ch);
    if (itch==chain_container_.end())
      itch=chain_container_.insert(std::make_pair(ch,Chain(line_format,line,*this))).first;
    return itch->second;
  }
  
  Chain& get_or_create_chain(char id){
    typename Chain_container::iterator itch=chain_container_.find(id);
    if (itch==chain_container_.end())
      itch=chain_container_.insert(std::make_pair(id,Chain(id,*this))).first;
    return itch->second;    
  }

  const Chain& get_chain(char id) const {
    typename Chain_container::iterator itch=chain_container_.find(id);
    assert (itch!=chain_container_.end());
    return itch->second;    
  }

  const typename System::Residue& get_residue(char ch_id,int ressn,char insc=' ') const {
    return get_chain(ch_id).get_residue(ressn,insc);
  }
  
  const typename System::Atom& get_atom(char ch_id,int ressn,char insc,unsigned atom_sn) const {
    return get_chain(ch_id).get_residue(ressn,insc).get_atom(atom_sn);
  }
  
  size_t number_of_chains() const {
    return chain_container_.size();
  }
  
  size_t number_of_residues() const {
    size_t total=0;
    for (typename Chain_container::const_iterator it=chain_container_.begin();it!=chain_container_.end();++it)
      total+=it->second.number_of_residues();
    return total;
  }

  size_t number_of_atoms() const {
    size_t total=0;
    for (typename Chain_container::const_iterator it=chain_container_.begin();it!=chain_container_.end();++it)
      total+=it->second.number_of_atoms();
    return total;
  }  
  //iterators
  //---------
  //functions needed by generic iterator in iterators.h 
  static const Chain& dereference (typename Chain_container::const_iterator it) {return it->second;}
  static Chain& dereference (typename Chain_container::iterator it) {return it->second;}
  
  /** \ingroup grp_iters*/
  typedef internal::Chains_iterator_from_model<Self,true>  Chains_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Chains_iterator_from_model<Self,false> Chains_iterator;
  
  /** \ingroup grp_iters*/
  Chains_iterator chains_begin() {return Chains_iterator(chain_container_.begin());}
  /** \ingroup grp_iters*/
  Chains_iterator chains_end()   {return Chains_iterator(chain_container_.end());}

  /** \ingroup grp_iters*/
  Chains_const_iterator chains_begin() const {return Chains_const_iterator(chain_container_.begin());}
  /** \ingroup grp_iters*/
  Chains_const_iterator chains_end()   const {return Chains_const_iterator(chain_container_.end());}  
  
  /** \ingroup grp_iters*/
  typedef internal::Residues_iterator_from_model<Self,true>  Residues_const_iterator;
  /** \ingroup grp_iters*/  
  typedef internal::Residues_iterator_from_model<Self,false> Residues_iterator;

  /** \ingroup grp_iters*/  
  Residues_iterator residues_begin() {return Residues_iterator(*this);}
  /** \ingroup grp_iters*/
  Residues_iterator residues_end()   {return Residues_iterator(*this,true);}

  /** \ingroup grp_iters*/
  Residues_const_iterator residues_begin() const {return Residues_const_iterator(*this);}
  /** \ingroup grp_iters*/
  Residues_const_iterator residues_end()   const {return Residues_const_iterator(*this,true);}  
  
  
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_model<Self,true>  Atoms_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_model<Self,false> Atoms_iterator;
  
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_begin() {return Atoms_iterator(*this);}
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_end()   {return Atoms_iterator(*this,true);}

  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_begin() const {return Atoms_const_iterator(*this);}
  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_end()   const {return Atoms_const_iterator(*this,true);}
  //---------
  
  
  //DECLARE_AND_ACCESS(chain_container,Chain_container)
  DECLARE_AND_ACCESS(model_number,int)

};


/** A class representing a chain.
  * \tparam System is a system (like ESBTL::Molecular_system for example).
  */
template<class System>
class Molecular_chain{
public:
  typedef typename System::Residue                  Residue;
  typedef typename Residue::Atom                    Atom;
  typedef std::map<std::pair<int,char>,Residue>     Residue_container;     
  typedef typename System::Model                    Model;
private:
  typedef Molecular_chain<System>                   Self;
  const Model& model_;
  Residue_container residue_container_;
public:
  
  template<class Line_format>  
  Molecular_chain(const Line_format& line_format,const std::string& line, const Model& mod):
    model_(mod),
    chain_identifier_(line_format.get_chain_identifier(line)) {}

  Molecular_chain(char id, const Model& mod):model_(mod),chain_identifier_(id) {}  
      
  const Model& model() const {return model_;}

  template<class Line_format>
  Residue& get_or_create_residue(const Line_format& line_format,const std::string& line){
    int ressn=line_format.get_residue_sequence_number(line);
    char insc=line_format.get_insertion_code(line);
    typename Residue_container::iterator it=residue_container_.find(std::make_pair(ressn,insc));
    if (it==residue_container_.end())
      it=residue_container_.insert( std::make_pair(std::make_pair(ressn,insc),Residue(line_format,line,*this) ) ).first;
    return it->second;
  }
  
  Residue& get_or_create_residue(const std::string& resname,int ressn,char insc){
    typename Residue_container::iterator it=residue_container_.find(std::make_pair(ressn,insc));
    if (it==residue_container_.end())
      it=residue_container_.insert( std::make_pair(std::make_pair(ressn,insc),Residue(resname,ressn,insc,*this) ) ).first;
    return it->second;
  }
  
  const Residue& get_residue(int ressn,char insc=' '){
    typename Residue_container::iterator it=residue_container_.find(std::make_pair(ressn,insc));
    assert (it!=residue_container_.end());
    return it->second;    
  }
  
  const typename System::Atom get_atom(int ressn,char insc,unsigned atom_sn) const {
    return get_residue(ressn,insc).get_atom(atom_sn);
  }
  
  size_t number_of_residues() const {
    return residue_container_.size();
  }

  size_t number_of_atoms() const {
    size_t total=0;
    for (typename Residue_container::const_iterator it=residue_container_.begin();it!=residue_container_.end();++it)
      total+=it->second.number_of_atoms();
    return total;
  }  
  
  //iterators
  //---------
  //functions needed by generic iterator in iterators.h 
  static const Residue& dereference (typename Residue_container::const_iterator it) {return it->second;}
  static Residue& dereference (typename Residue_container::iterator it) {return it->second;}
  
  /** \ingroup grp_iters*/
  typedef internal::Residues_iterator_from_chain<Self,true>  Residues_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Residues_iterator_from_chain<Self,false> Residues_iterator;
  
  /** \ingroup grp_iters*/
  Residues_iterator residues_begin() {return Residues_iterator(residue_container_.begin());}
  /** \ingroup grp_iters*/
  Residues_iterator residues_end()   {return Residues_iterator(residue_container_.end());}  

  /** \ingroup grp_iters*/
  Residues_const_iterator residues_begin() const {return Residues_const_iterator(residue_container_.begin());}
  /** \ingroup grp_iters*/
  Residues_const_iterator residues_end()   const {return Residues_const_iterator(residue_container_.end());}
  
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_chain<Self,true>  Atoms_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_chain<Self,false> Atoms_iterator;
  
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_begin() {return Atoms_iterator(*this);}
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_end()   {return Atoms_iterator(*this,true);}  

  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_begin() const {return Atoms_const_iterator(*this);}
  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_end()   const {return Atoms_const_iterator(*this,true);}  
  //-------
  
  DECLARE_AND_ACCESS(chain_identifier,char)
};

/** A class representing a residue.
  * \tparam System is a system (like ESBTL::Molecular_system for example).
  */
template<class System>
class Molecular_residue{
public:
  typedef typename System::Atom          Atom;
  typedef std::map<unsigned,Atom>        Atom_container;   
  typedef typename System::Chain         Chain;
private:
  typedef Molecular_residue<System>      Self;
  const Chain& chain_;
protected:
  Atom_container atom_container_;
public:
  const Chain& chain() const {return chain_;}
  char chain_identifier() const {return chain_.chain_identifier();}
  
  template<class Line_format>
  Molecular_residue(const Line_format& line_format, const std::string& line,const Chain& ch):
    chain_(ch),
    residue_name_(line_format.get_residue_name(line)),
    residue_sequence_number_(line_format.get_residue_sequence_number(line)),
    insertion_code_(line_format.get_insertion_code(line)){}

  Molecular_residue(const std::string& resname,int index,char insc,const Chain& ch):
    chain_(ch),residue_name_(resname),residue_sequence_number_(index),insertion_code_(insc){}
      
  template<class Line_format>    
  void add_atom(const Line_format& line_format, const std::string& line)
  {
    unsigned sn=line_format.get_atom_serial_number(line);
    assert (atom_container_.find(sn)==atom_container_.end() || !"Two atoms with same serial numbers.");
    atom_container_.insert(std::make_pair(sn,Atom(line_format,line,*this)));
  }

  size_t number_of_atoms() const {
    return atom_container_.size();
  }
  
  const Atom& get_atom(unsigned sn) const {
    typename Atom_container::const_iterator  it=atom_container_.find(sn);
    assert(atom_container_.find(sn)!=atom_container_.end());
    return it->second;
  }
  
  //iterators
  //---------
  //functions needed by generic iterator in iterators.h 
  static const Atom& dereference (typename Atom_container::const_iterator it) {return it->second;}
  static Atom& dereference (typename Atom_container::iterator it) {return it->second;}
  
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_residue<Self,true>  Atoms_const_iterator;
  /** \ingroup grp_iters*/
  typedef internal::Atoms_iterator_from_residue<Self,false> Atoms_iterator;
  
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_begin() {return Atoms_iterator(atom_container_.begin(),atom_container_.end());}
  /** \ingroup grp_iters*/
  Atoms_iterator atoms_end()   {return Atoms_iterator(atom_container_.end());}  

  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_begin() const {return Atoms_const_iterator(atom_container_.begin(),atom_container_.end());}
  /** \ingroup grp_iters*/
  Atoms_const_iterator atoms_end()   const {return Atoms_const_iterator(atom_container_.end());}
  //-------  
  
  DECLARE_AND_ACCESS(residue_name,std::string)
  DECLARE_AND_ACCESS(residue_sequence_number,int)
  DECLARE_AND_ACCESS(insertion_code,char)
};


/** A class representing an atom.
  * \tparam System_ is a system (like ESBTL::Molecular_system for example).
  * \tparam Point is a point type with coordinates const access methods x(), y() and z(). 
  */
template<class System_,class Point>
class Molecular_atom: public Point{
public:
  typedef System_                                      System;
  typedef Point                                        Point_3;
  typedef typename System::Residue                     Residue;

private:
  const Residue* residue_;
public:
  
  template <class Line_format,class Residue_type>
  Molecular_atom(const Line_format& line_format, const std::string& line,const Residue_type& res):
    Point(line_format.get_x(line),line_format.get_y(line),line_format.get_z(line)),
    residue_(static_cast<const Residue*>(&res)),
    is_hetatm_(line_format.is_hetatm()),
    atom_serial_number_(line_format.get_atom_serial_number(line)),
    atom_name_(line_format.get_atom_name(line)),
    alternate_location_(line_format.get_alternate_location(line)),
    occupancy_(line_format.get_occupancy(line)),
    temperature_factor_(line_format.get_temperature_factor(line)),
    element_(line_format.get_element(line)),
    charge_(line_format.get_charge(line))
  {
    //IN SOME FILES ... ELEMENT_SYMBOL IS CORRUPTED BY PREVIOUS FIELDS
    //~ if (element_==" ")
      //~ element_=line_format.get_element_using_atom_name(line);
    //~ assert(element_==line_format.get_element_using_atom_name(line));
  }
  
  Molecular_atom():Point(),residue_(NULL){}
  //TODO should use the number type of the Point 
  Molecular_atom(double x,double y,double z):Point(x,y,z),residue_(NULL){}

  //System functions
  int system_index() const {return residue_->chain().model().system().index();};
  
  //Chain functions
  char chain_identifier() const {return residue_->chain().chain_identifier();}
  //Residue functions
  const Residue& residue() const {return *residue_;}
  const std::string& residue_name() const {return residue_->residue_name();}
  int residue_sequence_number() const {return residue_->residue_sequence_number();}
  char insertion_code() const {return residue_->insertion_code();}
  
  //member variables
  DECLARE_AND_ACCESS(is_hetatm,bool)
  DECLARE_AND_ACCESS(atom_serial_number,int)
  DECLARE_AND_ACCESS(atom_name,std::string)
  DECLARE_AND_ACCESS(alternate_location,char)
  DECLARE_AND_ACCESS(occupancy,double)
  DECLARE_AND_ACCESS(temperature_factor,double)
  DECLARE_AND_ACCESS(element,std::string)
  DECLARE_AND_ACCESS(charge,int)
};


/** A default all-atom system items gathering wrappers
  * (needed to break the cyclic dependency) 
  * that allow to define a model, chain, residue and atom type
  * from the types Molecular_model, Molecular_chain, Molecular_residue and Molecular_atom respectively.
  */
struct Default_system_items{
  template <class System,class Point>
  struct Model_wrapper{
    typedef Molecular_model<System>         Type;
  };
  
  template <class System,class Point>
  struct Chain_wrapper{
    typedef Molecular_chain<System>         Type;
  };
  
  template <class System,class Point>
  struct Residue_wrapper{
    typedef Molecular_residue<System>       Type;
  };
  
  template <class System,class Point>
  struct Atom_wrapper{
    typedef Molecular_atom<System,Point>  Type;
  };  
};


//global access functions
//for chains
  template <class System>
  char get_chain_identifier(const Molecular_chain<System>& c) {return c.chain_identifier();}
//for residues
  template <class System>
  const std::string& get_residue_name(const Molecular_residue<System>& r) {return r.residue_name();}
  template <class System>
  int get_residue_sequence_number(const Molecular_residue<System>& r) {return r.residue_sequence_number();}
  template <class System>
  char get_insertion_code(const Molecular_residue<System>& r) {return r.insertion_code();}
//for atoms
  template <class System,class Point>
  bool get_is_hetatm(const Molecular_atom<System,Point>& a) {return a.is_hetatm() ;}
  template <class System,class Point>
  int get_atom_serial_number(const Molecular_atom<System,Point>& a) {return a.atom_serial_number() ;}
  template <class System,class Point>
  const std::string& get_atom_name(const Molecular_atom<System,Point>& a) {return a.atom_name() ;}
  template <class System,class Point>
  char get_alternate_location(const Molecular_atom<System,Point>& a) {return a.alternate_location() ;}
  template <class System,class Point>
  double get_occupancy(const Molecular_atom<System,Point>& a) {return a.occupancy() ;}
  template <class System,class Point>
  double get_temperature_factor(const Molecular_atom<System,Point>& a) {return a.temperature_factor() ;}
  template <class System,class Point>
  const std::string& get_element(const Molecular_atom<System,Point>& a) {return a.element() ;}
  template <class System,class Point>
  int get_charge(const Molecular_atom<System,Point>& a) {return a.charge() ;}
  template <class System,class Point>
  char get_chain_identifier(const Molecular_atom<System,Point>& a) {return a.chain_identifier() ;}
  template <class System,class Point>
  const std::string& get_residue_name(const Molecular_atom<System,Point>& a) {return a.residue_name() ;}
  template <class System,class Point>
  int get_residue_sequence_number(const Molecular_atom<System,Point>& a) {return a.residue_sequence_number() ;}
  template <class System,class Point>
  char get_insertion_code(const Molecular_atom<System,Point>& a) {return a.insertion_code() ;}
  template <class System,class Point>
  char get_x(const Molecular_atom<System,Point>& a) {return a.x() ;}
  template <class System,class Point>
  char get_y(const Molecular_atom<System,Point>& a) {return a.y() ;}
  template <class System,class Point>
  char get_z(const Molecular_atom<System,Point>& a) {return a.z() ;}
  
//~ //function to insert a set water molecule in as a system
//~ template<class Input_iterator,class System>
//~ void insert_atoms(Input_iterator first,Input_iterator last,System& system,int modelid=1,char chainid='Z',std::string resname="HOH",int starting_res_index=1){
  //~ typename System::Model& model=system.get_or_create_model(modelid);
  //~ typename System::Chain& chain=model.get_or_create_chain(chainid);
  //~ int res_index=starting_res_index-1;
  //~ char insc=' ';
  //~ for (Input_iterator it=first; it!=last;++it){
    //~ typename System::Residue& residue=chain.get_or_create_residue(resname,++res_index,insc);
    //~ residue.add_coarse_atom(*it,0);
  //~ }
//~ }
  

}// namespace ESBTL

#endif //ESBTL_MOLECULAR_SYSTEM_H
