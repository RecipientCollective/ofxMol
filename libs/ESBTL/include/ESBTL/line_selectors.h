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



#ifndef ESBTL_LINE_SELECTORS_H
#define ESBTL_LINE_SELECTORS_H

#include <ESBTL/global_functions.h>
#include <ESBTL/PDB.h>
#include <map>

namespace ESBTL{

  /** \defgroup linesel Line selectors
    * This module gathers all objects provided to select lines of a molecular data file to create systems.
    * Such an object must be a class of the following type:
    * \code
    * struct Line_selector_concept{
    *   //returns the index (1...n) of the system in which that line matches: 
    *   // RMK stands for remark (including MODEL), while DISCARD indicate that the line should be discarded.
    *   template <class Line_format,class Occupancy_handler>
    *   int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy)
    *   { ... }
    * 
    *   unsigned max_nb_systems() const {return 1;}
    * };
    * \endcode
    *  In function \b keep the template parameters are:
    *  - \c Line_format is a helper class to read a molecular data file (such as ESBTL::PDB::Line_format).
    *  - \c Occupancy_handler is a class following the concept of \ref occpol. To benefit from
    *  the occupancy policy, you should use the function Occupancy_handler::add_or_postpone(line_format,line,index_of_system)
    *  to return the index of the system the line must inserted in.
    */
  

/** This class is a simple line selector: all atoms and hetero-atoms are in the one system.
  *\ingroup linesel
  */
class PDB_line_selector{
  public:
  int discarded;
  
  PDB_line_selector():discarded(0){}
  
  //return the system pattern that matches the line: 0 remark, 1...n coordinates, -1 to discard
  template <class Line_format,class Occupancy_handler>
  int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy){
    if (line_format.record_type()==PDB::ATOM || line_format.record_type()==PDB::HETATM)
      return occupancy.add_or_postpone(line_format,line,1);
    
    if (line_format.record_type()==PDB::MODEL)
      return RMK;
    ++discarded;
    return DISCARD;
  }
  
  unsigned max_nb_systems() const {return 1;}
  
};

/** This class is a line selector defining two systems:
  * - the first one contains all heavy atoms that do not belong to a water molecule.
  * - the second one contains all containing heavy atoms of water molecules.
  *\ingroup linesel
  */
class PDB_line_selector_two_systems{
  public:
  
  //return the system pattern that matches the line: 0 remark, 1...n coordinates, -1 to discard
  template <class Line_format,class Occupancy_handler>
  int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy){
    
    if (line_format.record_type()==PDB::ATOM || line_format.record_type()==PDB::HETATM){
      if (is_hydrogen(std::make_pair(line_format,line)) )  return DISCARD;
      if ( is_water(std::make_pair(line_format,line)) ) return occupancy.add_or_postpone(line_format,line,2);
      return occupancy.add_or_postpone(line_format,line,1);
    }

    if (line_format.record_type()==PDB::MODEL)
      return RMK;
    return DISCARD;
  }
  
  unsigned max_nb_systems() const {return 2;}
  
};


/** This class is a line selector defining one system per group of chains provided.
  * \ingroup linesel
  */
class PDB_line_selector_chain{
  typedef std::map<char,int> System_indices;
  System_indices system_index_map_;
  bool keep_hydrogen_;
  bool keep_water_;
  bool keep_remaining_chains_;
  double Water_Bfactorlim_;
  int max_chain_index;
  
  unsigned nb_atm_htm_seen_;
  
  public:
  int discarded;

  const unsigned& nb_atm_htm_seen() const {
    return nb_atm_htm_seen_;
  }
  
  //TODO enable_if only for string
  //all characters of one string belong to the same system
  /**
    * Constructor from a range of string. 
    * A string \e AB indicates that the system is composed of chains \e A and \e B
    * \tparam Str_iterator is an iterator over std::string objects.
    * \param begin is the first iterator of the range of strings.
    * \param end is the past-end iterator of the range of strings.
    * \param keep_water indicate whether water molecules should be kept.
    * \param keep_remaining_chains indicates whether the last system should contains unselected chains.
    * \param Water_Bfactorlim water molecules with a temperature factor higher that this value are discarded.
    * \param keep_hydrogen indicates whether hydrogens atoms must be selected.
    */
  template<class Str_iterator>
  PDB_line_selector_chain
  (const Str_iterator& begin,const Str_iterator& end,
   bool keep_water=true,bool keep_remaining_chains=false,
   double Water_Bfactorlim=200,bool keep_hydrogen=false):
    keep_hydrogen_(keep_hydrogen),keep_water_(keep_water),keep_remaining_chains_(keep_remaining_chains),
    Water_Bfactorlim_(Water_Bfactorlim),max_chain_index(0),nb_atm_htm_seen_(0),discarded(0)
  {
    for (Str_iterator it=begin;it!=end;++it){
      const std::string& str=*it;
      ++max_chain_index;
      for (int i=0;i!= (int) str.size();++i)
        system_index_map_.insert(std::make_pair(str[i],max_chain_index));
    }
  }
  
  
  //return the system pattern that matches the line: 0 remark, 1...n coordinates, -1 to discard
  template <class Line_format,class Occupancy_handler>
  int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy){
    if (line_format.record_type()==PDB::ATOM || line_format.record_type()==PDB::HETATM){
      
      System_indices::iterator it=system_index_map_.find(line_format.get_chain_identifier(line));
      ++nb_atm_htm_seen_;
      
      //Filter hydrogen
      if (!keep_hydrogen_ && is_hydrogen(std::make_pair(line_format,line)) )
        return DISCARD;
      

      //Filter water
      if ( is_water(std::make_pair(line_format,line)) ){
        //Filter Temperature factor too high
        if ( !keep_water_ || line_format.get_temperature_factor(line) > Water_Bfactorlim_)
          return DISCARD;
        else
          return occupancy.add_or_postpone(line_format,line,max_chain_index+1);
      }

      if ( it==system_index_map_.end() ){
        if (keep_remaining_chains_)
          return occupancy.add_or_postpone(line_format,line,max_chain_index+(keep_water_?2:1));
        return DISCARD;
      }
      return occupancy.add_or_postpone(line_format,line,it->second);
    }

    if (line_format.record_type()==PDB::MODEL)
      return RMK;
    ++discarded;
    return DISCARD;
  }
  
  unsigned max_nb_systems() const {return max_chain_index+(keep_water_?1:0)+(keep_remaining_chains_?1:0); }
};


//generic line selector
//---------------------------------------
/** 
  *\ingroup linesel
  */

#ifndef __GXX_EXPERIMENTAL_CXX0X__

/** \cond */
#define GENERIC_LINE_SELECTOR_COMMON_PART \
  public: \
  \
  template <class Line_format,class Occupancy_handler> \
  int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy) const{ \
  \
    if (line_format.record_type()==PDB::ATOM || line_format.record_type()==PDB::HETATM){ \
      int ret=operator()(std::make_pair(line_format,line)); \
      if (ret!=-1) \
        return occupancy.add_or_postpone(line_format,line,ret);\
      return DISCARD;  \
    }\
  \
    if (line_format.record_type()==PDB::MODEL)\
      return RMK;\
    return DISCARD;\
  }\
  \
  static unsigned max_nb_systems() { return nb_system;}

template <class T1,class T2=void,class T3=void,class T4=void,class T5=void,class T6=void,class T7=void, class T8=void,class T9=void,class T10=void>
class Generic_line_selector;


template <class T1>
class Generic_line_selector<T1>
{
  T1 t1_;
  static const unsigned nb_system=1;
public:
  Generic_line_selector(){};  
  Generic_line_selector(const T1& t1):t1_(t1){}
protected:    
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    return -1;
  }
  
  GENERIC_LINE_SELECTOR_COMMON_PART
};

template <class T1,class T2>
class Generic_line_selector<T1,T2>
{
  T1 t1_;
  T2 t2_;
public:
  Generic_line_selector(){};  
  Generic_line_selector(const T1& t1,const T2& t2):t1_(t1),t2_(t2){}
protected:    
  static const unsigned nb_system=2;  
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    if (t2_(p))    
      return 2;
    return -1;
  }
    
  GENERIC_LINE_SELECTOR_COMMON_PART
};


template <class T1,class T2, class T3>
class Generic_line_selector<T1,T2,T3>: public Generic_line_selector<T2,T3>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3):t1_(t1),Generic_line_selector<T2,T3>(t2,t3){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3>::nb_system + 1;  
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

template <class T1,class T2, class T3,class T4>
class Generic_line_selector<T1,T2,T3,T4>: public Generic_line_selector<T2,T3,T4>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4):
    t1_(t1),Generic_line_selector<T2,T3,T4>(t2,t3,t4){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

template <class T1,class T2, class T3,class T4,class T5>
class Generic_line_selector<T1,T2,T3,T4,T5>: public Generic_line_selector<T2,T3,T4,T5>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5>(t2,t3,t4,t5){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

template <class T1,class T2, class T3,class T4,class T5, class T6>
class Generic_line_selector<T1,T2,T3,T4,T5,T6>: public Generic_line_selector<T2,T3,T4,T5,T6>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5,T6>(t2,t3,t4,t5,t6){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5,T6>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5,T6>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};


template <class T1,class T2, class T3,class T4,class T5, class T6,class T7>
class Generic_line_selector<T1,T2,T3,T4,T5,T6,T7>: public Generic_line_selector<T2,T3,T4,T5,T6,T7>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5,T6,T7>(t2,t3,t4,t5,t6,t7){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5,T6,T7>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5,T6,T7>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

template <class T1,class T2, class T3,class T4,class T5, class T6,class T7,class T8>
class Generic_line_selector<T1,T2,T3,T4,T5,T6,T7,T8>: public Generic_line_selector<T2,T3,T4,T5,T6,T7,T8>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5,T6,T7,T8>(t2,t3,t4,t5,t6,t7,t8){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

template <class T1,class T2, class T3,class T4,class T5, class T6,class T7,class T8, class T9>
class Generic_line_selector<T1,T2,T3,T4,T5,T6,T7,T8,T9>: public Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9>(t2,t3,t4,t5,t6,t7,t8,t9){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};

/** \endcond */

/**
  * This class is a line selector using atom selectors to defines the systems. Each template parameter
  * is a function object following the concept of \ref atomsel and define one system per template parameter.
  * This class can be provided one to ten template parameters. 
  * If the code is compiled using the c++0x standard, the number of template parameters is not limited (use -std=c++0x with gcc). 
  * \tparam Ti is a function object following the concept \ref atomsel. 
  */
template <class T1,class T2, class T3,class T4,class T5, class T6,class T7,class T8, class T9,class T10>
class Generic_line_selector: public Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9,T10>
{
  T1 t1_;
public:
  Generic_line_selector(){};
  Generic_line_selector(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9,const T10& t10):
    t1_(t1),Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9,T10>(t2,t3,t4,t5,t6,t7,t8,t9,t10){}
protected:    
  static const unsigned nb_system=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9,T10>::nb_system + 1;
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const 
  {
    if (t1_(p))
      return 1;
    int ret=Generic_line_selector<T2,T3,T4,T5,T6,T7,T8,T9,T10>::operator()(p);
    if (ret!=-1)
      return 1+ret;
    return -1;
  }

  GENERIC_LINE_SELECTOR_COMMON_PART  
};


#else
/** \cond */
namespace internal{
template <class ... T>
struct Generic_line_selector_base;

template <class T0,class ... T>
struct Generic_line_selector_base<T0,T...>: public Generic_line_selector_base<T...>{
  static const unsigned nb_system=Generic_line_selector_base<T...>::nb_system + 1;
  T0 t0_;

  Generic_line_selector_base(){}
  Generic_line_selector_base(const T0& t0,const T& ... ts):Generic_line_selector_base<T...>(ts...),t0_(t0){}
  
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const
  {
    if (t0_(p))
      return 1;
    int ret=Generic_line_selector_base<T...>::operator()(p);
    if (ret!=-1) return ret+1;
    return -1;
  }
};

template <class T0>
struct Generic_line_selector_base<T0>{
  static const unsigned nb_system=1;
  T0 t0_;
  
  Generic_line_selector_base() {}
  Generic_line_selector_base(const T0& t0):t0_(t0) {}
  
  template <class Line_format>
  int operator()(const std::pair<Line_format,std::string>& p) const {
    if (t0_(p))
      return 1;
    return -1;
  }
};

} //namespace internal


template <class ... T>
class Generic_line_selector: public internal::Generic_line_selector_base<T ... >{
  public:
  
  Generic_line_selector(){}
  Generic_line_selector(const T& ... ts):internal::Generic_line_selector_base<T...>(ts...){}
  
  //return the system pattern that matches the line: 0 remark, 1...n coordinates, -1 to discard
  template <class Line_format,class Occupancy_handler>
  int keep(const Line_format& line_format,const std::string& line,Occupancy_handler& occupancy){
    
    if (line_format.record_type()==PDB::ATOM || line_format.record_type()==PDB::HETATM){
      int ret=internal::Generic_line_selector_base<T ... >::operator()(std::make_pair(line_format,line));
      if (ret!=-1)
        return occupancy.add_or_postpone(line_format,line,ret);
      return DISCARD;  
    }

    if (line_format.record_type()==PDB::MODEL)
      return RMK;
    return DISCARD;
  }
  
  static unsigned max_nb_systems() {
    return internal::Generic_line_selector_base<T ... >::nb_system;
  }
  
};
/** \endcond */
#endif












}// namespace ESBTL 

#endif //ESBTL_LINE_SELECTORS_H
