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



#ifndef ESBTL_OCCUPANCY_HANDLERS_H
#define ESBTL_OCCUPANCY_HANDLERS_H

#include <set>
#include <iostream>

namespace ESBTL{

 /** \defgroup occpol Occupancy policies 
   * This module gathers all methods provided to handle atoms with no alternate location
   * identification and occupancy factor different from 1.
   * An occupancy policy must be a class of the following type: 
   * \code
   * template <class Line_format>
   * struct Occupancy_policy_concept{
   *   //check whether an atom can directly be added to a system. 
   *   //It returns -1 if not, and the number of the system otherwise.
   *   //Line_format must provide  similar functionalities to ESBTL::PDB::Line_format.
   *   //line is a line of a molecular data file representing an atom.
   *   //system_index is the index of the system the line should be added (provided by the selector).
   *   int add_or_postpone(const Line_format& line_format,const std::string& line,int system_index) const{ ... }
   *   
   *   //This function is called by the builder when the whole file have been read.
   *   //Potential atom that have not already been inserted can be. The integer returned is
   *   //the number of atoms added.
   *   template <class Builder> int finalize(Builder&) const { ... }
   * };
   * \endcode
   */
  
  /** This object provides no occupancy policy. This results in a crash of the application
    * when an atom with no alternate location and with occupancy!=1 is encountered.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */
  template <class Line_format>
  struct No_occupancy_policy{
    int add_or_postpone(const Line_format& line_format,const std::string& line,int system_index) const{
      double occupancy_value=line_format.get_occupancy(line);
      if (occupancy_value==1. || line_format.get_alternate_location(line)!=' ')
        return system_index;
      std::cerr << "Fatal error: Occupancy factor < 1. No occupancy selection policy is provided to handle line\n";
      std::cerr << "<|" << line << "|>\n";
      exit(EXIT_FAILURE);
    }
    
    template <class Builder>
    int finalize(Builder&) const {return 0;}
  };


  /** This object provides a no-restriction occupancy policy: All atoms are accepted.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */  
  template <class Line_format>
  struct Accept_all_occupancy_policy{
    int add_or_postpone(const Line_format& ,const std::string& ,int system_index) const {
      return system_index;
    }
    
    template <class Builder>
    int finalize(Builder&) const {return 0;}
  };

  /** This object provides a restrictive occupancy policy: All atoms are discarded.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */  
  template <class Line_format>
  struct Accept_none_occupancy_policy{
    int add_or_postpone(const Line_format& line_format,const std::string& line,int system_index) const {
      double occupancy_value=line_format.get_occupancy(line);
      if (occupancy_value==1. || line_format.get_alternate_location(line)!=' ')
        return system_index;      
      return DISCARD;
    }
    
    template <class Builder>
    int finalize(Builder&) const {return 0;}
  };

  /** \cond */
  namespace internal{
    
    template <class Line_format,bool take_max>
    class Min_or_max_occupancy_policy{
    protected:
      typedef std::map<double,std::pair<int,std::list<std::string> > > Info;
      Info lines;
      double min_or_max_occupancy;
      
      Min_or_max_occupancy_policy():min_or_max_occupancy(take_max?0.:1.){}
      
      void add_line(double occupancy_value,const std::string& line,int system_index){
        typename Info::iterator it=lines.find(occupancy_value);
        if (it==lines.end())
          it=lines.insert(
              std::make_pair(occupancy_value,
                             std::make_pair(system_index,std::list<std::string>())) 
          ).first;
        it->second.second.push_back(line);
      }
      
    public:
      int add_or_postpone(const Line_format& line_format,const std::string& line,int system_index){
        double occupancy_value=line_format.get_occupancy(line);
        if (occupancy_value==1. || line_format.get_alternate_location(line)!=' ')
          return system_index;
        
        bool do_take_it=take_max?
                          occupancy_value >= min_or_max_occupancy:
                          occupancy_value <= min_or_max_occupancy;
          
        if (do_take_it){
          add_line(occupancy_value,line,system_index);
          min_or_max_occupancy=occupancy_value;
        }
        
        return DISCARD;
      }
      
      template <class Builder>
      int finalize(Builder& builder) const {
        if (lines.empty())
          return 0;
        if (min_or_max_occupancy==0.5){
          std::cerr << "Fatal error: The atom occupancy selection policy is ambiguous is that case (0.5)\n";
          exit(EXIT_FAILURE);
        }
        int lines_added=0;
        Info::const_iterator it=lines.find(min_or_max_occupancy);
        for (std::list<std::string>::const_iterator cur=it->second.second.begin();
             cur!=it->second.second.end(); ++cur)
        {
          ++lines_added;
          builder.interpret_line(Line_format(*cur),*cur,it->second.first);
        }
        
        return lines_added;
      }
    };
    
  }//namespace internal
  /** \endcond */

  /** This object provides an occupancy policy. An atom with an occupancy != 1 and with no alternate location
    * is kept if its occupancy value corresponds to the highest found while reading a molecular structure.
    * Note that if a line have been discarded by the selector, then its occupancy value have not been taken
    * into account for computing the maximum. 
    * In case of ambiguous maximum (0.5), the program will crash with an error message.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */    
  template<class Line_format>
  class Max_occupancy_policy: public internal::Min_or_max_occupancy_policy<Line_format,true>{};


  /** This object provides an occupancy policy. An atom with an occupancy != 1 and with no alternate location
    * is kept if its occupancy value corresponds to the lowest found while reading a molecular structure.
    * Note that if a line have been discarded by the selector, then its occupancy value have not been taken
    * into account for computing the minimum.
    * In case of ambiguous minimum (0.5), the program will crash with an error message.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */
  template<class Line_format>
  class Min_occupancy_policy: public internal::Min_or_max_occupancy_policy<Line_format,false>{};
  
  /** Keep only a set of atoms given by their atom serial number
    * \ingroup occpol
    */

  /** This object provides an occupancy policy to keep atoms with an occupancy != 1 and with no alternate location
    * using their atom serial number.
    * @tparam Line_format is a helper class to read molecular data files such as ESBTL::PDB::Line_format.
    * \ingroup occpol
    */    
  template <class Line_format>
  class Atom_list_occupancy_policy{
    std::set<int> selected_atoms;
  public:  
    /** Constructor. The serial numbers of atoms to be kept are given by
      * an iterator range.
      * @tparam Iterator is an iterator over unsigned integers.
      * @param begin is the first iterator over the range of atom serial numbers.
      * @param end is the past-end iterator over the range of atom serial numbers.
      */
    template <class Iterator>
    Atom_list_occupancy_policy(Iterator begin,Iterator end){
      selected_atoms.insert(begin,end);
    }
  
    int add_or_postpone(const Line_format& line_format,const std::string& line,int system_index) const {
      double occupancy_value=line_format.get_occupancy(line);
      if (occupancy_value==1. || line_format.get_alternate_location(line)!=' ')
        return system_index;
      
      if (selected_atoms.find(line_format.get_atom_serial_number(line))!=selected_atoms.end())
        return system_index;
      
      return DISCARD;
    }
    
    template <class Builder>
    int finalize(Builder&) const {
      return 0;
    }
  };
} //namespace ESBTL


#endif //ESBTL_OCCUPANCY_HANDLERS_H
