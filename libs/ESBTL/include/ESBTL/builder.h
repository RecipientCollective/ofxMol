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



#ifndef ESBTL_BUILDER_H
#define ESBTL_BUILDER_H
#include <iostream>

namespace ESBTL {

/**
 * Class responsible for building a system.
 * @tparam System is a system as ESBTL::Molecular_system.
 */
template<class System>
class All_atom_system_builder{
  //TODO a selector inside the builder
public:
  /**
   * The container in which the system must be stored. 
   */
  typedef std::vector <System> System_container;
private:
  int current_model;//model_number : we expect to read one model after another
  System_container& systems_;
  bool prebuilt_systems;
  unsigned max_systems;
public:
  //I think this constructor is useless
  //All_atom_system_builder(System_container& systems):current_model(1),systems(systems),prebuilt_systems(false),max_systems(0){}
  /**
   * Constructor.
   * @param systems is a reference to the container in which created systems are inserted.
   * @param max_systems_ is an upper bound on the number of system to be created. Setting correctly this number
   *  is very important and the program will crash if not well set (in the current implementation it is used to set the size 
   *  of the system container.)
   */
  All_atom_system_builder(System_container& systems,unsigned max_systems_):
    current_model(1),systems_(systems),prebuilt_systems(true),max_systems(max_systems_){
    if (systems_.capacity()<max_systems)
      systems_.reserve(max_systems);
    //TODO warning cannot find a better method?
    for (unsigned i=0;i<max_systems;++i) systems_.push_back(System(i+1));
  }
  
//doxygen should not read this
/// \cond   
  template<class Line_format>
  void print_line(const Line_format& line_format,const std::string& line){
    std::cout << "[" << line_format.get_record_name(line) << "]" ;
    std::cout << "[" << line_format.get_atom_serial_number(line) << "]" ;
    std::cout << "[" << line_format.get_atom_name(line) << "]" ;
    std::cout << "[" << line_format.get_alternate_location(line) << "]" ;
    std::cout << "[" << line_format.get_residue_name(line) << "]" ;
    std::cout << "[" << line_format.get_chain_identifier(line) << "]" ;
    std::cout << "[" << line_format.get_residue_sequence_number(line) << "]" ;
    std::cout << "[" << line_format.get_insertion_code(line) << "]" ;
    std::cout << "[" << line_format.get_x(line) << "]" ;
    std::cout << "[" << line_format.get_y(line) << "]" ;
    std::cout << "[" << line_format.get_z(line) << "]" ;
    std::cout << "[" << line_format.get_occupancy(line) << "]" ;
    std::cout << "[" << line_format.get_temperature_factor(line) << "]" ;
    std::cout << "[" << line_format.get_element(line) << "]" ;
    std::cout << "[" << line_format.get_charge(line) << "]" ;
  }    
/// \endcond
  
  /**
   * Extract information from a line in a PDB file and add them to a system.
   * @tparam Line_format is a class providing functions to extract fields from a line of a PDB file (ESBTL::PDB::Line_format for example.).
   * @param line_format provides the functions to extract the fields.
   * @param line is a line of a PDB file.
   * @param system_info is the number of the system described by this line.
   */
  template<class Line_format>
  void interpret_line(const Line_format& line_format,const std::string& line,int system_info){
    if (system_info==RMK){
      if (line_format.record_type()==PDB::MODEL)
        current_model=line_format.get_model_number(line);
    }
    else{
      //ensure that references remains valid. You must use the function reserve for this purpose.
      assert(systems_.capacity()>=(unsigned) system_info);
      if (!prebuilt_systems && systems_.size()<(unsigned) system_info)
        systems_.push_back(System(system_info));
      systems_[system_info-1].interpret_line(line_format,line,current_model);
    }
  }

  
  
  /**
   * Finalize the creation of the systems incrementally built.
   * @param altloc indicate the alternate location identification used for systems.
   */
  void create_systems(char altloc){
    for (typename System_container::iterator it=systems_.begin();it!=systems_.end();++it)
      it->set_altloc(altloc);
  }
  
};

} //namespace ESBTL

#endif // ESBTL_BUILDER_H
