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



#ifndef ESBTL_GLOBAL_FUNCTIONS_H
#define ESBTL_GLOBAL_FUNCTIONS_H

#include <fstream>

namespace ESBTL {

  
/**
  * Check if an atom belong to the backbone of a protein. The check is done
  * on the atom name. An atom is considered to belong to the backbone if its atom name
  * is either \c N, \c C, \c CA, \c OXT or \c O.
    * @tparam Atom must represent an atom, global function ESBTL::get_atom_name and ESBTL::get_is_hetatm, 
    * both taking an object of this type must exits.
    * @param atom represents an atom.
  */
template <class Atom>
bool is_backbone(const Atom& atom){
  if (get_is_hetatm(atom))
    return false;
  std::string atomname=get_atom_name(atom);
  return atomname=="N"   ||
         atomname=="CA"  ||
         atomname=="C"   ||
         atomname=="O"   ||
         atomname=="OXT";
}

/**
  * Check if an atom belong to a protein side-chain or is carbon alpha atom. The check is done
  * on the atom name. An atom is considered to belong to the backbone if its atom name
  * is either \c N, \c C, \c CA, \c OXT or \c O. Carbon alpha atom name is \c CA.
    * @tparam Atom must represent an atom, global function ESBTL::get_atom_name and ESBTL::get_is_hetatm, 
    * both taking an object of this type must exits.
    * @param atom represents an atom.
  */
template <class Atom>
bool is_side_chain_or_CA(const Atom& atom){
  if (atom.is_hetatm())
    return false;
  return !is_backbone(atom) || get_atom_name(atom)=="CA";
}

/**
  * Check if an atom is an hydrogen.
  * @tparam Atom must represent an atom, global function ESBTL::get_element 
  * taking an object of this type must exits.
  * @param atom represents an atom.
  */
template <class Atom>
bool is_hydrogen(const Atom& atom){
  return get_element(atom) == std::string("H");
}

/**
  * Check if an atom belongs to a water molecule. An atom is considered to
  * belong to a water molecule if its residue name is either \c HOH, \c SOL or \c WAT.
  * @tparam Atom must represent an atom, global function ESBTL::get_residue_name 
  * taking an object of this type must exits.
  * @param atom represents an atom.
  */
template <class Atom>
bool is_water(const Atom& atom){
  std::string resname=get_residue_name(atom);
  return resname=="HOH" || resname=="SOL" || resname=="WAT";
}


/** 
 * Create a python script that can be run using PyMol and that load
 * a cgo objects representing a set of (pseudo-) atoms as a set of balls.
 * \tparam Atom_iterator iterator over a set of atoms.
 * \tparam Classifier_with_radius is a classifier that attributes a radius to an atom.
 * \tparam Color_of is a classifier that attributes a color to an atom (such as ESBTL::Color_of_atom).
 * \param filename is the name of the script to be created.
 * \param first is the iterator of the first atom of the set.
 * \param last is the past-end iterator of the last atom of the set.
 * \param classifier attributes a radius to an atom.
 * \param color_of attributes a color to an atom.
 */
template<class Atom_iterator,class Classifier_with_radius,class Color_of>
void write_to_cgo(const std::string& filename,Atom_iterator first, Atom_iterator last,
                  const Classifier_with_radius& classifier,const Color_of& color_of)
{
  std::ofstream st(filename.c_str());
  st << "from pymol.cgo import *\n";
  st << "from pymol import cmd\n";  
  st << "obj=[\n";  
  
  for (Atom_iterator it=first;it!=last;++it){
    st << "COLOR," << color_of(*it);
    st << ",SPHERE, " <<it->x() << ", " << it->y() << ", " << it->z() << ", ";
    st << get_radius(classifier.get_properties(*it)) << ",\n";
  }
  
  st << "]\ncmd.load_cgo(obj,\'cgo01\')\n";
}

}//namespace ESBTL

#endif //ESBTL_GLOBAL_FUNCTIONS_H
