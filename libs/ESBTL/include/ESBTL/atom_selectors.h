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


#include <set>

#ifndef ESBTL_ATOM_SELECTORS_H
#define ESBTL_ATOM_SELECTORS_H

/** \defgroup atomsel Atom selection
 * This module gathers function objects to make a selection using an atoms.
 * Such a function object must be a class of the following type: 
 * \code
 * class Atom_selection_concept{
 *   template <class Atom>
 *   bool operator() (const Atom& atom) const
 *   { ... }  
 * }
 * \endcode
 * Within such function objects, it is advised to use global function (prefixed by \c get_)
 * so as to be used with an atom class or a molecular file line.
 */

namespace ESBTL{

/** Function object to select an atom using its residue name.
  * \ingroup atomsel
  */
struct Select_by_resname{
  std::string name;
  /** Default constructor */
  Select_by_resname(){}; //needed by selected_atom_iterator::end and default constructor
  /**
    * Constructor.
    * @param str is the name of the residues to be selected.
    */
  Select_by_resname(const std::string& str):name(str){}
  /**
    * Checks if an atoms match a criteria.
    * @tparam Atom must represent an atom, a global function ESBTL::get_residue_name taking an element of this type must exits.
    * @param atom represents an atom.
    */
  template <class Atom>
  bool operator()(const Atom& atom) const {
    return get_residue_name(atom)==name;
  }
};


/** Function object to select an atom using its atom name.
  * \ingroup atomsel
  */
struct Select_by_atmname{
  std::string name;
  /** Default constructor */
  Select_by_atmname(){}; //needed by selected_atom_iterator::end and default constructor

  /**
    * Constructor.
    * @param str is the name of the atoms to be selected.
    */
  Select_by_atmname(const std::string& str):name(str){}
    
  /**
    * Checks if an atoms match a criteria.
    * @tparam Atom must represent an atom, a global function ESBTL::get_atom_name taking an element of this type must exits.
    * @param atom represents an atom.
    */
  template <class Atom>
  bool operator()(const Atom& atom) const {
    return get_atom_name(atom)==name;
  }
};

/** Function object to select an atom using its chemical name.
  * \ingroup atomsel
  */
struct Select_by_element{
  std::string name;
/** Default constructor */
  Select_by_element(){}; //needed by selected_atom_iterator::end and default constructor
  /**
    * Constructor.
    * @param str is the chemical name of the atoms to be selected.
    */
  Select_by_element(const std::string& str):name(str){}
  /**
    * Checks if an atoms match a criteria.
    * @tparam Atom must represent an atom, a global function ESBTL::get_element taking an element of this type must exits.
    * @param atom represents an atom.
    */
  template <class Atom>
  bool operator()(const Atom& atom) const {
    return get_element(atom)==name;
  }
};

/** Function object to select an atom using its chain identifier.
  * \ingroup atomsel
  */
class Select_by_chainids{
  std::set<char> chains;
public:
  /** Default constructor */  
  Select_by_chainids(){}; //needed by selected_atom_iterator::end and default constructor
  /**
    * Constructor.
    * @param chs is the chains identifiers of the atoms to be selected. If several chains must
    * be specified, they must be concatenated (for example selecting atoms of chains A or B corresponds to \e str="AB")
    */
  Select_by_chainids(const std::string& chs){
    for (unsigned i=0;i!=chs.size();++i)
      chains.insert(chs[i]);
  }
  
  /**
    * Checks if an atoms match a criteria.
    * @tparam Atom must represent an atom, a global function ESBTL::get_chain_identifier taking an element of this type must exits.
    * @param atom represents an atom.
    */  
  template <class Atom>
  bool operator()(const Atom& atom) const {
    return chains.find( get_chain_identifier( atom ) ) != chains.end();
  }
};

} //namespace ESBTL

#endif //ESBTL_ATOM_SELECTORS_H
