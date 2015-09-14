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


/** \mainpage ESBTL, a PDB parser and data structure for the structural and geometric analysis of biological macromolecules
  *
  * \b Authors: Julie Bernauer, Frédéric Cazals and Sébastien Loriot
  * 
  * This is the reference documentation for ESBTL. This page describes quickly the overall architecture of
  * the library. For other information (installation, FAQ,...) please have a look at the ESBTL website:
  * http://esbtl.sf.net.
  *
  *
  * \section sec_system System
  * A system provides access to the hierarchical data structure.
  * \n
  * It may contain one or several models.
  * Each model is made of chains. Chains are then made of residues, and residues are made of atoms.
  * \n
  * ESBTL provides a default instantiation for this data model when including the header <ESBTL/default.h>.
  * The types ESBTL::Default_system, ESBTL::Default_system::Model, ESBTL::Default_system::Chain, 
  * ESBTL::Default_system::Residue and ESBTL::Default_system::Atom provide access to all hierarchy levels.
  *
  * \section sec_readfile Reading a file
  * Reading a PDB file is performed in four stages:
  *  - A system type and a container for the systems have to be defined.
  *  - A line selector (see \ref linesel) is used to define what will the systems be made of.
  *  - A builder is used to fill the system container (see ESBTL::All_atom_system_builder for example).
  *  - An occupancy policy has to be chosen (see \ref occpol).
  * 
  * Example:
  *  \code
  *  #include <ESBTL/default.h>
  *  //Create one system with all atoms.
  *  ESBTL::PDB_line_selector sel;
  *  std::vector systems;
  *  //Build the system from the PDB file.
  *  ESBTL::All_atom_system_builder builder(systems,sel.max_nb_systems());
  *  ESBTL::read_a_pdb_file(filename,sel,builder,Accept_none_occupancy_policy()); 
  * \endcode
  *
  * \section iters Iterators
  *
  * Once a system has been built, iterators are provided to access the hierarchy information
  * from any higher level.
  * For 'father' standing for either system, model, chain or residue, and 'son'
  * standing for either model, chain, residue or atom, the following iterators
  * Father::Sons_const_iterator and Father::Sons_iterator, 
  * and the following functions Father::sons_begin, Father::sons_end are provided.
  * Note that if 'father' is system, 'son' can only be model as the
  * hierarchy must be respected.
  * \n
  * See \ref grp_iters for a list of a list of iteration facilities provided by ESBTL.
  *
  * \section assos_prop Associating properties to an object
  * The class ESBTL::Generic_classifier provides the association of a property to an object (for example a radius to an atom).
  * See \ref atomsel for the property concept and properties provided by ESBTL.
  * 
  * \section sec_atom_sel Making a selection of atoms.
  * Selections can be done at two different stages. Either while reading
  * a file thanks to a line selector (each selection will be stored within a different system),
  * or after having built the systems.
  * \n
  * A selection at the parsing stage is made using a line selector: see \ref linesel for details
  * on the object concept and line selectors provided by ESBTL.
  * \n
  * The class ESBTL::Selected_atom_iterator defines an iterator type over a subset of atoms of a model.
  * The subset is defined using a function object. See \ref atomsel for details on the object concept and
  * atom selection function objects provided by ESBTL.
  *
  */

#ifndef ESBTL_DEFAULT_H
#define ESBTL_DEFAULT_H

#include <ESBTL/constants.h>
#include <ESBTL/xyz_utils.h>
#include <ESBTL/molecular_system.h>
#include <ESBTL/PDB.h>
#include <ESBTL/line_selectors.h>
#include <ESBTL/builder.h>
#include <ESBTL/line_reader.h>
#include <ESBTL/occupancy_handlers.h>
#include <ESBTL/coarse_grain.h>

namespace ESBTL {
  /** A default declaration of an all-atom system.*/
  typedef Molecular_system<Default_system_items,Point_3 > Default_system;
  /** A default declaration of an all-atom and coarse grain system.*/
  typedef Molecular_system<System_items_with_coarse_grain,Point_3 > Default_system_with_coarse_grain;
} //namespace ESBTL

#endif //ESBTL_DEFAULT_H
