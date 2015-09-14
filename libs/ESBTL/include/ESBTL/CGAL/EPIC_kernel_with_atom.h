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



#ifndef ESBTL_CGAL_EPIC_KERNEL_WITH_ATOM
#define ESBTL_CGAL_EPIC_KERNEL_WITH_ATOM

//include files for convenience
#include <ESBTL/constants.h>
#include <ESBTL/xyz_utils.h>
#include <ESBTL/line_selectors.h>
#include <ESBTL/builder.h>
#include <ESBTL/line_reader.h>
#include <ESBTL/PDB.h>
#include <ESBTL/occupancy_handlers.h>
//include files needed
#include <ESBTL/molecular_system.h>
#include <ESBTL/coarse_grain.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>


namespace ESBTL{

namespace CGAL{
  
//~ template <class K>
//~ class Construct_atom_3 : public ::CGAL::CartesianKernelFunctors::Construct_point_3<K>{
//~ public:
  //~ #ifndef CGAL_CFG_MATCHING_BUG_6 
    //~ using ::CGAL::CartesianKernelFunctors::Construct_point_3<K>::operator();
  //~ #else 

  //~ #endif

  //~ template <class Items,class Point_base_3>
  //~ typename 
  //~ K::Point_3 operator()(const Molecular_system<Items, Point_base_3>& atom) const {
    //~ return typename K::Point_3(0,0,0);
  //~ }

//~ };


//Kernel for atoms
// K_ is the new kernel, and K_Base is the old kernel
template < typename K_, typename K_Base >
class Cartesian_kernel_base
  : public K_Base::template Base<K_>::Type
{
public:
  typedef K_                                    Kernel;
  typedef typename Molecular_system<Default_system_items,typename K_Base::Point_3 >::Atom  Point_3;
 
  //~ typedef Construct_atom_3<K_>                          Construct_point_3;
  //~ Construct_point_3 construct_point_3_object() const {return Construct_point_3();}

  template < typename Kernel2 >
  struct Base { typedef Cartesian_kernel_base<Kernel2, K_Base>  Type; };
};


template < typename FT_ >
struct Kernel_with_atom
  : public ::CGAL::Type_equality_wrapper<
                Cartesian_kernel_base<Kernel_with_atom<FT_>, ::CGAL::Simple_cartesian<FT_> >,
                Kernel_with_atom<FT_> >
{};
  
  
  
//Kernel for coarse grain
// K_ is the new kernel, and K_Base is the old kernel
template < typename K_, typename K_Base >
class Cartesian_kernel_base_cg
  : public K_Base::template Base<K_>::Type
{
  typedef typename Molecular_system<System_items_with_coarse_grain,typename K_Base::Point_3 >::Atom Atom;
public:
  typedef K_                                          Kernel;
  typedef Coarse_atom<Atom,typename K_Base::Point_3>  Point_3;

  template < typename Kernel2 >
  struct Base { typedef Cartesian_kernel_base_cg<Kernel2, K_Base>  Type; };
};


template < typename FT_ >
struct Kernel_with_coarse_atom
  : public ::CGAL::Type_equality_wrapper<
                Cartesian_kernel_base_cg<Kernel_with_coarse_atom<FT_>, ::CGAL::Simple_cartesian<FT_> >,
                Kernel_with_coarse_atom<FT_> >
{};
  
/** Short cut type to define an all atom system with atom type being a point of CGAL kernel EPIC_kernel_with_atom*/
typedef Molecular_system<Default_system_items, ::CGAL::Simple_cartesian<double>::Point_3 >             Default_system;
/** Short cut type to defined an exact predicates inexact constructions CGAL kernel with an atom as point type.
  *
  *See http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Kernel_23_ref/Class_Exact_predicates_inexact_constructions_kernel.html*/
typedef ::CGAL::Filtered_kernel< Kernel_with_atom<double> >                                            EPIC_kernel_with_atom;

/** Short cut type to define an coarse grain system with coarse grain atom type being a point of CGAL kernel EPIC_kernel_with_atom*/  
typedef Molecular_system<System_items_with_coarse_grain, ::CGAL::Simple_cartesian<double>::Point_3 >   System_with_coarse_grain;
  /** Short cut type to defined an exact predicates inexact constructions CGAL kernel with an atom as point type.
    *
    *See http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Kernel_23_ref/Class_Exact_predicates_inexact_constructions_kernel.html*/
typedef ::CGAL::Filtered_kernel< Kernel_with_coarse_atom<double> >                                     EPIC_kernel_with_coarse_atom;

}  } //namespace ESBTL::CGAL


#endif //ESBTL_CGAL_EPIC_KERNEL_WITH_ATOM

