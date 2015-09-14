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



#ifndef ESBTL_COARSE_CLASSIFIER_H
#define ESBTL_COARSE_CLASSIFIER_H


#include <ESBTL/atom_classifier.h>

namespace ESBTL{


/**
 * A property class associating a radius to an coarse atom.
 * @tparam NT is the number type used for the radius.
 * @tparam Coarse_atom is the pseudo-atom type.
  * \ingroup prop_classif
  */  
template <class NT,class Coarse_atom>
class Radius_of_coarse_atom{
private:
  typedef Radius_of_coarse_atom<NT,Coarse_atom>  Self;
  NT radius_;
public:
  typedef Coarse_atom           Query_type;
  typedef std::string           Key_type;
  typedef NT                    Value_type;
  

  Radius_of_coarse_atom(const NT& radius):radius_(radius){}
  NT value() const {return radius_;}
  
    
  template <class Dictionary,class Vector_properties>
  static unsigned default_loader(Dictionary& dict,Vector_properties& vect)
  {
    /** \cond */
    #include <ESBTL/properties/default_coarse_radii.h>
    /** \endcond */
  }
  
  /** The name of the residue is used as key */
  static std::string make_key(const Coarse_atom& atom)
  {
    return atom.residue().residue_name();
  }
  
  static
  int& index_of_default(){
    static int index_of_default=-1;
    return index_of_default;
  }  
};

/** Extract the radius from a property of type ESBTL::Radius_of_coarse_atom */
template <class NT,class Coarse_atom>
NT get_radius(const Radius_of_coarse_atom<NT,Coarse_atom>& property){
  return property.value();
}

/** Function object associating a color to an atom, according to
  * its residue name.
  * -yellow (1,1,0): ALA, CYS, GLY, PRO, SER, THR.
  * -green (0,0.8,0): VAL, LEU, ILE, MET, MSE, PHE, TYR, TRP.
  * -blue (0,0,0.93): HIS, LYS, ARG.
  * -purple (0.54,0.04,0.31): ASN, GLN.
  * -red (0.77,0,0): GLU, ASP.
  * -white (0.5,0.5,0.5): any other residue.
  * \tparam Atom is an atom type.
  */
template<class Atom>
struct Color_of_atom{
  /** returns the color associated to \c atm.
    * \return a string color (red corresponds to "1,0,0")
    */
  std::string operator() (const Atom& atm) const{
    std::string name=atm.residue().residue_name();
    if ( name == "ALA" || name == "CYS" || name == "GLY" ||
          name == "PRO" || name == "SER" || name == "THR" )
      return std::string("1,1,0"); // yellow: normal
    
    if ( name == "VAL" || name == "LEU" || name == "ILE" ||
           name == "MET" || name == "MSE" || name == "PHE" ||
           name == "TYR" || name == "TRP" )
        return std::string("0,0.80,0");// green: hydrophobic
    
    if ( name == "HIS" || name == "LYS" || name == "ARG")
      return std::string("0,0,0.93");//blue: positive charge
    
    if ( name == "ASN" || name == "GLN" )
      return std::string("0.54,0.04,0.31 ");//purple: +/- charge
    
    if ( name == "GLU" || name == "ASP" )
        return std::string("0.77,0,0 ");//red: negative charge
    return std::string("0.5,0.5,0.5 ");//unknown residue
  }    
};

/** 
  * A property class associating a color to an atom using its residue name.
  * Default color codes are the same as ESBTL::Color_of_atom function object.
  * \tparam Atom is an atom type.
  * \ingroup prop_classif
  */
template <class Atom>
class Color_of_residues{
  std::string color_;
public:  
  static
  int& index_of_default(){
    static int index_of_default=-1;
    return index_of_default;
  }
  
  
  //~ typedef std::pair<unsigned,unsigned> Key_type;
  typedef std::string Key_type;
  typedef Atom Query_type;

  Color_of_residues(const std::string& color):color_(color){}
    
  template <class Dictionary,class Vector_properties>
  static unsigned default_loader(Dictionary& dict,Vector_properties& vect)
  {
    dict["ALA"]=0;
    dict["CYS"]=0;
    dict["GLY"]=0;
    dict["PRO"]=0;
    dict["SER"]=0;
    dict["THR"]=0;
    vect.push_back(std::string("1,1,0")); // yellow: normal
    dict["VAL"]=1;
    dict["LEU"]=1;
    dict["ILE"]=1;
    dict["MET"]=1;
    dict["MSE"]=1;
    dict["PHE"]=1;
    dict["TYR"]=1;
    dict["TRP"]=1;
    vect.push_back(std::string("0,0.80,0"));// green: hydrophobic
    dict["HIS"]=2;
    dict["LYS"]=2;
    dict["ARG"]=2;
    vect.push_back(std::string("0,0,0.93"));//blue: positive charge
    dict["ASN"]=3;
    dict["GLN"]=3;
    vect.push_back(std::string("0.54,0.04,0.31 "));//purple: +/- charge
    dict["GLU"]=4;
    dict["ASP"]=4;
    vect.push_back(std::string("0.77,0,0 "));//red: negative charge
    vect.push_back(std::string("0.5,0.5,0.5 "));//unknown residue
    index_of_default()=4;
    return 5;
  }
  
  static Key_type make_key(const Query_type& atom)
  {
    return get_residue_name(atom);
  }
};

}//namespace ESBTL

#endif //ESBTL_COARSE_CLASSIFIER_H
