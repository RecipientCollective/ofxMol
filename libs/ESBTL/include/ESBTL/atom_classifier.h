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


#ifndef ESBTL_ATOM_CLASSIFIER_H
#define ESBTL_ATOM_CLASSIFIER_H

#include <boost/unordered_map.hpp>
#include<boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>
#include <fstream>


namespace ESBTL{

/**
  * \defgroup prop_classif Property class
  * Concept of Property class to be given to ESBTL::Generic_classifier
  *
  * \code
  *  struct Base_property_concept{
  *    typedef ...                   Query_type; // the object type to which properties are associated
  *    typedef ...                   Key_type;   // the hash key type
  *    
  *    //Constructor using a string stream (from file for example)
  *    //Defines a new property from a string. index indicates its position in the property vector.
  *    Base_property_concept(std::stringstream& ss,unsigned index){ ... }
  *
  *    //Create a key from a Query_type object
  *    static Key_type make_key(const Query_type& query){ ... }
  *
  *    //Static methods indicating what are the default properties associated  
  *    //  -Dictionary an associative container, mapping a Query_type to the index of a property.
  *    //  -Vector_properties is a vector containing all the properties.
  *    //  -The integer returned is the size of the property vector.
  *    template <class Dictionary,class Vector_properties>
  *    static unsigned default_loader(Dictionary& dict,Vector_properties& vect){ ... }
  *
  *    //Static methods indicating how to add to read a classification from a string stream (from a file for example) 
  *    //  Dictionary is a container associating an Key_type object to an unsigned integer (index of the property).
  *    //  The unsigned integer return is the index of the property the object is associated to.
  *    template <class Dictionary>
  *    static unsigned add_classification(std::stringstream& ss,Dictionary& dict){ ... }
  *
  *    //Static methods to handle extra information read from a string stream (from a file for example).  
  *    static void handle_extra(std::stringstream& ss){ ... }
  *
  *    //This function returns the index returned by default when no entry is found in the dictionary.
  *    //If the integer returned is -1, then this indicates that no default properties are expected to be provided
  *    //and will result in an error message.
  *    static int& index_of_default(){ ... }
  *  };
  * \endcode
  */

/**
 * A property class associating a radius to an atom.
 * @tparam NT is the number type used for the radius.
 * @tparam Atom is the atom type.
 * \ingroup prop_classif 
 */
template <class NT,class Atom>
class Radius_of_atom{
private:
  typedef Radius_of_atom<NT,Atom>  Self;
  NT radius_;
  unsigned index_;

public:
  /** The type on which the queries are made, */
  typedef Atom                  Query_type;
  /** The type of the hash key*/
  typedef std::string           Key_type;
  /** The number type of the property stored, here the radius.*/
  typedef NT                    Value_type;
  

  /**
   *Constructor
   * @param radius is the radius associated to this property class.
   * @param index is the index of the property.
   */
  //TODO check why index is useful
  Radius_of_atom(const NT& radius,const unsigned& index):radius_(radius),index_(index){}
  NT value() const {return radius_;}

  /**
   *Constructor
   * @param ss is a string stream containing information about the property.
   * ss must contains the radius associated to index (ex: 1.4).
   * @param index is the index of the property.
   */
  //TODO check why index is useful  
  Radius_of_atom(std::stringstream& ss,unsigned index):index_(index){
    ss >> radius_; 
  }
  
  /**
   * Function filling default radius of atoms. Current implementation uses radii
   * from Tsai J, Taylor R, Chothia C, Gerstein M. J Mol Biol. 1999 Jul 2;290(1):253-66.
   * @tparam Dictionary an associative container, mapping a Query_type to the index of a property.
   * @tparam Vector_properties is a vector storing the properties.
   * @param dict is the container mapping atoms to the index of a property.
   * @param vect is a vector storing properties.
   */    
  template <class Dictionary,class Vector_properties>
  static unsigned default_loader(Dictionary& dict,Vector_properties& vect)
  {
    /** \cond */
    #include <ESBTL/properties/Tsai_jmb_99_radii.h>
    /** \endcond */
  }
  
  /** 
   * Function defining a unique identifier of an atom type.
   * @param atom is an atom.
   */
  static std::string make_key(const Atom& atom)
  {
    return atom.residue_name()+atom.atom_name();
  }

 /**
   * Function adding the classification of an atom type.
   * @tparam Dictionary an associative container, mapping a Query_type to the index of a property.
   * @param ss is a string stream containing the classification of an atom type.
   * ss must contains in this order the residue name, the atom name and the index of the property (ex: ALA N 5)
   * @param dict is the container mapping atoms to the index of a property.
   */      
  template <class Dictionary>
  static 
  unsigned add_classification(std::stringstream& ss,Dictionary& dict){
    std::string resname,atmname;
    ss >> resname;
    ss >> atmname;
    unsigned prop;
    ss >> prop;
    dict[resname+atmname]=prop;
    return prop;
  }
  
  /** 
   * Unused function in that class (but definition needed by the classifier).
   * @param ss a string stream.
   */
  static void handle_extra(std::stringstream& ss){}
    
  /** 
   * Function indicating the index of the property assigned to an unknown atom type.
   * By default this value is -1, which indicate to the classifier that unknown atoms are not expected.
   */
  static
  int& index_of_default(){
    static int index_of_default=-1;
    return index_of_default;
  }
  
};

/** Extract the radius from a property of type ESBTL::Radius_of_atom */
template <class NT,class Atom>
NT get_radius(const Radius_of_atom<NT,Atom>& property){
  return property.value();
}


 /**
  * Property class associating a name and a radius to an atom.
  * \ingroup prop_classif
  */
template <class NT,class Atom>
class Name_and_radius_of_atom{
private:
  typedef Name_and_radius_of_atom<NT,Atom>  Self;
  /** Property index of oxygen of a water molecule.*/
  // this value must be updated according to file default_atom_properties.h
  static 
  unsigned& water_index(){
    static unsigned index=9;
    return index;
  }

  /** Residue name of water molecules */ 
  static
  std::string water_name(){
    static std::string water_name="Owat";
    return water_name;
  }

public:
  
  typedef Atom                  Query_type;
  typedef std::string           Key_type;
  

  Name_and_radius_of_atom(const std::string& name_,const NT& radius_,const unsigned& index_):name(name_),index(index_),radius(radius_)
  {}
  
  std::string name;
  unsigned index;
  NT radius;
  
  bool is_water() const {
    assert(index!=water_index() || name==water_name());
    return index==water_index();
  }
    
  template <class Dictionary,class Vector_properties>
  static unsigned default_loader(Dictionary& dict,Vector_properties& vect)
  {
    /** \cond */
    #include <ESBTL/properties/default_atom_properties.h>
    /** \endcond */
  }
  
  static std::string make_key(const Atom& atom)
  {
    return atom.residue_name()+atom.atom_name();
  }
  
  static
  int& index_of_default(){
    static int index_of_default=-1;
    return index_of_default;
  }
  
};

 /**
  * Property class associating an unordered name to a pair of atom types, 
  * each identified by an index (that of a property for example).
  * \ingroup prop_classif 
  */
class Name_of_pair{
private:
  typedef Name_of_pair Self;


public:
  

  static
  int& index_of_default(){
    static int index_of_default=-1;
    return index_of_default;
  }
  
  
  //~ typedef std::pair<unsigned,unsigned> Key_type;
  typedef std::string Key_type;
  typedef std::pair<unsigned,unsigned> Query_type;

  Name_of_pair(const std::string& name_,const unsigned& index_):name(name_),index(index_){}
  std::string name;
  unsigned index;
    
  template <class Dictionary,class Vector_properties>
  static unsigned default_loader(Dictionary& dict,Vector_properties& vect)
  {
    /** \cond */
    #include <ESBTL/properties/default_pair_properties.h>
    /** \endcond */
  }
  
  static Key_type make_key(const Query_type& indices)
  {
    unsigned i0=indices.first;
    unsigned i1=indices.second;
    if (i0<i1)
      return boost::lexical_cast<std::string>(i0)+
             boost::lexical_cast<std::string>(i1);
      //~ return indices;
    //~ return std::make_pair(i1,i0);
    return boost::lexical_cast<std::string>(i1)+
             boost::lexical_cast<std::string>(i0);
  }
};



/**
  * An object that help to associate properties to objects.
  * This can be used for example to to associate a radius or a color to an (pseudo-)atom type.
  * \tparam Properties_ must follow the concept of \ref prop_classif.
  */
template<class Properties_>
struct Generic_classifier{
  typedef Properties_ Properties;
  typedef typename Properties_::Key_type                      Key_type;
  typedef typename Properties_::Query_type                    Query_type;
  typedef typename std::vector<Properties_>::iterator         Properties_iterator;
  typedef typename std::vector<Properties_>::const_iterator   Properties_const_iterator;
private:
  typedef boost::unordered_map<Key_type, unsigned>      Internal_map_type;
  typedef std::vector<Properties_>                      Internal_vector_type;
  
  Internal_map_type hmap_;
  Internal_vector_type properties_;
  unsigned max_index_;
  
public:
  
  /** Default constructor. The method Properties_::default_loader is called to fill in
    * the dictionary and the property vector.
    */
  Generic_classifier(){
    max_index_=Properties_::default_loader(hmap_,properties_);
    assert(max_index_==properties_.size());
  }
  
  /** Constructor that uses a file to complete the default properties. Not yet implemented */
  Generic_classifier(std::string filename,int){
    /* may need a map from name to property position in the vector */
    /* This map can be filled by a static function of Properties_  */
  }
  
  /** Constructor using a file, default is not loaded.
    * The file must follow the following format:
    * \code
    * EXTRA
    * ...
    * END
    * CLASSIFICATION
    * ... i
    * END
    * PROPERTIES
    *  0  ...
    * ... ...
    *  n  ...
    * END
    * DEFAULT
    * i
    * END
    * \endcode
    * 
    * There are four different sections providing different kind of information:
    *   - EXTRA: For each line of the file, the function Properties_::handle_extra is called.
    *   - CLASSIFICATION: This section contains information to associate an object to an index.
    *     For each line of the file, the function Properties_::add_classification is called.
    *   - PROPERTIES: This section contains the properties. It must start by the index of the property
    *     followed by the information required to define the property. 
    *     For each line of the file, a property is constructed using the constructor with a string stream and an unsigned integer.
    *     The property index is first extracted and the rest of the line is put into a string stream.
    *   - DEFAULT: This section defines the index of the default property (if needed).
    *
    * Note that none of these sections are mandatory.
    */  
  Generic_classifier(const std::string& filename){
    enum Read_state {PROPERTIES, CLASSIFICATION, END, EXTRA,DEFAULT};
    
    std::ifstream input(filename.c_str());
    
    if (!input){
      std::cerr << "Cannot open file " << filename << std::endl;
      exit( EXIT_FAILURE);
    }
    
    std::string buffer;
    Read_state state=END;

    std::map<unsigned,Properties> properties_map;
    unsigned max_property=0;
    
    
    while(input){
      std::getline(input,buffer); //read a complete line
      
      if (!input) {
        if (input.eof()) break;
        std::cerr << "Fatal error: Error while reading file "<< filename <<"\n";
        exit( EXIT_FAILURE );
      }
      
      boost::trim(buffer); //remove white space at the beginning and at the end of the line
      if (buffer.empty() || buffer[0] == '#') continue; //handle comments and empty lines
      
      switch (state){
        case END:
          if (buffer == "PROPERTIES")     {state=PROPERTIES;     break;}
          if (buffer == "CLASSIFICATION") {state=CLASSIFICATION; break;}
          if (buffer == "EXTRA")          {state=EXTRA;          break;}
          if (buffer == "DEFAULT")          {state=DEFAULT;          break;}
          std::cerr << "Fatal error: Error while reading file "<< filename;
          std::cerr << ", unexpected <|" << buffer << "|> found.\n";
          exit( EXIT_FAILURE );
        case CLASSIFICATION:
        {
          if (buffer=="END")              {state=END;            break;}
          std::stringstream ss (std::stringstream::in | std::stringstream::out);
          ss << buffer;
          //inside a function
          unsigned index=Properties::add_classification(ss,hmap_);
          if (max_property <  index) max_property=index;
          break;
        }
        case PROPERTIES:
        {
          if (buffer=="END")              {state=END;            break;}
          std::stringstream ss (std::stringstream::in | std::stringstream::out);
          ss << buffer;
          unsigned index;
          ss >> index;
          //inside a function
          assert(properties_map.find(index)==properties_map.end());
          properties_map.insert(std::make_pair(index,Properties(ss,index)));
          break;
        }        
        case DEFAULT:
        {
          if (buffer=="END")              {state=END;            break;}
          unsigned i;
          std::stringstream ss (std::stringstream::in | std::stringstream::out);
          ss << buffer;
          ss >> i;
          if ( max_property<i ) max_property=i;
          Properties::index_of_default()=i;
          break;
        }
        case EXTRA:
        {
          if (buffer=="END")              {state=END;            break;}
          //call some function
          std::stringstream ss (std::stringstream::in | std::stringstream::out);
          ss << buffer;          
          Properties::handle_extra(ss);
          break;
        }
      };
    }
    
    
    unsigned i=0;
    for (typename std::map<unsigned,Properties>::const_iterator it=properties_map.begin();it!=properties_map.end();++it){
      if (i++ != it->first){
        std::cerr << "Fatal error: properties number " << i-1 << " could not be found in file " << filename <<"\n";
        exit( EXIT_FAILURE );
      }
      properties_.push_back(it->second);
    }  
    
    if (i <= max_property){
      std::cerr << "Fatal error: Up to " << max_property+1 << " properties are used, but only " << i << " are declared.\n";
      exit( EXIT_FAILURE );    
    }
    
    max_index_=i;
    
  }
  
  /** Returns a property given its index
    * \param i is the index of the desired property (starting from 0).
    */
  const Properties& get_properties(unsigned i) const{
    assert(i <= max_index_);
    return properties_[i];
  }
    
  /** Returns a property associated to a Query_type object.
    * \param query is the object a property is looking for in the dictionary.
    */
  const Properties& get_properties(const Query_type& query) const{
    typename Internal_map_type::const_iterator it=hmap_.find( Properties_::make_key(query) );
    if ( it==hmap_.end() ){
      if (Properties::index_of_default()!=-1)
        return properties_[Properties::index_of_default()];
      std::cerr << "Fatal error: Could not find an entry for " << Properties_::make_key(query);
      std::cerr << " and no default have been defined.\n";
      exit( EXIT_FAILURE );
    }
    return properties_[it->second];
  }
  
  /** returns the number of properties.*/
  unsigned number_of_properties() const {return max_index_;}
  
  /** Returns the first iterator on properties.*/
  Properties_iterator properties_begin(){ return properties_.begin();}
  /** Returns the first iterator on properties, const version.*/
  Properties_const_iterator properties_begin() const { return properties_.begin();}
  /** Returns the past-end iterator on properties.*/
  Properties_iterator properties_end(){ return properties_.end(); }
  /** Returns the past-end iterator on properties, const version.*/
  Properties_const_iterator properties_end() const { return properties_.end(); }
  
};

//TODO better specify this class: value and Value_type are not in the concept.
/** 
  * Function object providing the squared radius of an atom thanks to a classifier.
  * \tparam Atom_classifier is a classifier associating a radius to an atom.
  */
template <class Atom_classifier>
struct Weight_of_atoms{
  
  /** Constructor */
  Weight_of_atoms(const Atom_classifier* atom_classifier):atom_classifier_(atom_classifier){}
    
  /**
    * Returns the squared radius of atom. 
    * \param atom the atom used as query.
    */
  typename Atom_classifier::Properties::Value_type
  operator()(const typename Atom_classifier::Query_type& atom) const {
    typename Atom_classifier::Properties::Value_type value=get_radius(atom_classifier_->get_properties(atom));
    return value*value;
  }
private:
  const Atom_classifier* atom_classifier_;
};

}//namespace ESBTL

#endif //ESBTL_ATOM_CLASSIFIER_H
