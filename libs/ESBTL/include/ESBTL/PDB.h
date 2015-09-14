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




#ifndef ESBTL_PDB_H
#define ESBTL_PDB_H

#include<boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sstream>
#include <iostream>
#include <typeinfo>
#include <ESBTL/constants.h>



#define RECOVER_FIELD(FNAME,TYPE,FROM,TO,DEF) \
  TYPE get_##FNAME(const std::string& line) const { \
    if (line.length() < TO){\
      if (! Mandatory_fields::FNAME) return TYPE(DEF);\
      else{\
          std::cerr << "Fatal error: Cannot extract field \'" << #FNAME << "\' in \n";\
          std::cerr << "<|" << line << "|>\n";\
          exit (EXIT_FAILURE);\
      }\
    }\
    return PDB::extract_field<TYPE>(line,FROM,TO,DEF,#FNAME,Mandatory_fields::FNAME); \
  }
  

namespace ESBTL{
  namespace PDB{
  /** Enumeration of the type of a line read in a PDB file.*/
  enum Record_type{ATOM=0,HETATM,MODEL,ENDMDL,TER,END,ANISOU,CONECT,MASTER,UNKNOWN};
  
  /** \cond */
  template<class T>
  T extract_field(const std::string& line,unsigned from, unsigned to,T default_value,const char* name,bool is_mandatory=false){
    std::string s=line.substr(from,to-from+1);
    boost::trim(s);
    if (s.length()==0){
      if (!is_mandatory)
        return default_value;
      else{
        std::cerr << "Fatal error: field \'"<< name << "\' is empty in \n";
        std::cerr  << "<|" << line << "|>\n";
        exit (EXIT_FAILURE);
      }
    }
      try{
        return boost::lexical_cast<T>(s);
      }
      catch(boost::bad_lexical_cast &)
      {
        std::cerr << "Fatal error: Cannot convert \'"<< s << "\' to " << typeid(T).name() << " in line\n";
        std::cerr  << "<|" << line << "|>\n";
        exit (EXIT_FAILURE);
      }
  }
  /** \endcond */
  
  /** Default class to specify to ESBTL::PDB::Line_format which PDB fields of a coordinate line are mandatory.*/
  struct Mandatory_fields_default{
    static const bool record_name=true;
    static const bool atom_serial_number=true;
    static const bool atom_name=true;
    static const bool alternate_location=false;
    static const bool residue_name=true;
    static const bool chain_identifier=false;
    static const bool residue_sequence_number=true;
    static const bool insertion_code=false;
    static const bool x=true;
    static const bool y=true;
    static const bool z=true;
    static const bool occupancy=true;
    static const bool temperature_factor=true;
    static const bool element=true;
    static const bool charge_str=false;
    static const bool model_number=true;
  };
 
  /** 
    * Helper class handle to extract fields from line of a PDB file.
    * \tparam Mandatory_fields is a class indicating which fields are mandatory (see ESBTL::PDB::Mandatory_fields_default for requirements).
    */
  //TODO put the string line inside ?
  template <class Mandatory_fields=Mandatory_fields_default>
  class Line_format{
  private:  
    PDB::Record_type type;

  public:
    //Easy to do with a MACRO
    /** Constructor.
      *\param line is a line of a PDB file.
      */
    Line_format(const std::string& line){
      type=PDB::UNKNOWN;
      if (strncmp(line.c_str(),"ATOM",4)==0)    type=PDB::ATOM;
      if (strncmp(line.c_str(),"HETATM",6)==0)  type=PDB::HETATM;
      if (strncmp(line.c_str(),"MODEL",5)==0)   type=PDB::MODEL;
      if (strncmp(line.c_str(),"ENDMDL",6)==0)  type=PDB::ENDMDL;
      if (strncmp(line.c_str(),"TER",3)==0)     type=PDB::TER;
      if (strncmp(line.c_str(),"END",3)==0)     type=PDB::END;
      if (strncmp(line.c_str(),"ANISOU",6)==0)  type=PDB::ANISOU;
      if (strncmp(line.c_str(),"CONECT",6)==0)  type=PDB::CONECT;
      if (strncmp(line.c_str(),"MASTER",6)==0)  type=PDB::MASTER;
    }
    
    /** Indicates whether the line read is a coordinate line of an hetero-atom.*/
    bool is_hetatm() const {return (type==PDB::HETATM);}
    
    //ATOM and HETATM fields
    RECOVER_FIELD(record_name,std::string,0,5," ")
    RECOVER_FIELD(atom_serial_number,int,6,10,-1)
    RECOVER_FIELD(atom_name,std::string,12,15," ")
    RECOVER_FIELD(alternate_location,char,16,16,' ')
    RECOVER_FIELD(residue_name,std::string,17,19," ")
    RECOVER_FIELD(chain_identifier,char,21,21,' ')
    RECOVER_FIELD(residue_sequence_number,int,22,25,-1)
    RECOVER_FIELD(insertion_code,char,26,26,' ')
    RECOVER_FIELD(x,double,30,37,NO_FLOAT)
    RECOVER_FIELD(y,double,38,45,NO_FLOAT)
    RECOVER_FIELD(z,double,46,53,NO_FLOAT)
    RECOVER_FIELD(occupancy,double,54,59,NO_FLOAT)
    RECOVER_FIELD(temperature_factor,double,60,65,NO_FLOAT)
    RECOVER_FIELD(element,std::string,76,77," ")
    /** extract the field charge as a string. */
    RECOVER_FIELD(charge_str,std::string,78,79," ")
    
    /** extract the field charge as an integer. */
    int get_charge(const std::string& line) const {
      std::string charge_str=get_charge_str(line);
      if (charge_str==" ")
        return NO_CHARGE;
      else{
        int ch=1;
        if (charge_str.size()==2){
          if (charge_str[1]=='-') ch*=-1;
          else{
            if (charge_str[1]!='+'){
              std::cerr << "Fatal error: Cannot extract field \'charge\' in \n";
              std::cerr << "<|" << line << "|>\n";\
              exit(EXIT_FAILURE);
            }
          }
        }
        try{
          return ch*boost::lexical_cast<int>(charge_str[0]);
        }
        catch(boost::bad_lexical_cast &){
          std::cerr << "Fatal error: Cannot convert \'"<< charge_str[0] << "\' to " << typeid(int).name() << " in line\n";
          std::cerr  << "<|" << line << "|>\n";
          exit (EXIT_FAILURE);          
        }
      }
    }
    
    //~ //element from atom name
    //~ RECOVER_FIELD(element_using_atom_name,std::string,12,13," ")
    
    //MODEL fields
    RECOVER_FIELD(model_number,int,10,13,-1)
    
    
    PDB::Record_type record_type() const{
      return type;
    }
    
  };  
  

  template<class PDB_Atom>
  std::string 
  get_atom_pdb_format(const PDB_Atom& atom)
  {
    std::stringstream buf(std::stringstream::out);
    buf << ( !atom.is_hetatm() ? "ATOM  ":"HETATM" );
    buf <<  boost::format("%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s")
            % atom.atom_serial_number()
            % atom.atom_name().c_str()
            % atom.alternate_location()
            % atom.residue_name().c_str()
            % atom.chain_identifier()
            % atom.residue_sequence_number()
            % atom.insertion_code()
            % atom.x()
            % atom.y()
            % atom.z()
            % atom.occupancy()
            % atom.temperature_factor()
            % atom.element().c_str();
    if (atom.charge()==NO_CHARGE)   buf << "  ";
    else buf << boost::format("% 2i") % atom.charge();
    return buf.str();
  }

  template<class PDB_Atom>
  std::string 
  get_atom_pdb_reduced_format(const PDB_Atom& atom)
  {
    std::stringstream buf(std::stringstream::out);
    buf <<  boost::format("%5d %4s%c%3s %c%4d%c")
            % atom.atom_serial_number()
            % atom.atom_name().c_str()
            % atom.alternate_location()
            % atom.residue_name().c_str()
            % atom.chain_identifier()
            % atom.residue_sequence_number()
            % atom.insertion_code();
    return buf.str();
  }
  
} //namespace PDB

//global access functions
  template <class Mandatory_fields>
  bool get_is_hetatm(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.is_hetatm(p.second);}
  template <class Mandatory_fields>
  int get_atom_serial_number(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_atom_serial_number(p.second) ;}
  template <class Mandatory_fields>
  std::string get_atom_name(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_atom_name(p.second) ;}
  template <class Mandatory_fields>
  char get_alternate_location(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_alternate_location(p.second) ;}
  template <class Mandatory_fields>
  double get_occupancy(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_occupancy(p.second) ;}
  template <class Mandatory_fields>
  double get_temperature_factor(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_temperature_factor(p.second) ;}
  template <class Mandatory_fields>
  std::string get_element(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_element(p.second) ;}
  template <class Mandatory_fields>
  int get_charge(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_charge(p.second) ;}
  template <class Mandatory_fields>
  char get_chain_identifier(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_chain_identifier(p.second) ;}
  template <class Mandatory_fields>
  std::string get_residue_name(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_residue_name(p.second) ;}
  template <class Mandatory_fields>
  int get_residue_sequence_number(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_residue_sequence_number(p.second) ;}
  template <class Mandatory_fields>
  char get_insertion_code(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_insertion_code(p.second) ;}
  template <class Mandatory_fields>
  char get_x(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_x(p.second) ;}
  template <class Mandatory_fields>
  char get_y(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_y(p.second) ;}
  template <class Mandatory_fields>
  char get_z(const std::pair<PDB::Line_format<Mandatory_fields>,std::string>& p) {return p.first.get_z(p.second) ;}
} //namespace ESBTL

#endif //ESBTL_PDB_H
