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



#ifndef ESBTL_LINE_READER_H
#define ESBTL_LINE_READER_H

#include <ESBTL/internal/compressed_ifstream.h>

namespace ESBTL{
  
/** 
  * Class responsible for reading the lines of a file and provide
  * information to a builder to create a system.
  * \tparam Line_format is a helper class that can extract fields from a molecular data file (like ESBTL::PDB::Line_format for example).
  * \tparam Line_selector is a line selector object that must follow the concept of \ref linesel.
  * \tparam Builder is a class able to build a molecular system (like ESBTL::All_atom_system_builder for example).
  */
template <class Line_format,class Line_selector,class Builder>
class Line_reader{
  Line_selector& line_selector;
  Builder& builder;

  template <class Occupancy_handler,class Ifstream>
  bool read_stream(Ifstream& input,Occupancy_handler occupancy,char default_altloc=char(' ')){
    int nblines=0;

    while ( ! input.eof() ){
      std::string line;
      getline(input,line);
      if (line.empty()) continue;
      
      
      Line_format line_format(line);
      
      int system_index=line_selector.keep(line_format,line,occupancy);
      
      if (system_index!=DISCARD){
        //TODO warning Put inside selector: we could think of a different policy
        char altloc=line_format.get_alternate_location(line);
        if (altloc!=' '){
          if (default_altloc==' ')
            default_altloc=altloc;
          else
            if (altloc!=default_altloc) continue;
        }
        
        builder.interpret_line(line_format,line,system_index);
        ++nblines;        
      }
    }
    
    nblines+=occupancy.finalize(builder);
    builder.create_systems(default_altloc);
    
    #ifndef NDEBUG
    std::cout << "(ESBTL-DEBUG) Lines read " << nblines << std::endl;    
    #endif
    return true;
  }
  
public:
  /** Constructor.*/
  Line_reader(Line_selector& line_selector, Builder& builder):line_selector(line_selector),builder(builder){}
  
  
  //template parameter Occupancy_handler tells what to do with atoms with occupancy !=1 (when no altloc present)
  //TODO : think of the same think for altloc: when should have that the sum of the occupancy is 1 when considering these atoms!!!!!!
  //       This would be also a way to handle differently altloc and to avoid the update at the end to be sure all altloc selected are the same 
  //occupancy is not const for the moment since it may be used as internal line storage
  /** Reads the line of a file and give instruction to the builder to construct molecular system(s).
    * \tparam mode indicate the encoding of the file to be read. 
    * To read a compressed file, you must include the file <ESBTL/compressed_ifstream.h>,
    * and the boost_iostreams headers and library must be installed on your system.
    * \tparam Occupancy_handler is an occupancy policy that must be a model of the concept \ref occpol.
    * \param filename is the name of the file to be read.
    * \param occupancy is the occupancy policy
    * \param default_altloc is the alternate location identification that must be selected.
    * If the white space character is specified, the first alternate location identification
    * read will be used as default.
    *
    * A version of this function with reading mode ASCII is also provided.
    */
  template <Reading_mode mode,class Occupancy_handler>
  bool read(const std::string& filename,Occupancy_handler occupancy,char default_altloc=' '){
    //binary flags prevent from interpretation of some symbols \r\n for example
    std::ifstream input( filename.c_str(),std::ios_base::binary );

    if (! input){
      std::cerr << "Problem while trying to open file " << filename << ".\nPlease check that the file exists and that you have the right to read it."  << std::endl;
      return false;   
    }
    
    internal::Ifstream_compress<mode> stream;
    return read_stream(stream.get(input),occupancy,default_altloc);
  }

  template <class Occupancy_handler>
  bool read(const std::string& filename,Occupancy_handler occupancy,char default_altloc=' '){
    return read<ASCII>(filename,occupancy,default_altloc);
  }
  
  
};

/**
 * Short cuts to read a PDB file. See the documentation of ESBTL::Line_reader::read for more details.
 */
template<Reading_mode mode,class Line_selector,class Builder,class Occupancy_handler>
inline
bool read_a_pdb_file(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc=' '){
  return Line_reader<PDB::Line_format<>,Line_selector,Builder>(sel,builder).template read<mode>(filename,occupancy,altloc);
}

template<class Line_selector,class Builder,class Occupancy_handler>
inline
bool read_a_pdb_file(const std::string& filename,Line_selector& sel,Builder& builder,const Occupancy_handler& occupancy,char altloc=' '){
  return read_a_pdb_file<ASCII>(filename,sel,builder,occupancy,altloc);
}



} //namespace ESBTL

#endif // ESBTL_LINE_READER_H
