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



#ifndef ESBTL_INTERNAL_COMPRESSED_IFSTREAM_H
#define ESBTL_INTERNAL_COMPRESSED_IFSTREAM_H

#include <iostream>
#include <fstream>

namespace ESBTL{
/** Modes available to read a PDB file:
  * - ASCII stands for a non-compressed standard text file.
  * - BZIP2 stands for a compressed file using bzip2 algorithm (see http://www.bzip.org/).
  * - GZIP stands for a compressed file using gzip algorithm (see http://www.gzip.org/). 
  */
enum Reading_mode {ASCII,BZIP2,GZIP};
}//namespace ESBTL

namespace ESBTL{
/** \cond */
namespace internal{
  
  template <Reading_mode mode>
  class Ifstream_compress;

  template<>
  class Ifstream_compress<ASCII>
  {
    public:
    std::ifstream& get(std::ifstream& input) const {return input;}
  };

} //namespace internal
/** \endcond */
} //namespace ESBTL

#endif //ESBTL_INTERNAL_COMPRESSED_IFSTREAM_H
