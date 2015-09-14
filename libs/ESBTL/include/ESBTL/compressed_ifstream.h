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



#ifndef ESBTL_COMPRESSED_IFSTREAM_H
#define ESBTL_COMPRESSED_IFSTREAM_H

#include <ESBTL/internal/compressed_ifstream.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zlib.hpp>

/** \cond */
namespace ESBTL{
  namespace internal{
  template<>
  class Ifstream_compress<BZIP2>
  {
    boost::iostreams::filtering_stream<boost::iostreams::input> input_;
  public:
    boost::iostreams::filtering_stream<boost::iostreams::input>&
    get(std::ifstream& input){
      input_.push(boost::iostreams::bzip2_decompressor());
      input_.push(input);
      return input_;
    }
  };
  
  template<>
  class Ifstream_compress<GZIP>
  {
    boost::iostreams::filtering_stream<boost::iostreams::input> input_;
  public:
    boost::iostreams::filtering_stream<boost::iostreams::input>&
    get(std::ifstream& input){
      input_.push(boost::iostreams::gzip_decompressor());
      input_.push(input);
      return input_;
    }
  };
  
} } //ESBTL::internal
/** \endcond */
#endif //ESBTL_COMPRESSED_IFSTREAM_H
