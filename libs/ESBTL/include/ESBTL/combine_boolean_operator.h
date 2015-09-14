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



#ifndef ESBTL_COMBINE_BOOLEAN_OPERATOR_H
#define ESBTL_COMBINE_BOOLEAN_OPERATOR_H


namespace ESBTL{

/**
  * Class defining a function object returning the opposite result of a given function object.
  * @tparam S is a function object with an operator() taking one argument and returning a Boolean (such as \ref atomsel).
  */
template <class S>
struct Not_functor:public S{
  Not_functor(){}
  Not_functor(const S& s):S(s){}
  
  template <class T>
  bool operator()(const T& t) const {
    return ! S::operator()(t);
  }
};
  
#ifndef __GXX_EXPERIMENTAL_CXX0X__
/** \cond */
template <class S1=void,class S2=void,class S3=void,class S4=void,class S5=void,class S6=void,class S7=void,class S8=void,class S9=void,class S10=void>
class And_functors;


template <class S1,class S2>
class And_functors<S1,S2>: public S2{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2):s1_(s1),S2(s2){}

  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && S2::operator()(i);
  }  
};

template <class S1,class S2,class S3>
class And_functors<S1,S2,S3>: public And_functors<S2,S3>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3):s1_(s1),And_functors<S2,S3>(s2,s3){}
    
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4>
class And_functors<S1,S2,S3,S4>: public And_functors<S2,S3,S4>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4):s1_(s1),And_functors<S2,S3,S4>(s2,s3,s4){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5>
class And_functors<S1,S2,S3,S4,S5>: public And_functors<S2,S3,S4,S5>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5):s1_(s1),And_functors<S2,S3,S4,S5>(s2,s3,s4,s5){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6>
class And_functors<S1,S2,S3,S4,S5,S6>: public And_functors<S2,S3,S4,S5,S6>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6):s1_(s1),And_functors<S2,S3,S4,S5,S6>(s2,s3,s4,s5,s6){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5,S6>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6,class S7>
class And_functors<S1,S2,S3,S4,S5,S6,S7>: public And_functors<S2,S3,S4,S5,S6,S7>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7):s1_(s1),And_functors<S2,S3,S4,S5,S6,S7>(s2,s3,s4,s5,s6,s7){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5,S6,S7>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8>
class And_functors<S1,S2,S3,S4,S5,S6,S7,S8>: public And_functors<S2,S3,S4,S5,S6,S7,S8>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8):s1_(s1),And_functors<S2,S3,S4,S5,S6,S7,S8>(s2,s3,s4,s5,s6,s7,s8){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5,S6,S7,S8>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8,class S9>
class And_functors<S1,S2,S3,S4,S5,S6,S7,S8,S9>: public And_functors<S2,S3,S4,S5,S6,S7,S8,S9>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8,const S9& s9):s1_(s1),And_functors<S2,S3,S4,S5,S6,S7,S8,S9>(s2,s3,s4,s5,s6,s7,s8,s9){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5,S6,S7,S8,S9>::operator()(i);
  }  
};
/** \endcond */

/**
  * Class defining a function object returning the logical \b AND between the results of several given function objects.
  * The current implementation provides a version from two to ten parameters. If the code is compiled using
  * the c++0x standard, the number of parameters is not limited (use \c -std=c++0x with \c gcc).
  * @tparam Si are function objects with an operator() taking one argument and returning a Boolean (such as \ref atomsel).
  */
template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8,class S9,class S10>
class And_functors: public And_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8,const S9& s9,const S10& s10):And_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>::operator()(i);
  }  
};

/** \cond */
template <class S1=void,class S2=void,class S3=void,class S4=void,class S5=void,class S6=void,class S7=void,class S8=void,class S9=void,class S10=void>
class Or_functors;


template <class S1,class S2>
class Or_functors<S1,S2>: public S2{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2):s1_(s1),S2(s2){}

  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || S2::operator()(i);
  }  
};

template <class S1,class S2,class S3>
class Or_functors<S1,S2,S3>: public Or_functors<S2,S3>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3):s1_(s1),Or_functors<S2,S3>(s2,s3){}
    
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4>
class Or_functors<S1,S2,S3,S4>: public Or_functors<S2,S3,S4>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4):s1_(s1),Or_functors<S2,S3,S4>(s2,s3,s4){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5>
class Or_functors<S1,S2,S3,S4,S5>: public Or_functors<S2,S3,S4,S5>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5):s1_(s1),Or_functors<S2,S3,S4,S5>(s2,s3,s4,s5){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6>
class Or_functors<S1,S2,S3,S4,S5,S6>: public Or_functors<S2,S3,S4,S5,S6>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6):s1_(s1),Or_functors<S2,S3,S4,S5,S6>(s2,s3,s4,s5,s6){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5,S6>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6,class S7>
class Or_functors<S1,S2,S3,S4,S5,S6,S7>: public Or_functors<S2,S3,S4,S5,S6,S7>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7):s1_(s1),Or_functors<S2,S3,S4,S5,S6,S7>(s2,s3,s4,s5,s6,s7){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5,S6,S7>::operator()(i);
  }  
};

template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8>
class Or_functors<S1,S2,S3,S4,S5,S6,S7,S8>: public Or_functors<S2,S3,S4,S5,S6,S7,S8>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8):s1_(s1),Or_functors<S2,S3,S4,S5,S6,S7,S8>(s2,s3,s4,s5,s6,s7,s8){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5,S6,S7,S8>::operator()(i);
  }  
};


template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8,class S9>
class Or_functors<S1,S2,S3,S4,S5,S6,S7,S8,S9>: public Or_functors<S2,S3,S4,S5,S6,S7,S8,S9>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8,const S9& s9):s1_(s1),Or_functors<S2,S3,S4,S5,S6,S7,S8,S9>(s2,s3,s4,s5,s6,s7,s8,s9){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5,S6,S7,S8,S9>::operator()(i);
  }  
};
/** \endcond */

/**
  * Class defining a function object returning the logical \b OR between the results of several given function objects.
  * The current implementation provides a version from two to ten parameters. If the code is compiled using
  * the c++0x standard, the number of parameters is not limited (use \c -std=c++0x with \c gcc).
  * @tparam Si are function objects with an operator() taking one argument and returning a Boolean (such as \ref atomsel).
  */
template <class S1,class S2,class S3,class S4,class S5,class S6,class S7,class S8,class S9,class S10>
class Or_functors: public Or_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2,const S3& s3,const S4& s4,const S5& s5,const S6& s6,const S7& s7,const S8& s8,const S9& s9,const S10& s10):s1_(s1),Or_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>(s2,s3,s4,s5,s6,s7,s8,s9,s10){}    
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S2,S3,S4,S5,S6,S7,S8,S9,S10>::operator()(i);
  }  
};


#else
/** \cond */
template <class ... S>
class And_functors;

template <class S1,class ... S>
class And_functors<S1,S... >: And_functors<S...>{
  S1 s1_;
public:
  And_functors(){}
  And_functors(const S1& s1,S... s):s1_(s1),And_functors<S...>(s...) {}
    
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && And_functors<S...>::operator()(i);
  }
};

template <class S1,class S2>
struct And_functors<S1,S2>: public S2{
  S1 s1_;
public:  
  And_functors(){}
  And_functors(const S1& s1,const S2& s2):s1_(s1),S2(s2){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) && S2::operator()(i);
  }  
};

template <class ... S>
struct Or_functors;

template <class S1,class ... S>
struct Or_functors<S1,S... >: Or_functors<S...>{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,S... s):S1(s1),Or_functors<S...>(s...) {}
    
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || Or_functors<S...>::operator()(i);
  }
};

template <class S1,class S2>
struct Or_functors<S1,S2>: S2{
  S1 s1_;
public:
  Or_functors(){}
  Or_functors(const S1& s1,const S2& s2):s1_(s1),S2(s2){}
  
  template <class T>
  bool operator()(const T& i) const{
    return s1_(i) || S2::operator()(i);
  }
};
/** \endcond */
#endif

} //namespace ESBTL

#endif //ESBTL_COMBINE_BOOLEAN_OPERATOR_H

