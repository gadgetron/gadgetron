/** \file vector_td_io.h
    \brief Basic iostream "communication" using the vector_td class
*/

#pragma once

#include "vector_td.h"

#include <cmath>
#include <iostream>
#include <algorithm>

namespace Gadgetron{

  template<class T, unsigned int D> ::std::ostream& operator<<(std::ostream& os, const vector_td<T,D>& vec) {
    os <<'[' ;
    for (int i = 0; i < D-1; i++) os << vec[i] << ", ";
    return os << vec[D-1] <<']';
  }

  template<class T, unsigned int D> std::istream& operator>>(std::istream& is, vector_td<T,D>& vec) {
    char tmp;
    is.get(tmp);
    if (tmp != '['){
      is.setstate(std::ios::failbit);
      return is;
    }

    for (int i = 0; i < D-1; i++){
      T val;
      tmp = ' ';
      is >> val;
      vec[i]=val;
      while (tmp == ' ') is.get(tmp);
      if (tmp != ','){
	is.setstate(std::ios::failbit);
	return is;
      }
    }
    tmp = ' ';
    is >> vec[D-1];
    while (tmp == ' ') is.get(tmp);
    if (tmp != ']'){
      is.setstate(std::ios::failbit);
      return is;
    }
    return is;
  }
}
