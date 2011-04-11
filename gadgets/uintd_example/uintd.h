#ifndef __UINTD_H
#define __UINTD_H

template <int LENGTH> struct uintd
{
  unsigned int& operator[] (int i);
  unsigned int d[LENGTH];
};

//template <int ASIZE> void nidx_to_co(unsigned int& idx, uintd<ASIZE>& dims, uintd<ASIZE>& co);
//template <int ASIZE> unsigned int nco_to_idx(uintd<ASIZE>& dims, uintd<ASIZE>& co);
//template <int ASIZE> unsigned int nco_to_idx(uintd<ASIZE>& dims, uintd<ASIZE>& co, uintd<ASIZE>& order);

#endif //__UINTD_H
