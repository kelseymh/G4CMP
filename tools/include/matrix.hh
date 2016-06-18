//  matrix.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef _matrix_h
#define _matrix_h

#include <cstddef>

/* Two-dimension matrix class implemented analogously to std::vector;
   does not include all of the required STL features */

namespace G4CMP {
template <class T>
class matrix {
private:
  size_t nn,mm;
  T **v;
 public:
  typedef T value_type; 			// make T available externally
  typedef T& reference;

  matrix() : nn(0), mm(0), v(NULL) {;}
  matrix(size_t n, size_t m);				// Zero based array
  matrix(size_t n, size_t m, const T &a);		// Initialize to constant
  matrix(size_t n, size_t m, const T *a);		// Initialize to array
  matrix(const matrix &rhs);			// Copy constructor

  inline void clear() { if (v) { delete[] v[0]; delete[] v; v=NULL; } }
  ~matrix() { clear(); }
  
  matrix& operator=(const matrix &rhs);		// Assignment
  matrix& operator=(const T *a);		// Copy array into matrix

  inline void resize(size_t n, size_t m);	// Reallocated, no filling
  inline void resize(size_t n, size_t m, const T& a);	// Reallocated, filling

  inline T& at(size_t i, size_t j);		// double subscripting
  inline const T& at(size_t i, size_t j) const;

  inline T* operator[](size_t i);	// double subscripting: pointer to row i
  inline const T* operator[](size_t i) const;

  inline size_t size() const { return nn*mm; }
  inline size_t nrows() const { return nn; }
  inline size_t ncols() const { return mm; }

private:
  inline void fill(const T &a=T(0));
};

}

#include "matrix.icc"

#endif	/* matrix.h */
