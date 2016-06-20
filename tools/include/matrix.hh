//  matrix.h
//  Created by Daniel Palken in 2014 for G4CMP

#ifndef _matrix_h
#define _matrix_h

#include <vector>
using std::vector;

/* Two-dimension matrix class implemented analogously to std::vector;
   does not include all of the required STL features */

namespace G4CMP {
template <class T>
class matrix {
private:
  size_t nn, mm;
  vector<vector<T> > v;

 public:
  typedef T value_type; 			// make T available externally
  typedef T& reference;

  matrix() : nn(0), mm(0), v(0) {;}
  matrix(size_t n, size_t m);			// Zero based array
  matrix(size_t n, size_t m, const T &a);	// Initialize to constant
  matrix(size_t n, size_t m, const T *a);	// Initialize to array
  matrix(const matrix &rhs);			// Copy constructor

  ~matrix() {;}
  
  matrix& operator=(const matrix &rhs);		// Assignment
  matrix& operator=(const T *a);		// Copy array into matrix

  void resize(size_t n, size_t m, const T& a=T(0));
  inline void clear();

  inline T& at(size_t i, size_t j);		// double subscripting
  inline const T& at(size_t i, size_t j) const;

  // Support double subscripting to individual elements
  inline vector<T>& operator[](size_t i) { return v[i]; }
  inline const vector<T>& operator[](size_t i) const { return v[i]; }

  // Appending functions
  void vert_cat(const matrix<T>& rhs);
  void vert_cat( matrix<T>&& rhs);
  void horiz_cat(const matrix<T>& rhs);

  inline size_t size() const { return nn*mm; }
  inline size_t rows() const { return nn; }
  inline size_t columns() const { return mm; }
};

}

#include "matrix.icc"

#endif	/* matrix.h */
