//  matrix.h
//  Created by Daniel Palken in 2014 for G4CMP
//
//  20170525  M. Kelsey -- Add move semantics for copy ctor and assignment
//  20170728  M. Kelsey -- Replace n,m with nrow,ncol to avoid conflicts w/m,mm

#ifndef G4CMPMatrix_h
#define G4CMPMatrix_h

#include <cstddef>
#include <vector>
using std::vector;

/* Two-dimension matrix class implemented analogously to std::vector;
   does not include all of the required STL features */

namespace G4CMP {
template <class T>
class matrix {
private:
  size_t nrow, ncol;
  vector<vector<T> > v;

public:
  typedef T value_type; 			// make T available externally
  typedef T& reference;

  matrix() : nrow(0), ncol(0), v(0) {;}
  matrix(size_t nr, size_t nc);			// Zero based array
  matrix(size_t nr, size_t nc, const T &a);	// Initialize to constant
  matrix(size_t nr, size_t nc, const T *a);	// Initialize to array
  matrix(const matrix &rhs);			// Copy constructor
  matrix(matrix &&rhs);				// Move constructor

  ~matrix() {;}
  
  matrix& operator=(const matrix &rhs);		// Assignment
  matrix& operator=(matrix &&rhs);		// Move assignment
  matrix& operator=(const T *a);		// Copy array into matrix

  void resize(size_t nr, size_t nc, const T& a=T(0));
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

  inline size_t size() const { return nrow*ncol; }
  inline size_t rows() const { return nrow; }
  inline size_t columns() const { return ncol; }
};

}

#include "G4CMPMatrix.icc"

#endif	/* G4CMPMatrix_h */
