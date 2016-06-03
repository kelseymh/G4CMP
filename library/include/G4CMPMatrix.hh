#ifndef G4CMPMatrix_hh
#define G4CMPMatrix_hh 1

#include <vector>
#include <cstdlib>
#include <ostream>

// Forward declarations needed for friend operators
template <class T> class G4CMPMatrix;

template <class T>
bool operator==(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
bool operator!=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
bool operator<(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
bool operator<=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
bool operator>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
bool operator>=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

// This class is mostly just a wrapper around std::vector to mimic being 2D
template <class T> class G4CMPMatrix {
public:
  G4CMPMatrix();
  G4CMPMatrix(size_t rows, size_t ncols, const T& val = T(0));
  G4CMPMatrix(const std::vector<T>& vec, size_t ncols);
  explicit G4CMPMatrix(const std::vector<T>& vec);

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

  // Compound Assignment
  G4CMPMatrix& operator+=(const G4CMPMatrix& rhs);
  G4CMPMatrix& operator-=(const G4CMPMatrix& rhs);
  G4CMPMatrix& operator+=(const T& rhs);
  G4CMPMatrix& operator-=(const T& rhs);
  G4CMPMatrix& operator*=(const T& rhs);
  G4CMPMatrix& operator/=(const T& rhs);

  // Comparison
  friend bool operator== <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator!= <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator<  <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator<= <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator>  <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator>= <T>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

  // Access
  T& operator()(size_t i, size_t j);
  const T& operator()(size_t i, size_t j) const;
  T& operator()(size_t i);
  const T& operator()(size_t i) const;
  T& operator[](size_t i);
  const T& operator[](size_t i) const;

  // Modify
  void push_back(const std::vector<T>& vec);
  void push_back(std::vector<T>&& vec);
  void concat(const G4CMPMatrix<T>& rhs);
  void concat(G4CMPMatrix<T>&& rhs);

  // Information
  size_t size() const;
  size_t columns() const;
  size_t rows() const;

private:
  std::vector<T> data;
  size_t ncols;
};

// Negate
template <class T>
G4CMPMatrix<T>& operator-(G4CMPMatrix<T>& lhs);

// Matrix-Matrix math:
template <class T>
G4CMPMatrix<T> operator+(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator-(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator*(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

// Matrix-Number math:
template <class T>
G4CMPMatrix<T> operator+(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator+(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T> operator-(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator-(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T> operator*(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator*(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T> operator/(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator/(const T& lhs, G4CMPMatrix<T> rhs);

// Printing
template <class T>
std::ostream& operator<<(std::ostream& out, const G4CMPMatrix<T>& rhs);

#include "G4CMPMatrix.icc"
#endif
