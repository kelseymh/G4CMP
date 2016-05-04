#ifndef G4CMPMatrix_hh
#define G4CMPMatrix_hh 1

#include <vector>
#include <cstdlib>

// This class is mostly just a wrapper around std::vector to mimic being 2D
template <class T> class G4CMPMatrix {
public:
  G4CMPMatrix();
  G4CMPMatrix(size_t rows, size_t columns, const T& val = T(0));
  G4CMPMatrix(const std::vector<T>& vec, size_t columns);
  explicit G4CMPMatrix(const std::vector<T>& vec);

  using iterator = std::vector<T>::iterator;
  iterator begin() const;
  iterator end() const;

  // Assignment
  G4CMPMatrix<T>& operator=(const G4CMPMatrix<T>& rhs);
  G4CMPMatrix<T>& operator=(G4CMPMatrix<T>&& rhs);

  // Compound Assignment
  G4CMPMatrix<T>& operator+=(const G4CMPMatrix<T>& rhs);
  G4CMPMatrix<T>& operator-=(const G4CMPMatrix<T>& rhs);
  G4CMPMatrix<T>& operator+=(const T& rhs);
  G4CMPMatrix<T>& operator-=(const T& rhs);
  G4CMPMatrix<T>& operator*=(const T& rhs);
  G4CMPMatrix<T>& operator/=(const T& rhs);

  // Comparison
  friend bool operator==(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator!=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator<(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator<=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator>(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
  friend bool operator>=(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

  // Access
  T& operator()(size_t i, size_t j) const;
  const T& operator()(size_t i, size_t j) const;
  T& operator()(size_t i) const;
  const T& operator()(size_t i) const;
  T& operator[](size_t i) const;
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
  size_t columns;
};

// Negate
template <class T>
G4CMPMatrix<T>& operator-(G4CMPMatrix<T>& lhs);

// Matrix-Matrix math:
template <class T>
G4CMPMatrix<T>& operator+(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
const G4CMPMatrix<T>& operator+(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T>& operator-(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
const G4CMPMatrix<T>& operator-(G4CMPMatrix<T> lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T>& operator*(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
const G4CMPMatrix<T>& operator*(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

// Matrix-Number math:
template <class T>
G4CMPMatrix<T>& operator+(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T>& operator+(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
const G4CMPMatrix<T>& operator+(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
const G4CMPMatrix<T>& operator+(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T>& operator-(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T>& operator-(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
const G4CMPMatrix<T>& operator-(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
const G4CMPMatrix<T>& operator-(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T>& operator*(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T>& operator*(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
const G4CMPMatrix<T>& operator*(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
const G4CMPMatrix<T>& operator*(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
G4CMPMatrix<T>& operator/(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
G4CMPMatrix<T>& operator/(const T& lhs, G4CMPMatrix<T> rhs);
template <class T>
const G4CMPMatrix<T>& operator/(G4CMPMatrix<T> lhs, const T& rhs);
template <class T>
const G4CMPMatrix<T>& operator/(const T& lhs, G4CMPMatrix<T> rhs);

#include "G4CMPMatrix.icc"
#endif
