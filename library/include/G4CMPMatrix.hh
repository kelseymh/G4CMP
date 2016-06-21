#ifndef G4CMPMatrix_hh
#define G4CMPMatrix_hh 1

#include <ostream>
#include <vector>

// This class is mostly just a wrapper around std::vector to mimic being 2D
template <class T> class G4CMPMatrix {
public:
  G4CMPMatrix();
  G4CMPMatrix(size_t rows, size_t ncols, const T& val = T(0));
  G4CMPMatrix(const std::vector<T>& vec, size_t ncols);
  explicit G4CMPMatrix(const std::vector<T>& vec);
  G4CMPMatrix(const G4CMPMatrix<T> &rhs);

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

  // Assignment
  G4CMPMatrix& operator=(const G4CMPMatrix<T> &rhs);
  G4CMPMatrix& operator=(G4CMPMatrix<T>&& rhs);
  G4CMPMatrix& operator=(const std::vector<T>& rhs);
  G4CMPMatrix& operator=(const T& a);		// Filling/zeroing
  void clear() { operator=(T(0)); }

  // Compound Assignment
  G4CMPMatrix& operator+=(const G4CMPMatrix<T>& rhs);
  G4CMPMatrix& operator-=(const G4CMPMatrix<T>& rhs);
  G4CMPMatrix& operator+=(const T& rhs);
  G4CMPMatrix& operator-=(const T& rhs);
  G4CMPMatrix& operator*=(const T& rhs);
  G4CMPMatrix& operator/=(const T& rhs);

  // Comparison
  bool operator==(const G4CMPMatrix<T>& rhs) const;
  bool operator!=(const G4CMPMatrix<T>& rhs) const;
  bool operator< (const G4CMPMatrix<T>& rhs) const;
  bool operator<=(const G4CMPMatrix<T>& rhs) const;
  bool operator> (const G4CMPMatrix<T>& rhs) const;
  bool operator>=(const G4CMPMatrix<T>& rhs) const;

  // Access
  T& operator()(size_t i, size_t j);
  const T& operator()(size_t i, size_t j) const;
  T& operator()(size_t i);
  const T& operator()(size_t i) const;
  T& operator[](size_t i);
  const T& operator[](size_t i) const;

  // Access for use with pointers
  inline T& at(size_t i, size_t j) { return operator()(i,j); }
  inline const T& at(size_t i, size_t j) const { return operator()(i,j); }

  // Modify
  void push_back(const std::vector<T>& vec);
  void push_back(std::vector<T>&& vec);
  void vert_cat(const G4CMPMatrix<T>& rhs);
  void vert_cat(G4CMPMatrix<T>&& rhs);
  void horiz_cat(const G4CMPMatrix<T>& rhs);
  // No point in having a move version of horiz_cat

  // Information
  inline size_t size() const { return data.size(); }
  inline size_t columns() const { return ncols; }
  inline size_t rows() const { return (ncols ? data.size()/ncols : 0); }

private:
  std::vector<T> data;
  size_t ncols;
};

// Negate
template <class T>
G4CMPMatrix<T>& operator-(G4CMPMatrix<T>& lhs);

// Matrix-Matrix math:
template <class T>
G4CMPMatrix<T> operator+(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator+(G4CMPMatrix<T>&& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> operator-(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator-(G4CMPMatrix<T>&& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> operator*(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

// Matrix-Number math:
template <class T>
G4CMPMatrix<T> operator+(const G4CMPMatrix<T>& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator+(G4CMPMatrix<T>&& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator+(const T& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator+(const T& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> operator-(const G4CMPMatrix<T>& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator-(G4CMPMatrix<T>&& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator-(const T& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator-(const T& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> operator*(const G4CMPMatrix<T>& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator*(G4CMPMatrix<T>&& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator*(const T& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator*(const T& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> operator/(const G4CMPMatrix<T>& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator/(G4CMPMatrix<T>&& lhs, const T& rhs);
template <class T>
G4CMPMatrix<T> operator/(const T& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> operator/(const T& lhs, G4CMPMatrix<T>&& rhs);

// Matrix-Matrix operations:
template <class T>
G4CMPMatrix<T> vert_cat(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);
template <class T>
G4CMPMatrix<T> vert_cat(G4CMPMatrix<T>&& lhs, G4CMPMatrix<T>&& rhs);
template <class T>
G4CMPMatrix<T> horiz_cat(const G4CMPMatrix<T>& lhs, const G4CMPMatrix<T>& rhs);

// Printing
template <class T>
std::ostream& operator<<(std::ostream& out, const G4CMPMatrix<T>& rhs);

#include "G4CMPMatrix.icc"
#endif
