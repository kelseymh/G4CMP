#ifndef G4CMPBlockData_hh
#define G4CMPBlockData_hh 1

#include <ostream>
#include <vector>

// This class is mostly just a wrapper around std::vector to mimic being 2D
template <class T> class G4CMPBlockData {
public:
  G4CMPBlockData();
  G4CMPBlockData(size_t nrows, size_t ncols, const T& val = T(0));
  G4CMPBlockData(const std::vector<T>& vec, size_t ncols);
  explicit G4CMPBlockData(const std::vector<T>& vec);
  G4CMPBlockData(const G4CMPBlockData<T> &rhs);

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;

  // Assignment
  G4CMPBlockData& operator=(const G4CMPBlockData<T> &rhs);
  G4CMPBlockData& operator=(G4CMPBlockData<T>&& rhs);
  G4CMPBlockData& operator=(const std::vector<T>& rhs);
  G4CMPBlockData& operator=(const T& a);		// Filling/zeroing
  void clear() { operator=(T(0)); }

  // Compound Assignment
  G4CMPBlockData& operator+=(const G4CMPBlockData<T>& rhs);
  G4CMPBlockData& operator-=(const G4CMPBlockData<T>& rhs);
  G4CMPBlockData& operator+=(const T& rhs);
  G4CMPBlockData& operator-=(const T& rhs);
  G4CMPBlockData& operator*=(const T& rhs);
  G4CMPBlockData& operator/=(const T& rhs);

  // Comparison
  bool operator==(const G4CMPBlockData<T>& rhs) const;
  bool operator!=(const G4CMPBlockData<T>& rhs) const;
  bool operator< (const G4CMPBlockData<T>& rhs) const;
  bool operator<=(const G4CMPBlockData<T>& rhs) const;
  bool operator> (const G4CMPBlockData<T>& rhs) const;
  bool operator>=(const G4CMPBlockData<T>& rhs) const;

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
  void vert_cat(const G4CMPBlockData<T>& rhs);
  void vert_cat(G4CMPBlockData<T>&& rhs);
  void horiz_cat(const G4CMPBlockData<T>& rhs);
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
G4CMPBlockData<T>& operator-(G4CMPBlockData<T>& lhs);

// Matrix-Matrix math:
template <class T>
G4CMPBlockData<T> operator+(const G4CMPBlockData<T>& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator+(G4CMPBlockData<T>&& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> operator-(const G4CMPBlockData<T>& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator-(G4CMPBlockData<T>&& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> operator*(const G4CMPBlockData<T>& lhs, const G4CMPBlockData<T>& rhs);

// Matrix-Number math:
template <class T>
G4CMPBlockData<T> operator+(const G4CMPBlockData<T>& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator+(G4CMPBlockData<T>&& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator+(const T& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator+(const T& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> operator-(const G4CMPBlockData<T>& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator-(G4CMPBlockData<T>&& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator-(const T& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator-(const T& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> operator*(const G4CMPBlockData<T>& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator*(G4CMPBlockData<T>&& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator*(const T& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator*(const T& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> operator/(const G4CMPBlockData<T>& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator/(G4CMPBlockData<T>&& lhs, const T& rhs);
template <class T>
G4CMPBlockData<T> operator/(const T& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> operator/(const T& lhs, G4CMPBlockData<T>&& rhs);

// Matrix-Matrix operations:
template <class T>
G4CMPBlockData<T> vert_cat(const G4CMPBlockData<T>& lhs, const G4CMPBlockData<T>& rhs);
template <class T>
G4CMPBlockData<T> vert_cat(G4CMPBlockData<T>&& lhs, G4CMPBlockData<T>&& rhs);
template <class T>
G4CMPBlockData<T> horiz_cat(const G4CMPBlockData<T>& lhs, const G4CMPBlockData<T>& rhs);

// Printing
template <class T>
std::ostream& operator<<(std::ostream& out, const G4CMPBlockData<T>& rhs);

#include "G4CMPBlockData.icc"
#endif
