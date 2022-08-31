#include "G4CMPBlockData.hh"
#include <algorithm>
#include <iostream>
#include <cstdlib>

using std::cout; 
using std::endl;

int main() {
  int a = rand()%5;
  int b = a;
  do {
    b = rand()%5;
  } while (b == a); // I don't want to test on square matrices

  G4CMPBlockData<int> mat1(a, b);
  cout << "Constructing first matrix. Should have " << a << " rows and "
       << b << " columns: \n G4CMPBlockData::rows() = " << mat1.rows()
       << "\n G4CMPBlockData::columns() = " << mat1.columns() << endl;

  cout << "Filling matrix with random data: " << endl;
  std::generate(mat1.begin(), mat1.end(), 
               []{int c = rand()%10; cout << c << endl; return c;});

  cout << "Display matrix values in row major order: " << endl;
  for (int v : mat1) cout << v << endl;

  cout << "Print matrix: " << endl;
  cout << mat1 << endl;
  
  cout << "Make and print another matrix to work with: " << endl;
  G4CMPBlockData<int> mat2(a, b);
  std::generate(mat2.begin(), mat2.end(), []{return rand()%10;});
  cout << mat2 << endl;

  cout << "Test mat1 + mat2:" << endl;
  cout << mat1 + mat2 << endl;
  
  cout << "Test mat1 - mat2:" << endl;
  cout << mat1 - mat2 << endl;

  cout << "Make and print another matrix to work with: " << endl;
  G4CMPBlockData<int> mat3(b, a);
  std::generate(mat3.begin(), mat3.end(), []{return rand()%10;});
  cout << mat3 << endl;

  cout << "Test mat1 * mat3:" << endl;
  cout << mat1 * mat3 << endl;

  cout << "Test mat1 + 2:" << endl;
  cout << mat1 + 2 << endl;

  cout << "Test 2 + mat1:" << endl;
  cout << 2 + mat1 << endl;

  cout << "Test mat1 - 2:" << endl;
  cout << mat1 - 2 << endl;

  cout << "Test 2 - mat1:" << endl;
  cout << 2 - mat1 << endl;

  cout << "Test mat1 * 2:" << endl;
  cout << mat1 * 2 << endl;

  cout << "Test 2 * mat1:" << endl;
  cout << 2 * mat1 << endl;

  cout << "Test vert_cat(mat1, mat2):" << endl;
  cout << vert_cat(mat1, mat2) << endl;

  cout << "Test horiz_cat(mat1, mat2):" << endl;
  cout << horiz_cat(mat1, mat2) << endl;

  return 0;
}
