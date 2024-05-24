#include <iostream>
#include <math.h>
#include <vector>
using namespace std;
/* finding next square of 2*/
int nextPowerOf2(int k) { return pow(2, int(ceil(log2(k)))); }

//! addition and subtraction
void add(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C,
         int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
}
void sub(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C,
         int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      C[i][j] = A[i][j] - B[i][j];
    }
  }
}

void Strassen_algorithm(vector<vector<int>> &A, vector<vector<int>> &B,
                        vector<vector<int>> &C, int size) {
  if (size == 1) {
    C[0][0] = A[0][0] * B[0][0];
    return;
  } else {
    int newSize = size / 2;
    vector<vector<int>> A11(newSize, vector<int>(newSize)),
        A12(newSize, vector<int>(newSize)), A21(newSize, vector<int>(newSize)),
        A22(newSize, vector<int>(newSize)), B11(newSize, vector<int>(newSize)),
        B12(newSize, vector<int>(newSize)), B21(newSize, vector<int>(newSize)),
        B22(newSize, vector<int>(newSize)), c11(newSize, vector<int>(newSize)),
        c12(newSize, vector<int>(newSize)), c21(newSize, vector<int>(newSize)),
        c22(newSize, vector<int>(newSize)), R(newSize, vector<int>(newSize)),
        T(newSize, vector<int>(newSize)), Q(newSize, vector<int>(newSize)),
        S(newSize, vector<int>(newSize)), P(newSize, vector<int>(newSize)),
        V(newSize, vector<int>(newSize)), U(newSize, vector<int>(newSize)),
        fResult(newSize, vector<int>(newSize)),
        sResult(newSize, vector<int>(newSize));

    int i, j;

    for (i = 0; i < newSize; i++) {
      for (j = 0; j < newSize; j++) {
        A11[i][j] = A[i][j];
        A12[i][j] = A[i][j + newSize];
        A21[i][j] = A[i + newSize][j];
        A22[i][j] = A[i + newSize][j + newSize];

        B11[i][j] = B[i][j];
        B12[i][j] = B[i][j + newSize];
        B21[i][j] = B[i + newSize][j];
        B22[i][j] = B[i + newSize][j + newSize];
      }
    }
    // P = (A11 + A22) * (B11 + B22)
    add(A11, A22, fResult, newSize);
    add(B11, B22, sResult, newSize);
    Strassen_algorithm(fResult, sResult, P, newSize);

    // Q = (A21 + A22) * B11
    add(A21, A22, fResult, newSize);
    Strassen_algorithm(fResult, B11, Q, newSize);

    // R = A11 * (B12 - B22)
    sub(B12, B22, sResult, newSize);
    Strassen_algorithm(A11, sResult, R, newSize);

    // S = A22 * (B21 - B11)
    sub(B21, B11, sResult, newSize);
    Strassen_algorithm(A22, sResult, S, newSize);

    // T = (A11 + A12) * B22
    add(A11, A12, fResult, newSize);
    Strassen_algorithm(fResult, B22, T, newSize);

    // U = (A21 - A11) * (B11 + B12)
    sub(A21, A11, fResult, newSize);
    add(B11, B12, sResult, newSize);
    Strassen_algorithm(fResult, sResult, U, newSize);

    // V = (A12 - A22) * (B21 + B22)
    sub(A12, A22, fResult, newSize);
    add(B21, B22, sResult, newSize);
    Strassen_algorithm(fResult, sResult, V, newSize);

    // c11 = P + S - T + V
    add(P, S, fResult, newSize);
    sub(fResult, T, sResult, newSize);
    add(sResult, V, c11, newSize);

    // c12 = R + T
    add(R, T, c12, newSize);

    // c21 = Q + S
    add(Q, S, c21, newSize);

    // c22 = P + R - Q + U
    add(P, R, fResult, newSize);
    sub(fResult, Q, sResult, newSize);
    add(sResult, U, c22, newSize);
    // Grouping the results obtained in a single matrix:
    for (i = 0; i < newSize; i++) {
      for (j = 0; j < newSize; j++) {
        C[i][j] = c11[i][j];
        C[i][j + newSize] = c12[i][j];
        C[i + newSize][j] = c21[i][j];
        C[i + newSize][j + newSize] = c22[i][j];
      }
    }
  }
}

void multiplicar(vector<vector<int>> &A, vector<vector<int>> &B,
                 vector<vector<int>> &C) {
  int n = A.size();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = 0;
      for (int k = 0; k < n; k++) {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

int main() {

  // vector<vector<int>> A = { {17, 1, 25, 29}, {1, 22, 26, 1}, {19, 23, 1, 31},
  // {20, 24, 28, 32} }; vector<vector<int>> B = { {2, 5, 9, 13}, {3, 6, 1, 14},
  // {4, 1, 11, 15}, {1, 8, 12, 1} };

  vector<vector<int>> A = {
      {2, 3, 4, 1}, {5, 6, 1, 8}, {9, 1, 11, 12}, {13, 14, 15, 1}};
  vector<vector<int>> B = {
      {17, 1, 19, 20}, {1, 22, 23, 24}, {25, 26, 1, 28}, {29, 1, 31, 32}};
  int size = A.size();
  vector<vector<int>> C(size, vector<int>(size));

  // Llamar a la función Strassen_algorithm
  Strassen_algorithm(A, B, C, size);

  vector<vector<int>> D = {
      {2, 3, 4, 1}, {5, 6, 1, 8}, {9, 1, 11, 12}, {13, 14, 15, 1}};
  vector<vector<int>> E = {
      {17, 1, 19, 20}, {1, 22, 23, 24}, {25, 26, 1, 28}, {29, 1, 31, 32}};

  vector<vector<int>> F(A.size(), vector<int>(A.size()));
  multiplicar(D, E, F);

  cout << "Strassen" << endl;

  // Imprimir el resultado
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      cout << C[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl << "Normalito de toda la vida" << endl;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      cout << F[i][j] << " ";
    }
    cout << endl;
  }
  return 0;
}

// strassen.cpp Mostrando strassen.cpp.
