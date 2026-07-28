#include "Vector.h"
#include "TridiagonalMatrix.h"
#include "SpecialMatrix.h"
#include "ComputationalExperiment.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
using namespace std;

int main() {
    ofstream file("Lab_1.txt");
    cout << "=== TEST TRIDIAGONAL MATRIX ===" << endl;
    file << "=== TEST TRIDIAGONAL MATRIX ===" << endl;

    Vector A(4);
    A[1] = 0;
    A[2] = 2;
    A[3] = 3;
    A[4] = 4;
    A.print();

    Vector B(4);
    B[1] = 3;
    B[2] = 4;
    B[3] = 5;
    B[4] = 6;

    Vector C(4);
    C[1] = 5;
    C[2] = 2;
    C[3] = 3;
    C[4] = 0;

    cout << "matr: " << endl;
    file << "matr: " << endl;
    TridiagonalMatrix matr(C, B, A);
    matr.print_vectors();
    matr.print_to_file(file);

    Vector X(4);
    X[1] = 10;
    X[2] = 1;
    X[3] = -1;
    X[4] = 0;  

    cout << "X: " << endl;
    file << "X: " << endl;
    X.print();
    X.print_to_file(file);

    Vector d = matr * X;
    cout << "d: " << endl;
    file << "d: " << endl;
    d.print();
    d.print_to_file(file);
    
    cout << "L, M coefficients: " << endl;
    file << "L, M coefficients: " << endl;

    Vector X_computed = matr.Thomas_algorithm(d, file, true);
    cout << "X_computed: " << endl;
    file << "X_computed: " << endl;
    X_computed.print();
    X_computed.print_to_file(file);

    matr.test_Thomas_algorithm(file);

    cout << "-----------------------------------------------------------------------------" << endl;
    file << "-----------------------------------------------------------------------------" << endl;
    cout << "=== TEST SPECIAL MATRIX ===" << endl;
    file << "=== TEST SPECIAL MATRIX ===" << endl;

    Vector R(4);
    R[1] = 1;
    R[2] = 2;
    R[3] = 3;
    R[4] = 4;

    SpecialMatrix my_matrix(C, B, A, R);
    cout << "my_matrix: " << endl;
    file << "my_matrix: " << endl;
    my_matrix.print_vectors();
    my_matrix.print_to_file(file);

    Vector my_d = my_matrix * X;
    cout << "my_d: " << endl;
    file << "my_d: " << endl;
    my_d.print();
    my_d.print_to_file(file);

    cout << "L, M coefficients: " << endl;
    file << "L, M coefficients: " << endl;


    Vector my_X_computed = my_matrix.Thomas_algorithm(my_d, file, true);
    cout << "my_X_computed: " << endl;
    file << "my_X_computed: " << endl;
    my_X_computed.print();
    my_X_computed.print_to_file(file);

    Vector M(4);
    M[1] = 0;
    M[2] = 3;
    M[3] = 2;
    M[4] = 4;
    M.print();

    Vector N(4);
    N[1] = 5;
    N[2] = 10;
    N[3] = 3;
    N[4] = 1;

    Vector K(4);
    K[1] = 2;
    K[2] = 6;
    K[3] = 5;
    K[4] = 0;

    Vector L(4);
    L[1] = 1;
    L[2] = 3;
    L[3] = 4;
    L[4] = 1;

    SpecialMatrix my_matrix_(K, N, M, L);
    cout << "my_matrix_: " << endl;
    file << "my_matrix_: " << endl;
    my_matrix_.print_vectors();
    my_matrix_.print_to_file(file);

    Vector my_d_ = my_matrix_ * X;
    cout << "my_d_: " << endl;
    file << "my_d_: " << endl;
    my_d_.print();
    my_d_.print_to_file(file);

    cout << "L, M coefficients: " << endl;
    file << "L, M coefficients: " << endl;


    Vector my_X_computed_ = my_matrix_.Thomas_algorithm(my_d_, file, true);
    cout << "my_X_computed_: " << endl;
    file << "my_X_computed_: " << endl;
    my_X_computed_.print();
    my_X_computed_.print_to_file(file);

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "=== TEST LARGE SPECIAL MATRIX ===" << endl;
    my_matrix.test_Thomas_algorithm();

    ComputationalExperiment experiment;
    experiment.run_experiment();
}