#pragma once
#include "Vector.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
using namespace std;

class TridiagonalMatrix {
private:
    int size;
    Vector higher;
    Vector main;
    Vector lower;

    void validate_indices(int i, int j) const {
        if (i < 1 || i > size || j < 1 || j > size ) {
                    throw out_of_range("Matrix indices out of range");
        }
    }

public:
    // Конструктор
    TridiagonalMatrix(int size_): size(size_), lower(size_), main(size_), higher(size_){
        lower[1] = 0;
        higher[size] = 0;
        if (size_ < 2) throw invalid_argument("Matrix size must be at least 2");
    } 

    TridiagonalMatrix(Vector higher_, Vector main_, Vector lower_): higher(higher_), lower(lower_), main(main_) {
        size = main.get_size();
        lower[1] = 0;
        higher[size] = 0;
    }
    TridiagonalMatrix(const TridiagonalMatrix& other) = default;

    int get_size() const { return size; }
    const Vector& get_higher() const { return higher;  }
    const Vector& get_main() const { return main; }
    const Vector& get_lower() const { return lower; }


    double get_upper(int i) const {
        if ((i < 1) || (i > higher.get_size()))
            throw std::out_of_range("Index out of range");
            return higher[i];
    }

    double get_main(int i) const {
        if ((i < 1) || (i > main.get_size()))
            throw std::out_of_range("Index out of range");
        return main[i];
    }

    double get_lower(int i) const {
        if ((i < 1) || (i > lower.get_size()))
            throw std::out_of_range("Index out of range");
        return lower[i];
    }


    void set_upper(int i, double value) {
        if ((i < 1) || (i > higher.get_size()))
            throw std::out_of_range("Index out of range");
        higher[i] = value;
    }

    void set_main(int i, double value) {
        if ((i < 1) || (i > main.get_size()))
            throw std::out_of_range("Index out of range");
        main[i] = value;
    }

    void set_lower(int i, double value){
        if ((i < 1) || (i > lower.get_size()))
            throw std::out_of_range("Index out of range");
        lower[i] = value;
    }

    // Операции с векторами
    Vector operator*(const Vector& vec) const {
        if (vec.get_size() != size) throw invalid_argument("Vector size must match matrix size");
    
        Vector result(size);
        int first = 1;
        int last = size;
    
        // первая строка
        result[first] = main[first] * vec[first];

        if (size > 1) result[first] += higher[first] * vec[first + 1];
    
        // промежуточные строки
        for (int i = first + 1; i < last; ++i) {
            result[i] = lower[i] * vec[i - 1] + main[i] * vec[i] + higher[i] * vec[i + 1];
        }
    
        // последняя строка
        result[last] = lower[last] * vec[last - 1] + main[last] * vec[last];
    
        return result;
    }

    // Алгоритм Томаса
    Vector Thomas_algorithm(const Vector& d, ofstream& file, bool verbose = false) const {
        if (d.get_size() != size) throw invalid_argument("Vector size must match matrix size");

        Vector L(size);
        Vector M(size + 1);
        Vector x(size);

        L[2] = higher[1] / main[1];
        M[2] = d[1] / main[1];

        if (verbose) {
            cout << "L[" << 2 << "] = " << L[2] << endl;
            cout << "M[" << 2 << "] = " << M[2] << endl;
            file << "L[" << 2 << "] = " << L[2] << endl;
            file << "M[" << 2 << "] = " << M[2] << endl;
        }

        for (int i = 2; i <= size - 1; ++i) {
            double denominator = main[i] - lower[i] * L[i];
            if (fabs(denominator) < 1e-12) {
                throw runtime_error("Zero denominator in Thomas algorithm");
            }
            L[i + 1] = higher[i] / denominator;
            M[i + 1] = (d[i] - lower[i] * M[i]) / denominator;

            if (verbose) {
                cout << "L[" << i + 1 << "] = " << L[i + 1] << endl;
                cout << "M[" << i + 1 << "] = " << M[i + 1] << endl;
                file << "L[" << i + 1 << "] = " << L[i + 1] << endl;
                file << "M[" << i + 1 << "] = " << M[i + 1] << endl;
            }
        }

        double denominator = main[size] - lower[size] * L[size];
        if (fabs(denominator) < 1e-12) {
            throw runtime_error("Zero denominator in Thomas algorithm");
        }
        M[size + 1] = (d[size] - lower[size] * M[size]) / denominator;

        if (verbose) {
            cout << "M[" << size + 1 << "] = " << M[size + 1] << endl;
            file << "M[" << size + 1 << "] = " << M[size + 1] << endl;
        }

        x[size] = M[size + 1];

        if (verbose) {
            cout << "x[" << size << "] = " << x[size] << endl;
            file << "x[" << size << "] = " << x[size] << endl;
        }
        for (int i = size - 1; i >= 1; --i) {
            x[i] = M[i + 1] - L[i + 1] * x[i + 1];
            if (verbose) {
                cout << "x[" << i << "] = " << x[i] << endl;
                file << "x[" << i << "] = " << x[i] << endl;
            }
        }

        return x;
    }


    // Тестирование
    void test_Thomas_algorithm(ofstream& file) {
        Vector l(5), m(5), h(5);
        m[1] = 4;
        m[2] = 6;
        m[3] = 5;
        m[4] = 11;
        m[5] = 8;

        l[1] = 0;
        l[2] = 2;
        l[3] = 5;
        l[4] = 2;
        l[5] = 5;

        h[1] = 2;
        h[2] = 1;
        h[3] = 2;
        h[4] = 2;
        h[5] = 0;

        Vector d(5);
        d[1] = 2;
        d[2] = -2;
        d[3] = 3;
        d[4] = -5;
        d[5] = 3;

        cout << "Matrix: " << endl;
        TridiagonalMatrix matrix(h, m, l);
        matrix.print_vectors();
        cout << endl;

        cout << "Vector d: " << endl;
        d.print();

        cout << "Coefficients of Thomas algorithm: " << endl;
        Vector solution = matrix.Thomas_algorithm(d, file, true);

        cout << "Solution vector X: " << endl;
        solution.print();
        cout << endl;

        Vector expected(5);

        expected[1] = 1;
        expected[2] = -1;
        expected[3] = 2;
        expected[4] = -1;
        expected[5] = 1;

        cout << "Expected solution: " << endl;
        expected.print();
        cout << endl;

        cout << "Error: " << endl;
        double max_error = 0.0;
        for (int i = 1; i <= 5; i++) {
            double error = fabs(solution[i] - expected[i]);
            cout << "x[" << i << "]: error = " << error << endl;
            if (error > max_error) max_error = error;
        }
        cout << "Max error: " << max_error << endl;


        cout << "Verification: AX = d" << endl;
        Vector Ax = matrix * solution;
        cout << "A*x = ";
        Ax.print();
        cout << endl;

        cout << "Residual (A*x - d): " << endl;
        Vector residual = Ax - d;
        cout << "Residual vector: ";
        residual.print();
        cout << endl;
        cout << "Residual norm (Ax-d): " << residual.norm() << endl;

    }

    // Заполнение и вывод
    void fill_random(int min_val, int max_val) {
        main.fill_random(min_val, max_val);
        lower.fill_random(min_val, max_val);
        higher.fill_random(min_val, max_val);
    }
    //void print();
    void print_vectors() {
        cout << "Main diagonal: ";
        main.print();
    
        cout << "Lower diagonal: ";
        lower.print();
    
        cout << "Higher diagonal: ";
        higher.print();
    };

    void print_to_file(ofstream& file) const {
        file << "Main diagonal: ";
        main.print_to_file(file);

        file << "Lower diagonal: ";
        lower.print_to_file(file);

        file << "Higher diagonal: ";
        higher.print_to_file(file);
    }
};