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
class SpecialMatrix {
private:
    int size;
    Vector higher;
    Vector main;
    Vector lower;
    Vector last_row;

    void validate_indices(int i, int j) const {
        	if (i < 1 || i> size || j < 1 || j > size) {
    		throw out_of_range("Matrix indices out of range");
    	}
    }

public:
    // Конструктор
    SpecialMatrix(int size_) :  size(size_), lower(size_), main(size_), higher(size_), last_row(size_) {
        lower[1] = 0;
        higher[size] = 0;
        last_row[size - 1] = lower[size];
        last_row[size] = main[size];
	    if (size_ < 2) throw invalid_argument("Matrix size must be at least 2");
    };

    SpecialMatrix(Vector higher_, Vector main_, Vector lower_, Vector last_row_) : higher(higher_), lower(lower_), main(main_), last_row(last_row_) {
        size = main.get_size();
        lower[1] = 0;
        higher[size] = 0;
        last_row[size - 1] = lower[size];
        last_row[size] = main[size];
    }

    SpecialMatrix(const SpecialMatrix& other) = default;

    // Основные методы доступа
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

    double get_last_row(int i) const {
        if ((i < 1) || (i > last_row.get_size()))
            throw std::out_of_range("Index out of range");
        return last_row[i];
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

    void set_lower(int i, double value) {
        if ((i < 1) || (i > lower.get_size()))
            throw std::out_of_range("Index out of range");
        lower[i] = value;
    }

    void set_last_row(int i, double value) {
        if ((i < 1) || (i > last_row.get_size()))
            throw std::out_of_range("Index out of range");
        last_row[i] = value;
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
        for (int i = first; i <= size; i++) {
            result[last] += last_row[i]*vec[i];
        }
        return result;
    }
    

    // Алгоритм Томаса
    Vector Thomas_algorithm(const Vector& d, ofstream& file, bool verbose = false) const {
        if (d.get_size() != size) throw invalid_argument("Vector size must match matrix size");

        Vector d_reverse(d.get_size());
        for (int i = 1; i <= d.get_size(); i++) {
            d_reverse[d.get_size() - i + 1] = d[i];
        }

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

        double com_L = 1;
        double sum_denom = 0;
        int sign = 1;
        for (int i = size; i >= 1; i--) {
            sum_denom += sign * last_row[i] * com_L;
            com_L *= L[i];
            sign *= (-1);
        }
        double sum_numer = 0;
        double com_form = M[size];
        for (int i = size - 1; i >= 1; i--) {
            sum_numer += last_row[i] * com_form;
            com_form = M[i] - L[i] * com_form;
        }
        if (fabs(sum_denom) < 1e-12) {
            throw runtime_error("Zero denominator in Thomas algorithm");
        }
        M[size + 1] = (d[size] - sum_numer) / sum_denom;
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
    };

    // Тестирование
    void test_Thomas_algorithm() {
        ofstream file("Special_large_matrix.txt");

        int size = 12;

        cout << "=== LARGE SPECIAL MATRIX TEST (size " << size << ") ===" << endl;
        file << "=== LARGE SPECIAL MATRIX TEST (size " << size << ") ===" << endl;

        cout << "matrix: " << endl;
        file << "matrix: " << endl;
        SpecialMatrix matrix(size);
        matrix.fill_random(-10, 10);
        matrix.print_vectors();
        matrix.print_to_file(file);
        cout << endl;
        file << endl;

        matrix.check_compatibility_conditions(file);

        cout << "X: " << endl;
        file << "X: " << endl;
        Vector exact_solution(size);  
        exact_solution.fill_random(-5, 5);
        exact_solution.print();
        exact_solution.print_to_file(file);

        Vector d = matrix * exact_solution; 
        matrix.check_equivalence_step_by_step(exact_solution, d, file);

        cout << "Coefficients of Thomas algorithm L, M: " << endl;
        file << "Coefficients of Thomas algorithm L, M: " << endl;
        Vector X_computed = matrix.Thomas_algorithm(d, file, true);

        cout << "X_computed: " << endl;
        file << "X_computed: " << endl;
        X_computed.print();
        X_computed.print_to_file(file);
        cout << endl;

        matrix.check_solution_accuracy_verbose(X_computed, d, file); 

        double error_norm = (exact_solution - X_computed).norm(); 
        cout << "error_norm: " << error_norm << endl;
        file << "error_norm: " << error_norm << endl;
    }
    // Заполнение и вывод
    void fill_random(int min_val, int max_val) {
        main.fill_random(min_val, max_val);
        lower.fill_random(min_val, max_val);
        higher.fill_random(min_val, max_val);
        last_row.fill_random(min_val, max_val);
    };
    void print_vectors() {
        cout << "Main diagonal: ";
        main.print();

        cout << "Lower diagonal: ";
        lower.print();

        cout << "Higher diagonal: ";
        higher.print();

        cout << "Last row: ";
        last_row.print();
    }
    void print_to_file(ofstream& file) const {
        file << "Main diagonal: ";
        main.print_to_file(file);

        file << "Lower diagonal: ";
        lower.print_to_file(file);

        file << "Higher diagonal: ";
        higher.print_to_file(file);

        file << "Last row: ";
        last_row.print_to_file(file);
    }

    const Vector& get_main() const { return main; }
    const Vector& get_lower() const { return lower; }
    const Vector& get_higher() const { return higher; }
    const Vector& get_last_row() const { return last_row; }
    int get_size() const { return size; }

    double check_solution_accuracy(const Vector& x, const Vector& d) const {
        Vector Ax = operator*(x);
        Vector residual = Ax - d;
        return residual.norm();
    }

    double check_solution_accuracy_verbose(const Vector& x, const Vector& d, ofstream& file) const {
        Vector Ax = operator*(x);
        Vector residual = Ax - d;
        double residual_norm = residual.norm();

        cout << "Solution verification: " << endl;
        file << "Solution verification: " << endl;
        cout << "A*x =  ";
        Ax.print();
        Ax.print_to_file(file);

        cout << "d = ";
        file << "d = ";
        d.print();
        d.print_to_file(file);

        cout << "Residual A*x - d = ";
        file << "Residual A*x - d = ";
        residual.print();
        residual.print_to_file(file);

        cout << "Residual norm: " << residual_norm << endl;
        file << "Residual norm: " << residual_norm << endl;

        return residual_norm;
    }
    void check_compatibility_conditions(ofstream& file) const {

        bool condition_1 = (lower[size] - last_row[size - 1] < 1e-12);
        bool condition_2 = (main[size] - last_row[size] < 1e-12);

        if (condition_1 && condition_2) {
            cout << "All conditions are satisfied" << endl;
            file << "All conditions are satisfied" << endl;
        }
        else {
            cout << "Conditions aren't satisfied" << endl;
            file << "Conditions aren't satisfied" << endl;
        }
    }
    void check_equivalence_step_by_step(const Vector& exact_x, const Vector& d, ofstream& file) const {
        cout << "=== Equivalence Check ===" << endl;
        file << "=== Equivalence Check ===" << endl;

        // 1. Проверка исходной системы
        Vector Ax_original = operator*(exact_x);
        double original_error = (Ax_original - d).norm();

        cout << "1. Original system check:" << endl;
        cout << "   ||A * exact_x - d|| = " << original_error << endl;
        file << "1. Original system check:" << endl;
        file << "   ||A * exact_x - d|| = " << original_error << endl;

        // 2. Проверка после решения
        Vector computed_x = Thomas_algorithm(d, file, false); // без verbose
        Vector Ax_computed = operator*(computed_x);
        double computed_error = (Ax_computed - d).norm();

        cout << "2. After Thomas algorithm:" << endl;
        cout << "   ||A * computed_x - d|| = " << computed_error << endl;
        file << "2. After Thomas algorithm:" << endl;
        file << "   ||A * computed_x - d|| = " << computed_error << endl;

        // 3. Проверка эквивалентности
        double transformation_equivalence = fabs(original_error - computed_error);
        cout << "3. Transformation equivalence:" << endl;
        cout << "   |Original error - Computed error| = " << transformation_equivalence << endl;
        file << "3. Transformation equivalence:" << endl;
        file << "   |Original error - Computed error| = " << transformation_equivalence << endl;

        if (transformation_equivalence < 1e-10) {
            cout << "Transformations are equivalent" << endl;
            file << "Transformations are equivalent" << endl;
        }
        else {
            cout << "Transformations may not be equivalent" << endl;
            file << "Transformations may not be equivalent" << endl;
        }
    }

};