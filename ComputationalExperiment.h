#pragma once
#include "Vector.h"
#include "TridiagonalMatrix.h"
#include "SpecialMatrix.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <iomanip>

class ComputationalExperiment {
private:
    // Методы создания матриц
    TridiagonalMatrix create_well_conditioned_matrix(int size, double min_val, double max_val);
    TridiagonalMatrix create_ill_conditioned_matrix(int size, double min_val, double max_val);
    SpecialMatrix create_well_conditioned_special_matrix(int size, double min_val, double max_val);
    SpecialMatrix create_ill_conditioned_special_matrix(int size, double min_val, double max_val);

    // Вспомогательные методы
    double compute_relative_error(const Vector& computed, const Vector& exact);
    bool check_diagonal_dominance(const SpecialMatrix& A);

public:
    // Основные методы экспериментов
    void run_experiment();
    void analyze_dimension_dependence(ofstream& file);
    void analyze_range_dependence(ofstream& file);
    void compute_average_errors(const vector<double>& tridiagonal_well_errors, const vector<double>& tridiagonal_ill_errors, const vector<double>& special_well_errors, const vector<double>& special_ill_errors,
        double& avg_well_tridiagonal, double& avg_ill_tridiagonal, double& avg_well_special, double& avg_ill_special, ofstream& file);
};