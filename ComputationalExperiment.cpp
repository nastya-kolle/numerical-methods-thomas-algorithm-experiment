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
#include <fstream>
#include <iomanip>
using namespace std;

// Function to check diagonal dominance (condition 1 from the task)
bool ComputationalExperiment::check_diagonal_dominance(const SpecialMatrix& A) {
    int n = A.get_size();
    for (int i = 1; i <= n; i++) {
        double main_diag = fabs(A.get_main(i));
        double lower_diag = (i > 1) ? fabs(A.get_lower(i)) : 0.0;
        double upper_diag = (i < n) ? fabs(A.get_upper(i)) : 0.0;

        // Check condition: |b_i| ≥ |a_i| + |c_i|
        if (main_diag < lower_diag + upper_diag) {
            return false;
        }
    }
    return true;
}

TridiagonalMatrix ComputationalExperiment::create_well_conditioned_matrix(int size, double min_val, double max_val) {
    TridiagonalMatrix matrix(size);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min_val, max_val);

    for (int i = 1; i <= size; i++) {
        double lower_val = dist(gen);  // a_i
        double upper_val = dist(gen);  // c_i 

        // |b_i| ≥ |a_i| + |c_i|
        double min_main = fabs(lower_val) + fabs(upper_val);
        double max_main = min_main + (max_val - min_val); // Preserve range

        uniform_real_distribution<double> main_dist(min_main, max_main);
        double main_val = main_dist(gen);

        matrix.set_main(i, main_val);
        matrix.set_upper(i, upper_val);
        matrix.set_lower(i, lower_val);
    }
    return matrix;
}

TridiagonalMatrix ComputationalExperiment::create_ill_conditioned_matrix(int size, double min_val, double max_val) {
    TridiagonalMatrix matrix(size);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min_val, max_val);

    for (int i = 1; i <= size; i++) {
        double main_val = dist(gen)* 1e-1; 
        double upper_val = dist(gen);
        double lower_val = dist(gen);
        matrix.set_main(i, main_val);
        matrix.set_upper(i, upper_val);
        matrix.set_lower(i, lower_val);
    }
    return matrix;
}

SpecialMatrix ComputationalExperiment::create_well_conditioned_special_matrix(int size, double min_val, double max_val) {
    SpecialMatrix matrix(size);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min_val, max_val);

    for (int i = 1; i <= size; i++) {
        double lower_val = dist(gen);  // a_i
        double upper_val = dist(gen);  // c_i
        double last_row_val = dist(gen);

        // condition |b_i| ≥ |a_i| + |c_i|
        double min_main = fabs(lower_val) + fabs(upper_val);
        double max_main = min_main + (max_val - min_val);

        uniform_real_distribution<double> main_dist(min_main, max_main);
        double main_val = main_dist(gen);

        matrix.set_main(i, main_val);
        matrix.set_upper(i, upper_val);
        matrix.set_lower(i, lower_val);
        matrix.set_last_row(i, last_row_val);
    }

    // last row
    double last_lower = matrix.get_lower(size);
    double last_main = matrix.get_main(size);
 
    matrix.set_main(size, last_main + fabs(last_lower) + 1.0);

    if (!check_diagonal_dominance(matrix)) {
        throw runtime_error("Well-conditioned matrix does not satisfy diagonal dominance condition");
    }

    return matrix;
}

SpecialMatrix ComputationalExperiment::create_ill_conditioned_special_matrix(int size, double min_val, double max_val) {
    SpecialMatrix matrix(size);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(min_val, max_val);

    for (int i = 1; i <= size; i++) {
        double main_val, lower_val, upper_val, last_row_val;

        if (i == size) {
            main_val = dist(gen) * 1e-2;
            lower_val = dist(gen);
            upper_val = 0; 
            last_row_val = dist(gen);
        }
        else {
            main_val = dist(gen) * 1e-1;
            lower_val = dist(gen);
            upper_val = dist(gen);
            last_row_val = dist(gen);
        }

        matrix.set_main(i, main_val);
        matrix.set_upper(i, upper_val);
        matrix.set_lower(i, lower_val);
        matrix.set_last_row(i, last_row_val);
    }

    if (check_diagonal_dominance(matrix)) {
        throw runtime_error("Ill-conditioned matrix unexpectedly satisfies diagonal dominance condition");
    }

    return matrix;
}

double ComputationalExperiment::compute_relative_error(const Vector& computed, const Vector& exact) {
    double error = 0.0;
    int count = 0;
    for (int i = 1; i <= computed.get_size(); i++) {
        if (fabs(exact[i]) > 1e-12) {
            error += fabs(computed[i] - exact[i]) / fabs(exact[i]);
            count++;
        }
    }
    if (count > 0) {
        return error / count;
    }
    else {
        return 0.0;
    }
}

void ComputationalExperiment::run_experiment() {
    vector<int> dimensions = { 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 };
    vector<double> ranges = { 1e-1, 1e1, 1e2, 1e3 };
    const int num_trials = 10;
    vector<double> tridiagonal_well_errors;
    vector<double> tridiagonal_ill_errors; 
    vector<double> special_well_errors; 
    vector<double> special_ill_errors;

    ofstream file("Computational_Experiment_Results.txt");

    // Table for TridiagonalMatrix (Part 1)
    file << "Table 1 - Tridiagonal Matrices (Part 1)\n\n";
    file << "| Condition Type | System Dimension | Matrix Elements Range | Average Relative Error |\n";

    // Console output for Part 1
    cout << "|----------------------------------------------------------------|" << endl;
    cout << "|                    PART 1: TRIDIAGONAL MATRIX EXPERIMENT       |" << endl;
    cout << "|----------------------------------------------------------------|" << endl;
    cout << "|  Condition        Dimension    Range         Avg Relative Error|" << endl;
    cout << "|----------------------------------------------------------------|" << endl;

    // Test well-conditioned tridiagonal systems
    for (auto size : dimensions) {
        for (auto range : ranges) {
            double total_error = 0.0;
            int valid_trials = 0;

            for (int trial = 0; trial < num_trials; trial++) {
                try {
                    TridiagonalMatrix A = create_well_conditioned_matrix(size, -range, range);
                    Vector exact_solution(size);
                    exact_solution.fill_random(-range, range);
                    Vector b = A * exact_solution;
                    Vector computed_solution = A.Thomas_algorithm(b, file, false);
                    total_error += compute_relative_error(computed_solution, exact_solution);
                    valid_trials++;
                }
                catch (const exception& e) {
                    continue;
                }
            }

            if (valid_trials >= num_trials) {
                double avg_error = total_error / valid_trials;
                tridiagonal_well_errors.push_back(avg_error);
                file << "| Well-conditioned SLAE | " << size << " | " << scientific << setprecision(0) << range << " | " << scientific << setprecision(5) << avg_error << " |\n";

                // Console output
                cout << "|  Well-conditioned  " << setw(8) << size << "    " << setw(8) << scientific << setprecision(0) << range
                    << "    " << scientific << setprecision(5) << avg_error << " |" << endl;
            }
        }
    }

    // Test ill-conditioned tridiagonal systems
    for (auto size : dimensions) {
        for (auto range : ranges) {
            double total_error = 0.0;
            int valid_trials = 0;

            for (int trial = 0; trial < num_trials; trial++) {
                try {
                    TridiagonalMatrix A = create_ill_conditioned_matrix(size, -range, range);
                    Vector exact_solution(size);
                    exact_solution.fill_random(-range, range);
                    Vector b = A * exact_solution;
                    Vector computed_solution = A.Thomas_algorithm(b, file, false);
                    total_error += compute_relative_error(computed_solution, exact_solution);
                    valid_trials++;
                }
                catch (const exception& e) {
                    continue;
                }
            }

            if (valid_trials >= num_trials) {
                double avg_error = total_error / valid_trials;
                tridiagonal_ill_errors.push_back(avg_error);
                file << "| Ill-conditioned SLAE | " << size << " | " << scientific << setprecision(0) << range << " | " << scientific << setprecision(5) << avg_error << " |\n";

                // Console output
                cout << "|  Ill-conditioned   " << setw(8) << size << "    " << setw(8) << scientific << setprecision(0) << range
                    << "    " << scientific << setprecision(5) << avg_error << " |" << endl;
            }
        }
    }
    cout << "|----------------------------------------------------------------|" << endl;

    // Table for SpecialMatrix (Part 2)
    file << "\n\nTable 2 - Special Matrices (Part 2)\n\n";
    file << "| Condition Type | System Dimension | Matrix Elements Range | Average Relative Error |\n";

    // Console output for Part 2
    cout << "\n|----------------------------------------------------------------|" << endl;
    cout << "|                    PART 2: SPECIAL MATRIX EXPERIMENT           |" << endl;
    cout << "|----------------------------------------------------------------|" << endl;
    cout << "|  Condition        Dimension    Range         Avg Relative Error|" << endl;
    cout << "|----------------------------------------------------------------|" << endl;

    // Test well-conditioned special systems
    for (auto size : dimensions) {
        for (auto range : ranges) {
            double total_error = 0.0;
            int valid_trials = 0;

            for (int trial = 0; trial < num_trials; trial++) {
                try {
                    SpecialMatrix A = create_well_conditioned_special_matrix(size, -range, range);
                    Vector exact_solution(size);
                    exact_solution.fill_random(-range, range);
                    Vector b = A * exact_solution;
                    Vector computed_solution = A.Thomas_algorithm(b, file, false);
                    total_error += compute_relative_error(computed_solution, exact_solution);
                    valid_trials++;
                }
                catch (const exception& e) {
                    continue;
                }
            }

            if (valid_trials >= num_trials) {
                double avg_error = total_error / valid_trials;
                special_well_errors.push_back(avg_error);
                file << "| Well-conditioned SLAE | " << size << " | " << scientific << setprecision(0) << range << " | " << scientific << setprecision(5) << avg_error << " |\n";

                // Console output
                cout << "|  Well-conditioned  " << setw(8) << size << "    " << setw(8) << scientific << setprecision(0) << range
                    << "    " << scientific << setprecision(5) << avg_error << " |" << endl;
            }
        }
    }

    // Test ill-conditioned special systems
    for (auto size : dimensions) {
        for (auto range : ranges) {
            double total_error = 0.0;
            int valid_trials = 0;

            for (int trial = 0; trial < num_trials; trial++) {
                try {
                    SpecialMatrix A = create_ill_conditioned_special_matrix(size, -range, range);
                    Vector exact_solution(size);
                    exact_solution.fill_random(-range, range);
                    Vector b = A * exact_solution;
                    Vector computed_solution = A.Thomas_algorithm(b, file, false);
                    total_error += compute_relative_error(computed_solution, exact_solution);
                    valid_trials++;
                }
                catch (const exception& e) {
                    continue;
                }
            }

            if (valid_trials >= num_trials) {
                double avg_error = total_error / valid_trials;
                special_ill_errors.push_back(avg_error);
                file << "| Ill-conditioned SLAE | " << size << " | " << scientific << setprecision(0) << range << " | " << scientific << setprecision(5) << avg_error << " |\n";

                // Console output
                cout << "|  Ill-conditioned   " << setw(8) << size << "    " << setw(8) << scientific << setprecision(0) << range
                    << "    " << scientific << setprecision(5) << avg_error << "|" << endl;
            }
        }
    }
    analyze_dimension_dependence(file);
    analyze_range_dependence(file);
    double avg_well_tridiagonal = 0.0;
    double avg_ill_tridiagonal = 0.0;
    double avg_well_special = 0.0;
    double avg_ill_special = 0.0;
    compute_average_errors(tridiagonal_well_errors, tridiagonal_ill_errors, special_well_errors, special_ill_errors, avg_well_tridiagonal, avg_ill_tridiagonal, avg_well_special, avg_ill_special, file);
    file << "Average error of well-condioned tridiagonal matrix: " << scientific << setprecision(6) << avg_well_tridiagonal << "\n";
    file << "Average error of ill-condioned tridiagonal matrix: " << scientific << setprecision(6) << avg_ill_tridiagonal << "\n";
    file << "Average error of well-condioned special matrix: " << scientific << setprecision(6) << avg_well_special << "\n";
    file << "Average error of ill-condioned special matrix: " << scientific << setprecision(6) << avg_ill_special << "\n";
    cout << "Average error of well-condioned tridiagonal matrix: " << scientific << setprecision(6) << avg_well_tridiagonal << "\n";
    cout << "Average error of ill-condioned tridiagonal matrix: " << scientific << setprecision(6) << avg_ill_tridiagonal << "\n";
    cout << "Average error of well-condioned special matrix: " << scientific << setprecision(6) << avg_well_special << "\n";
    cout << "Average error of ill-condioned special matrix: " << scientific << setprecision(6) << avg_ill_special << "\n";
    cout << "|----------------------------------------------------------------|" << endl;
    cout << "\n Experimental results have been saved to 'Computational_Experiment_Results.txt'" << endl;
}

void ComputationalExperiment::analyze_dimension_dependence(ofstream& file) {
    vector<int> dimensions = { 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 };
    const int num_trials = 10;
    double fixed_range = 1.0;

    cout << "\n DIMENSION DEPENDENCE ANALYSIS" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << " Dimension    Well-Conditioned Error    Ill-Conditioned Error" << endl;
    cout << "----------------------------------------------------------------" << endl;
    file << "\n DIMENSION DEPENDENCE ANALYSIS" << endl;
    file << "----------------------------------------------------------------" << endl;
    file << " Dimension    Well-Conditioned Error    Ill-Conditioned Error" << endl;
    file << "----------------------------------------------------------------" << endl;

    for (auto size : dimensions) {
        double total_error_well = 0.0;
        double total_error_ill = 0.0;
        int valid_trials_well = 0;
        int valid_trials_ill = 0;

        for (int trial = 0; trial < num_trials; trial++) {
            try {
                SpecialMatrix A_well = create_well_conditioned_special_matrix(size, -fixed_range, fixed_range);
                Vector exact_well(size);
                exact_well.fill_random(-fixed_range, fixed_range);
                Vector b_well = A_well * exact_well;
                Vector computed_well = A_well.Thomas_algorithm(b_well, file, false);
                total_error_well += compute_relative_error(computed_well, exact_well);
                valid_trials_well++;
            }
            catch (...) {
                // Skip failed trials
            }

            try {
                SpecialMatrix A_ill = create_ill_conditioned_special_matrix(size, -fixed_range, fixed_range);
                Vector exact_ill(size);
                exact_ill.fill_random(-fixed_range, fixed_range);
                Vector b_ill = A_ill * exact_ill;
                Vector computed_ill = A_ill.Thomas_algorithm(b_ill, file, false);
                total_error_ill += compute_relative_error(computed_ill, exact_ill);
                valid_trials_ill++;
            }
            catch (...) {
                // Skip failed trials
            }
        }

        if (valid_trials_well > 0 && valid_trials_ill > 0) {
            double avg_well = total_error_well / valid_trials_well;
            double avg_ill = total_error_ill / valid_trials_ill;

            cout << " " << setw(8) << size << "      " << scientific << setprecision(5) << avg_well;
            cout << "          " << scientific << setprecision(5) << avg_ill << endl;
            file << " " << setw(8) << size << "      " << scientific << setprecision(5) << avg_well;
            file << "          " << scientific << setprecision(5) << avg_ill << endl;
        }
    }
}

void ComputationalExperiment::analyze_range_dependence(ofstream& file) {
    vector<double> ranges = { 1e-1, 1e1, 1e2, 1e3 };
    const int num_trials = 10;
    int fixed_size = 100;

    cout << "\n RANGE DEPENDENCE ANALYSIS" << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << " Range      Well-Conditioned Error    Ill-Conditioned Error" << endl;
    cout << "----------------------------------------------------------------" << endl;
    file << "\n RANGE DEPENDENCE ANALYSIS" << endl;
    file << "----------------------------------------------------------------" << endl;
    file << " Range      Well-Conditioned Error    Ill-Conditioned Error" << endl;
    file << "----------------------------------------------------------------" << endl;

    for (auto range : ranges) {
        double total_error_well = 0.0;
        double total_error_ill = 0.0;
        int valid_trials_well = 0;
        int valid_trials_ill = 0;

        for (int trial = 0; trial < num_trials; trial++) {
            try {
                SpecialMatrix A_well = create_well_conditioned_special_matrix(fixed_size, -range, range);
                Vector exact_well(fixed_size);
                exact_well.fill_random(-range, range);
                Vector b_well = A_well * exact_well;
                Vector computed_well = A_well.Thomas_algorithm(b_well, file, false);
                total_error_well += compute_relative_error(computed_well, exact_well);
                valid_trials_well++;
            }
            catch (...) {
                // Skip failed trials
            }

            try {
                SpecialMatrix A_ill = create_ill_conditioned_special_matrix(fixed_size, -range, range);
                Vector exact_ill(fixed_size);
                exact_ill.fill_random(-range, range);
                Vector b_ill = A_ill * exact_ill;
                Vector computed_ill = A_ill.Thomas_algorithm(b_ill, file, false);
                total_error_ill += compute_relative_error(computed_ill, exact_ill);
                valid_trials_ill++;
            }
            catch (...) {
                // Skip failed trials
            }
        }

        if (valid_trials_well > 0 && valid_trials_ill > 0) {
            double avg_well = total_error_well / valid_trials_well;
            double avg_ill = total_error_ill / valid_trials_ill;

            cout << " " << scientific << setprecision(0) << range << "      " << scientific << setprecision(5) << avg_well;
            cout << "          " << scientific << setprecision(5) << avg_ill << endl;
            file << " " << scientific << setprecision(0) << range << "      " << scientific << setprecision(5) << avg_well;
            file << "          " << scientific << setprecision(5) << avg_ill << endl;
        }
    }
}
// Функция считает среднее арифметическое двух наборов ошибок:
// - well_errors — ошибки для хорошо обусловленных матриц
// - ill_errors  — ошибки для плохо обусловленных матриц
// Возвращает их как два отдельных значения через параметры (ссылки)
void ComputationalExperiment::compute_average_errors(const vector<double>& tridiagonal_well_errors, const vector<double>& tridiagonal_ill_errors, const vector<double>& special_well_errors, const vector<double>& special_ill_errors,
    double& avg_well_tridiagonal, double& avg_ill_tridiagonal, double& avg_well_special, double& avg_ill_special, ofstream& file) {
    cout << "\n AVERAGE ERRORS" << endl;
    cout << "----------------------------------------------------------------" << endl;
    file << "\n AVERAGE ERRORS" << endl;
    file << "----------------------------------------------------------------" << endl;

    double sum_well_tridiagonal = 0.0;
    double sum_ill_tridiagonal = 0.0;
    double sum_well_special = 0.0;
    double sum_ill_special = 0.0;
    for (double val : tridiagonal_well_errors)
        sum_well_tridiagonal += val;
    for (double val : tridiagonal_ill_errors)
        sum_ill_tridiagonal += val;
    for (double val : special_well_errors)
        sum_well_special += val;
    for (double val : special_ill_errors)
        sum_ill_special += val;

    if (!tridiagonal_well_errors.empty())
        avg_well_tridiagonal = sum_well_tridiagonal / tridiagonal_well_errors.size();
    else
        avg_well_tridiagonal = 0.0;

    if (!tridiagonal_ill_errors.empty())
        avg_ill_tridiagonal = sum_ill_tridiagonal / tridiagonal_ill_errors.size();
    else
        avg_ill_tridiagonal = 0.0;

    if (!special_well_errors.empty())
        avg_well_special = sum_well_special / special_well_errors.size();
    else
        avg_well_special = 0.0;

    if (!special_ill_errors.empty())
        avg_ill_special = sum_ill_special / special_ill_errors.size();
    else
        avg_ill_special = 0.0;
}