#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <string>

class Vector {
private:
    double* data;
    int size;

public:
    // Конструкторы и деструктор
    
    Vector(int size_ = 2);
    Vector(const Vector& other);
    ~Vector();

    // Операторы доступа


    double& operator[](int index);
    const double& operator[](int index) const;

    // Операторы
    Vector operator-(const Vector& other) const;
    Vector& operator=(const Vector& other);

    // Методы
    Vector subtract(const Vector& other) const;
    double norm() const;
    void fill_random(double min_val, double max_val);
    void print() const;
    void copy_from(const Vector& other);
    void print_to_file(std::ofstream& file) const;

    // Геттеры
    int get_size() const;
};