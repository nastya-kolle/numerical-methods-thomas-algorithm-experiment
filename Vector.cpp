#include "Vector.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;

Vector::Vector(int size_) : size(size_) {
	data = new double[size];
	for (int i = 0; i < size; i++) {
		data[i] = 0.0;
	}
}

Vector::Vector(const Vector& other) : size(other.size) {
	data = new double[size];
	for (int i = 0; i < size ; i++) {
		data[i] = other.data[i];
	}
}

Vector::~Vector() {
	delete[] data;
}

double& Vector::operator[](int index) {
	if (index < 1 || index > size) {
		throw out_of_range("Vector index out of range");
	}
	return data[index - 1];
}

const double& Vector::operator[](int index) const {
	if (index < 1 || index > size) {
		throw out_of_range("Vector index out of range");
	}
	return data[index - 1];
}

Vector Vector::operator-(const Vector& other) const {
	return subtract(other);
}

Vector& Vector::operator=(const Vector& other) {
	if (this != &other) {
		delete[] data;
		size = other.size;
		data = new double[size];
		for (int i = 0; i < size; i++) {
			data[i] = other.data[i];
		}
	}
	return *this;
}

Vector Vector::subtract(const Vector& other) const {
	if (size != other.size) {
		throw invalid_argument("Vectors must have the same size and start index");
	}
	Vector result(size);
	for (int i = 0; i < size; i++) {
		result.data[i] = data[i] - other.data[i];
	}
	return result;
}

double Vector::norm() const {
	double max = 0.0;
	for (int i = 0; i < size; i++) {
		if (abs(data[i]) > abs(max)) max = abs(data[i]);
	}
	return max;
}

void Vector::fill_random(double min_val, double max_val) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(min_val, max_val);

	for (int i = 0; i < size; i++) data[i] = dist(gen);
}

int Vector::get_size() const { return size; }

void Vector::print() const {
	cout << "[";
	for (int i = 0; i < size; i++) {
		cout << data[i];
		if (i < size - 1) cout << ", ";
	}
	cout << "]" << endl;
}

void Vector::copy_from(const Vector& other) {
	if (size != other.size) {
		throw invalid_argument("Vector sizes don't match");
	}
	data = other.data;
}

void Vector::print_to_file(ofstream& file) const {
	file << "[";
	for (int i = 0; i < size; i++) {
		file << data[i];
		if (i < size - 1) file << ", ";
	}
	file << "]" << endl;
}