#pragma once

#include <iostream>
#include <iomanip>

#include "Vector.h"

// класс матрицы
class Matrix {
	double **values;
	int n;
	int m;

public:
	Matrix(int n); // конструктор квадратной матрицы из размера
	Matrix(int n, int m); // конструктор из размеров

	Matrix(const Matrix& matrix); // конструктор копирования

	Matrix& operator=(const Matrix& matrix); // оператор присваивания

	int rows() const; // получение числа строк
	int cols() const; // получение числа столбцов

	double operator()(int i, int j) const; // получение элемента по индексу
	double& operator()(int i, int j); // получение элемента по индексу

	Vector operator*(const Vector &vector) const; // умножение матрицы на вектор
	Matrix operator*(const Matrix &matrix) const; // умножение матрицы на матрицу
	
	void Read(); // заполнение матрицы с клавиатуры
	
	Matrix Transpose() const; // получение транспонированной матрицы

	~Matrix(); // деструктор (освобождение памяти)

	friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix); // вывод в поток
};