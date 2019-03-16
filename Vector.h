#pragma once

#include <iostream>
#include <cmath>

//класс n-мерного вектора
class Vector {
	int n; // размерность вектора
	double *values; // значения вектора

public:
	Vector(int n); // конструктор из размера
	Vector(const Vector& vector); // конструктор копирования

	Vector& operator=(const Vector& vector); // оператор присваивания

	int size() const; // получение размера вектора
	int NonZeroIndex() const; // поиск первого ненулевого элемента

	double& operator[](int index); // получение значения элемента по индексу
	double operator[](int index) const; // получение значения элемента по индексу

	friend std::ostream& operator<<(std::ostream& os, const Vector& vector); // вывод в поток

	~Vector(); // деструктор, освобождение памяти
};