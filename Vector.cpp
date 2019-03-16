#include "Vector.h"

// конструктор из размера
Vector::Vector(int n) {
	this->n = n;

	values = new double[n]; // выделяем память под значения

	// обнуляем значения
	for (int i = 0; i < n; i++)
		values[i] = 0;
}

// конструктор копирования
Vector::Vector(const Vector& vector) {
	n = vector.n; // копируем размерность

	values = new double[n]; // выделяем память

	// копируем значения
	for (int i = 0; i < n; i++)
		values[i] = vector.values[i];
}

// оператор присваивания
Vector& Vector::operator=(const Vector& vector) {
	if (this == &vector)
		return *this;

	delete[] values; // удаляем старый массив

	n = vector.n; // копируем размерность

	values = new double[n]; // выделяем память под новый массив

	// копируем значения
	for (int i = 0; i < n; i++)
		values[i] = vector.values[i];

	return *this; // возвращаем объект
}

// получение размера вектора
int Vector::size() const {
	return n;
}

// поиск первого ненулевого элемента
int Vector::NonZeroIndex() const {
	for (int i = 0; i < n; i++)
		if (fabs(values[i]) > 0)
			return i;

	return -1; // не нашли
}

// получение значения элемента по индексу
double& Vector::operator[](int index) {
	return values[index];
}

// получение значения элемента по индексу
double Vector::operator[](int index) const {
	return values[index];
}

// вывод в поток
std::ostream& operator<<(std::ostream& os, const Vector& vector) {
	for (int i = 0; i < vector.n; i++)
		os << vector.values[i] << " ";

	return os;
}

// деструктор, освобождение памяти
Vector::~Vector() {
	delete[] values; // освобождаем память, занимаемую массивом значений
}