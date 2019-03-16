#include <iostream>
#include <iomanip>
#include <cmath>

#include "Matrix.h"
#include "Vector.h"

using namespace std;

// считывание матрицы
Matrix ReadMatrix() {
	int n; // размер матрицы
	cout << "Enter size of matrix: ";
	cin >> n; // считываем размер матрицы

	Matrix matrix(n); // создаём матрицу для поиска собственных значений

	cout << "Enter values of matrix:" << endl;
	matrix.Read(); // считываем значения матрицы

	return matrix; // возвращаем матрицу
}

// степенной метод (находит максимальное собственное значение)
void PowerMethod(const Matrix& A,  double eps) {
	int n = A.rows();
	Vector x0(n); // случайное начальное приближение

	for (int i = 0; i < n; i++)
		x0[i] = rand() * 1.0 / RAND_MAX;

	Vector x1 = A * x0;

	int index = x0.NonZeroIndex();
	double value1;
	double value2 = x1[index] / x0[index];
	int k = 0;

	do {
		x0 = x1;
		x1 = A * x0; // находим новое приближение вектора

		index = x0.NonZeroIndex(); // ищем индекс ненулевого элемента
		value1 = value2;
		value2 = x1[index] / x0[index]; // находим новое приблжение числа

		k++; // увеличиваем номер итерации
	} while (fabs(value2 - value1) > eps);

	cout << "Power method" << endl;

	// если вектор нулевой
	if (index == -1) {
		cout << "Unable to find max eigenvalue" << endl; // сообщаем, что нет возможности найти
	}
	else {
		cout << "Max eigenvalue: " << setprecision(15) << value2 << endl;
		cout << "Iterations: " << k << endl;
	}
	
	cout << endl;
}

int main() {
	Matrix A = ReadMatrix(); // считываем матрицу

	cout << endl << "Entered matrix: " << endl << A;

	double eps; // точность
	cout << endl << "Enter eps: ";
	cin >> eps; // считываем точность

	PowerMethod(A, eps); // находим наибольшее собственное число
}	