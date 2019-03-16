#include "Matrix.h"

// конструктор rквадратной матрицы из размера
Matrix::Matrix(int n) {
	this->n = n; // сохраняем размер
	this->m = n;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[n];

		for (int j = 0; j < n; j++)
			values[i][j] = 0; // обнуляем все элементы
	}
}

// конструктор из размеров
Matrix::Matrix(int n, int m) {
	this->n = n; // сохраняем размер
	this->m = m;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[m];

		for (int j = 0; j < m; j++)
			values[i][j] = 0; // обнуляем все элементы
	}
}

// конструктор копирования
Matrix::Matrix(const Matrix& matrix) {
	n = matrix.n;
	m = matrix.m;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[m];

		for (int j = 0; j < m; j++)
			values[i][j] = matrix.values[i][j]; // копируем значения матрицы
	}
}

// оператор присваивания
Matrix& Matrix::operator=(const Matrix& matrix) {
	if (this == &matrix)
		return *this;

	for (int i = 0; i < n; i++)
		delete[] values[i];

	delete[] values;

	n = matrix.n;
	m = matrix.m;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[m];

		for (int j = 0; j < m; j++)
			values[i][j] = matrix.values[i][j]; // копируем значения матрицы
	}

	return *this;
}

// получение числа строк
int Matrix::rows() const {
	return n;
}

// получение числа столбцов
int Matrix::cols() const {
	return m;
}

// получение элемента по индексу
double Matrix::operator()(int i, int j) const {
	return values[i][j];
}

// получение элемента по индексу
double& Matrix::operator()(int i, int j) {
	return values[i][j];
}

// умножение матрицы на вектор
Vector Matrix::operator*(const Vector &vector) const {
	Vector result(n);

	for (int i = 0; i < n; i++) {
		double sum = 0;

		for (int j = 0; j < m; j++)
			sum += values[i][j] * vector[j];

		result[i] = sum;
	}

	return result;
}

// умножение матрицы на матрицу
Matrix Matrix::operator*(const Matrix &matrix) const {
	Matrix result(n, matrix.m);

	for (int i = 0; i < result.n; i++) {
		for (int j = 0; j < result.m; j++) {
			double sum = 0;

			for (int k = 0; k < m; k++)
				sum += values[i][k] * matrix.values[k][j];

			result.values[i][j] = sum;
		}
	}

	return result;
}

// заполнение матрицы с клавиатуры
void Matrix::Read() {
	for (int i = 0; i < n; i++) {
		std::cout << "Enter row " << (i + 1) << ": ";

		for (int j = 0; j < m; j++)
			std::cin >> values[i][j];
	}
}

// получение транспонированной матрицы
Matrix Matrix::Transpose() const {
	Matrix T(m, n); // создаём матрицу

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			T.values[i][j] = values[j][i]; // копируем транспонированные элементы

	return T; // возвращаем матрицу
}

Matrix Matrix::Ortonormalize() const {
	Matrix Q(n); // ортонормированная мтарица
	Matrix B(n); // ортогональная матрица

	double beta_1 = 0; // вычисляем первый коэффициент

	for (int j = 0; j < n; j++)
		beta_1 = beta_1 + values[0][j] * values[0][j]; // считаем сумму квадратов

	beta_1 = sqrt(beta_1); // извлекаем квадратный корень

	for (int j = 0; j < n; j++) {
		B.values[0][j] = values[0][j]; // запоминаем исходый вектор
		Q.values[0][j] = values[0][j] / beta_1; // запоминаем нормированный вектор
	}

	// проводим процесс ортогонализации Грама-Шмидта по столбцам
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B.values[i][j] = values[i][j];

			for (int k = 0; k <= i - 1; k++) {
				double s = 0;
				double mod = 0;

				for (int m = 0; m < n; m++) {
					s += values[i][m] * B.values[k][m];
					mod += B.values[k][m] * B.values[k][m];
				}

				B.values[i][j] -= s / mod * B.values[k][j];
			}
		}

		// аналогично находим нормирующий множитель
		double beta_i = 0;
		for (int j = 0; j < n; j++)
			beta_i += B.values[i][j] * B.values[i][j];

		beta_i = sqrt(beta_i);

		for (int j = 0; j < n; j++)
			Q.values[i][j] = B.values[i][j] / beta_i; // нормируем вектор
	}

	return Q; // возвращаем матрицу
}

// создание матрицы поворота
Matrix Matrix::RotationMatrix(int n, int i, int j, double theta) {
	Matrix R(n);

	// на главной диагонали единицы
	for (int k = 0; k < n; k++)
		R.values[k][k] = 1;

	// а в элементах j0 k0 косинусы и синусы угла theta
	R.values[i][i] = cos(theta);
	R.values[j][j] = cos(theta);

	R.values[i][j] = -sin(theta);
	R.values[j][i] = sin(theta);

	return R; // возвращаем матрицу поворота
}

// деструктор (освобождение памяти)
Matrix::~Matrix() {
	for (int i = 0; i < n; i++)
		delete[] values[i];

	delete[] values;
}

// вывод в поток
std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
	for (int i = 0; i < matrix.n; i++) {
		for (int j = 0; j < matrix.m; j++)
			os << std::setw(6) << matrix.values[i][j] << " ";

		os << std::endl;
	}

	return os;
}