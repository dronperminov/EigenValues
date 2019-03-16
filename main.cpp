#include <iostream>
#include <iomanip>
#include <cmath>

#include "Matrix.h"
#include "Vector.h"

using namespace std;

const int maxIterations = 10000; // максимальное число итераций

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

// метод вращений Якоби
void RotationJakobiMethod(Matrix A, double eps) {
	int n = A.rows();
	int iteration = 0; // номер итерации

	while (iteration < maxIterations) {
		double max = 0;

		// находим максимальный по модулю элемент в верхней треугольной части матрицы
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (fabs(A(i, j)) > max)
					max = fabs(A(i, j));

		// если он меньше заданной точности
		if (max < eps)
			break; // выходим

		// выполняем повороты
		for (int i = 0; i < n - 1; i++) {
			for (int j = i + 1; j < n; j++) {
				double theta;

				// если элементы равны
				if (A(i, i) == A(j, j)) {
					theta = M_PI / 4; // то угол равен пи / 4
				}
				else {
					// иначе вычисляем параметры
					double tau = (A(i, i) - A(j, j)) / (2 * A(i, j)); 
					int sign = tau > 0 ? 1 : tau < 0 ? -1 : 0;
					double t = sign / (fabs(tau) + sqrt(1 + tau * tau));

					theta = atan(t); // находим угол
				}

				Matrix R = Matrix::RotationMatrix(n, i, j, theta); // создаём матрицу поворота

				A = R.Transpose() * A * R; // выполняем преобразование матрицы
			}
		}

		iteration++;
	}

	cout << "Rotation Jakobi method:" << endl;

	if (iteration == maxIterations) {
		cout << "Unable to find eigenvalues" << endl;
	}
	else {
		cout << "Eigenvalues: ";

		for (int i = 0; i < n; i++) {
			if (fabs(A(i, i)) < eps)
				cout << "0 ";
			else
				cout << setprecision(15) << A(i, i) << " ";
		}

		cout << endl;
		cout << "Iterations: " << iteration << endl;
	}

	cout << endl;
}

// LU метод
void LUDecompositionMethod(Matrix A, double eps) {
	int n = A.rows();

	Matrix L(n);
	Matrix U(n);

	int iteration = 0;

	while (iteration < maxIterations) {
		// выполняем LU разложение матрицы
		for (int j = 0; j < n; j++) {
			U(0, j) = A(0, j);
			L(j, 0) = A(j, 0) / U(0, 0);
		}

		for (int i = 1; i < n; i++) {
			for (int j = i; j < n; j++) {
				double sum = 0;

				for (int k = 0; k < i; k++)
					sum += L(i, k) * U(k, j);

				U(i, j) = A(i, j) - sum;
			
				sum = 0;

				for (int k = 0; k < i; k++)
					sum += L(j, k) * U(k, i);

				L(j, i) = (A(j, i) - sum) / U(i, i);
			}
		}

		// ищем максимальный по модулю элемент под главной диагональю матрицы L
		double max = 0;

		for (int i = 0; i < n; i++)
			for (int j = 0; j < i; j++)
				if (fabs(L(i, j)) > max)
					max = fabs(L(i, j));

		// если он меньше eps
		if (max < eps)
			break; // выходим

		A = U * L; // иначе задаём новую матрицу
		iteration++; // увеличиваем число итераций
	}

	cout << "LU decomposition method:" << endl;
	
	if (iteration == maxIterations) {
		cout << "Unable to find eigenvalues" << endl;
	}
	else {
		cout << "Eigenvalues: ";

		for (int i = 0; i < n; i++) {
			if (fabs(U(i, i)) < eps)
				cout << "0 ";
			else
				cout << setprecision(15) << U(i, i) << " ";
		}

		cout << endl;
		cout << "Iterations: " << iteration << endl;
	}

	cout << endl;
}

// LR метод
void LRDecompositionMethod(Matrix A, double eps) {
	int n = A.rows();

	int iteration = 0; // номер итерации

	Matrix L(n);
	Matrix R(n);

	while (iteration < maxIterations) {
		// выполняем LR разложение
		for (int i = 0; i < n; i++) {
			L(i, i) = 1; // на диагонали у матрицы L единицы
			R(0, i) = A(0, i);
		}

		for (int i = 1; i < n; i++)
			L(i, 0) = A(i, 0) / R(0, 0);

		for (int i = 1; i < n; i++) {
			for (int k = i; k < n; k++) {
				R(i, k) = A(i, k);

				for (int j = 0; j < i; j++)
					R(i, k) -= L(i, j) * R(j, k);
			}

			for (int k = 1; k < i + 1 && i < n - 1; k++) {
				L(i + 1, k) = A(i + 1, k);

				for (int j = 0; j < k; j++)
					L(i + 1, k) -= L(i + 1, j) * R(j, k);

				L(i + 1, k) /= R(k, k);
			}
		}

		// ищем максимальный по модулю элемент под главной диагональю матрицы L
		double max = 0;

		for (int i = 0; i < n; i++)
			for (int j = 0; j < i; j++)
				if (fabs(L(i, j)) > max)
					max = fabs(L(i, j));

		// если он меньше eps
		if (max < eps)
			break; // выходим
		
		A = R * L; // новя матрица есть R * L
		iteration++;
	}

	cout << "LR decomposition method:" << endl;
	
	if (iteration == maxIterations) {
		cout << "Unable to find eigenvalues" << endl;
	}
	else {
		cout << "Eigenvalues: ";

		for (int i = 0; i < n; i++) {
			if (fabs(R(i, i)) < eps)
				cout << "0 ";
			else
				cout << setprecision(15) << R(i, i) << " ";
		}

		cout << endl;
		cout << "Iterations: " << iteration << endl;
	}

	cout << endl;
}

// QR метод
void QRDecompositionMethod(Matrix A, double eps) {
	int n = A.rows();

	int iteration = 0; // номер итерации

	while (iteration < maxIterations) {
		Matrix Q = A.Ortonormalize();
		Matrix R = Q.Transpose() * A;

		double max = 0;

		// ищем максимальный внедиагональный элемент матрицы Q
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (fabs(Q(i, j)) > max)
					max = fabs(Q(i, j));

		if (max < eps)
			break;

		A = R * Q;
		iteration++;
	}

	cout << "QR decomposition method:" << endl;
	
	if (iteration == maxIterations) {
		cout << "Unable to find eigenvalues" << endl;
	}
	else {
		cout << "Eigenvalues: ";

		for (int i = 0; i < n; i++) {
			if (fabs(A(i, i)) < eps)
				cout << "0 ";
			else
				cout << setprecision(15) << A(i, i) << " ";
		}

		cout << endl;
		cout << "Iterations: " << iteration << endl;
	}

	cout << endl;
}

// разложение Холецкого
void CholeskyDecompositionMethod(Matrix A, double eps) {
	int n = A.rows();

	int iteration = 0;

	Matrix L(n);

	// проверка матрицы на симметричность
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			if (A(i, j) != A(j, i))
				iteration = maxIterations;

	while (iteration < maxIterations) {
		for (int i = 0; i < n; i++) {
			double sum = A(i, i);

			for (int k = 0; k < i; k++)
				sum -= L(i, k) * L(i, k);

			// корень из отрицательных чисел извлечь нельзя
			if (sum < 0) {
				iteration = maxIterations - 1;
				break;
			}

			L(i, i) = sqrt(sum);

			for (int j = i; j < n; j++) {
				sum = 0;

				for (int k = 0; k < i; k++)
					sum += L(i, k) * L(j, k);

				L(j, i) = (A(j, i) - sum) / L(i, i);
			}
		}

		double max = 0;

		// ищем максимальное по модулю значение вне диагонали
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (fabs(A(i, j)) > max)
					max = fabs(A(i, j));

		if (max < eps)
			break;

		A = L.Transpose() * L; // находим матрицу следующей итерации
		iteration++;
	}

	cout << "Holetskii decomposition method:" << endl;
	
	if (iteration == maxIterations) {
		cout << "Unable to find eigenvalues" << endl;
	}
	else {
		cout << "Eigenvalues: ";

		for (int i = 0; i < n; i++) {
			if (fabs(A(i, i)) < eps)
				cout << "0 ";
			else
				cout << setprecision(15) << A(i, i) << " ";
		}

		cout << endl;
		cout << "Iterations: " << iteration << endl;
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
	RotationJakobiMethod(A, eps); // находим все собственные значения по методу вращений Якоби
	LUDecompositionMethod(A, eps); // находим все собственные значения по методу LU разложения
	LRDecompositionMethod(A, eps); // находим все собственные значения по методу LR разложения
	QRDecompositionMethod(A, eps); // находим все собственные значения по методу QR разложения
	CholeskyDecompositionMethod(A, eps); // находим все собственные значения по методу разложения Холецкого
}