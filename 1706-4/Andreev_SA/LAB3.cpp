
#include <iostream>
#include <time.h>
#include <mpi.h>

using namespace std;

double* Str_alg(double* matrix1, double* matrix2, double N, double threshold);

double* CreateMatrix(double N)
{
	double *matrix = new double[N*N];
	return matrix;
}

void PrdoubleMatrix(double* matrix, double N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << matrix[i*N + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void GenerateRandomMatrix(double* matrix1, double* matrix2, double N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matrix1[i*N + j] = rand() % 10;
			matrix2[i*N + j] = rand() % 10;
		}
	}
}

double* defaultMult(double* matrix1, double* matrix2, double N)
{
	double* tmp = new double[N*N];
	double sum;
	for (double i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) {
			sum = 0;
			for (int k = 0; k < N; k++)
			{
				sum += matrix1[i*N + k] * matrix2[k*N + j];
			}
			tmp[i*N + j] = sum;
		}
	}
	return tmp;
}

double* Add(double* matrix1, double* matrix2, double N)
{
	double* tmp = new double[N*N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tmp[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j];
		}
	}
	return tmp;
}

double* Add(double* matrix1, double* matrix2, double* matrix3, double* matrix4, double N)
{
	double* tmp = new double[N*N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tmp[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j] + matrix3[i*N + j] + matrix4[i*N + j];
		}
	}
	return tmp;
}

double* Sub(double* matrix1, double* matrix2, double N)
{
	double* tmp = new double[N*N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tmp[i*N + j] = matrix1[i*N + j] - matrix2[i*N + j];
		}
	}
	return tmp;
}

double* Sub(double* matrix1, double* matrix2, double* matrix3, double* matrix4, double N)
{
	double* tmp = new double[N*N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tmp[i*N + j] = matrix1[i*N + j] + matrix2[i*N + j] + matrix3[i*N + j] - matrix4[i*N + j];
		}
	}
	return tmp;
}

double* Str_alg_Pp_MPI(double* matr_A, double* matr_B, double N, double thr, int ProcSize, int ProcRank)
{
	MPI_Status Status;
	double* matr_Rez_Str_PP = NULL;
	//ProcSize =Numbers(ProcSize);
	if (N <= thr)
	{
		if (ProcRank == 0)
		{
			matr_Rez_Str_PP = defaultMult(matr_A, matr_B, N);
		}
		return matr_Rez_Str_PP;
	}
	double i, j;
	if (ProcRank == 0)
	{
		matr_Rez_Str_PP = new double[N*N];
		N = N / 2;

		double* TMP1 = NULL; double* TMP2 = NULL; double* TMP3 = NULL; double* TMP4 = NULL; double* TMP5 = NULL;
		double* TMP6 = NULL; double* TMP7 = NULL; double* TMP8 = NULL; double* TMP9 = NULL; double* TMP10 = NULL;

		double* A[4]; double* B[4]; double* C[4]; double* P[7];

		for (int i = 0; i < 4; i++)
		{
			A[i] = new double[N*N];
			B[i] = new double[N*N];
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A[0][i*N + j] = matr_A[2 * i*N + j];
				A[1][i*N + j] = matr_A[2 * i*N + j + N];
				A[2][i*N + j] = matr_A[2 * i*N + j + 2 * N*N];
				A[3][i*N + j] = matr_A[2 * i*N + j + 2 * N*N + N];

				B[0][i*N + j] = matr_B[2 * i*N + j];
				B[1][i*N + j] = matr_B[2 * i*N + j + N];
				B[2][i*N + j] = matr_B[2 * i*N + j + 2 * N*N];
				B[3][i*N + j] = matr_B[2 * i*N + j + 2 * N*N + N];
			}
		}

		TMP1 = Add(A[0], A[3], N);
		TMP2 = Add(B[0], B[3], N);
		TMP3 = Add(A[2], A[3], N);
		TMP4 = Sub(B[1], B[3], N);
		TMP5 = Sub(B[2], B[0], N);
		TMP6 = Add(A[0], A[1], N);
		TMP7 = Sub(A[2], A[0], N);
		TMP8 = Add(B[0], B[1], N);
		TMP9 = Sub(A[1], A[3], N);
		TMP10 = Add(B[2], B[3], N);

		P[0] = Str_alg(TMP1, TMP2, N, thr);

		if (ProcSize == 7)
		{
			MPI_Send(TMP3, N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(B[0], N*N, MPI_double, 1, 1, MPI_COMM_WORLD);
			MPI_Send(A[0], N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(TMP4, N*N, MPI_double, 2, 1, MPI_COMM_WORLD);
			MPI_Send(A[3], N*N, MPI_double, 3, 0, MPI_COMM_WORLD);
			MPI_Send(TMP5, N*N, MPI_double, 3, 1, MPI_COMM_WORLD);
			MPI_Send(TMP6, N*N, MPI_double, 4, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 4, 1, MPI_COMM_WORLD);
			MPI_Send(TMP7, N*N, MPI_double, 5, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 5, 1, MPI_COMM_WORLD);
			MPI_Send(TMP9, N*N, MPI_double, 6, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 6, 1, MPI_COMM_WORLD);

			for (int i = 1; i < 7; i++)
			{
				P[i] = new double[N*N];
			}

			for (int i = 1; i < ProcSize; i++)
			{
				MPI_Recv(P[i], N*N, MPI_double, i, i, MPI_COMM_WORLD, &Status);
			}

		}
		else if (ProcSize == 6)
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);

			MPI_Send(A[0], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP4, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);
			MPI_Send(A[3], N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(TMP5, N*N, MPI_double, 2, 1, MPI_COMM_WORLD);
			MPI_Send(TMP6, N*N, MPI_double, 3, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 3, 1, MPI_COMM_WORLD);
			MPI_Send(TMP7, N*N, MPI_double, 4, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 4, 1, MPI_COMM_WORLD);
			MPI_Send(TMP9, N*N, MPI_double, 5, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 5, 1, MPI_COMM_WORLD);

			for (int i = 2; i < 7; i++)
			{
				P[i] = new double[N*N];
			}

			for (int i = 2, j = 1; i < ProcSize; i++, j++)
			{
				MPI_Recv(P[3], N*N, MPI_double, j, j, MPI_COMM_WORLD, &Status);
			}
		}
		else if (ProcSize == 5)
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);

			for (double i = 2; i < 7; i++)
			{
				P[i] = new double[N*N];
			}

			MPI_Send(A[0], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP4, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(A[3], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP5, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP6, N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 2, 1, MPI_COMM_WORLD);
			MPI_Send(TMP7, N*N, MPI_double, 3, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 3, 1, MPI_COMM_WORLD);
			MPI_Send(TMP9, N*N, MPI_double, 4, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 4, 1, MPI_COMM_WORLD);

			for (double i = 2, j = 1; i < 6; i++, j++)
			{
				MPI_Recv(P[i], N*N, MPI_double, j, j, MPI_COMM_WORLD, &Status);
			}
		}
		else if (ProcSize == 4)
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);

			for (double i = 2; i < 7; i++)
			{
				P[i] = new double[N*N];
			}

			MPI_Send(A[0], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP4, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(A[3], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP5, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP6, N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 2, 1, MPI_COMM_WORLD);

			MPI_Send(TMP7, N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 2, 1, MPI_COMM_WORLD);

			MPI_Send(TMP9, N*N, MPI_double, 3, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 3, 1, MPI_COMM_WORLD);

			for (int i = 2, j = 1; i < ProcSize; i++, j++)
			{
				MPI_Recv(P[i], N*N, MPI_double, j, j, MPI_COMM_WORLD, &Status);
			}
			
		}
		else if (ProcSize == 3)
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);
			P[2] = Str_alg(A[0], TMP4, N, thr);

			for (double i = 3; i < 7; i++)
			{
				P[i] = new double[N*N];
			}
			MPI_Send(A[3], N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP5, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP6, N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP7, N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 2, 1, MPI_COMM_WORLD);

			MPI_Send(TMP9, N*N, MPI_double, 2, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 2, 1, MPI_COMM_WORLD);
			
			for (int i = 3,j=1; i < ProcSize; i++,j++)
			{
				MPI_Recv(P[i], N*N, MPI_double, j, j, MPI_COMM_WORLD, &Status);
			}
		}
		else if (ProcSize == 2)
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);
			P[2] = Str_alg(A[0], TMP4, N, thr);
			P[3] = Str_alg(A[3], TMP5, N, thr);

			for (double i = 4; i < 7; i++)
			{
				P[i] = new double[N*N];
			}
			MPI_Send(TMP6, N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(B[3], N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP7, N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP8, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			MPI_Send(TMP9, N*N, MPI_double, 1, 0, MPI_COMM_WORLD);
			MPI_Send(TMP10, N*N, MPI_double, 1, 1, MPI_COMM_WORLD);

			for (int i = 4, j = 1; i < ProcSize; i++, j++)
			{
				MPI_Recv(P[i], N*N, MPI_double, j, j, MPI_COMM_WORLD, &Status);
			}
		}
		else
		{
			P[1] = Str_alg(TMP3, B[0], N, thr);
			P[2] = Str_alg(A[0], TMP4, N, thr);
			P[3] = Str_alg(A[3], TMP5, N, thr);
			P[4] = Str_alg(TMP6, B[3], N, thr);
			P[5] = Str_alg(TMP7, TMP8, N, thr);
			P[6] = Str_alg(TMP9, TMP10, N, thr);
		}

		C[0] = Sub(P[0], P[3], P[6], P[4], N);
		C[1] = Add(P[2], P[4], N);
		C[2] = Add(P[1], P[3], N);
		C[3] = Sub(P[0], P[2], P[5], P[1], N);

		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++) {
				matr_Rez_Str_PP[i * 2 * N + j] = C[0][i*N + j];
				matr_Rez_Str_PP[i * 2 * N + j + N] = C[1][i*N + j];
				matr_Rez_Str_PP[i * 2 * N + j + 2 * N*N] = C[2][i*N + j];
				matr_Rez_Str_PP[i * 2 * N + j + 2 * N*N + N] = C[3][i*N + j];
			}
		}
	}
	else if (ProcRank > 0 && ProcRank < 7)
	{
		double* tmp1; double* tmp2; double* tmp3;

		N = N / 2;

		tmp1 = new double[N*N];
		tmp2 = new double[N*N];

		double end = 7 / ProcSize;
		double size_work = 7 % ProcSize;

		if (ProcRank < size_work)
		{
			end++;
		}

		for (double i = 0; i < end; i++)
		{
			MPI_Recv(tmp1, N*N, MPI_double, 0, 0, MPI_COMM_WORLD, &Status);
			MPI_Recv(tmp2, N*N, MPI_double, 0, 1, MPI_COMM_WORLD, &Status);
			tmp3 = Str_alg(tmp1, tmp2,N, thr);
			MPI_Send(tmp3, N*N, MPI_double, 0, ProcRank, MPI_COMM_WORLD);
			delete[] tmp3;
		}
	}

	return matr_Rez_Str_PP;
}

double Numbers(double& size);

int thr = 1;

int main(int argc, char** argv)
{
	double *matr_A = NULL;
	double *matr_B = NULL;
	double *matr_Rez_Str = NULL;
	double *matr_Rez_Str_PP = NULL;
	int ProcSize, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	double num_th = 1;
	double k = 2;
	double N = (double)pow(2.0, k);
	if (ProcRank == 0)
	{
		cout << "Size of matrix: " << N << " x " << N << endl;
		matr_A = new double[N*N];
		matr_B = new double[N*N];
		GenerateRandomMatrix(matr_A, matr_B, N);
	}

	if (ProcRank == 0)
	{
		//double *matr_Rez_Check = NULL;
		cout << "Matrix A: " << endl;
		PrdoubleMatrix(matr_A, N);
		cout << "Matrix B: " << endl;
		PrdoubleMatrix(matr_B, N);
		cout << "Matrix C: " << endl;
		PrdoubleMatrix(matr_Rez_Str_PP, N);

		delete[] matr_Rez_Str_PP;
		delete[] matr_B;
		delete[] matr_A;

	}
	MPI_Finalize();
	return 0;
}

double* Str_alg(double* matrix1, double* matrix2, double N, double threshold)
{
	double* Rez;

	if (N <= threshold)
		Rez = defaultMult(matrix1, matrix2, N);
	else
	{
		Rez = new double[N*N];
		N = N / 2;

		double* A[4]; double* B[4]; double* C[4]; double* P[7];

		double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
		double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

		/* Выделение памяти под вспомогательные матрицы */
		for (int i = 0; i < 4; i++)
		{
			A[i] = new double[N*N];
			B[i] = new double[N*N];
		}

		/* Разбиение матрицы на 4 блока */
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A[0][i*N + j] = matrix1[2 * i*N + j];
				A[1][i*N + j] = matrix1[2 * i*N + j + N];
				A[2][i*N + j] = matrix1[2 * i*N + j + 2 * N*N];
				A[3][i*N + j] = matrix1[2 * i*N + j + 2 * N*N + N];

				B[0][i*N + j] = matrix2[2 * i*N + j];
				B[1][i*N + j] = matrix2[2 * i*N + j + N];
				B[2][i*N + j] = matrix2[2 * i*N + j + 2 * N*N];
				B[3][i*N + j] = matrix2[2 * i*N + j + 2 * N*N + N];
			}
		}

		TMP1 = Add(A[0], A[3], N);
		TMP2 = Add(B[0], B[3], N);
		TMP3 = Add(A[2], A[3], N);
		TMP4 = Sub(B[1], B[3], N);
		TMP5 = Sub(B[2], B[0], N);
		TMP6 = Add(A[0], A[1], N);
		TMP7 = Sub(A[2], A[0], N);
		TMP8 = Add(B[0], B[1], N);
		TMP9 = Sub(A[1], A[3], N);
		TMP10 = Add(B[2], B[3], N);

		P[0] = Str_alg(TMP1, TMP2, N, threshold); // (A11 + A22)*(B11 + B22)
		P[1] = Str_alg(TMP3, B[0], N, threshold); // (A21 + A22)*B11
		P[2] = Str_alg(A[0], TMP4, N, threshold); // A11*(B12 - B22)
		P[3] = Str_alg(A[3], TMP5, N, threshold); // A22*(B21 - B11)
		P[4] = Str_alg(TMP6, B[3], N, threshold); // (A11 + A12)*B22
		P[5] = Str_alg(TMP7, TMP8, N, threshold); // (A21 - A11)*(B11 + B12)
		P[6] = Str_alg(TMP9, TMP10, N, threshold); // (A12 - A22)*(B21 + B22)

		C[0] = Sub(P[0], P[3], P[6], P[4], N); // P1 + P4 - P5 + P7
		C[1] = Add(P[2], P[4], N); // P3 + P5
		C[2] = Add(P[1], P[3], N); // P2 + P4
		C[3] = Sub(P[0], P[2], P[5], P[1], N); // P1 - P2 + P3 + P6

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++) {
				Rez[i * 2 * N + j] = C[0][i*N + j];
				Rez[i * 2 * N + j + N] = C[1][i*N + j];
				Rez[i * 2 * N + j + 2 * N*N] = C[2][i*N + j];
				Rez[i * 2 * N + j + 2 * N*N + N] = C[3][i*N + j];
			}
		}
	}

	return Rez;
}

double Numbers(double& size)
{
	size = 7;
	return size;
}