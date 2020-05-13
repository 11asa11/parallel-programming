#include <iomanip>
#include <iostream>
#include <time.h>
#include <omp.h>

using namespace std;

void PrintMatrix(double* matrix, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::cout.width(6);
            //cout << matrix[i*N + j] << " ";
            std::cout << matrix[i * N + j] << " " << std::setw(6);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void GenerateRandomMatrix(double* matrix1, double* matrix2, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            matrix1[i * N + j] = ((rand() % 101 - 50) / 10.0);
            matrix2[i * N + j] = ((rand() % 101 - 50) / 10.0);
            /*matrix1[i*N + j] = (double)rand() / (double)RAND_MAX* 5;
            matrix2[i*N + j] = (double)rand() / (double)RAND_MAX* 5 ;*/
        }
    }
}

double* CreateMatrix(int N)
{
    double* matrix = new double[N * N];
    return matrix;
}

double* defaultMult(double* matrix1, double* matrix2, int N)
{
    double* tmp = CreateMatrix(N);
    double sum;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++) {
            sum = 0;
            for (int k = 0; k < N; k++)
            {
                sum += matrix1[i * N + k] * matrix2[k * N + j];
            }
            tmp[i * N + j] = sum;
        }
    }
    return tmp;
}

double* Add(double* matrix1, double* matrix2, int N)
{
    double* tmp = CreateMatrix(N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            tmp[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j];
        }
    }
    return tmp;
}


double* Add(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N)
{
    double* tmp = CreateMatrix(N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            tmp[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j] + matrix3[i * N + j] + matrix4[i * N + j];
        }
    }
    return tmp;
}


double* Sub(double* matrix1, double* matrix2, int N)
{
    double* tmp = CreateMatrix(N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            tmp[i * N + j] = matrix1[i * N + j] - matrix2[i * N + j];
        }
    }
    return tmp;
}


double* Sub(double* matrix1, double* matrix2, double* matrix3, double* matrix4, int N)
{
    double* tmp = CreateMatrix(N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            tmp[i * N + j] = matrix1[i * N + j] + matrix2[i * N + j] + matrix3[i * N + j] - matrix4[i * N + j];
        }
    }
    return tmp;
}


bool matComparison(double* matrix1, double* matrix2, int N)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (matrix1[i * N + j] != matrix2[i * N + j])
            {
                return false;
                //std::cout << "Not Equal" << std::endl;
            }
        }
    }
    return true;
}


double* Str_alg(double* matrix1, double* matrix2, int N, int threshold);

double* StrassenParallel(double* matrix1, double* matrix2, int N, int threshold)
{
    double* Rez;

    if (N <= threshold)
        Rez = defaultMult(matrix1, matrix2, N);
    else
    {
        Rez = CreateMatrix(N);
        N = N / 2;

        double* A[4]; double* B[4]; double* C[4]; double* P[7];

        double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
        double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

        for (int i = 0; i < 4; i++)
        {
            A[i] = CreateMatrix(N);
            B[i] = CreateMatrix(N);
        }

        int i, j;
#pragma omp parallel //shared(ThreadsNum)
        {
            //std::cout << "Num: " << omp_get_num_threads() << std::endl;
#pragma omp for private(i,j) schedule(static) 
// schedule(static) — статическое планирование. При использовании такой опции итерации цикла будут поровну (приблизительно) поделены между потоками
// schedule(static, 10) — блочно-циклическое распределение итераций. Каждый поток получает заданное число итераций в начале цикла 
            for (i = 0; i < N; i++)
            {
                for (j = 0; j < N; j++)
                {
                    A[0][i * N + j] = matrix1[2 * i * N + j];
                    A[1][i * N + j] = matrix1[2 * i * N + j + N];
                    A[2][i * N + j] = matrix1[2 * i * N + j + 2 * N * N];
                    A[3][i * N + j] = matrix1[2 * i * N + j + 2 * N * N + N];

                    B[0][i * N + j] = matrix2[2 * i * N + j];
                    B[1][i * N + j] = matrix2[2 * i * N + j + N];
                    B[2][i * N + j] = matrix2[2 * i * N + j + 2 * N * N];
                    B[3][i * N + j] = matrix2[2 * i * N + j + 2 * N * N + N];
                }
            }
#pragma omp sections
            {
#pragma omp section
                {
                    TMP1 = Add(A[0], A[3], N);  // (A11 + A22)
                }
#pragma omp section
                {
                    TMP2 = Add(B[0], B[3], N);  // (B11 + B22)
                }
#pragma omp section
                {
                    TMP3 = Add(A[2], A[3], N); // (A21 + A22)
                }
#pragma omp section
                {
                    TMP4 = Sub(B[1], B[3], N); // (B12 - B22)
                }
#pragma omp section
                {
                    TMP5 = Sub(B[2], B[0], N); // (B21 - B11)
                }
#pragma omp section
                {
                    TMP6 = Add(A[0], A[1], N); // (A11 + A12)
                }
#pragma omp section
                {
                    TMP7 = Sub(A[2], A[0], N); // (A21 - A11)
                }
#pragma omp section
                {
                    TMP8 = Add(B[0], B[1], N); // (B11 + B12)
                }
#pragma omp section
                {
                    TMP9 = Sub(A[1], A[3], N); // (A12 - A22)
                }
#pragma omp section
                {
                    TMP10 = Add(B[2], B[3], N); // (B21 + B22)
                }
            }

#pragma omp sections
            {
#pragma omp section
                {
                    P[0] = Str_alg(TMP1, TMP2, N, threshold); // (A11 + A22)*(B11 + B22)
                }
#pragma omp section
                {
                    P[1] = Str_alg(TMP3, B[0], N, threshold); // (A21 + A22)*B11
                }
#pragma omp section
                {
                    P[2] = Str_alg(A[0], TMP4, N, threshold); // A11*(B12 - B22)
                }
#pragma omp section
                {
                    P[3] = Str_alg(A[3], TMP5, N, threshold); // A22*(B21 - B11)
                }
#pragma omp section
                {
                    P[4] = Str_alg(TMP6, B[3], N, threshold); // (A11 + A12)*B22
                }
#pragma omp section
                {
                    P[5] = Str_alg(TMP7, TMP8, N, threshold); // (A21 - A11)*(B11 + B12)
                }
#pragma omp section
                {
                    P[6] = Str_alg(TMP9, TMP10, N, threshold); // (A12 - A22)*(B21 + B22)
                }
            }


#pragma omp sections
            {
#pragma omp section
                {
                    C[0] = Sub(P[0], P[3], P[6], P[4], N); // P1 + P4 - P5 + P7
                }
#pragma omp section
                {
                    C[1] = Add(P[2], P[4], N); // P3 + P5
                }
#pragma omp section
                {
                    C[2] = Add(P[1], P[3], N); // P2 + P4
                }
#pragma omp section
                {
                    C[3] = Sub(P[0], P[2], P[5], P[1], N); // P1 - P2 + P3 + P6
                }
            }

#pragma omp for private(i,j) schedule(static)
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    Rez[i * 2 * N + j] = C[0][i * N + j];
                    Rez[i * 2 * N + j + N] = C[1][i * N + j];
                    Rez[i * 2 * N + j + 2 * N * N] = C[2][i * N + j];
                    Rez[i * 2 * N + j + 2 * N * N + N] = C[3][i * N + j];
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            delete[] A[i];
            delete[] B[i];
            delete[] C[i];
        }

        for (int i = 0; i < 7; i++) {
            delete[] P[i];
        }

        delete[]TMP1; delete[]TMP2; delete[]TMP3; delete[]TMP4; delete[]TMP5;
        delete[]TMP6; delete[]TMP7; delete[]TMP8; delete[]TMP9; delete[]TMP10;
    }

    //std::cout << "Threads: "<<omp_get_num_threads() << std::endl;

    return Rez;
}


int main()
{
    int degree;
    int Size;
    int ThreadsNum;

    double startDefaultMult = 0;
    double endDefaultMult = 0;
    double startStras = 0;
    double endStras = 0;
    double startStrasParallel = 0;
    double endStrasParallel = 0;

    std::cout << "Degree of matrix size: ";
    std::cin >> degree;
    Size = (int)pow(2, degree);
    std::cout << "Size of matrix: " << Size << std::endl << std::endl;

    std::cout << "Number of threads: ";
    std::cin >> ThreadsNum;
    std::cout << std::endl;
    omp_set_num_threads(ThreadsNum);

    double* matA = nullptr;
    double* matB = nullptr;
    double* matResDefaultMult = nullptr;
    double* matRes = nullptr;
    double* matResParallel = nullptr;

    matA = CreateMatrix(Size);
    matB = CreateMatrix(Size);
    matRes = CreateMatrix(Size);

    GenerateRandomMatrix(matA, matB, Size);

    startStras = omp_get_wtime();
    matRes = Str_alg(matA, matB, Size, 64);
    endStras = omp_get_wtime();

    startStrasParallel = omp_get_wtime();
    matResParallel = StrassenParallel(matA, matB, Size, 64);
    endStrasParallel = omp_get_wtime();

    if (matComparison(matRes, matResParallel, Size) != true) { std::cout << "Mats are not equal" << std::endl << std::endl; }
    else { std::cout << "Mats are equal" << std::endl << std::endl; }

    std::cout << "Strassen alg: " << (endStras - startStras) << std::endl;
    std::cout << "Strassen alg parallel: " << (endStrasParallel - startStrasParallel) << std::endl;


    delete[] matA;
    delete[] matB;
    delete[] matRes;
    delete[] matResParallel;

    return 0;
}

double* Str_alg(double* matrix1, double* matrix2, int N, int threshold)
{
    double* Rez;

    if (N <= threshold)
        Rez = defaultMult(matrix1, matrix2, N);
    else
    {
        Rez = CreateMatrix(N);
        N = N / 2;

        double* A[4]; double* B[4]; double* C[4]; double* P[7];

        double* TMP1; double* TMP2; double* TMP3; double* TMP4; double* TMP5;
        double* TMP6; double* TMP7; double* TMP8; double* TMP9; double* TMP10;

        for (int i = 0; i < 4; i++)
        {
            A[i] = CreateMatrix(N);
            B[i] = CreateMatrix(N);
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                A[0][i * N + j] = matrix1[2 * i * N + j];
                A[1][i * N + j] = matrix1[2 * i * N + j + N];
                A[2][i * N + j] = matrix1[2 * i * N + j + 2 * N * N];
                A[3][i * N + j] = matrix1[2 * i * N + j + 2 * N * N + N];

                B[0][i * N + j] = matrix2[2 * i * N + j];
                B[1][i * N + j] = matrix2[2 * i * N + j + N];
                B[2][i * N + j] = matrix2[2 * i * N + j + 2 * N * N];
                B[3][i * N + j] = matrix2[2 * i * N + j + 2 * N * N + N];
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
            for (int j = 0; j < N; j++) {
                Rez[i * 2 * N + j] = C[0][i * N + j];
                Rez[i * 2 * N + j + N] = C[1][i * N + j];
                Rez[i * 2 * N + j + 2 * N * N] = C[2][i * N + j];
                Rez[i * 2 * N + j + 2 * N * N + N] = C[3][i * N + j];
            }

        for (int i = 0; i < 4; i++) {
            delete[] A[i];
            delete[] B[i];
            delete[] C[i];
        }

        for (int i = 0; i < 7; i++) {
            delete[] P[i];
        }

        delete[]TMP1; delete[]TMP2; delete[]TMP3; delete[]TMP4; delete[]TMP5;
        delete[]TMP6; delete[]TMP7; delete[]TMP8; delete[]TMP9; delete[]TMP10;
    }

    return Rez;
}