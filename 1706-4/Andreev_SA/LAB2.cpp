#define MSMPI_NO_DEPRECATE_20
#include <iostream>
#include "mpi.h"
#include <typeinfo>
using namespace std;
#define ELEMS(x)  ( sizeof(x) / sizeof(x[0]) )

bool check(int num)
{
	if (num >= 2) { return true; }
	else { false; }
}

int getPart(int size, int part)
{
	int res = size - part;
	return res;
}

int rProc(int r, int root, int procNum)
{
	return (r - 1 + root) % procNum;
}

//template<typename T>
void treeScatter(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
	int ProcNum, ProcRank;
	MPI_Comm_size(comm, &ProcNum);
	MPI_Aint sendExtent;
	MPI_Aint recvExtent;
	MPI_Comm_rank(comm, &ProcRank);
	MPI_Type_extent(sendtype, &sendExtent);
	MPI_Type_extent(recvtype, &recvExtent);
	MPI_Status status;
	int r = ((ProcRank - root) + ProcNum) % ProcNum + 1;
	int recvProc = rProc(2 * r, root, ProcNum);
	for (int i = 2 * r; i <= ProcNum; i++)
	{
		for (int j = recvProc; j <= recvProc + 1; j++)
		{
			MPI_Send((char*)sendbuf + j * sendcount*sendExtent, sendcount, sendtype, j, 0, comm);
		}
	}

	if (ProcRank != root)
	{
		MPI_Recv(recvbuf, recvcount, recvtype, rProc(r / 2, root, ProcNum), 0, comm, MPI_STATUS_IGNORE);
	}
}

void printArr(int*& arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << arr[i] << endl;
	}
}

const int root = 6;

int n = 10;

int countProc = 2;

int plus = 0;

int main()
{
	setlocale(LC_ALL, "rus");
	int ProcNum, ProcRank;
	//int size = rand() % 15 + 13;
	//int* arr = new int[n];

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status status;
	int part, shift;
	part = n / ProcNum;
	int size;

	double* data = new double[7];
	double num;
	double k = 0.3;
	for (int i = 0; i < 7; i++)
	{
		data[i] = i+k;
		//cout << "DATA " << i << " = "<< data[i] << endl;;
	}
	part = 7 / ProcNum;

	treeScatter(data, part, MPI_DOUBLE, &num, part, MPI_DOUBLE, root, MPI_COMM_WORLD);
	//MPI_Scatter(data, part, MPI_INT, &num, part, MPI_INT, 6, MPI_COMM_WORLD);
	if (ProcRank != root)
	{
		cout << ProcRank << " recieved " << num << endl;
	}
	MPI_Finalize();
	//system("pause");
	return 0;
}