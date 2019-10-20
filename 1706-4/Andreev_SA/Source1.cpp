#include <iostream>
#include "mpi.h"
using namespace std;
int countOfPrintSymbol(const char *str, int size)
{
	int count = 0;
	for (int i = 0; i < size; ++i) {
		if ((' ' != str[i]) && ('\0' != str[i]) && ('\n' != str[i]) && ('\t' != str[i]) && ((int)str[i] >= 33) && ((int)str[i] <= 255))
			count++;
	}
	return count;
}

int main(int argc, char **argv)
{
	setlocale(LC_ALL, "rus");
	int ProcRank, ProcNum;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Status status;
	char* str = NULL; int count = 0, buffer;
	int part;
	int len;
	double start = 0;
	double finish,dTime;
	int sizeStr, length;
	if (ProcRank == 0)
	{
		cout << "Size str: ";
		cin >> sizeStr;
		str = new char[sizeStr];
		str[sizeStr - 1] = '\0';
		length = strlen(str);
		for (int i = 0; i < sizeStr - 1; i++)
		{
			str[i] = (char)(rand() % (35 - 32) + 32);
			//str[i] = (char)(rand() % (100 - 32) + 32);
			//cout<< "Symbol:" << str[i] << endl;
		}
		cout << endl;
		part = length / ProcNum;
		len = length % ProcNum;
		start = MPI_Wtime();
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(str + len + part * i, part, MPI_CHAR, i, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&part, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (ProcRank>0)
	{
		char* recvStr = new char[length];
		MPI_Recv(recvStr, length, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		count = countOfPrintSymbol(recvStr, length);
		MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	if (ProcRank == 0)
	{
		//count = countOfPrintSymbol(str, len + part);
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Recv(&buffer, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			count += buffer;
		}
		count += countOfPrintSymbol(str, part+len);
		finish = MPI_Wtime();
		dTime = finish - start;
	}
	if(ProcRank!=0) { cout << ProcRank << "process count:" << count << endl; }
	else {cout<< "Result:" << count; }
	MPI_Finalize();
	return 0;
}
