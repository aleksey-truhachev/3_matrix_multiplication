#include <iostream>
#include "math.h"
#include <vector>
using namespace std;

typedef double datatype;


# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#define MINI_DATASET
# endif
# if !defined(NI) && !defined(NJ) && !defined(NK) && !defined(NL) && !defined(NM)
# ifdef MINI_DATASET
#define NI 16
#define NJ 18
#define NK 20
#define NL 22
#define NM 24
# endif
# ifdef SMALL_DATASET
#define NI 40
#define NJ 50
#define NK 60
#define NL 70
#define NM 80
# endif
# ifdef MEDIUM_DATASET
#define NI 180
#define NJ 190
#define NK 200
#define NL 210
#define NM 220
# endif
# ifdef LARGE_DATASET
#define NI 800
#define NJ 900
#define NK 1000
#define NL 1100
#define NM 1200
# endif
# ifdef EXTRALARGE_DATASET
#define NI 1600
#define NJ 1800
#define NK 2000
#define NL 2200
#define NM 2400
# endif
#endif


void init(int ni, int nk, int nj, int nm, int nl, vector<datatype> &A, vector<datatype> &B, vector<datatype> &C, vector<datatype> &D)
{
	int i, j;

	for (i = 0; i < ni; i++)
		for (j = 0; j < nk; j++)
			A[i*nk + j] = (double)((i*j + 1) % ni) / (5 * ni);
	for (i = 0; i < nk; i++)
		for (j = 0; j < nj; j++)
			B[i*nj + j] = (double)((i*(j + 1) + 2) % nj) / (5 * nj);
	for (i = 0; i < nj; i++)
		for (j = 0; j < nm; j++)
			C[i*nm + j] = (double)(i*(j + 3) % nl) / (5 * nl);
	for (i = 0; i < nm; i++)
		for (j = 0; j < nl; j++)
			D[i*nl + j] = (double)((i*(j + 2) + 2) % nk) / (5 * nk);
}


void multiplication(vector<datatype> &X, vector<datatype> &Y, vector<datatype> &XY, int n1, int n2, int n3)
{
	// (n1 n2)*(n2 n3) = (n1 n3)
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n3; j++)
		{
			for (int k = 0; k < n2; k++)
			{
				XY[i*n3 + j] += X[i*n2 + k] * Y[k*n3 + j];
			}
		}
	}
}


void print_array(vector<datatype> &X, int n1, int n2)
{
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			cout << X[i*n2 + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Russian");

	int ni = NI;
	int nk = NK;
	int nj = NJ;
	int nm = NM;
	int nl = NL;

	vector<datatype> A(ni*nk);
	vector<datatype> B(nk*nj);
	vector<datatype> AB(ni*nj);
	vector<datatype> C(nj*nm);
	vector<datatype> D(nm*nl);
	vector<datatype> CD(nj*nl);
	vector<datatype> ABCD(ni*nl);

	init(ni, nk, nj, nm, nl, A, B, C, D); // initialization A, B, C, D 
	multiplication(A, B, AB, ni, nk, nj);
	multiplication(C, D, CD, nj, nm, nl);
	multiplication(AB, CD, ABCD, ni, nj, nl);
	print_array(ABCD, ni, nl);
	/*
	int ni = 3;
	int nk = 2;
	int nj = 4;
	vector<datatype> A = {1, 5, 7, 7, 6, 9};
	vector<datatype> B = {1, 2, 3, 4, 5, 6, 3, 8};
	vector<datatype> AB(ni*nj);
	multiplication(A, B, AB, ni, nk, nj);
	print_array(AB, ni, nj);*/
	printf("jgf");
	cin.get();
	return 0;
}