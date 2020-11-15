#include <iostream>
#include "math.h"
#include <vector>
#include "omp.h"
using namespace std;


typedef double datatype;


# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET) && !defined(EXTRAMINI_DATASET)
#define LARGE_DATASET
# endif
# if !defined(NI) && !defined(NJ) && !defined(NK) && !defined(NL) && !defined(NM)
# ifdef EXTRAMINI_DATASET
#define NI 3
#define NJ 3
#define NK 3
#define NL 2
#define NM 2
# endif
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


void init_test(int ni, int nk, int nj, int nm, int nl, vector<datatype> &A, vector<datatype> &B, vector<datatype> &C, vector<datatype> &D)
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

void transpose(vector<datatype> &X_transp, vector<datatype> &X, int n1, int n2)
{
	// X = (n1 n2)
	// X_transp = (n2 n1)
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			X_transp[j*n1 + i] = X[i*n2 + j];
		}
	}
}


void transpose_parall(vector<datatype> &X_transp, vector<datatype> &X, int n1, int n2)
{
	// X = (n1 n2)
	// X_transp = (n2 n1)
#pragma omp parallel for
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			X_transp[j*n1 + i] = X[i*n2 + j];
		}
	}
}


void multiplication_transp(vector<datatype> &X, vector<datatype> &Y_transp, vector<datatype> &XY, int n1, int n2, int n3)
{
	// X = (n1 n2)
	// Y = (n2 n3)
	// Y_transp = (n3 n2)
	// XY = (n1 n3)
	//перемножаем по стокам. XY[i][j] += X[i][k] * Y_transp[j][k]
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n3; j++)
		{
			for (int k = 0; k < n2; k++)
			{
				XY[i*n3 + j] += X[i*n2 + k] * Y_transp[j*n2 + k];
			}
		}
	}
}


void multiplication_transp_parall(vector<datatype> &X, vector<datatype> &Y_transp, vector<datatype> &XY, int n1, int n2, int n3)
{
	// X = (n1 n2)
	// Y = (n2 n3)
	// Y_transp = (n3 n2)
	// XY = (n1 n3)
	//перемножаем по стокам. XY[i][j] += X[i][k] * Y_transp[j][k]
#pragma omp parallel for
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n3; j++)
		{
			for (int k = 0; k < n2; k++)
			{
				XY[i*n3 + j] += X[i*n2 + k] * Y_transp[j*n2 + k];
			}
		}
	}
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
	vector<datatype> B_transp(nj*nk); // B transposed
	vector<datatype> AB(ni*nj);
	vector<datatype> C(nj*nm);
	vector<datatype> D(nm*nl);
	vector<datatype> D_transp(nl*nm); // D transposed
	vector<datatype> CD(nj*nl);
	vector<datatype> CD_transp(nl*nj); // CD transposed
	vector<datatype> ABCD(ni*nl);
	vector<datatype> ABCD1(ni*nl);

	double t1, t2;

	init(ni, nk, nj, nm, nl, A, B, C, D); // initialization A, B, C, D 
	t1 = omp_get_wtime();
	transpose(B_transp, B, nk, nj);
	multiplication_transp(A, B_transp, AB, ni, nk, nj);
	/*cout << "A" << endl;
	print_array(A, ni, nk);
	cout << "B" << endl;
	print_array(B, nk, nj);
	cout << "AB" << endl;
	print_array(AB, ni, nj);*/
	transpose(D_transp, D, nm, nl);
	multiplication_transp(C, D_transp, CD, nj, nm, nl);
	/*cout << "C" << endl;
	print_array(C, nj, nm);
	cout << "D" << endl;
	print_array(D, nm, nl);
	cout << "CD" << endl;
	print_array(CD, nj, nl);*/
	transpose(CD_transp, CD, nj, nl);
	multiplication_transp(AB, CD_transp, ABCD, ni, nj, nl);
	t2 = omp_get_wtime();
	cout << "С транспонированием: " << t2 - t1 << " секунд (посчитано верно)" << endl; //Выполнение для LARGE_DATASET за  2,9 - 3 секунды
	

	init(ni, nk, nj, nm, nl, A, B, C, D); // initialization A, B, C, D 
	t1 = omp_get_wtime();
	transpose_parall(B_transp, B, nk, nj);
	multiplication_transp_parall(A, B_transp, AB, ni, nk, nj);
	/*cout << "A" << endl;
	print_array(A, ni, nk);
	cout << "B" << endl;
	print_array(B, nk, nj);
	cout << "AB" << endl;
	print_array(AB, ni, nj);*/
	transpose_parall(D_transp, D, nm, nl);
	multiplication_transp_parall(C, D_transp, CD, nj, nm, nl);
	/*cout << "C" << endl;
	print_array(C, nj, nm);
	cout << "D" << endl;
	print_array(D, nm, nl);
	cout << "CD" << endl;
	print_array(CD, nj, nl);*/
	transpose_parall(CD_transp, CD, nj, nl);
	multiplication_transp_parall(AB, CD_transp, ABCD1, ni, nj, nl);
	t2 = omp_get_wtime();
	cout << "С параллельным транспонированием: " << t2 - t1 << " секунд" << endl; //Выполнение для LARGE_DATASET за  0,85  секунды

	int counter = 0;
	for (int i = 0; i < ABCD1.size(); i++)
	{
		if (abs(ABCD1[i] - 4*ABCD[i]) > 1.e-8)
		{
			counter++;
		}
	}
	cout << counter << endl;


	//init(ni, nk, nj, nm, nl, A, B, C, D); // initialization A, B, C, D 
	//t1 = omp_get_wtime();
	//multiplication(A, B, AB, ni, nk, nj);
	///*cout << "A" << endl;
	//print_array(A, ni, nk);
	//cout << "B" << endl;
	//print_array(B, nk, nj);
	//cout << "AB" << endl;
	//print_array(AB, ni, nj);*/
	//multiplication(C, D, CD, nj, nm, nl);
	///*cout << "C" << endl;
	//print_array(C, nj, nm);
	//cout << "D" << endl;
	//print_array(D, nm, nl);
	//cout << "CD" << endl;
	//print_array(CD, nj, nl);*/
	//multiplication(AB, CD, ABCD1, ni, nj, nl);
	//t2 = omp_get_wtime();
	////cout << t2 - t1 << endl; //Выполнение для LARGE_DATASET за 22,5 - 23,1 секунды
	//cout << "Без транспонирования: " << t2 - t1 << " секунд" << endl;
	////print_array(ABCD1, ni, nl);

	//int counter = 0;
	//for (int i = 0; i < ABCD1.size(); i++)
	//{
	//	if (abs(ABCD1[i] - 4*ABCD[i]) > 1.e-8)
	//	{
	//		counter++;
	//	}
	//}
	//cout << counter << endl;

	








	//int ni = 3;
	//int nk = 4;
	//int nj = 3;
	//int nm = 2;
	//int nl = 3;

	//vector<datatype> A = { 1,2,3,4,5,6,7,8,9,10,11,12 };
	//vector<datatype> B = { 12,11,10,9,8,7,6,5,4,3,2,1 };
	//vector<datatype> B_transp(nj*nk); // B transposed
	//vector<datatype> AB(ni*nj);
	//vector<datatype> C = {1,2,3,4, 5, 6};
	//vector<datatype> D = {6,5,4,3,2,1};
	//vector<datatype> D_transp(nl*nm); // D transposed
	//vector<datatype> CD(nj*nl);
	//vector<datatype> CD_transp(nl*nj); // CD transposed
	//vector<datatype> ABCD(ni*nl);

	//double t1, t2;

	//t1 = omp_get_wtime();
	//transpose(B_transp, B, nk, nj);
	//multiplication_transp(A, B_transp, AB, ni, nk, nj);
	//print_array(AB, ni, nj);
	//transpose(D_transp, D, nm, nl);
	//multiplication_transp(C, D_transp, CD, nj, nm, nl);
	//print_array(CD, nj, nl);
	//transpose(CD_transp, CD, nj, nl);
	//multiplication_transp(AB, CD_transp, ABCD, ni, nj, nl);
	//t2 = omp_get_wtime();
	//cout << "С транспонированием: " << t2 - t1 << " секунд" << endl; //Выполнение для LARGE_DATASET за  8.72 - 8.74  секунды
	//print_array(ABCD, ni, nl);


	////vector<datatype> ABCD1(ni*nl);
	////t1 = omp_get_wtime();
	////multiplication(A, B, AB, ni, nk, nj);
	////print_array(AB, ni, nj);
	////multiplication(C, D, CD, nj, nm, nl);
	////print_array(CD, nj, nl);
	////multiplication(AB, CD, ABCD1, ni, nj, nl);
	////t2 = omp_get_wtime();
	////cout << "Без транспонирования: " << t2 - t1 << " секунд" << endl;
	////print_array(ABCD1, ni, nl);

	///*int counter = 0;
	//for (int i = 0; i < ABCD1.size(); i++)
	//{
	//	if (abs(ABCD1[i] - ABCD[i]) > 1.e-7)
	//	{
	//		counter++;
	//	}
	//}
	//cout << counter << endl;*/

	cin.get();
	return 0;
}