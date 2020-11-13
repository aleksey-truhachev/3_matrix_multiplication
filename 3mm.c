/* Include benchmark-specific header. */
#include "3mm.h"

typedef datatype double;

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



int main(int argc, char** argv)
{
	int ni = NI;
	int nk = NK;
	int nj = NJ;
	int nm = NM;
	int nl = NL;

	vector<datatype> A[ni*nk];
	vector<datatype> B[nk*nj];
	vector<datatype> AB[ni*nj];
	vector<datatype> C[nj*nm];
	vector<datatype> D[nm*nl];
	vector<datatype> CD[nj*nl];
	vector<datatype> F[ni*nl];
	return 0;
}
