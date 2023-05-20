#include "PseudoDouble.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "PseudoDouble_iostream.h"

using std::cout;
using std::endl;

static const int MAX_MATRIX_SIZE=10;

double determinant(double[][MAX_MATRIX_SIZE], int);
void cofactors(double[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int);
void trans(double[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int);
PseudoDouble determinant(PseudoDouble[][MAX_MATRIX_SIZE], int);
void cofactors(PseudoDouble[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int);
void trans(PseudoDouble[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int);

uint64_t get_thread_time() {
	timespec tvalue;
	clock_gettime(CLOCK_REALTIME, &tvalue);
	uint64_t start=tvalue.tv_sec;
	start=start*1000000000+tvalue.tv_nsec;
	return start;
}

double determinant(double a[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int k) {
	double s = 1, det = 0, b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
	int i, j, m, n, c;
	if (k == 1) {
		return (a[0][0]);
	} else {
		det = 0;
		for (c = 0; c < k; c++) {
			m = 0;
			n = 0;
			for (i = 0; i < k; i++) {
				for (j = 0; j < k; j++) {
					b[i][j] = 0;
					if (i != 0 && j != c) {
						b[m][n] = a[i][j];
						if (n < (k - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinant(b, k - 1));
			s = -1 * s;
		}
	}
	return (det);
}
 
void cofactors(double num[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int f) {
	double b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], fac[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
	int p, q, m, n, i, j;
	for (q = 0; q < f; q++) {
		for (p = 0; p < f; p++) {
			m = 0;
			n = 0;
			for (i = 0; i < f; i++) {
				for (j = 0; j < f; j++) {
					b[i][j] = 0;
					if (i != q && j != p) {
						b[m][n] = num[i][j];
						if (n < (f - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = (((q+p)&1)?-1:1) * determinant(b, f - 1);
		}
	}
	trans(num, inv, fac, f);
}
 
void trans(double num[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double fac[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int r) {
	int i, j;
	double b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], d;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			b[i][j] = fac[j][i];
		}
	}
	d = determinant(num, r);
	//inv[i][j] = 0;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			inv[i][j] = b[i][j] / d;
		}
	}
}

PseudoDouble determinant(PseudoDouble a[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int k) {
	PseudoDouble s = 1, det = 0, b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
	int i, j, m, n, c;
	if (k == 1) {
		return (a[0][0]);
	} else {
		det = 0;
		for (c = 0; c < k; c++) {
			m = 0;
			n = 0;
			for (i = 0; i < k; i++) {
				for (j = 0; j < k; j++) {
					b[i][j] = 0;
					if (i != 0 && j != c) {
						b[m][n] = a[i][j];
						if (n < (k - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinant(b, k - 1));
			s = PseudoDouble(-1) * s;
		}
	}
	return (det);
}
 
void cofactors(PseudoDouble num[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int f) {
	PseudoDouble b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], fac[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
	int p, q, m, n, i, j;
	for (q = 0; q < f; q++) {
		for (p = 0; p < f; p++) {
			m = 0;
			n = 0;
			for (i = 0; i < f; i++) {
				for (j = 0; j < f; j++) {
					b[i][j] = 0;
					if (i != q && j != p) {
						b[m][n] = num[i][j];
						if (n < (f - 2))
						       n++; else {
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = PseudoDouble((((q+p)&1)?-1:1)) * determinant(b, f - 1);
		}
	}
	trans(num, inv, fac, f);
}
 
void trans(PseudoDouble num[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], PseudoDouble fac[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int r) {
	int i, j;
	PseudoDouble b[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], d;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			b[i][j] = fac[j][i];
		}
	}
	d = determinant(num, r);
	//inv[i][j] = 0;
	for (i = 0; i < r; i++) {
		for (j = 0; j < r; j++) {
			inv[i][j] = b[i][j] / d;
		}
	}
}

int main(int, char**) {
	srand(0);
	double a[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], inv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], aa[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], d;
	PseudoDouble pa[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], pinv[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], paa[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], pd;
	int i, j, n=MAX_MATRIX_SIZE;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			a[i][j]=(double)rand()/ RAND_MAX;
			aa[i][j]=a[i][j];
			pa[i][j]=a[i][j];
			paa[i][j]=a[i][j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf(" %2.12f",a[i][j]);
		}
		printf("\n");
	}
	d = determinant(a, n);
	printf("\nTHE DETERMINANT IS=%2f\n", d);
	if (d == 0)
		printf("\nMATRIX IS NOT INVERSIBLE\n");
	else {
		const static int loop=20;
		uint64_t t0=get_thread_time();
		for(int i=0;i<loop;i++) {
			printf("%d\n",i);
			cofactors(a, inv, n);
			cofactors(inv,aa, n);
		}
		uint64_t t1=get_thread_time();
		for(int i=0;i<loop;i++) {
			printf("%d\n",i);
			cofactors(pa, pinv, n);
			cofactors(pinv,paa, n);
		}
		uint64_t t2=get_thread_time();
		cout << "Matrix time double=" << (t1-t0)*0.000000001 << endl;
		cout << "Matrix time pseudo=" << (t2-t1)*0.000000001 << endl;
		printf("\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << ' ' << a[i][j];
			}
			printf("\n");
		}
		printf("\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << ' ' << a[i][j]-aa[i][j];
			}
			printf("\n");
		}
		printf("\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << ' ' << pa[i][j];
			}
			printf("\n");
		}
		printf("\n");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << ' ' << pa[i][j]-paa[i][j];
			}
			printf("\n");
		}
	}
	{
		double a=1,b=0.3;
		PseudoDouble pa=1,pb=b;
		uint64_t t0=get_thread_time();
		for(int i=0;i<1000000000;i++) {
			a=a+(b*b);
			b=1-b;
		}
		PseudoDouble dec(0.1);
		uint64_t t1=get_thread_time();
		for(int i=0;i<1000000000;i++) {
			pa=pa+pb*pb;
			pb=PD_ONE-pb;
		}
		uint64_t t2=get_thread_time();
		cout << "Converge time nodiv double=" << (t1-t0)*0.000000001 << endl;
		cout << "Converge time nodiv pseudo=" << (t2-t1)*0.000000001 << endl;
		cout << a << ':' << b << endl;
		cout << pa << ':' << pb << endl;
	}
	{
		double a=1,b=1;
		PseudoDouble pa=1,pb=1;
		uint64_t t0=get_thread_time();
		for(int i=0;i<1000000000;i++) {
			a=a+1/(b*b);
			b=b-a*a;
		}
		uint64_t t1=get_thread_time();
		for(int i=0;i<1000000000;i++) {
			pa=pa+PD_ONE/(pb*pb);
			pb=pa-pa*pa;
		}
		uint64_t t2=get_thread_time();
		cout << "Converge time double=" << (t1-t0)*0.000000001 << endl;
		cout << "Converge time pseudo=" << (t2-t1)*0.000000001 << endl;
		cout << a << ':' << b << endl;
		cout << pa << ':' << pb << endl;
	}
	return 0;
}
