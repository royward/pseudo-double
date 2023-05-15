#include "PseudoFloat.h"
#include "PseudoFloat_iostream.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <iomanip>

using namespace std;

const static double NEAR_EXACT14=0.99999999999999;
const static double NEAR_EXACT13=0.9999999999999;
const static double NEAR_EXACT12=0.999999999999;
const static double NEAR_EXACT11=0.99999999999;
const static double NEAR_EXACT10=0.9999999999;
const static double NEAR_EXACT9 =0.999999999;
const static double NEAR_EXACT8 =0.99999999;
const static double NEAR_EXACT7 =0.9999999;
const static double NEAR_EXACT6 =0.999999;
const static double NEAR_EXACT5 =0.99999;
const static double NEAR_EXACT4 =0.9999;
const static double NEAR_EXACT3 =0.999;

void debug_pf_output(pseudo_float x) {
	signed_pf_internal mant=x&EXP_MASK_INV;
	int64_t exponent=x&EXP_MASK;
	printf("m=");
	for(int32_t i=PSEUDO_FLOAT_TOTAL_BITS-1;i>=PSEUDO_FLOAT_EXP_BITS;i--) {
		printf("%ld",(mant>>i)&1);
	}
	printf("=%lf, e=",pow(2,-PSEUDO_FLOAT_TOTAL_BITS)*mant);
	for(int32_t i=PSEUDO_FLOAT_EXP_BITS-1;i>=0;i--) {
		printf("%ld",(exponent>>i)&1);
	}
	printf("=%ld, n=%lf",exponent-PSEUDO_FLOAT_EXP_BIAS,pow(2,exponent-PSEUDO_FLOAT_EXP_BIAS-PSEUDO_FLOAT_TOTAL_BITS)*mant);
}

bool compare(double d1, double d2, double exactness) {
	if(d1>=0) {
		return (d1*exactness<=d2 && d2*exactness<=d1);
	} else {
		return (d1*exactness>=d2 && d2*exactness>=d1);
	}
}

uint64_t atan_rev_64_internal(uint64_t x);

int main() {
	uint32_t count=0;
	uint32_t failures=0;
	srand(0);
	vector<double> list;
	for(int i=-20;i<20;i++) {
		list.push_back(ldexp(1.0,i));
		list.push_back(-ldexp(1.0,i));
		list.push_back(ldexp(3.0,i));
		list.push_back(-ldexp(3.0,i));
		list.push_back(i);
		list.push_back(i+0.5);
	}
	for(int i=0;i<100;i++) {
		double r=(double)rand()/ RAND_MAX;
		double f=r*1000000.0;
		list.push_back(f);
		list.push_back(-f);
		f=r/1000000.0;
		list.push_back(f);
		list.push_back(-f);
	}
	for(uint32_t i=0;i<list.size();i++) {
		double f1=list[i];
		PseudoFloat pf1=f1;
		for(uint32_t j=0;j<list.size();j++) {
			double f2=list[j];
			PseudoFloat pf2=f2;
			double ff;
			ff=pf1+pf2;
			count++;
			if(!compare(f1+f2,ff,NEAR_EXACT8)) {
				failures++;
				cout << "add  " << f1 << '+' << f2 << "==" << f1+f2 << "!=" << ff << endl;
			}
			ff=pf1-pf2;
			count++;
			if(!compare(f1-f2,ff,NEAR_EXACT8)) {
				failures++;
				cout << "sub  " << f1 << '-' << f2 << "==" << f1-f2 << "!=" << ff << endl;
			}
			ff=pf1*pf2;
			count++;
			if(!compare(f1*f2,ff,NEAR_EXACT13)) {
				failures++;
				cout << "mult " << f1 << '*' << f2 << "==" << f1*f2 << "!=" << ff << endl;
			}
			if(f2!=0) {
				count++;
				ff=pf1/pf2;
				if(!compare(f1/f2,ff,NEAR_EXACT13)) {
					failures++;
					cout << "div  " << f1 << '/' << f2 << "==" << f1/f2 << "!=" << ff << endl;
				}
			}
			count++;
			if((f1<f2) != (pf1<pf2)) {
				failures++;
				cout << "comp " << f1 << "<" << f2 << "==" << (int)(f1<f2) << "!=" << (int)(pf1<pf2) << endl;
			}
			count++;
			if((f1<=f2) != (pf1<=pf2)) {
				failures++;
				cout << "comp " << f1 << "<=" << f2 << "==" << (int)(f1<=f2) << "!=" << (int)(pf1<=pf2) << endl;
			}
			count++;
			if((f1>f2) != (pf1>pf2)) {
				failures++;
				cout << "comp " << f1 << ">" << f2 << "==" << (int)(f1>f2) << "!=" << (int)(pf1>pf2) << endl;
			}
			count++;
			if((f1>=f2) != (pf1>=pf2)) {
				failures++;
				cout << "comp " << f1 << ">=" << f2 << "==" << (int)(f1>=f2) << "!=" << (int)(pf1>=pf2) << endl;
			}
			count++;
			ff=atan2(pf1,pf2);
			if((!compare(atan2(f1,f2),ff,NEAR_EXACT9)) && (fabs(f2/f1)<1e9 || !compare(atan2(f1,f2),ff,NEAR_EXACT4)) && (f1>=0 || f2!=0)) {
				failures++;
				cout << "atan2(" << f1 << ',' << f2 << ")==" << atan2(f1,f2) << "!=" << ff << endl;
			}
		}
	}
	for(int i=-20;i<20;i++) {
		list.push_back(i+0.4999999);
		list.push_back(i+0.5000001);
	}
	for(uint32_t i=0;i<list.size();i++) {
		double f=list[i];
		PseudoFloat pf=f;
		double ff=pf;
		count++;
		if(!compare(f,ff,NEAR_EXACT13)) {
			failures++;
			cout << "conv " << f << ' ' << ff << endl;
		}
		ff=-pf;
		count++;
		if(!compare(-f,ff,NEAR_EXACT13)) {
			failures++;
			cout << "neg  " << -f << ' ' << ff << endl;
		}
		if(f>0) {
			count++;
			ff=inv_sqrt(pf);
			if(!compare(1.0/sqrt(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "inv_sqrt  " << 1.0/sqrt(f) << ' ' << ff << endl;
			}
		}
		if(f>0) {
			count++;
			ff=sqrt(pf);
			if(!compare(sqrt(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "sqrt  " << sqrt(f) << ' ' << ff << endl;
			}
		}
		ff=floor(pf);
		count++;
		if(!compare(floor(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "floor  " << floor(f) << ' ' << ff << endl;
		}
		ff=ceil(pf);
		count++;
		if(!compare(ceil(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "ceil  " << ceil(f) << ' ' << ff << endl;
		}
		ff=round(pf);
		count++;
		if(!compare(round(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "round  " << round(f) << ' ' << ff << endl;
		}
		if(!isinf(exp2(f)) && f!=-1024.0) {
			count++;
			ff=exp2(pf);
			if(!compare(exp2(f),ff,NEAR_EXACT12)) {
				failures++;
				cout << "exp2  " << setprecision(15) << f << ' ' << setprecision(15) << exp2(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(!isinf(exp(f))) {
			ff=exp(pf);
			count++;
			if(!compare(exp(f),ff,NEAR_EXACT11)) {
				failures++;
				cout << "exp  " << setprecision(15) << f << ' ' << setprecision(15) << exp(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(f>0) {
			ff=log2(pf);
			count++;
			if(!compare(log2(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log2  " << setprecision(15) << f << ' ' << setprecision(15) << log2(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(f>0) {
			ff=log(pf);
			count++;
			if(!compare(log(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log  " << setprecision(15) << f << ' ' << setprecision(15) << log(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(f>0) {
			ff=log10(pf);
			count++;
			if(!compare(log10(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log10  " << setprecision(15) << f << ' ' << setprecision(15) << log10(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(-10000.0<f && f<10000.0) {
			ff=sin(pf);
			count++;
			if(!compare(sin(f),ff,NEAR_EXACT10)) {
				failures++;
				cout << "sin  " << f << ' ' << sin(f) << ' ' << ff << endl;
			}
			ff=cos(pf);
			count++;
			if(!compare(cos(f),ff,NEAR_EXACT10)) {
				failures++;
				cout << "cos  " << f << ' ' << cos(f) << ' ' << ff << endl;
			}
		}
	}
	for(uint32_t i=0;i<list.size();i++) {
		double f1=list[i];
		PseudoFloat pf1=f1;
		for(uint32_t j=0;j<list.size();j++) {
			double f2=list[j];
			if(f1>0 && !isinf(pow(f1,f2)) && pow(f1,f2)>1e-300) {
				PseudoFloat pf2=f2;
				double ff;
				ff=pow(pf1,pf2);
				count++;
				if(!compare(pow(f1,f2),ff,NEAR_EXACT9)) {
					failures++;
					cout << "pow(" << f1 << ',' << f2 << ")==" << pow(f1,f2) << "!=" << ff << endl;
				}
			}
		}
	}
	for(int64_t i=-1000;i<1000;i++) {
		PseudoFloat pf(i);
		int64_t ii=pf;
		count++;
		if(i!=ii) {
			failures++;
			cout << "sint convert" << i << ' ' << ii << endl;
		}
	}
	for(uint64_t i=0;i<1000;i++) {
		PseudoFloat pf(i);
		uint64_t ii=pf;
		count++;
		if(i!=ii) {
			failures++;
			cout << "uint convert" << i << ' ' << ii << endl;
		}
	}
	cout << "Tests done, passed " << (count-failures) << '/' << count << endl;

// This can be modified to test indivual internal functions
// 	for(uint64_t i=0x4000000000000000ULL;i<=0xFF00000000000000ULL;i+=0x40000000000000ULL) {
// 		cout << hex << i << ' ' << inv_sqrt64_internal(i) << dec << endl;
// 	}
	return 0;
}