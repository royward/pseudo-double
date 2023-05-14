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
			if(!compare(f1+f2,ff,NEAR_EXACT8)) {
				cout << "add  " << f1 << '+' << f2 << "==" << f1+f2 << "!=" << ff << endl;
				PseudoFloat pff=pf1+pf2;
				debug_pf_output(pf1.get_internal());
				cout << endl;
				debug_pf_output(pf2.get_internal());
				cout << endl;
				debug_pf_output(pff.get_internal());
				cout << endl;
				__asm__("int3");
				double ff=pf1+pf2;
				cout << ff << endl;
			}
			ff=pf1-pf2;
			if(!compare(f1-f2,ff,NEAR_EXACT8)) {
				cout << "sub  " << f1 << '-' << f2 << "==" << f1-f2 << "!=" << ff << endl;
			}
			ff=pf1*pf2;
			if(!compare(f1*f2,ff,NEAR_EXACT13)) {
				cout << "mult " << f1 << '*' << f2 << "==" << f1*f2 << "!=" << ff << endl;
			}
			if(f2!=0) {
				ff=pf1/pf2;
				if(!compare(f1/f2,ff,NEAR_EXACT13)) {
					cout << "div  " << f1 << '/' << f2 << "==" << f1/f2 << "!=" << ff << endl;
				}
			}
			if((f1<f2) != (pf1<pf2)) {
				cout << "comp " << f1 << "<" << f2 << "==" << (int)(f1<f2) << "!=" << (int)(pf1<pf2) << endl;
			}
			if((f1<=f2) != (pf1<=pf2)) {
				cout << "comp " << f1 << "<=" << f2 << "==" << (int)(f1<=f2) << "!=" << (int)(pf1<=pf2) << endl;
			}
			if((f1>f2) != (pf1>pf2)) {
				cout << "comp " << f1 << ">" << f2 << "==" << (int)(f1>f2) << "!=" << (int)(pf1>pf2) << endl;
// 				debug_pf_output(pf1.get_internal());
// 				cout << endl;
// 				debug_pf_output(pf2.get_internal());
// 				cout << endl;
// 				__asm__("int3");
// 				bool t=(pf1>pf2);
// 				cout << t << endl;
			}
			if((f1>=f2) != (pf1>=pf2)) {
				cout << "comp " << f1 << ">=" << f2 << "==" << (int)(f1>=f2) << "!=" << (int)(pf1>=pf2) << endl;
			}
			ff=atan2(pf1,pf2);
			if((!compare(atan2(f1,f2),ff,NEAR_EXACT9)) && (fabs(f2/f1)<1e9 || !compare(atan2(f1,f2),ff,NEAR_EXACT4)) && (f1>=0 || f2!=0)) {
				cout << "atan2(" << f1 << ',' << f2 << ")==" << atan2(f1,f2) << "!=" << ff << endl;
			}
		}
	}
// 	for(int32_t i1=-2;i1<=2;i1++) {
// 		for(int32_t i2=-2;i2<=2;i2++) {
// 			PseudoFloat pf1=i1;
// 			PseudoFloat pf2=i2;
// 			double ff=atan2(pf1,pf2);
// 			if(!compare(atan2((double)i1,(double)i2),ff,NEAR_EXACT13)) {
// 				cout << "atan2(" << i1 << ',' << i2 << ")==" << atan2((double)i1,(double)i2) << "!=" << ff << endl;
// 				__asm__("int3");
// 				atan2(pf1,pf2);
// 			}
// 		}
// 	}
	for(int i=-20;i<20;i++) {
		list.push_back(i+0.4999999);
		list.push_back(i+0.5000001);
	}
	for(uint32_t i=0;i<list.size();i++) {
		double f=list[i];
		PseudoFloat pf=f;
		double ff=pf;
		if(!compare(f,ff,NEAR_EXACT13)) {
			cout << "conv " << f << ' ' << ff << endl;
		}
		ff=-pf;
		if(!compare(-f,ff,NEAR_EXACT13)) {
			cout << "neg  " << -f << ' ' << ff << endl;
		}
		if(f>0) {
			ff=inv_sqrt(pf);
			if(!compare(1.0/sqrt(f),ff,NEAR_EXACT13)) {
				cout << "inv_sqrt  " << 1.0/sqrt(f) << ' ' << ff << endl;
// 				__asm__("int3");
// 				ff=inv_sqrt(pf);
			}
		}
		if(f>0) {
			ff=sqrt(pf);
			if(!compare(sqrt(f),ff,NEAR_EXACT13)) {
				cout << "sqrt  " << sqrt(f) << ' ' << ff << endl;
// 				__asm__("int3");
// 				ff=sqrt(pf);
			}
		}
		ff=floor(pf);
		if(!compare(floor(f),ff,NEAR_EXACT13)) {
			cout << "floor  " << floor(f) << ' ' << ff << endl;
			__asm__("int3");
			ff=floor(pf);
		}
		ff=ceil(pf);
		if(!compare(ceil(f),ff,NEAR_EXACT13)) {
			cout << "ceil  " << ceil(f) << ' ' << ff << endl;
// 			__asm__("int3");
// 			ff=ceil(pf);
		}
		ff=round(pf);
		if(!compare(round(f),ff,NEAR_EXACT13)) {
			cout << "round  " << round(f) << ' ' << ff << endl;
// 				__asm__("int3");
// 				ff=round(pf);
		}
		//cout << f << ':' << pf.get_internal() << endl;
		if(!isinf(exp2(f)) && f!=-1024.0) {
			ff=exp2(pf);
			if(!compare(exp2(f),ff,NEAR_EXACT12)) {
				cout << "exp2  " << setprecision(15) << f << ' ' << setprecision(15) << exp2(f) << ' ' << setprecision(15) << ff << endl;
// 				__asm__("int3");
// 				pf=exp2(pf);
// 				f=pf;
			}
		}
		if(!isinf(exp(f))) {
			ff=exp(pf);
			if(!compare(exp(f),ff,NEAR_EXACT11)) {
				cout << "exp  " << setprecision(15) << f << ' ' << setprecision(15) << exp(f) << ' ' << setprecision(15) << ff << endl;
// 				__asm__("int3");
// 				pf=exp(pf);
// 				f=pf;
			}
		}
		if(f>0) {
			ff=log2(pf);
			if(!compare(log2(f),ff,NEAR_EXACT13)) {
				cout << "log2  " << setprecision(15) << f << ' ' << setprecision(15) << log2(f) << ' ' << setprecision(15) << ff << endl;
// 				__asm__("int3");
// 				pf=log2(pf);
// 				f=pf;
			}
		}
		if(f>0) {
			ff=log(pf);
			if(!compare(log(f),ff,NEAR_EXACT13)) {
				cout << "log  " << setprecision(15) << f << ' ' << setprecision(15) << log(f) << ' ' << setprecision(15) << ff << endl;
// 				__asm__("int3");
// 				pf=log(pf);
// 				f=pf;
			}
		}
		if(f>0) {
			ff=log10(pf);
			if(!compare(log10(f),ff,NEAR_EXACT13)) {
				cout << "log10  " << setprecision(15) << f << ' ' << setprecision(15) << log10(f) << ' ' << setprecision(15) << ff << endl;
// 				__asm__("int3");
// 				pf=log10(pf);
// 				f=pf;
			}
		}
		if(-10000.0<f && f<10000.0) {
			ff=sin(pf);
			if(!compare(sin(f),ff,NEAR_EXACT10)) {
				cout << "sin  " << f << ' ' << sin(f) << ' ' << ff << endl;
			}
			ff=cos(pf);
			if(!compare(cos(f),ff,NEAR_EXACT10)) {
				cout << "cos  " << f << ' ' << cos(f) << ' ' << ff << endl;
					__asm__("int3");
					ff=cos(pf);
			}
		}
		//pf_log2(pf.get_internal());
// 		ff=exp2(pf);
// 		cout << f << ' ' << exp2(f) << ' ' << ff << endl;
// 		exp2(pf);
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
				if(!compare(pow(f1,f2),ff,NEAR_EXACT9)) {
					cout << "pow(" << f1 << ',' << f2 << ")==" << pow(f1,f2) << "!=" << ff << endl;
// 					PseudoFloat pff=pow(pf1,pf2);
// 					debug_pf_output(pf1.get_internal());
// 					cout << endl;
// 					debug_pf_output(pf2.get_internal());
// 					cout << endl;
// 					debug_pf_output(pff.get_internal());
// 					cout << endl;
// 					__asm__("int3");
// 					double ff=pow(pf1,pf2);
// 					cout << ff << endl;
				}
			}
		}
	}
	for(int64_t i=-1000;i<1000;i++) {
		PseudoFloat pf(i);
		int64_t ii=pf;
		if(i!=ii) {
			cout << "sint convert" << i << ' ' << ii << endl;
		}
	}
	for(uint64_t i=0;i<1000;i++) {
		PseudoFloat pf(i);
		uint64_t ii=pf;
		if(i!=ii) {
			cout << "uint convert" << i << ' ' << ii << endl;
		}
	}
	cout << "Tests done" << endl;
// 	for(uint32_t i=0;i<21;i++) {
// 		cout << i << ' ' << ldexp(pf_atan2_rev(sint64_to_pf(i),sint64_to_pf(20)),-62) << endl;
// 	}
// 	for(uint32_t i=0;i<257;i++) {
// 		uint64_t x=0x3FFFFFFFFFFFFF00ULL+i;
// 		uint64_t p0=atan_rev_64_internal(x);
// 		cout << hex << x << ' ' << p0 << endl;
// 	}
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL*0.0),-62) << ' ' << atan(0.0)*4/M_PI << endl;
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL*0.2),-62) << ' ' << atan(0.2)*4/M_PI << endl;
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL*0.4),-62) << ' ' << atan(0.4)*4/M_PI << endl;
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL*0.6),-62) << ' ' << atan(0.6)*4/M_PI << endl;
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL*0.8),-62) << ' ' << atan(0.8)*4/M_PI << endl;
// 	cout << ldexp(atan_rev_64_internal(0x4000000000000000ULL),-62)     << ' ' << atan(1.0)*4/M_PI << endl;
// 	for(int i=-500;i<500;i++) {
// 		double d=i/100.0;
// 		PseudoFloat dd=d;
// 		cout << d << ' ' << sin(dd) << endl;
// 	}
// 	cout << '$' << endl;
// 	cout << PF_PI << endl;
// 	cout << PF_INV_TAU << endl;
// 	cout << PF_get_fixed2(PF_PI,10) << endl;
// 	cout << setprecision(12) << PF_create_fixed2(1234567,-10) << endl;
// 	cout << setprecision(12) << PF_get_fixed2(PF_create_fixed2(1234567,-10),-10) << endl;

// C/C++: uint64_t inv_sqrt64_internal(uint64_t x);
// calculate 1/sqrt(x)
// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
// 
// C/C++: uint64_t exp2_64_internal(uint64_t x);
// calculate 2^x
// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
// 


	for(uint64_t i=0x4000000000000000ULL;i<=0xFF00000000000000ULL;i+=0x40000000000000ULL) {
		cout << hex << i << ' ' << inv_sqrt64_internal(i) << dec << endl;
	}
	return 0;
}
