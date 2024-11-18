// BSD 3-Clause License
// 
// Copyright (c) 2023, Roy Ward
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "PseudoDouble.h"
#include "PseudoDouble_iostream.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <iomanip>

using std::cout;
using std::endl;
using std::isinf;
using std::setprecision;
using namespace pseudodouble;

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

void debug_pd_output(pseudo_double_i x) {
	signed_pd_internal mant=x&EXP_MASK_INV;
	int64_t exponent=x&EXP_MASK;
	printf("m=");
	for(int32_t i=PSEUDO_DOUBLE_TOTAL_BITS-1;i>=PSEUDO_DOUBLE_EXP_BITS;i--) {
#ifdef _MSC_VER
		printf("%lld", (mant >> i) & 1);
#else
		printf("%ld", (mant >> i) & 1);
#endif
	}
	printf("=%lf, e=",pow(2,-PSEUDO_DOUBLE_TOTAL_BITS)*mant);
	for(int32_t i=PSEUDO_DOUBLE_EXP_BITS-1;i>=0;i--) {
#ifdef _MSC_VER
		printf("%lld", (exponent >> i) & 1);
#else
		printf("%ld", (exponent >> i) & 1);
#endif
	}
#ifdef _MSC_VER
	printf("=%lld, n=%lf", exponent - PSEUDO_DOUBLE_EXP_BIAS, pow(2, exponent - PSEUDO_DOUBLE_EXP_BIAS - PSEUDO_DOUBLE_TOTAL_BITS) * mant);
#else
	printf("=%ld, n=%lf", exponent - PSEUDO_DOUBLE_EXP_BIAS, pow(2, exponent - PSEUDO_DOUBLE_EXP_BIAS - PSEUDO_DOUBLE_TOTAL_BITS) * mant);
#endif
}

bool compare(double d1, double d2, double exactness) {
	if(d1>=0) {
		return (d1*exactness<=d2 && d2*exactness<=d1);
	} else {
		return (d1*exactness>=d2 && d2*exactness>=d1);
	}
}

int main() {
	// PseudoDouble aa=0.99999999;
	// cout << inv_sqrt(aa) << endl;
	// PseudoDouble bb=0.25000001;
	// cout << inv_sqrt(bb) << endl;
	// asm("int3");
	uint32_t count=0;
	uint32_t failures=0;
	srand(0);
	std::vector<double> list;
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
		PseudoDouble pd1=f1;
		for(uint32_t j=0;j<list.size();j++) {
			double f2=list[j];
			PseudoDouble pd2=f2;
			double ff;
			ff=pd1+pd2;
			count++;
			if(!compare(f1+f2,ff,NEAR_EXACT8)) {
				failures++;
				cout << "add  " << f1 << '+' << f2 << "==" << f1+f2 << "!=" << ff << endl;
			}
			ff=pd1-pd2;
			count++;
			if(!compare(f1-f2,ff,NEAR_EXACT8)) {
				failures++;
				cout << "sub  " << f1 << '-' << f2 << "==" << f1-f2 << "!=" << ff << endl;
			}
			ff=pd1*pd2;
			count++;
			if(!compare(f1*f2,ff,NEAR_EXACT13)) {
				failures++;
				cout << "mult " << f1 << '*' << f2 << "==" << f1*f2 << "!=" << ff << endl;
			}
			count++;
			if((pseudodouble::max(pd1,pd2)==pd1) != (std::max(f1,f2)==f1)) {
				failures++;
				cout << "difference in max(" << f1 << ',' << f2 << ')' << endl;
			}
			count++;
			if((pseudodouble::min(pd1,pd2)==pd1) != (std::min(f1,f2)==f1)) {
				failures++;
				cout << "difference in min(" << f1 << ',' << f2 << ')' << endl;
			}
			if(f2!=0) {
				count++;
				ff=pd1/pd2;
				if(!compare(f1/f2,ff,NEAR_EXACT13)) {
					failures++;
					cout << "div  " << f1 << '/' << f2 << "==" << f1/f2 << "!=" << ff << endl;
				}
			}
			count++;
			if((f1<f2) != (pd1<pd2)) {
				failures++;
				cout << "comp " << f1 << "<" << f2 << "==" << (int)(f1<f2) << "!=" << (int)(pd1<pd2) << endl;
			}
			count++;
			if((f1<=f2) != (pd1<=pd2)) {
				failures++;
				cout << "comp " << f1 << "<=" << f2 << "==" << (int)(f1<=f2) << "!=" << (int)(pd1<=pd2) << endl;
			}
			count++;
			if((f1>f2) != (pd1>pd2)) {
				failures++;
				cout << "comp " << f1 << ">" << f2 << "==" << (int)(f1>f2) << "!=" << (int)(pd1>pd2) << endl;
			}
			count++;
			if((f1>=f2) != (pd1>=pd2)) {
				failures++;
				cout << "comp " << f1 << ">=" << f2 << "==" << (int)(f1>=f2) << "!=" << (int)(pd1>=pd2) << endl;
			}
			count++;
			ff=atan2(pd1,pd2);
			if((!compare(atan2(f1,f2),ff,NEAR_EXACT9)) && (fabs(f2/f1)<1e9 || !compare(atan2(f1,f2),ff,NEAR_EXACT3)) && (f1>=0 || f2!=0)) {
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
		PseudoDouble pd=f;
		double ff=pd;
		count++;
		if(!compare(f,ff,NEAR_EXACT13)) {
			failures++;
			cout << "conv " << f << ' ' << ff << endl;
		}
		ff=-pd;
		count++;
		if(!compare(-f,ff,NEAR_EXACT13)) {
			failures++;
			cout << "neg  " << -f << ' ' << ff << endl;
		}
		if(f>0) {
			count++;
			ff=inv_sqrt(pd);
			if(!compare(1.0/sqrt(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "inv_sqrt  " << 1.0/sqrt(f) << ' ' << ff << endl;
			}
		}
		if(f>0) {
			count++;
			ff=sqrt(pd);
			if(!compare(sqrt(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "sqrt  " << sqrt(f) << ' ' << ff << endl;
			}
		}
		ff=floor(pd);
		count++;
		if(!compare(floor(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "floor  " << floor(f) << ' ' << ff << endl;
		}
		ff=ceil(pd);
		count++;
		if(!compare(ceil(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "ceil  " << ceil(f) << ' ' << ff << endl;
		}
		ff=round(pd);
		count++;
		if(!compare(round(f),ff,NEAR_EXACT13)) {
			failures++;
			cout << "round  " << round(f) << ' ' << ff << endl;
		}
		if(!isinf(exp2(f))/* && f!=-1024.0*/) {
			if(f<128 && f>-128) {
				count++;
				ff=exp2(pd);
				if(!compare(exp2(f),ff,NEAR_EXACT12)) {
					failures++;
					cout << "exp2  " << setprecision(15) << f << ' ' << setprecision(15) << exp2(f) << ' ' << setprecision(15) << ff << endl;
				}
			}
		}
		if(!isinf(exp(f))) {
			if(f<96 && f>-96) {
				ff=exp(pd);
				count++;
				if(!compare(exp(f),ff,NEAR_EXACT11)) {
					failures++;
					cout << "exp  " << setprecision(15) << f << ' ' << setprecision(15) << exp(f) << ' ' << setprecision(15) << ff << endl;
				}
			}
		}
		if(f>0) {
			ff=log2(pd);
			count++;
			if(!compare(log2(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log2  " << setprecision(15) << f << ' ' << setprecision(15) << log2(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(f>0) {
			ff=log(pd);
			count++;
			if(!compare(log(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log  " << setprecision(15) << f << ' ' << setprecision(15) << log(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(f>0) {
			ff=log10(pd);
			count++;
			if(!compare(log10(f),ff,NEAR_EXACT13)) {
				failures++;
				cout << "log10  " << setprecision(15) << f << ' ' << setprecision(15) << log10(f) << ' ' << setprecision(15) << ff << endl;
			}
		}
		if(-10000.0<f && f<10000.0) {
			ff=sin(pd);
			count++;
			if(!compare(sin(f),ff,NEAR_EXACT9)) {
				failures++;
				cout << "sin  " << f << ' ' << sin(f) << ' ' << ff << endl;
			}
			ff=cos(pd);
			count++;
			if(!compare(cos(f),ff,NEAR_EXACT9)) {
				failures++;
				cout << "cos  " << f << ' ' << cos(f) << ' ' << ff << endl;
			}
		}
	}
	for(uint32_t i=0;i<list.size();i++) {
		double f1=list[i];
		PseudoDouble pd1=f1;
		for(uint32_t j=0;j<list.size();j++) {
			double f2=list[j];
			if(f1>0 && !isinf(pow(f1,f2)) && pow(f1,f2)>1e-35 && pow(f1,f2)<1e35) {
				PseudoDouble pd2=f2;
				double ff;
				ff=pow(pd1,pd2);
				count++;
				if(!compare(pow(f1,f2),ff,NEAR_EXACT9)) {
					failures++;
					cout << "pow(" << f1 << ',' << f2 << ")==" << pow(f1,f2) << "!=" << ff << endl;
				}
			}
		}
	}
	for(int64_t i=-1000;i<1000;i++) {
		PseudoDouble pd(i);
		int64_t ii=pd;
		count++;
		if(i!=ii) {
			failures++;
			cout << "sint convert" << i << ' ' << ii << endl;
		}
	}
	for(uint64_t i=0;i<1000;i++) {
		PseudoDouble pd(i);
		uint64_t ii=pd;
		count++;
		if(i!=ii) {
			failures++;
			cout << "uint convert" << i << ' ' << ii << endl;
		}
	}
	cout << "Tests done, passed " << (count-failures) << '/' << count << endl;


	// Code examples from README.md
	{
		double a=0.3;
		double b=-4.0;
		double c=6.0;
		double disc=sqrt(b*b-4.0*a*c);
		double sol1=(-b-disc)/(2.0*a);
		double sol2=(-b+disc)/(2.0*a);
		printf("C: Solution 1 = %lf\n",sol1);
		printf("C: Solution 2 = %lf\n",sol2);
	}
	{
		PseudoDouble a=PD_create_fixed10(3,-1); // 0.3
		PseudoDouble b=-4;
		PseudoDouble c=6;
		PseudoDouble disc=sqrt(b*b-PseudoDouble(4)*a*c);
		PseudoDouble sol1=(-b-disc)/(PseudoDouble(2)*a);
		PseudoDouble sol2=(-b+disc)/(PseudoDouble(2)*a);
		std::cout << "C++: Solution 1 = " << sol1 << std::endl;
		std::cout << "C++: Solution 2 = " << sol2 << std::endl;
	}
	{
		pseudo_double a=int64fixed10_to_pd(3,-1); // 0.3
		pseudo_double b=int64_to_pd(-4);
		pseudo_double c=int64_to_pd(6);
		pseudo_double disc=pd_sqrt(pd_sub(pd_mult(b,b),pd_mult(pd_mult(int64_to_pd(4),a),c)));
		pseudo_double sol1=pd_div(pd_sub(pd_neg(b),disc),pd_mult(int64_to_pd(2),a));
		pseudo_double sol2=pd_div(pd_add(pd_neg(b),disc),pd_mult(int64_to_pd(2),a));
		printf("C: Solution 1 = %lf\n",pd_to_double(sol1));
		printf("C: Solution 2 = %lf\n",pd_to_double(sol2));
	}
	{
		pseudo_double_i a=int64fixed10_to_pdi(3,-1); // 0.3
		pseudo_double_i b=int64_to_pdi(-4);
		pseudo_double_i c=int64_to_pdi(6);
		pseudo_double_i disc=pdi_sqrt(pdi_sub(pdi_mult(b,b),pdi_mult(pdi_mult(int64_to_pdi(4),a),c)));
		pseudo_double_i sol1=pdi_div(pdi_sub(pdi_neg(b),disc),pdi_mult(int64_to_pdi(2),a));
		pseudo_double_i sol2=pdi_div(pdi_add(pdi_neg(b),disc),pdi_mult(int64_to_pdi(2),a));
		printf("C (unsafe): Solution 1 = %lf\n",pdi_to_double(sol1));
		printf("C (unsafe): Solution 2 = %lf\n",pdi_to_double(sol2));
	}
// This can be modified to test indivual internal functions
// 	for(uint64_t i=0x4000000000000000ULL;i<=0xFF00000000000000ULL;i+=0x40000000000000ULL) {
// 		cout << hex << i << ' ' << inv_sqrt64_internal(i) << dec << endl;
// 	}
	return 0;
}
