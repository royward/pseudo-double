// ---- conversion to decimal strings ----
// |x| is scaled by a power of 10 to a 128-bit value (a floor) with a sticky flag for every bit dropped, giving 19
// significant digits plus 64 more fraction bits. Rounding to up to 19 digits is then exact (ties to even when the
// value is exact), except within about 2^-59 of a unit in the 19th digit of a rounding boundary. Digits asked for
// beyond the 19th are written as 0

typedef struct {
	uint64_t d19;                                                      // 19 digits, in [10^18,10^19)
	uint64_t frac;                                                     // the next 64 bits, as 0.64
	uint64_t sticky;                                                   // nonzero: more bits below frac
	int32_t exp10;                                                     // power of 10 of the first digit
} pd_decimal_t;

// |x| for a positive nonzero x: the floor of |x|*10^(18-exp10), as d19 then frac and sticky
static inline void pd_to_decimal(pseudo_double_i a, pd_decimal_t* r) {
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // |x| = ax*2^(ex-64), in [2^(ex-2),2^(ex-1))
	uint64_t ax=(uint64_t)(a&EXP_MASK_INV);
	int32_t d0=(int32_t)(((int64_t)(ex-2)*1292913986LL)>>32);          // floor((ex-2)*log10(2)): exp10 is d0 or d0+1
	for(;;) {
		// |x|*10^(17-d0) is then in [10^17,10^19), below 2^64 (if d0 were off, try again)
		uint64_t hi=ax<<1, lo=0, st=0;
		int32_t f=ex-129;                                              // |x| = (hi:lo)*2^f
		pd_scale10(&hi,&lo,&f,17-d0,&st);
		int32_t s=-f;                                                  // value = (hi:lo)/2^s
		if(s<64) {
			d0++;
			continue;
		}
		uint64_t dd,frac;
		if(s==64) {
			dd=hi;
			frac=lo;
		} else if(s<128) {
			dd=hi>>(s-64);
			frac=(hi<<(128-s))|(lo>>(s-64));
			st|=lo<<(128-s);
		} else {
			d0--;
			continue;
		}
		if(dd>=pd_pow10[19]) {
			d0++;
			continue;
		}
		if(dd<pd_pow10[17]) {
			d0--;
			continue;
		}
		if(dd<pd_pow10[18]) {                                          // 18 digits: one more from frac (exact)
			dd=dd*10+multu64hi(frac,10);
			frac*=10;
			r->exp10=d0;
		} else {
			r->exp10=d0+1;
		}
		r->d19=dd;
		r->frac=frac;
		r->sticky=st;
		return;
	}
}

// d19 rounded to n significant digits (1 <= n <= 19), to nearest with ties to even: returns q in [10^(n-1),10^n),
// with the power of 10 of its first digit in *exp10
static inline uint64_t pd_decimal_round(const pd_decimal_t* r, int n, int32_t* exp10) {
	uint64_t q,up;
	if(n>=19) {
		q=r->d19;
		up=r->frac>(1ULL<<63) || (r->frac==(1ULL<<63) && (r->sticky!=0 || (q&1)));
	} else {
		uint64_t p=pd_pow10[19-n], tail=r->d19%p, half=p>>1;
		q=r->d19/p;
		up=tail>half || (tail==half && ((r->frac|r->sticky)!=0 || (q&1)));
	}
	q+=up;
	*exp10=r->exp10;
	if(q==pd_pow10[n]) {                                               // rounded up through a power of 10
		q=pd_pow10[n-1];
		*exp10+=1;
	}
	return q;
}

// the fewest significant digits that read back (with string_to_pdi) as the positive x. For n digits the candidates
// are the two n-digit decimals either side of x: the one x rounds to if it reads back, else the other one (which
// happens when x is close to the middle, and at powers of 2 where the gap below is half the gap above). Some n
// digits read back for all n from some point on, so a binary search finds it. Precision p bits always reads back
// with 2+floor(p*log10(2)) digits
static inline int pd_decimal_n_reads_back(pseudo_double_i a, const pd_decimal_t* r, int n, uint64_t* q, int32_t* exp10) {
	uint64_t c=pd_decimal_round(r,n,exp10);
	if(decimal_to_pdi(0,c,false,*exp10-n+1,false)==a) {
		*q=c;
		return 1;
	}
	uint64_t f=r->d19/pd_pow10[19-n];                                  // the n digits below x
	uint64_t other=(c==f && *exp10==r->exp10)?f+1:f;                   // whichever of f and f+1 c isn't
	if(decimal_to_pdi(0,other,false,r->exp10-n+1,false)==a) {
		*exp10=r->exp10;
		if(other==pd_pow10[n]) {
			other=pd_pow10[n-1];
			*exp10+=1;
		}
		*q=other;
		return 1;
	}
	return 0;
}

static inline int pd_decimal_shortest(pseudo_double_i a, const pd_decimal_t* r, uint64_t* q, int32_t* exp10) {
	int lo=1, hi=2+(int)(((63-PSEUDO_DOUBLE_EXP_BITS)*1233)>>12);
	while(lo<hi) {
		int mid=(lo+hi)>>1;
		uint64_t c;
		int32_t e;
		if(pd_decimal_n_reads_back(a,r,mid,&c,&e)) {
			hi=mid;
		} else {
			lo=mid+1;
		}
	}
	if(!pd_decimal_n_reads_back(a,r,lo,q,exp10)) {
		*q=pd_decimal_round(r,lo,exp10);                               // (not reached)
	}
	return lo;
}

// the exact decimal digits of a positive integer x below 2^127, into digs (at least 40 chars): returns how many
static int pd_int_digits(pseudo_double_i a, char* digs) {
	int32_t e=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64;         // x = m*2^e, e <= 64
	uint64_t m=(uint64_t)(a&EXP_MASK_INV), hi, lo;
	if(e>=64) {
		hi=m<<(e-64);
		lo=0;
	} else if(e>0) {
		hi=m>>(64-e);
		lo=m<<e;
	} else {
		hi=0;
		lo=m>>-e;                                                      // exact: x is an integer
	}
	char t[45];
	int k=0;
	do {
		// divide (hi:lo) by 10^9 in 32-bit limbs, taking the remainder's 9 digits
		uint64_t w[4]={hi>>32,hi&0xFFFFFFFFULL,lo>>32,lo&0xFFFFFFFFULL}, r=0;
		for(int i=0;i<4;i++) {
			uint64_t cur=(r<<32)|w[i];
			w[i]=cur/1000000000ULL;
			r=cur%1000000000ULL;
		}
		hi=(w[0]<<32)|w[1];
		lo=(w[2]<<32)|w[3];
		for(int j=0;j<9;j++) {
			t[k++]=(char)('0'+r%10);
			r/=10;
		}
	} while((hi|lo)!=0);
	while(k>1 && t[k-1]=='0') {
		k--;                                                           // leading zeros
	}
	for(int i=0;i<k;i++) {
		digs[i]=t[k-1-i];
	}
	return k;
}

// whether a positive x is below 2^127 (and so at least 2^precision, an integer once it has 20 digits)
static inline int pd_below_2_127(pseudo_double_i a) {
	return (int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64<=64;
}

// ---- exact digits for long outputs (more than 19 significant digits) ----
// |x|*10^k is formed exactly as a big integer (|x| = m*2^e, so m*10^k shifted by e), then written out in decimal.
// This covers anything up to PD_BIG_LIMBS*32 bits; longer still (only for huge exponent widths or precisions)
// falls back to 19 correct digits followed by zeros

#ifndef PD_BIG_LIMBS
#define PD_BIG_LIMBS 160                                               // 5120 bits
#endif

typedef struct {
	uint32_t w[PD_BIG_LIMBS];                                          // little endian
	int n;                                                             // limbs in use
} pd_big_t;

static inline void pd_big_mul(pd_big_t* b, uint32_t m) {
	uint64_t c=0;
	for(int i=0;i<b->n;i++) {
		c+=(uint64_t)b->w[i]*m;
		b->w[i]=(uint32_t)c;
		c>>=32;
	}
	if(c!=0) {
		b->w[b->n++]=(uint32_t)c;
	}
}

static inline uint32_t pd_big_div(pd_big_t* b, uint32_t d) {             // returns the remainder
	uint64_t r=0;
	for(int i=b->n-1;i>=0;i--) {
		uint64_t cur=(r<<32)|b->w[i];
		b->w[i]=(uint32_t)(cur/d);
		r=cur%d;
	}
	while(b->n>0 && b->w[b->n-1]==0) {
		b->n--;
	}
	return (uint32_t)r;
}

// |x|*10^k for a positive x, as decimal digits (most significant first, no leading zeros, "0" for 0) of its floor,
// returning the count (0 if it won't fit). *round_info: 0 exact, 1 below a half, 2 exactly a half, 3 above a half
static int pd_exact_scaled(pseudo_double_i a, int32_t k, char* digs, int maxdigs, int* round_info) {
	int32_t e=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64;         // |x| = m*2^e
	uint64_t m=(uint64_t)(a&EXP_MASK_INV);
	if(k<0 || (int64_t)64+(int64_t)k*10/3+1+(e>0?e:0)>(int64_t)PD_BIG_LIMBS*32-32) {
		return 0;
	}
	pd_big_t b;
	b.w[0]=(uint32_t)m;
	b.w[1]=(uint32_t)(m>>32);
	b.n=2;
	for(;k>=9;k-=9) {
		pd_big_mul(&b,1000000000u);
	}
	if(k>0) {
		pd_big_mul(&b,(uint32_t)pd_pow10[k]);
	}
	*round_info=0;
	if(e>0) {                                                          // shift left
		int32_t ws=e>>5, bs=e&31;
		if(bs!=0) {
			pd_big_mul(&b,1u<<bs);
		}
		if(ws>0) {
			for(int i=b.n-1;i>=0;i--) {
				b.w[i+ws]=b.w[i];
			}
			for(int i=0;i<ws;i++) {
				b.w[i]=0;
			}
			b.n+=ws;
		}
	} else if(e<0) {                                                   // shift right, classifying what drops off
		int32_t s=-e;
		int half=0, rest=0;
		for(int32_t i=0;i<s;i++) {
			int32_t limb=i>>5;
			int bit=limb<b.n?(int)((b.w[limb]>>(i&31))&1):0;
			if(i==s-1) {
				half=bit;
			} else {
				rest|=bit;
			}
		}
		*round_info=half?(rest?3:2):(rest?1:0);
		int32_t ws=s>>5, bs=s&31;
		int nn=b.n-ws;
		if(nn<=0) {
			b.n=0;
		} else {
			for(int i=0;i<nn;i++) {
				uint64_t v=b.w[i+ws];
				if(bs!=0 && i+ws+1<b.n) {
					v|=(uint64_t)b.w[i+ws+1]<<32;
				}
				b.w[i]=(uint32_t)(v>>bs);
			}
			b.n=nn;
		}
		while(b.n>0 && b.w[b.n-1]==0) {
			b.n--;
		}
	}
	// to decimal, 9 digits at a time from the bottom
	char t[PD_BIG_LIMBS*10+10];
	int nt=0;
	do {
		uint32_t r=pd_big_div(&b,1000000000u);
		for(int j=0;j<9;j++) {
			t[nt++]=(char)('0'+r%10);
			r/=10;
		}
	} while(b.n>0);
	while(nt>1 && t[nt-1]=='0') {
		nt--;
	}
	if(nt>maxdigs) {
		return 0;
	}
	for(int i=0;i<nt;i++) {
		digs[i]=t[nt-1-i];
	}
	return nt;
}

// rounds the digit string digs[0..len-1] up by one in the last place: returns 1 if it carried out of the top
// (digs is then 1000...)
static inline int pd_digits_inc(char* digs, int len) {
	for(int i=len-1;i>=0;i--) {
		if(digs[i]!='9') {
			digs[i]++;
			return 0;
		}
		digs[i]='0';
	}
	digs[0]='1';
	return 1;
}

// |x| to exactly n significant digits (n > 19), correctly rounded: digs gets n digits, *ex the power of 10 of the
// first one (exp10 is the first digit's power before rounding). Returns 0 if it won't fit
static int pd_exact_sig(pseudo_double_i a, int n, int32_t exp10, char* digs, int maxdigs, int32_t* ex) {
	int32_t k=n-1-exp10;
	int info,len;
	if(k>=0) {
		len=pd_exact_scaled(a,k,digs,maxdigs,&info);                   // exactly n digits
		if(len!=n) {
			return 0;
		}
	} else {
		// at least 10^20, so an integer: all its digits, then round at n from the ones after
		len=pd_exact_scaled(a,0,digs,maxdigs,&info);
		if(len!=exp10+1) {
			return 0;
		}
		int rest=0;
		for(int i=n+1;i<len;i++) {
			rest|=digs[i]!='0';
		}
		info=digs[n]>'5'?3:(digs[n]=='5'?(rest?3:2):((digs[n]!='0' || rest)?1:0));
	}
	*ex=exp10;
	if(info==3 || (info==2 && ((digs[n-1]-'0')&1))) {
		*ex+=pd_digits_inc(digs,n);
	}
	return 1;
}

// ---- formatting ----

#ifndef PD_FMT_PLUS                                                    // (these belong in pseudo_double.h)
#define PD_FMT_PLUS 1    // '+': always write a sign
#define PD_FMT_SPACE 2   // ' ': a space where a plus sign would go
#define PD_FMT_ALT 4     // '#': always write the decimal point, and %g keeps its trailing zeros
#endif

typedef struct {
	char* buf;
	size_t size;
	size_t len;                                                        // what the whole result needs
} pd_out_t;

static inline void pd_out(pd_out_t* o, char c) {
	if(o->len+1<o->size) {
		o->buf[o->len]=c;
	}
	o->len++;
}

// the significant digits digs[0..n-1] (first at 10^x), as %f with p digits after the point; digits past n are 0
static void pd_out_fixed(pd_out_t* o, const char* digs, int n, int32_t x, int32_t p, int alt, char point) {
	if(x<0) {
		pd_out(o,'0');
	} else {
		for(int32_t i=0;i<=x;i++) {
			pd_out(o,i<n?digs[i]:'0');
		}
	}
	if(p>0 || alt) {
		pd_out(o,point);
	}
	for(int32_t j=1;j<=p;j++) {
		int32_t i=x+j;
		pd_out(o,(i>=0 && i<n)?digs[i]:'0');
	}
}

// the same as %e with p digits after the point
static void pd_out_exp(pd_out_t* o, const char* digs, int n, int32_t x, int32_t p, int alt, int upper, char point) {
	pd_out(o,digs[0]);
	if(p>0 || alt) {
		pd_out(o,point);
	}
	for(int32_t i=1;i<=p;i++) {
		pd_out(o,i<n?digs[i]:'0');
	}
	pd_out(o,upper?'E':'e');
	pd_out(o,x<0?'-':'+');
	uint32_t ux=x<0?0u-(uint32_t)x:(uint32_t)x;
	char t[12];
	int k=0;
	do {
		t[k++]=(char)('0'+ux%10);
		ux/=10;
	} while(ux!=0);
	if(k<2) {
		t[k++]='0';                                                    // at least two exponent digits, as printf
	}
	while(k>0) {
		pd_out(o,t[--k]);
	}
}

// the digits of q, exactly n of them (n <= 19)
static inline void pd_digits(char* digs, uint64_t q, int n) {
	for(int i=n-1;i>=0;i--) {
		digs[i]=(char)('0'+q%10);
		q/=10;
	}
}

// Formats x as printf's %.{precision}{conv} does, for conv 'e', 'f', 'g' (or 'E', 'F', 'G'), with PD_FMT_ flags,
// into buf (size bytes, NUL terminated when size > 0). A negative precision means the shortest digits that read
// back as x: 'e' and 'f' write all of them, and 'g' takes whichever of the two is shorter (fixed on a tie), as
// C++'s std::to_chars does. point is the decimal point character. Returns the length the whole result needs, as
// snprintf does. Correctly rounded (ties to even when exact) to any number of digits, except that outputs needing
// more than PD_BIG_LIMBS*32 bits of working (only possible with very wide exponents or huge precisions) get 19
// correct digits and then zeros
int pdi_format_point(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags, char point) {
	pd_out_t o={buf,size,0};
	int upper=(conv>='A' && conv<='Z');
	char c=(char)(conv|0x20);
	int alt=(flags&PD_FMT_ALT)!=0;
	int negative=x!=PD_NAN && ((signed_pd_internal)x)<0;
	if(negative) {
		pd_out(&o,'-');
	} else if(flags&PD_FMT_PLUS) {
		pd_out(&o,'+');
	} else if(flags&PD_FMT_SPACE) {
		pd_out(&o,' ');
	}
	if(x==PD_NAN) {
		const char* s=upper?"NAN":"nan";
		while(*s) {
			pd_out(&o,*s++);
		}
	} else {
		pseudo_double_i a=negative?pdi_neg(x):x;
		char digs[PD_BIG_LIMBS*10+10];
		int n=1;
		int32_t ex=0;                                                  // power of 10 of digs[0]
		pd_decimal_t r={0,0,0,0};
		if(x==0) {
			digs[0]='0';
		} else {
			pd_to_decimal(a,&r);
		}
		if(precision<0) {
			// shortest digits that read back as x
			if(x!=0) {
				uint64_t q;
				n=pd_decimal_shortest(a,&r,&q,&ex);
				pd_digits(digs,q,n);
			}
			int32_t fp=n-1-ex>0?n-1-ex:0;                              // digits needed after the point
			if(c=='g') {
				int32_t flen=(ex>=0?ex+1:1)+(fp>0?fp+1:0);
				int32_t elen=n+(n>1)+4;                                // digits, point, e, sign, two exponent digits
				for(int32_t ae=(ex<0?-ex:ex)/100;ae!=0;ae/=10) {
					elen++;                                            // more exponent digits
				}
				c=elen<flen?'e':'f';
			}
			if(c=='e') {
				pd_out_exp(&o,digs,n,ex,n-1,alt,upper,point);
			} else {
				if(ex+1>n && pd_below_2_127(a)) {
					n=pd_int_digits(a,digs);                           // an integer: its exact digits, as std::to_chars
					ex=n-1;
				}
				pd_out_fixed(&o,digs,n,ex,fp,alt,point);
			}
		} else if(c=='e') {
			int want=precision+1;
			if(x!=0) {
				n=want<19?want:19;
				if(want<=19 || !pd_exact_sig(a,n=want,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=want<19?want:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				}
			}
			pd_out_exp(&o,digs,n,ex,precision,alt,upper,point);
		} else if(c=='f') {
			if(x!=0) {
				int64_t want=(int64_t)r.exp10+1+precision;             // significant digits down to 10^-precision
				if(want>19 && want<(int64_t)sizeof digs && pd_exact_sig(a,(int)want,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=(int)want;                                       // more than 19 digits: exact
				} else if(want>19 && r.exp10>=19 && pd_below_2_127(a)) {
					n=pd_int_digits(a,digs);                           // an integer of 20 to 39 digits: exact
					ex=n-1;
				} else if(want>=1) {
					n=want<19?(int)want:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				} else {
					// below the last place: rounds to 0 or to one unit of 10^-precision (ties to even, so to 0)
					int up=want==0 && (r.d19>5*pd_pow10[18] || (r.d19==5*pd_pow10[18] && (r.frac|r.sticky)!=0));
					digs[0]=up?'1':'0';
					ex=up?-precision:0;
				}
			}
			pd_out_fixed(&o,digs,n,ex,precision,alt,point);
		} else {
			// 'g': precision significant digits (at least 1), then %e if the exponent is below -4 or at least the
			// precision, else %f; trailing zeros dropped unless '#'
			int pr=precision==0?1:precision;
			if(x!=0) {
				if(pr<=19 || !pd_exact_sig(a,n=pr,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=pr<19?pr:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				}
			}
			int last=n-1;                                              // last nonzero digit
			while(last>0 && digs[last]=='0') {
				last--;
			}
			if(ex<-4 || ex>=pr) {
				pd_out_exp(&o,digs,n,ex,alt?pr-1:last,alt,upper,point);
			} else {
				int32_t fp=alt?pr-1-ex:(last-ex>0?last-ex:0);
				pd_out_fixed(&o,digs,n,ex,fp,alt,point);
			}
		}
	}
	if(size>0) {
		buf[o.len<size?o.len:size-1]='\0';
	}
	return (int)o.len;
}

int pdi_format(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags) {
	return pdi_format_point(buf,size,x,conv,precision,flags,'.');
}

// the shortest string that reads back as x (string_to_pdi), fixed or scientific, whichever is shorter. 32 bytes
// is always enough
int pdi_to_string(char* buf, size_t size, pseudo_double_i x) {
	return pdi_format_point(buf,size,x,'g',-1,0,'.');
}
