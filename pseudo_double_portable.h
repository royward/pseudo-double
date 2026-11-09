// Portable 64-bit helpers for compilers without __int128 (e.g. MSVC on ARM64)
// Pure integer code: results are bit-identical to the __int128 / intrinsic versions
#include <stdint.h>

// Count leading zeros. x must be non-zero
static inline int pd_clz64_portable(uint64_t x) {
	int n=0;
	if(!(x>>32)) {n+=32; x<<=32;}
	if(!(x>>48)) {n+=16; x<<=16;}
	if(!(x>>56)) {n+=8;  x<<=8;}
	if(!(x>>60)) {n+=4;  x<<=4;}
	if(!(x>>62)) {n+=2;  x<<=2;}
	if(!(x>>63)) {n+=1;}
	return n;
}

// High 64 bits of the unsigned 128-bit product
static inline uint64_t pd_umulh_portable(uint64_t a, uint64_t b) {
	uint64_t a0=a&0xFFFFFFFFu, a1=a>>32, b0=b&0xFFFFFFFFu, b1=b>>32;
	uint64_t p00=a0*b0, p01=a0*b1, p10=a1*b0, p11=a1*b1;
	uint64_t mid=(p00>>32)+(p01&0xFFFFFFFFu)+(p10&0xFFFFFFFFu);
	return p11+(p01>>32)+(p10>>32)+(mid>>32);
}

// High 64 bits of the signed 128-bit product
static inline int64_t pd_mulh_portable(int64_t a, int64_t b) {
	uint64_t h=pd_umulh_portable((uint64_t)a,(uint64_t)b);
	h-=(a<0?(uint64_t)b:0)+(b<0?(uint64_t)a:0);
	return (int64_t)h;
}

// Unsigned (u1*2^64+u0)/v, requires u1<v. Knuth algorithm D with 32-bit digits
//(Hacker's Delight divlu). Returns the quotient, stores the remainder in *r
static inline uint64_t pd_divlu_portable(uint64_t u1, uint64_t u0, uint64_t v, uint64_t *r) {
	const uint64_t b=1ULL<<32;
	int s=pd_clz64_portable(v);
	v<<=s;
	uint64_t vn1=v>>32, vn0=v&0xFFFFFFFFu;
	uint64_t un32=(u1<<s)|(s?(u0>>(64-s)):0);
	uint64_t un10=u0<<s;
	uint64_t un1=un10>>32, un0=un10&0xFFFFFFFFu;
	uint64_t q1=un32/vn1, rhat=un32-q1*vn1;
	while(q1>=b || q1*vn0>b*rhat+un1) {q1--; rhat+=vn1; if(rhat>=b) break;}
	uint64_t un21=un32*b+un1-q1*v;      /* wraps mod 2^64 by design */
	uint64_t q0=un21/vn1; rhat=un21-q0*vn1;
	while(q0>=b || q0*vn0>b*rhat+un0) {q0--; rhat+=vn1; if(rhat>=b) break;}
	*r=(un21*b+un0-q0*v)>>s;
	return q1*b+q0;
}

// floor(x*2^64/y). Requires |x|<|y| so the quotient fits (true in pdi_div)
static inline int64_t pd_divs64hi_portable(int64_t x, int64_t y) {
	uint64_t ux=x<0?0-(uint64_t)x:(uint64_t)x;
	uint64_t uy=y<0?0-(uint64_t)y:(uint64_t)y;
	uint64_t r;
	uint64_t uq=pd_divlu_portable(ux,0,uy,&r);
	if((x^y)<0) {
		return (int64_t)(0-uq-(r!=0));   /* negative quotient: floor = -(uq+1) if inexact */
	}
	return (int64_t)uq;
}
