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
