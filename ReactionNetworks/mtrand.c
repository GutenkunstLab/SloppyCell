#include "mtrand.h"
/* 
A C-program for MT19937, with initialization improved
2002/1/26.
Coded by Takuji Nishimura and Makoto Matsumoto.

Before using, initialize the state by using 
init_genrand(seed) or init_by_array(init_key, key_length).

Copyright(C) 1997-2002, Makoto Matsumoto and Takuji Nishimura
All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the
following conditions are met: 
  1. Redistributions of source code must retain the above
     copyright notice, this list of conditions and the
     following disclaimer.

  2. Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the
     following disclaimer in the documentation and/or other
     materials provided with the distribution.

  3. The names of its contributors may not be used to endorse
     or promote products derived from this software without
     specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Any feedback is very welcome.
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

unsigned long __mtrand_mt[624]; /* the array for the state vector  */
int __mtrand_mti=624+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s) {
  __mtrand_mt[0]= s & 0xffffffffUL;
  for (__mtrand_mti=1; __mtrand_mti<624; __mtrand_mti++) {
    __mtrand_mt[__mtrand_mti]=(1812433253UL*(__mtrand_mt[__mtrand_mti-1]^(__mtrand_mt[__mtrand_mti-1]>>30))+__mtrand_mti); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    __mtrand_mt[__mtrand_mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
  unsigned long y;
  int kk;
  static unsigned long mag01[2]={0x0UL, 0x9908b0dfUL};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */
  
  if (__mtrand_mti >= 624) { /* generate N words at one time */
    if (__mtrand_mti==624+1) /* if init_genrand has not been called, */
      init_genrand(5489UL); /* use a default initial seed */

    for (kk=0;kk<624-397;kk++) {
      y = (__mtrand_mt[kk]&0x80000000UL)|(__mtrand_mt[kk+1]&0x7fffffffUL);
      __mtrand_mt[kk] = __mtrand_mt[kk+397] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<624-1;kk++) {
      y = (__mtrand_mt[kk]&0x80000000UL)|(__mtrand_mt[kk+1]&0x7fffffffUL);
      __mtrand_mt[kk] = __mtrand_mt[kk+(397-624)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (__mtrand_mt[624-1]&0x80000000UL)|(__mtrand_mt[0]&0x7fffffffUL);
    __mtrand_mt[624-1] = __mtrand_mt[397-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
    __mtrand_mti = 0;
  }
  
  y = __mtrand_mt[__mtrand_mti++];
    
  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);
  return y;
}

/* generates a random number on [0,1)-real-interval */
double genrand_real32(void) {
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}
