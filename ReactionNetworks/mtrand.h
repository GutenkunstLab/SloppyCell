#ifndef __MTRAND_HEADER__
#define __MTRAND_HEADER__

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void);

/* generates a random number on [0,1)-real-interval */
double genrand_real32(void);

#endif /* __MTRAND_HEADER__ */
