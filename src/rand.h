#ifndef _DEF_H_
#define _DEF_H_ 1

extern double seed;
extern double oldrand[55];
extern int jrand;

void randomize(void);
void warmup_random (double);
void advance_random (void);
double randomperc(void);
int rnd (int , int);
double rndreal (double , double);

#endif
