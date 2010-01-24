/*
 *
 * Short program to put out a random number based on the time.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

extern void   srand48();
extern long   lrand48();
extern int    gettimeofday();

int main()
{

  struct timeval tp;
  struct timezone tz;


  gettimeofday(&tp,&tz);

  srand48(tp.tv_usec);
  fprintf(stdout,"%ld",lrand48());

  exit(0);
}
