/******************************
 *
 * Global Variables and 
 * Function Declarations
 *
 ******************************/

#define UTILS

#ifndef _H_STDIO
   #include <stdio.h>
#endif

#ifndef _H_STDLIB
   #include <stdlib.h>
#endif

#ifndef _H_MATH
   #include <math.h>
#endif

#ifndef _H_MALLOC
 #include <malloc.h>
#endif

#ifndef _H_MEMORY
 #include <memory.h>
#endif

/*** These are needed because AIX is brain-dead ***/

#define TRUE         1
#define FALSE        0

#define NEW          0
#define OLD          1

extern int EventLevel;
#define Event(A) if((A)>EventLevel)

extern int DebugLevel;
#define DEBUG(A) if((A)<DebugLevel)

/***********
 * 
 * For Debug 
 *
 *************/

extern void *rcalloc(size_t,size_t,char *);
extern void *rrealloc(void *,size_t,size_t,char *);
extern char *PrintDNAType(int);
extern char *FindReactionType(int);

/**********
 *
 * Combinatorics
 *
 *
 **********/

extern double gammln(double);
extern double factln(int);
extern double bico(int,int);
extern void FillBicoTable();
extern unsigned long choose(int,int);
extern unsigned long factorial(int,int);
