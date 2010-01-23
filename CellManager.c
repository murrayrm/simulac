/**************************
 *
 * Subroutines to handle 
 * Cell Reactions such as volume and
 * position changes
 **************************/

/****************************/
/******* Includes ***********/
/****************************/

#include <assert.h>

#ifndef _H_STDIO
   #include <stdio.h>
#endif

#ifndef _H_STDLIB
   #include <stdlib.h>
#endif

#ifndef _H_MATH
   #include <math.h>
#endif

#ifndef DataStructures
   #include "DataStructures.h"
#endif

#ifndef UTILS
 #include "Util.h"
#endif

#ifndef MEMORY
 #include "Memory.h"
#endif

/******************************/
/******* Submission ***********/
/******************************/

void SubmitCellReactions()
{
  int i,j,nspecies,order;
  double prob;
  REACTION  *reaction;
  REACTDATA *rdata;
  void Balloon();
  void SubmitReaction();
  void MassAction();

  if(EColi->GrowthRate== 0.0) return;


  reaction= (REACTION *)   AllocReaction();
  reaction->Type=          Reaction_Type_ChangeCellVolume;
  reaction->ReactionFunc=  Balloon;
  reaction->Probability=   EColi->GrowthRate;    
  reaction->ReactionData=  NULL;
  SubmitReaction(reaction);

}

void Balloon(rdata)
void *rdata;
{
  int i;
  static long idum= -12;
  float bnldev();


  EColi->V += 1e-18;

  if((EColi->V/EColi->VI) >= 2.0){ /* Then Time to Divide */
    EColi->V /= 2.0;

    /* Binomial Partition of Chemical Species */

    for(i=0; i<NSpecies; i++)
      Concentration[i]=(int)bnldev(0.5,Concentration[i],&idum);
  }
    
}


#define PI 3.141592654

float bnldev(pp,n,idum)
float pp;
int n;
long *idum;
{
	double gammln();
	float ran1();
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg= (float) gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(idum)
long *idum;
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
