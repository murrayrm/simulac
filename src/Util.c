/*****************************
 * 
 * Utility Functions
 * 
 *****************************/

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

#ifndef DataStructures
 #include "DataStructures.h"
#endif

/*** These are needed because AIX is brain-dead ***/

extern int strcasecmp();
extern int strncasecmp();

#define TRUE         1
#define FALSE        0

#define NEW          0
#define OLD          1

int EventLevel=0;
#define Event(A) if((A)>EventLevel)

int DebugLevel=0;
#define DEBUG(A) if((A)<DebugLevel)

/***********
 * 
 * For Debug 
 *
 *************/

void *rcalloc(nelem,size,id)
size_t nelem;
size_t size;
char *id;
{
  void *ptr;

  if((ptr=calloc(nelem,size))==NULL){
    fprintf(stderr,"%s: Error Allocating Memory (%d,%d)\n",
	    progid, (int) nelem, (int) size);
    if(id) 
      fprintf(stderr,"%s: ID= %s\n", progid, id);
    exit(-1);
  }

  return(ptr);
}

void *rrealloc(happyptr,nelem,size,id)
void *happyptr;
size_t nelem;
size_t size;
char *id;
{
  void *ptr;

  if((ptr=calloc(nelem,size))==NULL){
    fprintf(stderr,"%s: Error reallocating memory\n",progid);
    if(id) 
      fprintf(stderr,"%s: ID= %s\n",progid,id);
    exit(-1);
  }

  
  /***### This -1 assumes we are only changing things by 1 ###***/

  memcpy(ptr,happyptr,(nelem-1)*size);
  
  return(ptr);
}

char *PrintDNAType(type)
int type;
{

  switch(type){
  case DNA_Type_Promotor:       return("Promotor"); 
  case DNA_Type_Coding:         return("Coding"); 
  case DNA_Type_NonCoding:      return("NonCoding"); 
  case DNA_Type_Terminator:     return("Terminator"); 
  case DNA_Type_AntiTerminator: return("AntiTerminator"); 
  default: return("*** Unknown ***"); 
  }
}

char *FindReactionType(type)
int type;
{
  switch(type){

  case Reaction_Type_Kinetic:
    return("Kinetic");
  case Reaction_Type_TransInit:
    return("Transcription Initiation");
  case Reaction_Type_MoveRNAP:         
    return("Move RNAP");
  case Reaction_Type_DNAAction:        
    return("DNA Action");
  case Reaction_Type_NextSegment:      
    return("Jump Segment");
  case Reaction_Type_EatmRNA:          
    return("mRNA Degradation");
  case Reaction_Type_MoveRibosome:     
    return("Move Ribosome");
  case Reaction_Type_ProduceProtein:   
    return("Produce Protein");
  case Reaction_Type_RNAP_RNAP:        
    return("RNAP Collision");
  case Reaction_Type_Recombine:        
    return("Recombination");
  case Reaction_Type_BindRibosome:     
    return("Bind Ribosome");
  case Reaction_Type_ProduceNewProtein:
    return("Produce New Protein");
  default:
    return("Unknown");
  }

}



/**********
 *
 * Combinatorics
 *
 *
 **********/

double gammln(xx)
double xx;
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double factln(n)
int n;
{
        static short flag=0;
	double gammln();
	static double a[101];
	
	if(flag==0){
	  int i;
	  for(i=0; i<101; i++)
	    a[i]=0.0;
	  flag=1;
	}

	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}


double bico(n,k)
int k,n;
{
  static short flag=0;
  double factln();
  static double table[1000][4];
  
  if(flag==0){
    int i,j;
    for(i=1; i<=1000; i++)
      for(j=1; j<=4; j++)
	table[i-1][j-1] = floor(0.5+exp(factln(i)-factln(j)-factln(i-j)));
    flag=1;
  }
  
  if(k==0)                return(1.0);
  if(n==0 || k>n)         return(0.0);    
  if(k==n)                return(1.0);
  
  if(n<=1000 && k<=4) return(table[n-1][k-1]);
    
  return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

void FillBicoTable()
{
  int i,j;
  double bico();

  for(i=1; i<=20; i++)
    for(j=1; j<=4; j++)
      (void) bico(i,j);  
}

unsigned long choose(nn,mm)
int mm,nn;
{
  unsigned long i;
  unsigned long factorial();

  if(mm==0) return(1);
  if(nn==0 || mm>nn) return(0);
  
  i= factorial(nn,nn-mm)/(factorial(mm,0));
  return(i);
}

unsigned long factorial(nn,t)
int nn,t;
{
  long i;
 
  i=1;
  while(nn>t){
    i *= nn;
    nn--;
  }

return(i);
}
