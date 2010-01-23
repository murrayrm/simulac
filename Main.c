/**************************************
 *
 * Driver for Main Stochastic Simulator
 *
 **************************************/

/*******************************/
/********* Includes ************/
/*******************************/

#ifndef _H_STDIO
#include <stdio.h>
#endif 

#ifndef _H_STDLIB
#include <stdlib.h>
#endif

#ifndef _H_STRING
#include <string.h>
#endif

#ifndef _H_MATH
#include <math.h>
#endif

#include <assert.h>

#ifndef DataStructures
   #include "DataStructures.h"
#endif

#ifndef UTILS
   #include "Util.h"
#endif

#ifndef MEMORY
  #include "Memory.h"
#endif

/****** Object Parameters (see DataStructures.h) ********/

double Rate_Of_Polymerase_Motion,Rate_Of_RNAP_Collision_Escape,Rate_Of_RNAP_Collision_Failure;
double Rate_Of_Ribosome_Binding,Rate_Of_Ribosome_Motion; 

/****** Program Globals (see DataStructures.h)*********/

int        NMechanisms=0;
char       **Mechanism;

int        NSpecies=2;
char     **SpeciesName;
int       *Concentration;

int        NOperators=0;
SHEADATA  *Operator;

int        NSequences=0;
DNA       *Sequence;

int        NTranscripts=0;
mRNA      *Transcript;

int        NMassAction=0;
int      **StoMat1;
int      **StoMat2;
double    *ReactionProbability;

int        NReactions=0;
REACTION  *Reaction;
double     TotalProbability=0.0;
double    *Probabilities;


CELL      *EColi;

double Time;
double MaximumTime;
double PrintTime,WriteTime;

char progid[80];




int main(argc,argv)
int argc;
char **argv;
{
  int i,j;
  int rcnt;
  double tau;
  REACTION *reaction=NULL;
  long SEED;
  
  void SetAckersState();
  void Polymerize();
  void SubmitKinetics();
  REACTION *SelectReaction();
  void ExecuteReaction();
  void FreeReactionQueue();
  void ParseOutline();
  void WriteSpeciesState();
  void FillBicoTable();

  strcpy(progid,argv[0]);
  
  if(argc!=5){
    fprintf(stderr,"usage: %s outline_file Maximum_Time Print_Time SEED\n",argv[0]);
    exit(-1);
  }

  /************
   *
   * Basic Initializations
   *
   ************/

  FillReactionBlock();
  FillRNAPBlock();
  FillRibosomeBlock();
  FillMRNAPBlock();
  FillMRibosomeBlock();

  SpeciesName= (char **) rcalloc((size_t) 2,(size_t) sizeof(char *),"main");
  SpeciesName[0]= (char *) rcalloc(5,sizeof(char),"main");
  strcpy(SpeciesName[0],"RNAP");
  SpeciesName[1]= (char *) rcalloc(9,sizeof(char),"main");
  strcpy(SpeciesName[1],"Ribosome");
  Concentration= (int *) rcalloc(2,sizeof(int),"main");
  Concentration[0]=Concentration[1]=0;

  StoMat1=    (int **) rcalloc(1,sizeof(int *),"main");
  StoMat1[0]= (int *)  rcalloc(2,sizeof(int),"main");
  StoMat2=    (int **) rcalloc(1,sizeof(int *),"main");
  StoMat2[0]= (int *)  rcalloc(2,sizeof(int),"main");
  ReactionProbability= (double *) rcalloc(1,sizeof(double),"main");  

  Rate_Of_Polymerase_Motion        =    30.;
  Rate_Of_RNAP_Collision_Escape    =     5.;
  Rate_Of_RNAP_Collision_Failure   =    30.;
  Rate_Of_Ribosome_Binding         =     0.002;
  Rate_Of_Ribosome_Motion          =   100.; 

  FillBicoTable();

  /****************
   *
   * Read in simulation model
   *
   *****************/
  
  ParseOutline(argv[1]);
  MaximumTime =  atof(argv[2]);
  PrintTime   =  atof(argv[3]);
  SEED        =  atol(argv[4]);
  srand48(SEED);

  DEBUG(20){
    fprintf(stderr,"@@@ NOperators  = %d\n", NOperators);
    fprintf(stderr,"@@@ NSequences  = %d\n", NSequences);
    fprintf(stderr,"@@@ NMassAction = %d\n", NMassAction);
    fprintf(stderr,"@@@ NSpecies    = %d\n", NSpecies);
  }

  /********************
   *
   * Main Loop 
   *
   ********************/

  fprintf(stdout,"!\t\tNR\tRPQ\t");
  for(i=0;i<NSpecies;i++)
    fprintf(stdout,"%6s\t",SpeciesName[i]);
  fprintf(stdout,"%6s\t","Volume");
  for(i=0;i<NOperators; i++)
    fprintf(stdout,"%6s\t",&Operator[i].Name[8]);
  fprintf(stdout,"\n");
  WriteSpeciesState(0.0,0,0.0);
  Time=0.0;
  WriteTime= Time+PrintTime;
  rcnt=0;
  SEED=0;
  do {
    /* Set Promotor States (Assumed Rapid-Equilibrium) */

    /*** Randomly set operator precedence ***/
    j= (int) ((double) NOperators*drand48());
    for(i=0; i<NOperators; i++){     
      SetAckersState(&Operator[j]);
      j++;
      if(j==NOperators) j=0;
    }
        
    /* Submit All Genetic Reactions, execute all current genetic Actions */
    
    Polymerize();
    
    /* Submit all non-genetic chemical reactions */
    
    SubmitKinetics();

    /* Submit Volume Changes */

    SubmitCellReactions();
    
    /* Reaction Clock */

    reaction= (REACTION *) SelectReaction(&tau);

    if(Time+tau > MaximumTime) break;
    while(Time+tau > WriteTime){
      WriteSpeciesState(WriteTime,rcnt,(rcnt> 0 ? (double) SEED/rcnt : 0.0));
      WriteTime= WriteTime+PrintTime;
      rcnt=0;
    }

    ExecuteReaction(reaction);
    rcnt++;
    SEED+=NReactions;
    /*    
    fprintf(stderr,"NR= %d\t",NReactions);    
    fprintf(stderr,"react= %d\trnap= %d\tribo= %d\tmrnap= %d\tmribo= %d\n",
	    react_mptr_full,
	    rnap_mptr_full,
	    ribo_mptr_full,
	    mrnap_mptr_full,
	    mribo_mptr_full
	    );
    */
    FreeReactionQueue();
    Time += tau;    

  } while(Time<=MaximumTime);
  
  while(Time<MaximumTime){
    WriteSpeciesState(WriteTime,rcnt,(rcnt> 0 ? (double) SEED/rcnt : 0.0));
    WriteTime= WriteTime+PrintTime;
    rcnt=0;
    Time += PrintTime;
  }

  WriteSpeciesState(WriteTime,rcnt,(rcnt> 0 ? (double) SEED/rcnt : 0.0));

  EmptyReactionBlock();
  EmptyRNAPBlock();
  EmptyRibosomeBlock();
  EmptyMRNAPBlock();
  EmptyMRibosomeBlock();

  exit(0);
}


void WriteSpeciesState(t,cnt,rpq)
double t,rpq;
int cnt;
{

  int i;

  fprintf(stdout,"%e\t%d\t%e\t",t,cnt,rpq);
  for(i=0; i<NSpecies; i++)
    fprintf(stdout,"%7d\t",Concentration[i]);

  fprintf(stdout,"%e\t",EColi->V/EColi->VI);

  for(i=0; i<NOperators; i++)
    fprintf(stdout,"%7d\t",Operator[i].CurrentState);

  fprintf(stdout,"\n");
  fflush(stdout);
}

