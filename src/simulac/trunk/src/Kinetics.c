/**************************
 *
 * Subroutines to handle 
 * homogeneous kinetics
 *
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

void SubmitKinetics()
{
  int i,j,nmolecs,order;
  double prob;
  REACTION  *reaction;
  REACTDATA *rdata;

  void SubmitReaction();
  void MassAction();
  double bico();
  unsigned long choose();


  for(i=0; i<NMassAction; i++){
    prob= ReactionProbability[i];

    nmolecs=0;
    for(j=0; j<NSpecies; j++){
      prob *= (double) bico(Concentration[j],StoMat1[i][j]);
      nmolecs += StoMat1[i][j];
    }
    
    /* Correction for Volume Changes*/
    
    if(prob!=0.0 && nmolecs!=1){ /* if First order then no correction */
      order= nmolecs-1; /* Note that if reaction is zeroth order then it is multiplied by an inverse volume
			 * This assumes that these reactions work to maintain molarity....
			 */
      prob *= pow(EColi->V0/EColi->V,(double) order);
    }
    /* Make Reaction */

    if(prob>1e-20){
      reaction= (REACTION *)   AllocReaction();
      reaction->Type=          Reaction_Type_Kinetic;
      reaction->ReactionFunc=  MassAction;
      reaction->Probability=   prob;    
      rdata= (REACTDATA *) rcalloc(1,sizeof(REACTDATA), "SubmitKinetics");
      rdata->Mu=i;
      reaction->ReactionData=  (void *) rdata;
      SubmitReaction(reaction);
    }
  }
}


void MassAction(rdata)
void *rdata;
{
  int i;
  REACTDATA *data;
  int mu;

  data= (REACTDATA *) rdata;
  mu= data->Mu;
  
   for(i=0; i<NSpecies; i++){
      Concentration[i] += StoMat2[mu][i] - StoMat1[mu][i];
   }
}





