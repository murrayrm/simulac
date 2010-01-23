/*******************
 *
 * These routines
 * manage the reaction queue
 * and form the heart of the reaction clock method
 * for stochastic simulation 
 *
 ******************/

/****************************/
/******* Includes ***********/
/****************************/


#ifndef _H_STDIO
   #include <stdio.h>
#endif

#ifndef _H_STDLIB
   #include <stdlib.h>
#endif

#ifndef _H_MATH
   #include <math.h>
#endif

#include <assert.h>

#ifndef DataStructures
   #include "DataStructures.h"
#endif

#ifndef MEMORY
 #include "Memory.h"
#endif

/****************************/
/** Reaction Queue Manager **/
/****************************/

void SubmitReaction(reaction)
REACTION *reaction;
{
  void FreeReactionData();

  if(reaction->Probability==0.0){
    if(reaction->ReactionData!=NULL)
      FreeReactionData(reaction->Type,reaction->ReactionData);
    FreeReaction(reaction);
    return;
  }

  if(NReactions==0){
    Reaction=reaction;
    reaction->LastReaction=NULL;
    reaction->NextReaction=NULL;
  } else {
    reaction->NextReaction=Reaction;
    Reaction->LastReaction=reaction;
    Reaction=reaction;
    Reaction->LastReaction=NULL;
  }

  TotalProbability += reaction->Probability;
  NReactions++;
}

/**********************
 *
 * Step 2 from Gillespie 
 *
 ***********************/

#define TINY 1e-16

REACTION *SelectReaction(tau)
double *tau;
{
  double r1,r2,sum;
  REACTION *reaction, *lastreaction;
  
  r1= drand48();
  r2= drand48()*TotalProbability;

  *tau= (double) (r1 > TINY ? -log(r1) : -log(TINY))/TotalProbability;
  
  lastreaction=reaction= Reaction;
  sum= reaction->Probability;

  while(sum<r2){
    reaction     = reaction->NextReaction;

    if(reaction==NULL){
      if((TotalProbability-sum)<1e-6) {break;}
      fprintf(stderr,"%s: SelectReaction() found inconsistent reaction probabilities. Premature end.\n",
	      progid);
      exit(-1);
    }

    sum += reaction->Probability;

    if((sum-TotalProbability)>1e-5){
      fprintf(stderr,"%s: SelectReaction() found inconsistent reaction probabilities (%e, %e). Too Much.\n",
	      progid,sum,TotalProbability);
      exit(-1);
    }
    lastreaction = reaction;
  }
  
  return(lastreaction);
}

#undef TINY 

void ExecuteReaction(reaction)
REACTION *reaction;
{

  (void) reaction->ReactionFunc( (void *) reaction->ReactionData);

}

void FreeReactionQueue()
{
  REACTION *rptr;
  void FreeReactionData();

  TotalProbability= 0.0;

  while(Reaction!=NULL){
    rptr= Reaction->NextReaction;
    if(Reaction->ReactionData!=NULL)
      FreeReactionData(Reaction->Type,Reaction->ReactionData);
    FreeReaction(Reaction);
    Reaction=rptr;
  }

  NReactions=0;

}


void FreeReactionData(type,rdata)
int type;
void *rdata;
{
  REACTDATA    *react;
  MOVERIBO     *ribo;
  MOVERNAP     *rnap;

  switch(type){
  case Reaction_Type_TransInit:
    /*** Don't free rdata, it points to DNA ***/
    return;
  case Reaction_Type_ChangeCellVolume:
    /*** This Reaction has no data ****/
    return;
  case Reaction_Type_MoveRNAP:         
  case Reaction_Type_NextSegment:
  case Reaction_Type_RNAP_RNAP:        
  case Reaction_Type_DNAAction:        
    rnap= (MOVERNAP *) rdata;
    FreeMRNAP(rnap);
    return;
  case Reaction_Type_Kinetic:
    react= (REACTDATA *) rdata;
    free(react);
    return;
  case Reaction_Type_EatmRNA:          
  case Reaction_Type_MoveRibosome:     
  case Reaction_Type_ProduceProtein:   
  case Reaction_Type_ProduceNewProtein:
  case Reaction_Type_BindRibosome:     
    ribo= (MOVERIBO *) rdata;
    FreeMRibosome(ribo);
    return;
  default:
    fprintf(stderr,"%s: Trying to free unknown reaction data type %d\n",progid,type);
    exit(-1);
  }
}






