/************************
 * 
 * These subroutines 
 * 
 * 1) Calculate the state of 
 *    promotors. 
 *
 * 2) Do the associated Housekeeping of
 *    adding and subtracting chemical species from
 *    the pool.
 *
 * 3) Submit reactions for the movement of RNAP
 *    from the promotor to the first linked piece of
 *    DNA;
 *
 * 4) Execute such reactions
 *
 ***************************/

#define PromotorDynamics

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

/****************************/
/**** Calculate States ******/
/****************************/

#define kcal_per_joule (0.001/4.184)
#define RT (8.314*310.15*kcal_per_joule)            /* (J/(mol K))*K = J/mol */   
#define Molec_to_Molar (1.0/(6.023e23*EColi->V))     /* 1 Mole/(6.023*10^23 (molecs) * 1.41e-15 (L)) */

void CalculateAckersProbabilities(data,prob)
SHEADATA *data;
double *prob;
{
  register int i,j;
  double tprob;
  int spec,cnt;
  double bico();
  double powm[11];

  tprob=0.0;

  for(i=0; i<11; i++) powm[i]=0.0;

  for(i=0; i<data->NConfigs; i++){
    /*** Preset to exp(-deltaG/RT) when loading operator ***/
    prob[i]= data->DeltaG[i];

    for(j=0; j<data->CList[i][0]; j++){
        spec=data->CList[i][1+2*j];
	cnt =data->CList[i][1+2*j+1];

	prob[i] *= (double) bico(Concentration[spec],cnt);
	if(prob[i]==0.0) break;
	prob[i] *= (powm[cnt] ? powm[cnt] :
		    (powm[cnt]=pow(Molec_to_Molar,cnt)));
    }
    tprob+= prob[i];
  }

  for(i=0; i<data->NConfigs; i++)
    prob[i] /= tprob;
}

int CalculateAckersState(n,prob)
int n;
double *prob;
{
  register int i,j;
  double rndm,running_prob;

  rndm= 1.0-drand48(); /* Interval Now (0,1] instead of [0,1) */
  running_prob=prob[0];

  /**** Roulette Wheel Selection ****/

  j=0;
  i=1;
  while(i<n && running_prob<rndm){
    if(prob[i]>1e-20) j=i;
    running_prob += prob[i];
    i++;
  }

  if(i==n && running_prob<rndm){
    fprintf(stderr,"%s: Total Running State Probability (%e) less than Roulette Selection (%e)\n",
	    progid,running_prob,rndm);
    exit(-1);
  }

return(j);
}

void SetAckersState(data)
SHEADATA *data;
{
  static double *prob=NULL;
  static int MaxConfig=0;
  int i,newstate,species;

  void CalculateAckersProbabilities();
  int  CalculateAckersState();


  /**** Release Current State *******/ 
  /**** (Rapid Equilibrium is assumed ****/

  for(i=0; i<data->CList[data->CurrentState][0]; i++)
    Concentration[data->CList[data->CurrentState][1+2*i]] += data->CList[data->CurrentState][1+2*i+1]; 

  if(prob==NULL){
    prob= (double *) rcalloc(data->NConfigs,sizeof(double),"SetSheaAckersState");
    MaxConfig=data->NConfigs;
  }else
    if(data->NConfigs>MaxConfig){
      prob= (double *) rrealloc((void *) prob,data->NConfigs,sizeof(double),"SetSheaAckersState");
      MaxConfig=data->NConfigs;
    }
  
  CalculateAckersProbabilities(data,prob);
  newstate=CalculateAckersState(data->NConfigs,prob);
    
  
  /**** Bind up new molecules ****/

  for(i=0; i<data->CList[newstate][0]; i++)
    Concentration[data->CList[newstate][1+2*i]] -= data->CList[newstate][1+2*i+1]; 
  

  
  data->CurrentState=newstate;
}

/*****************************************
 *
 * Promotor Action Functions
 *
 * These functions essentially submit RNAP movement reactions 
 * to the reaction queue given the current state of the promotor
 *
 * They accept a promotor state and a promotor activity structure.
 *
 * Next Set of Functions: Promotor Housekeeping functions.
 *
 ******************************************/

void PromotorAction(pfragment)
DNA *pfragment;
{
  int pstate;
  REACTION *reaction;
  PROMOTOR *promotor;
  DNA      *dna;
  int       type;

  void      SubmitReaction();
  char     *PrintDNAType();
  void      InitiateTranscription();

  if(pfragment->Type!= DNA_Type_Promotor){
    fprintf(stderr,"%s: Promotor Action passed an incorrect Data Type %s named %s.\n",progid,PrintDNAType(pfragment->Type),pfragment->Name);
    exit(-1);
  }

  promotor= (PROMOTOR *) pfragment->DNAStruct;

  pstate= Operator[promotor->Data].CurrentState;

  if(promotor->IsoRate[pstate]==0.0) return;

  /* Check here if there is a blocking RNAP ahead */

  if(promotor->TranscriptionDirection == LEFT) dna= pfragment->LeftSegment;
  else dna= pfragment->RightSegment;

  type= dna->Type;
  
  if(type != DNA_Type_Promotor){
    SEGMENT *dstruct;
    RNAP     *queue;
    
    dstruct= (SEGMENT *) dna->DNAStruct;
 
    queue= dstruct->RNAPQueue;
    while(queue!=NULL){

      if(queue->Direction == promotor->TranscriptionDirection &&
	 queue->CurrentPosition<=17)
	return;                        /* Reaction is Blocked */

      queue= queue->NextRNAP;
    }
  }else{
    PROMOTOR *dstruct;
    RNAP     *queue;

    dstruct= (PROMOTOR *) dna->DNAStruct;

    queue= dstruct->RNAPQueue;
    while(queue!=NULL){

      if(queue->Direction == promotor->TranscriptionDirection &&
	 queue->CurrentPosition<=17)
	return;                       /* Reaction is blocked */

      queue= queue->NextRNAP;
    }
  }

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_TransInit;
  reaction->ReactionData= (void *) pfragment;
  reaction->ReactionFunc= InitiateTranscription;
  reaction->Probability=  promotor->IsoRate[pstate];

  SubmitReaction(reaction);
}

/*****************************
 *
 * This executes a transcription
 * initiation event.
 *
 *****************************/

void InitiateTranscription(reactdata)
void *reactdata;
{
  DNA      *pfragment, *dna;
  PROMOTOR *promotor;
  RNAP     *tmprnap;
  int       type;
  MOVERNAP  rdata;
  
  char *PrintDNAType();
  void  ChangePromotorState();
  void  SimpleRNAPMover();

  /* First we recast the void reactdata into a recognizable form */

  pfragment= (DNA *) reactdata;

  type= pfragment->Type;

  if(type!= DNA_Type_Promotor){
    fprintf(stderr,"%s: InitiateTranscription() was passed bad reaction data (type= %s)\n",progid,PrintDNAType(type));
    exit(-1);
  }

  promotor= (PROMOTOR *) pfragment->DNAStruct;
#ifdef RMM_MODS
  /* Keep track of the number of RNA polymerases that have gone by */
  promotor->RNAPCount++;

  DEBUG(20) {
    fprintf(stderr, "Initiating transcription on %s, count = %d\n",
	    promotor->Name, promotor->RNAPCount);
  }
#endif

  if(promotor->TranscriptionDirection == LEFT) dna= pfragment->LeftSegment;
  else dna= pfragment->RightSegment;

  

  DEBUG(20)
    fprintf(stderr,"@@@@@ Transcription Initiation to gene %s\n", dna->Name);

  type= dna->Type;
  
  if(type != DNA_Type_Promotor){
    SEGMENT *dstruct;

    dstruct= (SEGMENT *) dna->DNAStruct;

    if(dstruct->RNAPQueue!=NULL){
      tmprnap=(RNAP *) AllocRNAP();
      dstruct->RNAPQueue->LastRNAP= tmprnap;
      tmprnap->NextRNAP= dstruct->RNAPQueue;
      tmprnap->LastRNAP= NULL;
      tmprnap->Direction=promotor->TranscriptionDirection;
      dstruct->RNAPQueue=tmprnap;
  
      if(dna->Direction == tmprnap->Direction)    /* Make Sure RNAP is at correct end of segment */
	dstruct->RNAPQueue->CurrentPosition=0;
      else
	dstruct->RNAPQueue->CurrentPosition=dna->Length+1; /* This will be decremented by SimpleRNAPMover */
      
      /*** Note that Transcript will be inititiated by SimpleRNAPMover() */
      dstruct->RNAPQueue->Transcript= NULL;      
      
      dstruct->RNAPQueue->NBound=0;
      dstruct->RNAPQueue->SpeciesIndex=NULL;
      ChangePromotorState(promotor);  
    }else {
      tmprnap=dstruct->RNAPQueue= (RNAP *) AllocRNAP();
      dstruct->RNAPQueue->LastRNAP= NULL;
      dstruct->RNAPQueue->NextRNAP= NULL;
      dstruct->RNAPQueue->Direction= promotor->TranscriptionDirection;
      if(dna->Direction==dstruct->RNAPQueue->Direction)
	dstruct->RNAPQueue->CurrentPosition=0;
      else 
	dstruct->RNAPQueue->CurrentPosition=dna->Length+1;
      /*** Note that Transcript will be inititiated by SimpleRNAPMover() */
      dstruct->RNAPQueue->Transcript= NULL;
      dstruct->RNAPQueue->NBound=0;
      dstruct->RNAPQueue->SpeciesIndex=NULL;
      ChangePromotorState(promotor);
    }
  }else{
    PROMOTOR *dstruct;

    dstruct= (PROMOTOR *) dna->DNAStruct;
    if(dstruct->RNAPQueue!=NULL){
      tmprnap=(RNAP *) AllocRNAP();
      dstruct->RNAPQueue->LastRNAP= tmprnap;
      tmprnap->NextRNAP= dstruct->RNAPQueue;
      tmprnap->LastRNAP= NULL;
      tmprnap->Direction=promotor->TranscriptionDirection;
      dstruct->RNAPQueue=tmprnap;
      if(dna->Direction==dstruct->RNAPQueue->Direction)
	dstruct->RNAPQueue->CurrentPosition=0;
      else 
	dstruct->RNAPQueue->CurrentPosition=dna->Length+1;
      dstruct->RNAPQueue->Transcript= NULL;
      dstruct->RNAPQueue->NBound=0;
      dstruct->RNAPQueue->SpeciesIndex=NULL;
      ChangePromotorState(promotor);
    }else {
      tmprnap=dstruct->RNAPQueue= (RNAP *) AllocRNAP();
      dstruct->RNAPQueue->LastRNAP= NULL;
      dstruct->RNAPQueue->NextRNAP= NULL;
      dstruct->RNAPQueue->Direction= promotor->TranscriptionDirection;  
      if(dna->Direction==dstruct->RNAPQueue->Direction)
	dstruct->RNAPQueue->CurrentPosition=0;
      else 
	dstruct->RNAPQueue->CurrentPosition=dna->Length+1;
      dstruct->RNAPQueue->Transcript= NULL;
      dstruct->RNAPQueue->NBound=0;
      dstruct->RNAPQueue->SpeciesIndex=NULL;
      ChangePromotorState(promotor);
    }
  }  
  

rdata.rnap1=tmprnap;
rdata.dna=dna;
SimpleRNAPMover(&rdata);
}  
  

void ChangePromotorState(promotor)
PROMOTOR *promotor;
{
  int ns,nc,cs,i,j,k,spec,cnt,nspec,splist[20],ctlist[20];
  int **clist;

  /* 
   * 
   * 
   * 
   * Since we will be recalculating the promotor 
   * state anyway, all we need to do is to put ourselves
   * into a new state identical to our own minus a polymerase.
   * This make sure we got no extra polymerases hanging around.
   * Its kludgy.
   *
   */

  nc   = Operator[promotor->Data].NConfigs;
  cs   = Operator[promotor->Data].CurrentState;
  clist= Operator[promotor->Data].CList;
  
  nspec=0;
  for(i=0; i<clist[cs][0]; i++){
    spec= clist[cs][1+2*i];
    cnt=  clist[cs][1+2*i+1];
    if(spec==0) /* Then its RNAP */
      cnt--;    /* Subtract 1 */

    if(!cnt) continue;
    nspec++;
    splist[nspec-1]=spec;
    ctlist[nspec-1]=cnt;
  }

  for(i=0; i<nc; i++){ /* Search the config table */

    if(clist[i][0]!=nspec) continue; /* definitely no match */

    for(j=0; j<nspec; j++){
      for(k=0; k<nspec; k++)
	if(clist[i][1+2*k]==splist[j] && clist[i][1+2*k+1]==ctlist[j])
	  break;
      if(k==nspec) break; /* No match */
    }
    
    if(j==nspec) break; /* We gotta match (see, 'cause k always broke early) */
  }
  
  if(i==nc){
    fprintf(stderr,"@@@@ No Matching Configuration!");
    exit(-1);
  }
  Operator[promotor->Data].CurrentState= i;
} 
    

