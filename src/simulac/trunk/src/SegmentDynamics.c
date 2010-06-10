/*******************
 *
 * These routines are for the 
 * creation, submission and execution of 
 * DNA reactions. These reactions include
 * (but are not limited to):
 * 
 * 1) RNAP movement and possibly concurrent mRNA elongation
 * 2) ribsome binding to and motion on elongating mRNA
 * 3) mRNA degradation by RNAse
 * 4) attenuation sites in DNA
 * 5) protein production by mRNA/ribosome complexes
 * 6) attenuation by convergent transcription
 * 7) mRNA inhibition by antisense mRNA
 * 8) RNAP modification (antitermination)
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
/**** RNAP Motion Routines ****/
/******************************/

/***********************
 *
 * This Routine cycles over
 * all Sequences to:
 *
 * 1) Move the RNAP's down the genes
 * 2) Move all Ribosomes down the bound transcripts
 * 3) Move all Ribosomes down the free  transcripts
 *
 ************************/


void Polymerize()
{
  int         i;
  DNA        *dna;
  RNAP       *queue;
  mRNA       *trans;
  SEGMENT    *seg;
  
  void MoveRNAPs();
  void MoveRibosomes();



  for(i=0; i<NSequences; i++){
    
    dna= &Sequence[i]; /* Sequences Start with Left-Most Member of Structure */
    
    while(dna!=NULL){
      /* Process RNAP actions for this segment */
      MoveRNAPs(dna);
      dna= dna->RightSegment;
    }

    dna= &Sequence[i]; /* Sequences Start with Left-Most Member of Structure */
    while(dna!=NULL){ 

      if(dna->Type == DNA_Type_Coding){ /*** Only Coding Segments with have polyribosomes attached ****/
	seg= (SEGMENT *) dna->DNAStruct;
	queue= (RNAP *) seg->RNAPQueue;
   
	while(queue!=NULL){

	  if(queue->Transcript!=NULL)
	    if(queue->Transcript->Type== mRNA_Type_Sense) /*** Assume AntiSense doesn't have RBS ? ****/
	      MoveRibosomes(queue->Transcript);

	  queue= queue->NextRNAP;
	}
      }
      dna= dna->RightSegment;
    }
  }

  /* Move ribosomes down the free transcripts */
  trans = Transcript;
      
  while (trans != NULL) {
    /* Save the next transcript to look at */
    mRNA *next = trans->NextTranscript;    

    /* Move the ribosomes (might free up memory for this transcript) */
    MoveRibosomes(trans);

    /* Move on to next transcript */
    trans = next;
  }
}
       
/************************
 *
 * This routine is for the simple
 * stepping of RNAP down a DNA 
 * chain. 
 * 
 * It takes into account:
 *
 *   1) Blocking by an RNAP just ahead of the propagating RNAP
 *   2) Collision with a counter-propagating RNAP on the anti-sense strand
 *   3) Possible Movement to next functional segment of DNA
 *
 *   (Elongation of mRNA is taken care of by the reaction execution routine)
 **************************/

void MoveRNAPs(dna)
DNA *dna;
{
  RNAP *queue,*rnap1,*rnap2;
  int   blockflag,eflag;
  SEGMENT  *seg;
  PROMOTOR *prom;

  void SubmitSimpleTranscription();
  void SubmitConvergentTranscription();
  void SubmitSimpleJumpSegment();
  void PromotorAction();

  /**** First we must look for transcription initiations *****/

  if(dna->Type== DNA_Type_Promotor){
    PromotorAction(dna);
  }
    
  if(dna->Type!= DNA_Type_Promotor){
    seg= (SEGMENT *) dna->DNAStruct;
    queue= (RNAP *) seg->RNAPQueue;
  }else{
    prom= (PROMOTOR *) dna->DNAStruct;
    queue= (RNAP *) prom->RNAPQueue;
  }

  if(queue==NULL) return;

  /**** Loop through RNAP's on this segment ************/

  rnap1=queue;
  while(rnap1!=NULL){
    
    /***** Check to see if we're at end of Segment *******/

    DEBUG(100)
    fprintf(stderr,"@@@@ RNAP: GENE= %s [%s =? %s], cp = %d lp= %d\n",
	    dna->Name,
	    (rnap1->Direction == LEFT ? "LEFT" : "RIGHT"),
	    (dna->Direction == LEFT ? "LEFT" : "RIGHT"),
	    rnap1->CurrentPosition,dna->Length);
    
    eflag=0;
    if(rnap1->Direction==dna->Direction){ 
      if(rnap1->CurrentPosition == dna->Length){

	/*################# NEED TO CHECK FOR COLLISIONS BETWEEN RNAPS ON DIFFERENT SEGMENTS ##############*/
	eflag=1;

	if(dna->Type != DNA_Type_Promotor){	  
	  seg= (SEGMENT *) dna->DNAStruct;
	  seg->SegmentFunc(dna,rnap1);
	  rnap1=rnap1->NextRNAP;
	  continue;
	} else {
	  /* Note, this is for RNAP on strand complementary to actual promotor */
	  SubmitSimpleJumpSegment(dna,rnap1); 
	  rnap1=rnap1->NextRNAP;
	  continue;
	}
      }
    } else 
      if(rnap1->CurrentPosition == 1) { /* Reached End of Wrong-Way Transcription */
	eflag=1;
	if(dna->Type != DNA_Type_Promotor){
	  seg= (SEGMENT *) dna->DNAStruct;
	  seg->SegmentFunc(dna,rnap1);              /** SegmentFunc also must detect direct of RNAP Motion **/
	  rnap1=rnap1->NextRNAP;
	  continue;	  
	} else {
	  /* Note, this is for RNAP on strand complementary to actual promotor */
	  SubmitSimpleJumpSegment(dna,rnap1); 
	  rnap1=rnap1->NextRNAP;
	  continue;	  
	}
      }


    if(eflag==0){ /*** We are in the middle of a segment *****/

      /***** Check if we are blocked by slow RNAP    *******/
      /***** Check for counter propagating RNAP      *******/

      blockflag=0;
      rnap2= queue;

      while(rnap2!=NULL){
	if(rnap2 == rnap1){
	  if(rnap2==rnap2->NextRNAP){
	    fprintf(stderr,"%s: Circular linkage of RNAP on fragment %s\n",progid,dna->Name);
	    exit(-1);
	  }
	  rnap2=rnap2->NextRNAP; 
	  continue;
	}
	
	if(abs(rnap2->CurrentPosition - rnap1->CurrentPosition) <= 17){ 
	  /*###### This assumes we will pick up rnap2 later ###########*/	  
	  /* If counterpropagating then collision */
	  if(rnap2->Direction != rnap1->Direction){
	    SubmitConvergentTranscription(dna,rnap1);
	    blockflag=1;
	    break;
	  } 
	  else 
	    if(rnap1->CurrentPosition<rnap2->CurrentPosition){
	      blockflag=1;
	      break;
	    }
	}
	rnap2= rnap2->NextRNAP;
      }
      
      if(blockflag==0){
	SubmitSimpleTranscription(dna,rnap1);
      }
    }
    rnap1= rnap1->NextRNAP;
  }

}

/*********************************/
/***** RNAP Reaction Submitters **/
/*********************************/

void SubmitSimpleTranscription(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  MOVERNAP *rdata;
  REACTION *reaction;
  void SimpleRNAPMover();
  void SubmitReaction();

  rdata= (MOVERNAP *) AllocMRNAP();

  rdata->dna=dna;
  rdata->rnap1=rnap;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_MoveRNAP;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= SimpleRNAPMover;
  reaction->Probability=  Rate_Of_Polymerase_Motion;

  SubmitReaction(reaction);  
}

void SubmitConvergentTranscription(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  MOVERNAP *rdata;
  REACTION *reaction;
  void SimpleRNAPMover();
  void RNAPFallsOff();
  void SubmitReaction();

  /* This submits two reactions:
   * One in which the RNAP falls off due to the collision.
   * The other in which the RNAP makes it through the collision
   * without mishap.
   */

/*######## Note that this doesn't allow BOTH RNAP's to fall off at once #########*/

  DEBUG(20)
  fprintf(stderr,"@@@@@@ Collision!!!! Time= %e : %s going Lorax is on %s\n",Time,
	  (rnap->Direction == LEFT ? "LEFT" : "RIGHT"),dna->Name);
  
  /**** Makes it through ********/

/****** NO ESCAPE ALLOWED!!!!! ******

  rdata= (MOVERNAP *) AllocMRNAP();

  rdata->dna=dna;
  rdata->rnap1=rnap;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_MoveRNAP;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= SimpleRNAPMover;
  reaction->Probability=  Rate_Of_RNAP_Collision_Escape;

  SubmitReaction(reaction);  

********/

  /**** Falls Off ********/

  rdata= (MOVERNAP *) AllocMRNAP();

  rdata->dna=dna;
  rdata->rnap1=rnap;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_RNAP_RNAP;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= RNAPFallsOff;
  reaction->Probability=  Rate_Of_RNAP_Collision_Failure;

  SubmitReaction(reaction);  
}

void SubmitSimpleJumpSegment(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  REACTION *reaction;
  MOVERNAP *rdata;
  void SimpleJumpSegment();
  void SubmitReaction();

  rdata= (MOVERNAP *) AllocMRNAP();

  rdata->dna=dna;
  rdata->rnap1=rnap;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_NextSegment;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= SimpleJumpSegment;
  reaction->Probability=  Rate_Of_Polymerase_Motion;

  SubmitReaction(reaction);  
}

void SubmitTermination(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  int i,antiterm=0;
  MOVERNAP *rdata;
  TERMDATA *tdata;
  REACTION *reaction=NULL;
  SEGMENT  *seg;
  void SubmitSimpleJumpSegment();
  void RNAPFallsOff();
  void SimpleJumpSegment();
  void SubmitReaction();

  /***** Assumes that Termination only works in ONE direction of transcription *******/


  if(rnap->Direction== dna->Direction){


    seg= (SEGMENT *) dna->DNAStruct;
    tdata= (TERMDATA *) seg->SegmentData;
    
    /**** Check Termination State of RNAP *****/
    
    if(tdata->SpeciesIndex!= -1){

      if(rnap->NBound>0 && rnap->SpeciesIndex!=NULL){
	for(i=0; i<rnap->NBound; i++)
	  if(rnap->SpeciesIndex[i]==tdata->SpeciesIndex) break;
    
	if(i==rnap->NBound) antiterm=0;
	else                antiterm=1;
      } else antiterm=0;
    }else antiterm=0;

    if(antiterm==0){
      // RNAP is not antiterminated - use BaseFallOffRate
      if(tdata->BaseFallOffRate !=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
	
	rdata->dna=dna;
	rdata->rnap1=rnap;
	
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= RNAPFallsOff;
	reaction->Probability=  tdata->BaseFallOffRate;
	
	SubmitReaction(reaction);  
      }
      
      if(tdata->BaseRNAPMotion!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
	
	rdata->dna=dna;
	rdata->rnap1=rnap;
	
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= SimpleJumpSegment;
	reaction->Probability=  tdata->BaseRNAPMotion;
	
	SubmitReaction(reaction);  
      }
    } else {
      
      // RNAP is antiterminated => use AntiTerminatedFallOffRate
      if(tdata->AntiTerminatedFallOffRate!=0.0){

	rdata= (MOVERNAP *) AllocMRNAP();
	
	rdata->dna=dna;
	rdata->rnap1=rnap;
	
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= RNAPFallsOff;
	reaction->Probability=  tdata->AntiTerminatedFallOffRate;
	
	SubmitReaction(reaction);  
      }
      
      if(tdata->AntiTerminatedRNAPMotion!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();

	rdata->dna=dna;
	rdata->rnap1=rnap;
	
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= SimpleJumpSegment;
	reaction->Probability=  tdata->AntiTerminatedRNAPMotion;
	
	SubmitReaction(reaction);  
      }
      
    }
  } else { /* Propagating a different direction */
    DEBUG(50)
    fprintf(stderr,"@@@@@@ Wrong Way Terminator:  %s\n",dna->Name);
    SubmitSimpleJumpSegment(dna,rnap);
  }
    
}

void SubmitAntiTermination(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  int          i,antiterm=0;
  MOVERNAP     *rdata;
  ANTITERMDATA *tdata;
  REACTION *reaction;
  SEGMENT  *seg;
  void SubmitSimpleJumpSegment();
  void AntiTerminateRNAP();
  void UnAntiTerminateRNAP();
  void SimpleJumpSegment();
  void SubmitReaction();

  /***** Assumes AntiTermination Only works in One direction *********/

  if(rnap->Direction==dna->Direction){

    seg= (SEGMENT *) dna->DNAStruct;
    tdata= (ANTITERMDATA *) seg->SegmentData;

    if(rnap->SpeciesIndex!=NULL){
      for(i=0; i<rnap->NBound; i++)
	if(rnap->SpeciesIndex[i]==tdata->SpeciesIndex) break;
      
      if(i==rnap->NBound) antiterm=0;
      else                antiterm=1;
    } else antiterm=0;

      if(antiterm==0){

      if(tdata->UnBoundRNAPMotion!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
      
	rdata->dna=dna;
	rdata->rnap1=rnap;
      
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= SimpleJumpSegment;
	reaction->Probability=  tdata->UnBoundRNAPMotion;
      

	SubmitReaction(reaction);  
      }

      if(tdata->BindingRate!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
      
	rdata->dna=dna;
	rdata->rnap1=rnap;
      
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= AntiTerminateRNAP;
	/* 
	   Note that this is a de Facto bimolecular reaction...it must be modified by a volume
	   element
	 */
	   
	reaction->Probability=  tdata->BindingRate*Concentration[tdata->SpeciesIndex]*(EColi->V0/EColi->V);
      
	SubmitReaction(reaction);  
      }

    } else {
      
      if(tdata->BoundRNAPMotion!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
      
	rdata->dna=dna;
	rdata->rnap1=rnap;
      
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= SimpleJumpSegment;
	reaction->Probability=  tdata->BoundRNAPMotion;
    
	SubmitReaction(reaction);  
      }
      
      if(tdata->UnBindingRate!=0.0){
	rdata= (MOVERNAP *) AllocMRNAP();
      
	rdata->dna=dna;
	rdata->rnap1=rnap;
      
	reaction= (REACTION *)  AllocReaction();
	reaction->Type=         Reaction_Type_DNAAction;
	reaction->ReactionData= (void *) rdata;
	reaction->ReactionFunc= UnAntiTerminateRNAP;
	reaction->Probability=  tdata->UnBindingRate;
      
	SubmitReaction(reaction);  
      }
    }
  } else { /* Counter Propagating */

    SubmitSimpleJumpSegment(dna,rnap);
  }
  
  
}

void SubmitProduceTranscript(dna,rnap)
DNA  *dna;
RNAP *rnap;
{
  int i;
  void SubmitSimpleJumpSegment();
  void SubmitReaction();

/**** This function both submits a movement of RNAP 
 **** to the next segment AND transfers the now complete 
 **** transcript to the UnboundTranscriptQueue
 ****/

  /*** Note that Transcript might already have been released 
   *** We don't do it twice!
   ***/

  if(rnap->Transcript!=NULL){
    NTranscripts++;

    /* Put RNAP Transcript at beginning of Queue */

    rnap->Transcript->CurrentLength++; /* Finish making Transcript */


    if(Transcript!=NULL){
      Transcript->LastTranscript=rnap->Transcript;
      rnap->Transcript->NextTranscript=Transcript;
      rnap->Transcript->LastTranscript=NULL;
      Transcript=rnap->Transcript;
    } else {  
      Transcript=rnap->Transcript;
      Transcript->NextTranscript=NULL;
      Transcript->LastTranscript=NULL;
    }
    Transcript->Rnap=NULL;
    /* Cleave Transcript from RNAP */

    rnap->Transcript=NULL;
  }

  /*
  for(i=0; i<NSequences; i++){
    DNA *sequence;
    sequence= &Sequence[i];
    while(sequence!=NULL){	
      if(sequence->DNAStruct==NULL)
	fprintf(stderr,"@@@ Time= %e. Null data in sequence %s\n",Time,sequence->Name);
	sequence= sequence->RightSegment;
    }
  }
  */

  SubmitSimpleJumpSegment(dna,rnap);  
}

/***********************************/
/****** RNAP Reactions *************/
/***********************************/

void SimpleRNAPMover(rdata)
void *rdata;
{

  MOVERNAP *data;
  RNAP *rnap;
  DNA  *dna;

  DEBUG(50)
  fprintf(stderr,"@@@ SimpleRNAPMover()\n");
  data= (MOVERNAP *) rdata;

  rnap= data->rnap1;
  dna=  data->dna;

  if(rnap->Direction==dna->Direction)
    rnap->CurrentPosition += 1;
  else
    rnap->CurrentPosition -= 1;
  
  if(dna->Type == DNA_Type_Coding){
    if(rnap->Transcript== NULL){
      rnap->Transcript= (mRNA *) rcalloc(1,sizeof(mRNA),"SimpleRNAPMover");
      rnap->Transcript->Gene=dna;
      rnap->Transcript->Rnap=rnap;
      rnap->Transcript->CurrentLength=2;
      rnap->Transcript->Type= (dna->Direction == rnap->Direction ? mRNA_Type_Sense : mRNA_Type_AntiSense);
    } else {
      rnap->Transcript->CurrentLength++;
    }
  }

}

void RNAPFallsOff(rdata)
void *rdata;
{
  int i;
  MOVERNAP *data;
  RNAP *rnap,*queue;
  DNA  *dna;
  SEGMENT  *seg;
  PROMOTOR *prom;



  data= (MOVERNAP *) rdata;

  rnap= data->rnap1;
  dna=  data->dna;

  DEBUG(50)
  fprintf(stderr,"@@@ RNAPFallsOff(): %s going Lorax falls off %s\n",
	  (rnap->Direction==LEFT ? "LEFT": "RIGHT"), dna->Name);

  Concentration[0] += 1; /* Return Polymerase to Pool */
  
  /* Free Memory */

  if(dna->Type!= DNA_Type_Promotor){
    seg= (SEGMENT *) dna->DNAStruct;

    if(rnap->LastRNAP==NULL){
      seg->RNAPQueue=rnap->NextRNAP;
      if(seg->RNAPQueue!=NULL) seg->RNAPQueue->LastRNAP=NULL;
    }else{
      queue=rnap->LastRNAP;      
      queue->NextRNAP= rnap->NextRNAP;
      if(rnap->NextRNAP!=NULL) rnap->NextRNAP->LastRNAP=queue;
    }
    
  } else {
    prom= (PROMOTOR *) dna->DNAStruct;
    if(rnap->LastRNAP==NULL){
      prom->RNAPQueue=rnap->NextRNAP;
      if(prom->RNAPQueue!=NULL) prom->RNAPQueue->LastRNAP=NULL;
    } else{
      queue=rnap->LastRNAP;      
      queue->NextRNAP= rnap->NextRNAP;
      if(rnap->NextRNAP!=NULL) rnap->NextRNAP->LastRNAP=queue;
    }  
  }

  if(rnap->SpeciesIndex!=NULL){
  
    /* Return Bound Modifiers to Pool */
    
    for(i=0; i<rnap->NBound; i++)
      Concentration[rnap->SpeciesIndex[i]]++;
    
    free(rnap->SpeciesIndex);
    rnap->SpeciesIndex = NULL;
  }

  if(rnap->Transcript!=NULL){
    RIBOSOME *rqueue;

    rqueue= (RIBOSOME *) rnap->Transcript->RiboQueue;
    while(rqueue!=NULL){
      
      /* Now free ribosome structure */
      if(rqueue->SpeciesIndex!=NULL){
	for(i=0; i<rqueue->NBound; i++)
	  Concentration[rqueue->SpeciesIndex[i]]++;
	free(rqueue->SpeciesIndex);
	rqueue->SpeciesIndex = NULL;
      }

      if(rqueue->NextRibosome==NULL){
	rnap->Transcript->RiboQueue= NULL;
	FreeRibosome(rqueue);
	break;
      } else {
	rnap->Transcript->RiboQueue= rqueue->NextRibosome;
	FreeRibosome(rqueue);
      }
      rqueue= rnap->Transcript->RiboQueue;
    }   
    free(rnap->Transcript);
    rnap->Transcript = NULL;
  }

  FreeRNAP(rnap);

}

void SimpleJumpSegment(rdata)
void *rdata;
{
  int i;
  MOVERNAP *data;
  RNAP *rnap,*queue;
  DNA  *dna;
  SEGMENT  *seg;
  PROMOTOR *prom;


  data= (MOVERNAP *) rdata;

  rnap= data->rnap1;

  /* Remove rnap from the original queue */
  
  dna= data->dna;

  DEBUG(20)
  fprintf(stderr,"@@@@ RNAP Jumping from %s\n",dna->Name);

  if(dna->Type!= DNA_Type_Promotor){
    seg= (SEGMENT *) dna->DNAStruct;
    if(rnap->LastRNAP==NULL){
      seg->RNAPQueue=rnap->NextRNAP;
      if(seg->RNAPQueue!=NULL) seg->RNAPQueue->LastRNAP=NULL;
    } else{
      queue=rnap->LastRNAP;   
      queue->NextRNAP= rnap->NextRNAP;
      if(rnap->NextRNAP!=NULL) rnap->NextRNAP->LastRNAP=queue;
    }
  } else {
    prom= (PROMOTOR *) dna->DNAStruct;
    if(rnap->LastRNAP==NULL){
      prom->RNAPQueue=rnap->NextRNAP;
      if(prom->RNAPQueue!=NULL) prom->RNAPQueue->LastRNAP=NULL;
    } else{
      queue=rnap->LastRNAP;      
      queue->NextRNAP= rnap->NextRNAP;
      if(rnap->NextRNAP!=NULL) rnap->NextRNAP->LastRNAP=queue;
    }
  }
  
  /**** Add to Next Queue *******/

  if(rnap->Direction == LEFT)
    dna=  data->dna->LeftSegment;
  else
    dna=  data->dna->RightSegment;
  
  /*** Add RNAP to next segment's RNAP queue *****/
  if(dna!=NULL){ /* If not at END of segment */

    DEBUG(20)
    fprintf(stderr,"@@@ New RNAP jumped to segment %s\n",dna->Name);
    
    if(dna->Type!= DNA_Type_Promotor){
      seg= (SEGMENT *) dna->DNAStruct;
      if(seg->RNAPQueue != NULL){
	queue= (RNAP *) seg->RNAPQueue;
	rnap->NextRNAP=queue;
	queue->LastRNAP=rnap;
	seg->RNAPQueue=rnap;
	rnap->LastRNAP=NULL;
      } else {
	rnap->LastRNAP=NULL;
	rnap->NextRNAP=NULL;
	seg->RNAPQueue= rnap;
      }
    } else {

      prom= (PROMOTOR *) dna->DNAStruct;

      if(prom->RNAPQueue != NULL){

	queue= (RNAP *) prom->RNAPQueue;

	rnap->NextRNAP=queue;
	queue->LastRNAP=rnap;
	prom->RNAPQueue=rnap;
	rnap->LastRNAP=NULL;
      } else {
	
	rnap->LastRNAP=NULL;
	rnap->NextRNAP=NULL;
	prom->RNAPQueue=rnap;
      }
    }
          
    /*** Reset Relative Position ****/
    
    if(rnap->Direction != dna->Direction)
      rnap->CurrentPosition= dna->Length;
    else 
      rnap->CurrentPosition= 1;

    /***** Gotta deal with Transcript ******/
    /***** Transcripts can't be carried
    ****** across segments in this version
    ****** If RNAP has transcript and it hasn't been
    ****** freed for some reason we have an error.
    ****** else if it needs a transcript, give one to it.
    ******/

    if(rnap->Transcript!=NULL){
      fprintf(stderr,"@@@@ RNAP Stuck to Transcript!\n");
      exit(-1);
    }
    
    data->dna=dna;
    data->rnap1=rnap;
    SimpleRNAPMover(data);

  } else {
    
    /*** free RNAP ****/
     
    DEBUG(20)
      fprintf(stderr,"@@@ Freeing RNAP\n");

     if(rnap->SpeciesIndex!=NULL){
       
       for(i=0; i<rnap->NBound; i++)
	Concentration[rnap->SpeciesIndex[i]]++;
       
       free(rnap->SpeciesIndex);
       rnap->SpeciesIndex = NULL;
     }
    
     if(rnap->Transcript!=NULL){
       RIBOSOME *rqueue;

       fprintf(stderr,"@@@ Removing a weird transcript\n");
       rqueue= (RIBOSOME *) rnap->Transcript->RiboQueue;
       while(rqueue!=NULL){
	 
	 /* Now free ribosome structure */
	 if(rqueue->SpeciesIndex!=NULL){
	   for(i=0; i<rqueue->NBound; i++)
	     Concentration[rqueue->SpeciesIndex[i]]++;
	   
	   if(rqueue->SpeciesIndex!=NULL) {
	     free(rqueue->SpeciesIndex);
	     rqueue->SpeciesIndex = NULL;
	   }

	 }
	 
	 if(rqueue->NextRibosome==NULL){
	   FreeRibosome(rqueue);
	   break;
	 } else {
	   rqueue=rqueue->NextRibosome;
	   FreeRibosome(rqueue->LastRibosome);
	 }
       }         
       fprintf(stderr,"%%%%%% Delete Transcript\n");
       free(rnap->Transcript);
       rnap->Transcript = NULL;
     }
     
     FreeRNAP(rnap);
     
   }

}
  
void AntiTerminateRNAP(rdata)
void *rdata;
{
  int i;
  MOVERNAP *data;
  RNAP *rnap;
  DNA  *dna;
  ANTITERMDATA *tdata;
  SEGMENT *seg;
  void *rrealloc();

  data= (MOVERNAP *) rdata;
  dna=  data->dna;
  rnap= data->rnap1;

  DEBUG(50)
    fprintf(stderr,"@@@ AntiTerminateRNAP()\n");
  seg= (SEGMENT *) dna->DNAStruct;
  tdata= (ANTITERMDATA *) seg->SegmentData;

  /**** Check Termination State of RNAP *****/

  if(rnap->SpeciesIndex!=NULL){
    for(i=0; i<rnap->NBound; i++)
      if(rnap->SpeciesIndex[i]==tdata->SpeciesIndex) break;

    if(i!=rnap->NBound){
      fprintf(stderr,"%s: AntiTerminateRNAP() passed an already antiterminated polymerase!.\n",progid);
      exit(-1);
    }
  }

  if(rnap->NBound==0){
    if(rnap->SpeciesIndex==NULL) rnap->SpeciesIndex= (int *) rcalloc(1,sizeof(int),"AntiterminateRNAP");
    else                         rnap->SpeciesIndex= (int *) rrealloc(rnap->SpeciesIndex,rnap->NBound+1,sizeof(int),"AntiterminateRNAP");
    
  } else {
    
    if(rnap->SpeciesIndex==NULL){
      fprintf(stderr,"%s: AntiTerminateRNAP() came across a polymerase with a memory problem.\n",progid);
      exit(-1);
    } else rnap->SpeciesIndex= (int *) rrealloc(rnap->SpeciesIndex,rnap->NBound+1,sizeof(int),"AntiterminateRNAP");

  }    
  
  Concentration[tdata->SpeciesIndex] -= 1;
  rnap->SpeciesIndex[rnap->NBound]=tdata->SpeciesIndex;
  rnap->NBound++;

}

void UnAntiTerminateRNAP(rdata)
void *rdata;
{
  int i,j;
  MOVERNAP *data;
  RNAP *rnap;
  DNA  *dna;
  ANTITERMDATA *tdata;
  SEGMENT *seg;
  void *rrealloc();

  DEBUG(50)
    fprintf(stderr,"@@@ UnAntiTerminateRNAP()\n");
  data= (MOVERNAP *) rdata;
  dna=  data->dna;
  rnap= data->rnap1;

  seg= (SEGMENT *) dna->DNAStruct;
  tdata= (ANTITERMDATA *) seg->SegmentData;

  /**** Check Termination State of RNAP *****/
  
  if(rnap->SpeciesIndex==NULL){
    fprintf(stderr,"%s: UnAntiTerminateRNAP() passed an RNAP with a memory problem.\n",progid);
    exit(-1);
  }
  
  for(i=0; i<rnap->NBound; i++)
    if(rnap->SpeciesIndex[i]==tdata->SpeciesIndex) break;
  
  if(i==rnap->NBound){
    fprintf(stderr,"%s: UnAntiTerminateRNAP() passed a non-antiterminated polymerase!.\n",progid);
    exit(-1);
  }
  
  
  if(rnap->NBound==1) {
    free(rnap->SpeciesIndex);
    rnap->SpeciesIndex = NULL;
  } else {
    /***** Compress Array ***/
    if (rnap->NBound>1) {
      for(j=i+1; j<rnap->NBound; j++)
	rnap->SpeciesIndex[j-1]=rnap->SpeciesIndex[j];    

      rnap->SpeciesIndex= (int *) rrealloc(rnap->SpeciesIndex,rnap->NBound-1,sizeof(int),"UnAntiterminateRNAP");
    } else {
      /* RMM: I think this code is unreachable... */
      free(rnap->SpeciesIndex);
      rnap->SpeciesIndex = NULL;
    }
  }
  
  Concentration[tdata->SpeciesIndex] += 1;
  rnap->NBound--;  

}

/***************************************/
/******** Ribsome Dynamics *************/
/***************************************/

/*****************************
 *
 * MoveRibosomes handle all aspects of:
 * 
 * 1) Ribosome Binding to the RBS of the transcript
 * 2) Ribosome Motion down transcript
 * 3) Blockage between transcripts
 * 4) Protein production at end of transcript
 * 5) RNAse chewing of RBS
 * 6) AntiSense Blockage of RBS
 *
 * Assumptions:
 *
 *    1) RBS is exposed when 14 or more nucleotides have be synthesized (14 is arbitrary)
 *    2) Ribosomes cover approximately 10 nucleotides.
 *    3) Translation Rate is approximately 60nt/sec (Cell Biology, p.104)
 *    4) Ribsome Binding Rate is the same for all transcripts
 *    5) mRNA Stabilitity is not the same for all transcripts
 *    6) Once Bound, ribosomes do NOT fall off except at the end.
 *   
 ********************************/

void MoveRibosomes(trans)
mRNA *trans;
{
  int FullLength,CurrLength,cp1;
  RIBOSOME *queue;

  void SubmitClearRBS();
  void SubmitMoveRibosome();
  void SubmitProduceProtein();
  
  if(trans->Type==mRNA_Type_AntiSense) return;
  if(trans->CurrentLength<20) return;

  FullLength= trans->Gene->Length;
  CurrLength= trans->CurrentLength;

  /***** First take care of possible binding events ******/

  /**** Remember, Riboqueues are ALWAYS in order! *****/
  queue= trans->RiboQueue;


  if(trans->Rnap==NULL && trans->RBSState==mRNA_RBS_Chewed && queue==NULL){ 
    /* Its a Free Transcript with no ribosomes bound */
    /* So free it up */
    
    if(trans->LastTranscript==NULL){
      Transcript=trans->NextTranscript;
      if(Transcript!=NULL) Transcript->LastTranscript=NULL;
    } else{
      trans->LastTranscript->NextTranscript=trans->NextTranscript;
      if(trans->NextTranscript!=NULL)
	trans->NextTranscript->LastTranscript=trans->LastTranscript;
    }
    NTranscripts--;
    free(trans);
    return;
  }
  
  
  if(trans->RBSState!=mRNA_RBS_Chewed){
    if(queue==NULL || queue->CurrentPosition>14) /**** RBS IS CLEAR!!!!!! *****/
      SubmitClearRBS(trans);
  } 
  
  /**** Move them puppies *****/

  while(queue!=NULL){
    
    cp1= queue->CurrentPosition;
    
    if(queue->NextRibosome==NULL){ /* Can't collide with next ribosome */
      
      /* Check if at end of transcript */

      if(CurrLength<FullLength){ /* Not a free transcript */
	if(cp1< (CurrLength-5)) /* Not at end of partial transcript */
	  SubmitMoveRibosome(trans,queue);
	/* else We are abutting polymerase */
      } else {                    /* We are a free transcript */
	if(cp1== FullLength)      /* We have reached end of transcript */
	  SubmitProduceProtein(trans,queue);
	else 
	  SubmitMoveRibosome(trans,queue);
      }
    } else { /* Same thing except we check for collisions here */

      if((queue->NextRibosome->CurrentPosition - queue->CurrentPosition)>=10){ /* Then No Collision imminent */

	if(CurrLength<FullLength){ /* Not a free transcript */
	  if(cp1< (CurrLength-5)) /* Not at end of partial transcript */
	    SubmitMoveRibosome(trans,queue);
	  /* else We are abutting polymerase */
	} else {                    /* We are a free transcript */
	  
	  if(cp1== FullLength)      /* We have reached end of transcript */
	    SubmitProduceProtein(trans,queue);
	  else 
	    SubmitMoveRibosome(trans,queue);
	}
	
      }
    }

    /* Next Ribosome */

    queue= queue->NextRibosome;
  }

}

/**************************************************************/
/*********** Ribosome Reaction Submission Dynamics *************/
/**************************************************************/

void SubmitClearRBS(trans)
mRNA    *trans;
{

  
  REACTION   *reaction;
  MOVERIBO   *rdata;
  CODINGDATA *SegData;
  SEGMENT    *seg;

  void SubmitReaction();
  void BindRibosome();
  void EatmRNA();

  rdata= (MOVERIBO *) AllocMRibosome();

  rdata->trans=trans;
  rdata->ribosome=NULL;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_BindRibosome;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= BindRibosome;
# ifdef RMM_MODS
  /* Grab the ribosome binding rate from the coding data */
  seg = (SEGMENT *) trans->Gene->DNAStruct;
  SegData = (CODINGDATA *) seg->SegmentData;
  reaction->Probability = SegData->RibosomeBindingRate * Concentration[1] * 
    (EColi->V0/EColi->V);
# else
  /* DeFacto Bimolecular: Volume Element Included */
  reaction->Probability=  Rate_Of_Ribosome_Binding*Concentration[1]*(EColi->V0/EColi->V);
# endif

  SubmitReaction(reaction);  

  seg= (SEGMENT *) trans->Gene->DNAStruct;
  SegData= (CODINGDATA *) seg->SegmentData;

  rdata= (MOVERIBO *) AllocMRibosome();

  rdata->trans=trans;
  rdata->ribosome=NULL;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_EatmRNA;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= EatmRNA;
  reaction->Probability=  SegData->mRNADegradationRate;

  SubmitReaction(reaction);  
}

void SubmitMoveRibosome(trans,ribosome)
mRNA    *trans;
RIBOSOME *ribosome;
{
  REACTION *reaction;
  MOVERIBO *rdata;
  void SubmitReaction();
  void SimpleRibosomeMover();


  rdata= (MOVERIBO *) AllocMRibosome();

  rdata->trans=trans;
  rdata->ribosome=ribosome;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_MoveRibosome;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= SimpleRibosomeMover;
  reaction->Probability=  Rate_Of_Ribosome_Motion;

  SubmitReaction(reaction); 
}

void SubmitProduceProtein(trans,ribosome)
mRNA    *trans;
RIBOSOME *ribosome;
{
  REACTION *reaction;
  MOVERIBO *rdata;
  void SubmitReaction();
  void ProduceNewProtein();

  rdata= (MOVERIBO *) AllocMRibosome();

  rdata->trans=trans;
  rdata->ribosome=ribosome;

  reaction= (REACTION *)  AllocReaction();
  reaction->Type=         Reaction_Type_ProduceNewProtein;
  reaction->ReactionData= (void *) rdata;
  reaction->ReactionFunc= ProduceNewProtein;
  reaction->Probability=  Rate_Of_Ribosome_Motion;

  SubmitReaction(reaction);  

}

/***************************************************/
/*********** Ribosome Reaction Dynamics *************/
/***************************************************/
  

void SimpleRibosomeMover(rdata)
void *rdata;
{

  RIBOSOME  *ribosome;
  MOVERIBO  *data;


  DEBUG(100)
    fprintf(stderr,"@@@ SimpleRiboMover()\n");
  data= (MOVERIBO *) rdata;

  ribosome= data->ribosome;

  ribosome->CurrentPosition += 1;

}

void BindRibosome(rdata)
void *rdata;
{
  mRNA      *trans;
  RIBOSOME  *ribosome;
  MOVERIBO  *data;
  
  DEBUG(50)
    fprintf(stderr,"@@@ SimpleBindRibo()\n");
  data= (MOVERIBO *) rdata;

  trans=    data->trans;
  ribosome= (RIBOSOME *) AllocRibosome();

  ribosome->CurrentPosition= 1;
  ribosome->Transcript=trans;
  ribosome->LastRibosome=NULL;
  if(trans->RiboQueue!=NULL) trans->RiboQueue->LastRibosome= ribosome;
  ribosome->NextRibosome=trans->RiboQueue;
  trans->RiboQueue=ribosome;

  Concentration[1]--;  

}

void EatmRNA(rdata)
void *rdata;
{
  mRNA      *trans;
  MOVERIBO  *data;
  

  data= (MOVERIBO *) rdata;

  trans=    data->trans;
  DEBUG(50)
    fprintf(stderr,"@@@ EatmRNA %s\n",trans->Gene->Name);
  trans->RBSState= mRNA_RBS_Chewed;

}

void ProduceNewProtein(rdata)
void *rdata;
{
  int         i,produced;
  mRNA       *trans;
  RIBOSOME   *ribosome;
  MOVERIBO   *data;
  CODINGDATA *SegData;
  SEGMENT    *seg;

  DEBUG(50)
    fprintf(stderr,"@@@ ProduceNewProtein()\n");
  data= (MOVERIBO *) rdata;

  trans=    data->trans;
  ribosome= data->ribosome;

  
  seg= (SEGMENT *) trans->Gene->DNAStruct;
  SegData= (CODINGDATA *) seg->SegmentData;

  /** Produce Protein **/

  produced= SegData->SpeciesIndex;
  Concentration[produced]++;

 /** Unbind Ribosome **/

  Concentration[1]++; /* Release Ribosome back to pool */

  if(ribosome->LastRibosome!=NULL)
    ribosome->LastRibosome->NextRibosome= NULL;
  else 
    trans->RiboQueue=NULL;

  if(ribosome->NextRibosome!=NULL){
    fprintf(stderr,"%s: Improperly terminated ribosome queue in ProduceNewProtein().\n",progid);
    exit(-1);
  }
  
  /* Now free ribosome structure */
  
  if(ribosome->SpeciesIndex!=NULL){
    for(i=0; i<ribosome->NBound; i++)
      Concentration[ribosome->SpeciesIndex[i]]++;
    
    free(ribosome->SpeciesIndex);
    ribosome->SpeciesIndex = NULL;
  }

  FreeRibosome(ribosome);
}
