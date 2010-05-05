/*******************
 * 
 * The memory management routines
 *
 ********************/

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

/* 
#ifndef MEMORY
 #include "Memory.h"
#endif
*/

/* Size of some memory preallocation blocks */
#define MEM_BLOCK_REACTION  500
#define MEM_BLOCK_RNAP      150      
#define MEM_BLOCK_RIBOSOME  150
#define MEM_BLOCK_MOVERNAP  150
#define MEM_BLOCK_MOVERIBO  150

int      react_mptr_full= -1;
REACTION *ReactionMemory[MEM_BLOCK_REACTION];

void FillReactionBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_REACTION; i++)
    ReactionMemory[i]= (REACTION *) 
      rcalloc(1,sizeof(REACTION),"FillReactionBlock");
}

REACTION *AllocReaction()
{
  react_mptr_full++;

  if(react_mptr_full>=MEM_BLOCK_REACTION) {
    /* Warn the user if he or she is interested */
    if (DebugLevel > 4)
      fprintf(stderr, "AllocReaction: allocating reaction %d\n",
	      react_mptr_full);

    return( (REACTION *) rcalloc(1,sizeof(REACTION),"AllocReaction"));
  }

  /* Take one of the reactions out of the block */
  REACTION *reaction = ReactionMemory[react_mptr_full];
  ReactionMemory[react_mptr_full] = NULL;

  return(reaction);
}

void FreeReaction(react)
REACTION *react;
{

  if(react_mptr_full>=MEM_BLOCK_REACTION) free(react);
  else
    /* Put the reaction back on the list for later use */
    ReactionMemory[react_mptr_full]=react;
  
    react_mptr_full--;
}

void EmptyReactionBlock()
{
  int i, cnt = 0;

  for(i=0; i<MEM_BLOCK_REACTION; i++)
    if(ReactionMemory[i]!=NULL) {
      free(ReactionMemory[i]);
      ++cnt;
    }

  if (DebugLevel > 3 && cnt > 0) 
    fprintf(stderr, "EmptyReactionBlock: freed %d reactions\n", cnt);
}


int      rnap_mptr_full= -1;
RNAP *RNAPMemory[MEM_BLOCK_RNAP];

void FillRNAPBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_RNAP; i++)
    RNAPMemory[i]= (RNAP *) rcalloc(1,sizeof(RNAP),"FillRNAPBlock");
}

RNAP *AllocRNAP()
{

  rnap_mptr_full++;

  if(rnap_mptr_full>=MEM_BLOCK_RNAP) {
    /* Warn the user if he or she is interested */
    if (DebugLevel > 4)
      fprintf(stderr, "AllocRNAP: allocating RNAP %d\n",
	      react_mptr_full);

    return( (RNAP *) rcalloc(1,sizeof(RNAP),"AllocRNAP"));
  }

  /* Take one of the RNAPs out of the block */
  RNAP *rnap = RNAPMemory[rnap_mptr_full];
  RNAPMemory[rnap_mptr_full] = NULL;

  return(rnap);
}

void FreeRNAP(rnap)
RNAP *rnap;
{


  if(rnap_mptr_full>=MEM_BLOCK_RNAP) free(rnap);
  else
    RNAPMemory[rnap_mptr_full]=rnap;

  rnap_mptr_full--;

}

void EmptyRNAPBlock()
{
  int i, cnt = 0;

  for(i=0; i<MEM_BLOCK_RNAP; i++)
    if(RNAPMemory[i]!=NULL) {
      free(RNAPMemory[i]);
      ++cnt;
    }

  if (DebugLevel > 3 && cnt > 0) 
    fprintf(stderr, "EmptyRNAPBlock: freed %d RNAPs\n", cnt);
}

int      ribo_mptr_full= -1;
RIBOSOME *RibosomeMemory[MEM_BLOCK_RIBOSOME];

void FillRibosomeBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_RIBOSOME; i++)
    RibosomeMemory[i]= (RIBOSOME *) rcalloc(1,sizeof(RIBOSOME),"FillRibosomeBlock");
}

RIBOSOME *AllocRibosome()
{

  ribo_mptr_full++;

  if(ribo_mptr_full>=MEM_BLOCK_RIBOSOME) {
    if (DebugLevel > 4)
      fprintf(stderr, "AllocRibosome: allocating ribosome %d\n",
	      react_mptr_full);

    return( (RIBOSOME *) rcalloc(1,sizeof(RIBOSOME),"AllocRibosome"));
  }

  /* Take one of the ribosomes out of the block */
  RIBOSOME *ribosome = RibosomeMemory[ribo_mptr_full];
  RibosomeMemory[ribo_mptr_full] = NULL;
  return(ribosome);
}

void FreeRibosome(ribo)
RIBOSOME *ribo;
{


  if(ribo_mptr_full>=MEM_BLOCK_RIBOSOME) free(ribo);
  else
    RibosomeMemory[ribo_mptr_full]=ribo;

  ribo_mptr_full--;

}

void EmptyRibosomeBlock()
{
  int i, cnt=0;

  for(i=0; i<MEM_BLOCK_RIBOSOME; i++)
    if(RibosomeMemory[i]!=NULL) {
      free(RibosomeMemory[i]);
      cnt++;
    }

  if (DebugLevel > 3 && cnt > 0) 
    fprintf(stderr, "EmptyMRibosomeBlock: freed %d MRibosomes\n", cnt);
}

int      mrnap_mptr_full= -1;
MOVERNAP *MRNAPMemory[MEM_BLOCK_MOVERNAP];

void FillMRNAPBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_MOVERNAP; i++)
    MRNAPMemory[i]= (MOVERNAP *) rcalloc(1,sizeof(MOVERNAP),"FillMRNAPBlock");
}
 
MOVERNAP *AllocMRNAP()
{

  mrnap_mptr_full++;

  if(mrnap_mptr_full>=MEM_BLOCK_MOVERNAP) {
    if (DebugLevel > 4)
      fprintf(stderr, "AllocMRNAP: allocating MRNAP %d\n",
	      react_mptr_full);

    return( (MOVERNAP *) rcalloc(1,sizeof(MOVERNAP),"AllocMRNAP"));
  }

  /* Take one of the MRNAPs off of the block */
  MOVERNAP *movernap = MRNAPMemory[mrnap_mptr_full];
  MRNAPMemory[mrnap_mptr_full] = NULL;

  return(movernap);
}

void FreeMRNAP(mrnap)
MOVERNAP *mrnap;
{

  if(mrnap_mptr_full>=MEM_BLOCK_MOVERNAP) free(mrnap);
  else
    MRNAPMemory[mrnap_mptr_full]=mrnap;

  mrnap_mptr_full--;

}

void EmptyMRNAPBlock()
{
  int i, cnt = 0;

  for(i=0; i<MEM_BLOCK_MOVERNAP; i++)
    if(MRNAPMemory[i]!=NULL) {
      free(MRNAPMemory[i]);
      ++cnt;
    }

  if (DebugLevel > 3 && cnt > 0) 
    fprintf(stderr, "EmptyMNAPBlock: freed %d MRNAPs\n", cnt);
}

int      mribo_mptr_full= -1;
MOVERIBO *MRibosomeMemory[MEM_BLOCK_MOVERIBO];

void FillMRibosomeBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_MOVERIBO; i++)
    MRibosomeMemory[i]= (MOVERIBO *) rcalloc(1,sizeof(MOVERIBO),"FillMRibosomeBlock");
}

MOVERIBO *AllocMRibosome()
{
  mribo_mptr_full++;

  if(mribo_mptr_full>=MEM_BLOCK_MOVERIBO) {
    if (DebugLevel > 4)
      fprintf(stderr, "AllocMRibosome: allocating MRibosome %d\n",
	      react_mptr_full);

    return( (MOVERIBO *) rcalloc(1,sizeof(MOVERIBO),"AllocMRibosome"));
  }

  /* Take one of the MOVERIBOs off of the block */
  MOVERIBO *moveribo = MRibosomeMemory[mribo_mptr_full];
  MRibosomeMemory[mribo_mptr_full] = NULL;
  return(moveribo);
}

void FreeMRibosome(mribo)
MOVERIBO *mribo;
{


  if(mribo_mptr_full>=MEM_BLOCK_MOVERIBO) free(mribo);
  else
    MRibosomeMemory[mribo_mptr_full]=mribo;

  mribo_mptr_full--; 
}

void EmptyMRibosomeBlock()
{
  int i, cnt = 0;

  for(i=0; i<MEM_BLOCK_MOVERIBO; i++) {
    if(MRibosomeMemory[i]!=NULL) {
      free(MRibosomeMemory[i]);
      ++cnt;
    }
  }

  if (DebugLevel > 3 && cnt > 0) 
    fprintf(stderr, "EmptyMRibosomeBlock: freed %d MRibosomes\n", cnt);
}
