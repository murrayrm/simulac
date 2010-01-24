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
    ReactionMemory[i]= (REACTION *) rcalloc(1,sizeof(REACTION),"FillReactionBlock");
}

REACTION *AllocReaction()
{


  react_mptr_full++;

  if(react_mptr_full>=MEM_BLOCK_REACTION) 
    return( (REACTION *) rcalloc(1,sizeof(REACTION),"AllocReaction"));

  return((REACTION *) ReactionMemory[react_mptr_full]);
}

void FreeReaction(react)
REACTION *react;
{

  if(react_mptr_full>=MEM_BLOCK_REACTION) free(react);
  else
    ReactionMemory[react_mptr_full]=react;
  
    react_mptr_full--;
}

void EmptyReactionBlock()
{
  int i;

  for(i=0; i<MEM_BLOCK_REACTION; i++)
    if(ReactionMemory[i]!=NULL)
      free(ReactionMemory[i]);
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

  if(rnap_mptr_full>=MEM_BLOCK_RNAP) 
    return( (RNAP *) rcalloc(1,sizeof(RNAP),"AllocRNAP"));

  return((RNAP *) RNAPMemory[rnap_mptr_full]);
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
  int i;

  for(i=0; i<MEM_BLOCK_RNAP; i++)
    if(RNAPMemory[i]!=NULL)
      free(RNAPMemory[i]);
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

  if(ribo_mptr_full>=MEM_BLOCK_RIBOSOME) 
    return( (RIBOSOME *) rcalloc(1,sizeof(RIBOSOME),"AllocRibosome"));

  return((RIBOSOME *) RibosomeMemory[ribo_mptr_full]);
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
  int i;

  for(i=0; i<MEM_BLOCK_RIBOSOME; i++)
    if(RibosomeMemory[i]!=NULL)
      free(RibosomeMemory[i]);
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

  if(mrnap_mptr_full>=MEM_BLOCK_MOVERNAP) 
    return( (MOVERNAP *) rcalloc(1,sizeof(MOVERNAP),"AllocMRNAP"));

  return((MOVERNAP *) MRNAPMemory[mrnap_mptr_full]);
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
  int i;

  for(i=0; i<MEM_BLOCK_MOVERNAP; i++)
    if(MRNAPMemory[i]!=NULL)
      free(MRNAPMemory[i]);
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
  if(mribo_mptr_full>=MEM_BLOCK_MOVERIBO) 
    return( (MOVERIBO *) rcalloc(1,sizeof(MOVERIBO),"AllocMRibosome"));

  return((MOVERIBO *) MRibosomeMemory[mribo_mptr_full]);
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
  int i;

  for(i=0; i<MEM_BLOCK_MOVERIBO; i++)
    if(MRibosomeMemory[i]!=NULL)
      free(MRibosomeMemory[i]);
}

