/**********************
 *
 * A little memory manager for 
 * Simulac.
 *
 **********************/
#define MEMORY

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

#define MEM_BLOCK_REACTION  500
#define MEM_BLOCK_RNAP      150      
#define MEM_BLOCK_RIBOSOME  150
#define MEM_BLOCK_MOVERNAP  150
#define MEM_BLOCK_MOVERIBO  150

extern int      react_mptr_full;
extern REACTION *ReactionMemory[MEM_BLOCK_REACTION];

extern void FillReactionBlock();
extern REACTION *AllocReaction();
extern void FreeReaction(REACTION *);
extern void EmptyReactionBlock();

extern int      rnap_mptr_full;
extern RNAP *RNAPMemory[MEM_BLOCK_RNAP];
extern void FillRNAPBlock();
extern RNAP *AllocRNAP();
extern void FreeRNAP(RNAP *);
extern void EmptyRNAPBlock();

extern int      ribo_mptr_full;
extern RIBOSOME *RibosomeMemory[MEM_BLOCK_RIBOSOME];

extern void FillRibosomeBlock();
extern RIBOSOME *AllocRibosome();
extern void FreeRibosome(RIBOSOME *);
extern void EmptyRibosomeBlock();

extern int      mrnap_mptr_full;
extern MOVERNAP *MRNAPMemory[MEM_BLOCK_MOVERNAP];

extern void FillMRNAPBlock();
extern MOVERNAP *AllocMRNAP();
void FreeMRNAP(MOVERNAP *);
void EmptyMRNAPBlock();

extern int      mribo_mptr_full;
extern MOVERIBO *MRibosomeMemory[MEM_BLOCK_MOVERIBO];
extern void FillMRibosomeBlock();
extern MOVERIBO *AllocMRibosome();
extern void FreeMRibosome(MOVERIBO *);
extern void EmptyMRibosomeBlock();






