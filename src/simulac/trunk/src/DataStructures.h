/*******************************
 *
 * Fundamental Data Structures
 *
 *          and 
 *
 * Global Chemical Variables 
 *
 *
 * 2/13/96- Added Volume Effects
 * This means that all probability constants in the input
 * files are multiplied by the appropriate volume factors
 * This rates are then divided by 
 *******************************/


#define DataStructures

#ifndef _H_STDLIB
  #include <stdlib.h>
#endif

/***********************
 *
 * Simulation Object Types
 *
 ***********************/

typedef struct rnap       RNAP;
typedef struct dna        DNA;
typedef struct segment    SEGMENT;
typedef struct transcript mRNA;
typedef struct promotor   PROMOTOR;
typedef struct ackersshea SHEADATA;
typedef struct reaction   REACTION;
typedef struct ribosome   RIBOSOME;
typedef struct cell       CELL;

#define LEFT  0
#define RIGHT 1

struct rnap {
  int        NBound;
  int        *SpeciesIndex;

  short      Direction;
  int        CurrentPosition;
  mRNA      *Transcript;
  RNAP      *LastRNAP;
  RNAP      *NextRNAP;
};

#define DNA_Type_Promotor       0
#define DNA_Type_Coding         1
#define DNA_Type_NonCoding      2
#define DNA_Type_Terminator     3
#define DNA_Type_AntiTerminator 4

struct dna {

  DNA         *LeftSegment;
  DNA         *RightSegment;

  CELL       *Bug;
  char       *Name;
  int        Length;
  short      Direction;

  short      Type;
  void       *DNAStruct;

};

struct promotor {
  int          Data;
  double      *IsoRate;   /* Isomerization rates for formation of OC for each config */
  short        TranscriptionDirection;
  RNAP        *RNAPQueue;
#ifdef RMM_MODS
  char        *Name;	/* Name of the promotor (for output file) */
  int RNAPCount;	/* Count the number of RNAP polymerases that go by  */
#endif
};

struct segment {
  RNAP*      RNAPQueue;
  void       *SegmentData;
  void      (*SegmentFunc)();
};

struct ribosome {
  int        NBound;
  int        *SpeciesIndex;


  int        CurrentPosition;
  mRNA      *Transcript;
  RIBOSOME  *LastRibosome;
  RIBOSOME  *NextRibosome;
};

#define mRNA_Type_Sense         0
#define mRNA_Type_AntiSense     1

#define mRNA_RBS_Intact         0
#define mRNA_RBS_Chewed         1
#define mRNA_RBS_AntiSensed     2

struct transcript {
  DNA        *Gene;
  RNAP       *Rnap;
  short       RBSState;
  int         CurrentLength;
  int         Type;
  RIBOSOME   *RiboQueue;
  mRNA       *LastTranscript;
  mRNA       *NextTranscript;
};

struct ackersshea {
  char        *Name;
  int         NSites;
  int         NConfigs;
  int       **CList;
  double     *DeltaG;
  int         CurrentState;  /* Index into configuration matrix */
};

#define Reaction_Type_Kinetic            0
#define Reaction_Type_TransInit          1
#define Reaction_Type_MoveRNAP           2
#define Reaction_Type_DNAAction          3
#define Reaction_Type_NextSegment        4
#define Reaction_Type_EatmRNA            5
#define Reaction_Type_MoveRibosome       6
#define Reaction_Type_ProduceProtein     7
#define Reaction_Type_RNAP_RNAP          8
#define Reaction_Type_Recombine          9
#define Reaction_Type_BindRibosome      10
#define Reaction_Type_ProduceNewProtein 11
#define Reaction_Type_ChangeCellVolume  12

struct reaction {
  int       Type;
  
  void      *ReactionData;
  void     (*ReactionFunc)();
  double     Probability;
  REACTION  *LastReaction;
  REACTION  *NextReaction;
};

struct cell {
  int Type;

  int id;
  int x,y;
  double GrowthRate;  /* In units of 10^-18 liters/sec */
  double VI;      /* Initial Volume */
  double V0;      /* Standard Volume by which all probability coefficient 
                     are normalized. Remember: 
                     c= prob coeff, k= rate coeff
		     k proportional to c*V so since we are defining
		     inputs with respect to a standard volume
		     all reaction probabilities get multiplied
		     by a factor V0/Volume to the appropriate order 
		     of the reaction:
		             zeroth and first order - no factor
			     second order - (V0/V)
			     third  order - (V0/V)^2
			     etc.
		   */
  double V;
};

/***************************
 *
 * Object Parameter Data Types
 *
 ***************************/


extern double Rate_Of_Polymerase_Motion,Rate_Of_RNAP_Collision_Escape,Rate_Of_RNAP_Collision_Failure;
extern double Rate_Of_Ribosome_Binding,Rate_Of_Ribosome_Motion; 

typedef struct movernap MOVERNAP;
typedef struct moveribo MOVERIBO;

typedef struct antitermdata ANTITERMDATA;
typedef struct termdata     TERMDATA;
typedef struct reactdata    REACTDATA;
typedef struct codingdata   CODINGDATA;

struct reactdata {

  int Mu;

};

struct movernap {
  DNA *dna;
  RNAP *rnap1;
  RNAP *rnap2;
};

struct moveribo {

  mRNA     *trans;
  RIBOSOME *ribosome;
};

struct codingdata {

  int SpeciesIndex;
  double mRNADegradationRate;
# ifdef RMM_MODS
  double RibosomeBindingRate;
# endif
};

struct antitermdata {

  int        SpeciesIndex;
  double     UnBoundRNAPMotion;
  double     BindingRate;
  double     BoundRNAPMotion;
  double     UnBindingRate;

};

struct termdata {

  int       SpeciesIndex;
  double    BaseFallOffRate;
  double    BaseRNAPMotion;
  double    AntiTerminatedFallOffRate;
  double    AntiTerminatedRNAPMotion;

};

/****** Globals *********/

extern int        NMechanisms;
extern char       **Mechanism;

extern int        NSpecies;
extern char     **SpeciesName;
extern int       *Concentration;

extern int        NOperators;
extern SHEADATA  *Operator;

extern int        NSequences;
extern DNA       *Sequence;

#ifdef RMM_MODS
/* Keep track of the promotors so that we can report on transcription rates */
extern int        NPromotors;
extern PROMOTOR   **Promotor;
#endif

extern int        NTranscripts;
extern mRNA      *Transcript;

extern int        NMassAction;
extern int      **StoMat1;
extern int      **StoMat2;
extern double    *ReactionProbability;

extern int       NReactions;
extern REACTION  *Reaction;
extern double     TotalProbability;
extern double    *Probabilities;


extern CELL      *EColi;

extern double Time;
extern double MaximumTime;
extern double PrintTime,WriteTime;

extern char progid[80];



