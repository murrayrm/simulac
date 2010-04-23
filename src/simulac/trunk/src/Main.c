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

/* Modifications by R. M. Murray, 1 Oct 09 */
#ifdef RMM_MODS
#include <time.h>
#include "param.h"
#include "cmdline.h"
struct gengetopt_args_info args_info;
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

#ifdef RMM_MODS
int        NPromotors=0;
PROMOTOR   **Promotor = NULL;
#endif

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

#ifdef RMM_MODS
  /* Command line argument parsing */
  if (cmdline_parser(argc, argv, &args_info) != 0) exit(1);

  /* See if the user wants some help with options */
  if (args_info.help_given) {
    cmdline_parser_print_help();
    exit(1);
  }

  /* Make sure we got enough input arguments */
  if (args_info.inputs_num < 3 || args_info.inputs_num > 4) {
    fprintf(stderr, "%s\n", gengetopt_args_info_usage);
    exit(1);
  }

  /* Perform processing of cmdline arguments that affect init/parsing */
  extern int global_MOI;
  if (args_info.moi_given) global_MOI = args_info.moi_arg;

  if (args_info.param_given) {
    /* Process the command line parameters and store them for later use */
    if (param_init(args_info.param_given, args_info.param_arg) < 0) 
      /* If we got a fatal error, stop all processing */
      exit(1);
  }

#else
  if (argc!=5) {
    fprintf(stderr, "usage: %s outline_file Maximum_Time Print_Time SEED\n",
	    argv[0]);
    exit(-1);
  }
#endif

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

#ifdef RMM_MODS
  /* Use the command line inputs to get basic parameters */
  ParseOutline(args_info.inputs[0]);
  MaximumTime =  atof(args_info.inputs[1]);
  PrintTime   =  atof(args_info.inputs[2]);
  SEED = (args_info.inputs_num > 3) ? atol(args_info.inputs[3]) : time(NULL);

  /* 
   * Process optional arguments 
   */

  /* Set cell volume */
  if (args_info.volume_given) {
    fprintf(stderr, "Resetting initial cell volume from %g", EColi->VI);
    EColi->V = (EColi->VI *= args_info.volume_arg);
    fprintf(stderr, " to %g\n", EColi->VI);
  }

  /* Set growth rate */
  if (args_info.growth_given) {
    fprintf(stderr, "Resetting cell growth rate from %g", EColi->GrowthRate);
    EColi->GrowthRate *= args_info.growth_arg;
    fprintf(stderr, " to %g\n", EColi->GrowthRate);
  }

  /* Determine if we should allow cell division */
  if (args_info.single_given) {
    fprintf(stderr, "Resetting cell volume to stop cell division\n");
    /* HACK: set cell division size to something huge */
    EColi->VI *= 1000;
  }

  /* Set parameter values */
  if (args_info.rate_given) {
    fprintf(stderr, "Processing %d rate parameters\n",
	    args_info.rate_given);

    /* Go through all initial condition commands */
    int i;
    for (i = 0; i < args_info.rate_given; ++i) {
      int index;
      float val;

      /* Parse the intial condition */
      fprintf(stderr, "  processing %s: ", args_info.rate_arg[i]);
      sscanf(args_info.rate_arg[i], "k%d=%g", &index, &val);
      fprintf(stderr, "reaction = %d, value = %g\n", index, val);

      /* Make sure the reaction number is OK */
      if (index < NMassAction) {
	ReactionProbability[index] = val;
      } else {
	fprintf(stderr, "  reaction number %d out of range (%d)\n",
		index, NMassAction); 
      }
    }
  }

  /* Set initial conditions */
  if (args_info.init_given) {
    fprintf(stderr, "Processing %d initial condition arguments\n",
	    args_info.init_given);

    /* Go through all initial condition commands */
    int i;
    for (i = 0; i < args_info.init_given; ++i) {
      char *name = (char *) calloc(1, 256);
      int val;
      assert(name != NULL);

      /* Parse the intial condition */
      fprintf(stderr, "  processing %s: ", args_info.init_arg[i]);
      sscanf(args_info.init_arg[i], "%[^=]=%d", name, &val);
      fprintf(stderr, "species = %s, value = %d\n", name, val);

      /* Try to find this argument in the list of species */
      int j;
      for (j = 0; j < NSpecies; ++j)
	if (strcmp(name, SpeciesName[j]) == 0) break;

      /* Release storage */
      free(name);

      if (j < NSpecies) {
	/* We found a match */
	fprintf(stderr, "  resetting initial concentration for species %s",
		SpeciesName[j]);
	fprintf(stderr, " from %d to %d\n", Concentration[j], val);
	Concentration[j] = val;
      } else {
	fprintf(stderr, "  couldn't find species %s\n", name);
      }
    }
  }
#else  
  /****************
   *
   * Read in simulation model
   *
   *****************/

  ParseOutline(argv[1]);
  MaximumTime =  atof(argv[2]);
  PrintTime   =  atof(argv[3]);
  SEED        =  atol(argv[4]);
#endif
  fprintf(stderr, "SEED = %ld\n", SEED);
  srand48(SEED);

  DEBUG(20){
    fprintf(stderr,"@@@ NOperators  = %d\n", NOperators);
    fprintf(stderr,"@@@ NSequences  = %d\n", NSequences);
    fprintf(stderr,"@@@ NMassAction = %d\n", NMassAction);
    fprintf(stderr,"@@@ NSpecies    = %d\n", NSpecies);
  }

  /* 
   * Print out headers and setup files
   *
   */
#ifdef RMM_MODS
  FILE *matlab_fp = NULL;

  /* Check to see if we should generate setup scripts */
  if (args_info.matlab_setup_given) {
    fprintf(stderr, "Generating MATLAB setup script '%s'\n",
	    args_info.matlab_setup_arg);

    /* Open up the file for writing (overwrite mode) */
    if ((matlab_fp = fopen(args_info.matlab_setup_arg, "w")) == NULL) { 
      perror(args_info.matlab_setup_arg); 
      exit(-1); 
    }

    /* Print out some header information */
    time_t curtime = time(NULL);
    fprintf(matlab_fp, "%% Simulac MATLAB setup file\n");
    fprintf(matlab_fp, "%% Run generated: %s", ctime(&curtime));
   
    /* Print out parameters governing the simulation */
    fprintf(matlab_fp, "\n%% Simulation parameters\n" );
    fprintf(matlab_fp, "sl_config_file = '%s';\n", args_info.inputs[0]);
    fprintf(matlab_fp, "sl_maxtime = %s;\n", args_info.inputs[1]);
    fprintf(matlab_fp, "sl_stepsize = %s;\n", args_info.inputs[2]);
    fprintf(matlab_fp, "sl_seed = %ld;\n", SEED);

    /* Print out information about the number of objects of each type */
    fprintf(matlab_fp, "\n%% System size\n" );
    fprintf(matlab_fp, "sl_n_species = %d;\n", NSpecies);
    fprintf(matlab_fp, "sl_n_operators = %d;\n", NOperators);
    fprintf(matlab_fp, "sl_n_promoters = %d;\n", NPromotors);

    /* Print the column indices for the output file */
    fprintf(matlab_fp, "\n%% Data indices\n" );
    int col = 1;
    fprintf(matlab_fp, "sl_time_index = %d;\n", col++);
    fprintf(matlab_fp, "sl_NR_index = %d;\n", col++);
    fprintf(matlab_fp, "sl_RPQ_index = %d;\n", col++);
    for (i = 0; i < NSpecies; ++i) 
      fprintf(matlab_fp, "sl_species_%s_index = %d;\n", SpeciesName[i], col++);
    fprintf(matlab_fp, "sl_volume_index = %d;\n", col++);
    for (i=0 ; i < NOperators; ++i)
      fprintf(matlab_fp, "sl_operator_%s_index = %d;\n",
	      Operator[i].Name, col++);
    if (args_info.pops_given) {
      for (i = 0; i < NPromotors; ++i)
	fprintf(matlab_fp, "sl_promoter_%s_index = %d;\n",
		Promotor[i]->Name, col++);
    }
    
    /* Close up the file */
    fclose(matlab_fp);
  }

  if (args_info.header_flag) {
# endif
    fprintf(stdout,"%% Time\tNR\tRPQ\t");
    for(i=0;i<NSpecies;i++)
      fprintf(stdout,"%6s\t",SpeciesName[i]);
    fprintf(stdout,"%6s\t","Volume");
    for(i=0;i<NOperators; i++)
      fprintf(stdout,"%6s\t",&Operator[i].Name[8]);
#ifdef RMM_MODS
    /* RNAP counts */
    if (args_info.pops_given) {
      for (i = 0; i < NPromotors; ++i)
	fprintf(stdout, "%6s-RNAP\t", Promotor[i]->Name);
    }
    fprintf(stdout,"\n");
  }
# else
  fprintf(stdout,"\n");
#endif

  /********************
   *
   * Main Loop 
   *
   ********************/

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

#ifdef RMM_MODS
  /* Use the reference cell size for normalization */
  fprintf(stdout,"%e\t",EColi->V/EColi->V0);
#else
  fprintf(stdout,"%e\t",EColi->V/EColi->VI);
#endif

  for(i=0; i<NOperators; i++)
    fprintf(stdout,"%7d\t",Operator[i].CurrentState);

#ifdef RMM_MODS
  if (args_info.pops_given) {
    /* Print out the number of RNApolymerases for each promoter */
    for (i = 0; i < NPromotors; ++i)
      fprintf(stdout, "%7d\t", Promotor[i]->RNAPCount);
  }
#endif

  fprintf(stdout,"\n");
  fflush(stdout);
}

