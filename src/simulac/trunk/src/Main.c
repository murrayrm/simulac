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
long SEED;

/* Configuration parameters */
char *SystemFile, *ConfigPath;
double PrintTime,WriteTime;
char progid[80];
FILE *ofp = NULL, *logfp = NULL;

/* Function declarations */
void generateSetupScript(char *, char *, char *, int);

int main(argc, argv)
int argc;
char **argv;
{
  int i,j;
  int rcnt;
  double tau;
  REACTION *reaction=NULL;
  
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

  /* Check to see if a configuration file was given */
  if (args_info.config_file_given) {
    struct cmdline_parser_params *params = cmdline_parser_params_init();
    cmdline_parser_config_file(args_info.config_file_arg, &args_info, params);
  }

  /* Make sure we got enough input arguments */
  if (args_info.system_file_given) {
    /* Assume that all data about simulation is on command line */
    ConfigPath = args_info.config_path_arg;
    SystemFile = args_info.system_file_arg;
    MaximumTime = args_info.maximum_time_arg;
    PrintTime = args_info.print_time_arg;
    SEED = (args_info.seed_given) ? args_info.seed_arg : time(NULL);

  } else if (args_info.inputs_num < 3 || args_info.inputs_num > 4) {
    /* Check to make sure that we have 3 or 4 input arguments */
    fprintf(stderr, "%s\n", gengetopt_args_info_usage);
    exit(1);
  } else {
    /* Process the arguments in the original manner */
    SystemFile = args_info.inputs[0];
    MaximumTime =  atof(args_info.inputs[1]);
    PrintTime   =  atof(args_info.inputs[2]);
    SEED = (args_info.inputs_num > 3) ? atol(args_info.inputs[3]) : time(NULL);
  }

  /* Perform processing of cmdline arguments that affect init/parsing */
  extern int global_MOI;
  if (args_info.moi_given) global_MOI = args_info.moi_arg;
  if (args_info.debug_given) DebugLevel = args_info.debug_arg;

  if (args_info.param_given) {
    /* Process the command line parameters and store them for later use */
    if (param_init(args_info.param_given, args_info.param_arg) < 0) 
      /* If we got a fatal error, stop all processing */
      exit(1);
  }

  /* Set the output file handle */
  ofp = stdout;
  if (args_info.output_file_given) {
    /* Open up a file for storing the simulation output */
    if ((ofp = fopen(args_info.output_file_arg, "w")) == NULL) {
      perror(args_info.output_file_arg);
      exit(1);
    }
  }

  /* Set the log file handle */
  logfp = stderr;
  if (args_info.log_file_given) {
    /* Open up a file for storing the simulation output */
    if ((logfp = fopen(args_info.log_file_arg, "w")) == NULL) {
      perror(args_info.log_file_arg);
      exit(1);
    }
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
  /* Parse the system description file */
  ParseOutline(SystemFile);

  /* 
   * Process optional arguments 
   */

  /* Set cell volume */
  if (args_info.volume_given) {
    if (DebugLevel > 1)
      fprintf(logfp, "Resetting initial cell volume from %g", EColi->VI);
    EColi->V = (EColi->VI *= args_info.volume_arg);
    if (DebugLevel > 1)
      fprintf(logfp, " to %g\n", EColi->VI);
  }

  /* Set growth rate */
  if (args_info.growth_given) {
    if (DebugLevel > 1)
      fprintf(logfp, "Resetting cell growth rate from %g", EColi->GrowthRate);
    EColi->GrowthRate *= args_info.growth_arg;
    if (DebugLevel > 1)
      fprintf(logfp, " to %g\n", EColi->GrowthRate);
  }

  /* Determine if we should allow cell division */
  if (args_info.single_given) {
    if (DebugLevel > 1)
      fprintf(logfp, "Resetting cell volume to stop cell division\n");
    /* HACK: set cell division size to something huge */
    EColi->VI *= 1000;
  }

  /* Set parameter values */
  if (args_info.rate_given) {
    fprintf(logfp, "Processing %d rate parameters\n",
	    args_info.rate_given);

    /* Go through all initial condition commands */
    int i;
    for (i = 0; i < args_info.rate_given; ++i) {
      int index;
      float val;

      /* Parse the intial condition */
      fprintf(logfp, "  processing %s: ", args_info.rate_arg[i]);
      sscanf(args_info.rate_arg[i], "k%d=%g", &index, &val);
      fprintf(logfp, "reaction = %d, value = %g\n", index, val);

      /* Make sure the reaction number is OK */
      if (index < NMassAction) {
	ReactionProbability[index] = val;
      } else {
	fprintf(logfp, "  reaction number %d out of range (%d)\n",
		index, NMassAction); 
      }
    }
  }

  /* Set initial conditions */
  if (args_info.init_given) {
    fprintf(logfp, "Processing %d initial condition arguments\n",
	    args_info.init_given);

    /* Go through all initial condition commands */
    int i;
    for (i = 0; i < args_info.init_given; ++i) {
      char *name = (char *) calloc(1, 256);
      int val;
      assert(name != NULL);

      /* Parse the intial condition */
      fprintf(logfp, "  processing %s: ", args_info.init_arg[i]);
      sscanf(args_info.init_arg[i], "%[^=]=%d", name, &val);
      fprintf(logfp, "species = %s, value = %d\n", name, val);

      /* Try to find this argument in the list of species */
      int j;
      for (j = 0; j < NSpecies; ++j)
	if (strcmp(name, SpeciesName[j]) == 0) break;

      /* Release storage */
      free(name);

      if (j < NSpecies) {
	/* We found a match */
	fprintf(logfp, "  resetting initial concentration for species %s",
		SpeciesName[j]);
	fprintf(logfp, " from %d to %d\n", Concentration[j], val);
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
  if (DebugLevel > 2)
    fprintf(logfp, "SEED = %ld\n", SEED);
  srand48(SEED);

  DEBUG(20){
    fprintf(logfp,"@@@ NOperators  = %d\n", NOperators);
    fprintf(logfp,"@@@ NSequences  = %d\n", NSequences);
    fprintf(logfp,"@@@ NMassAction = %d\n", NMassAction);
    fprintf(logfp,"@@@ NSpecies    = %d\n", NSpecies);
  }

  /* 
   * Print out headers and setup files
   *
   */
#ifdef RMM_MODS
  FILE *matlab_fp = NULL;

  /* Check to see if we should generate setup scripts */
  if (args_info.matlab_setup_given) {
    if (DebugLevel > 1) 
      fprintf(logfp, "Generating MATLAB setup script '%s'\n",
	      args_info.matlab_setup_arg);
    generateSetupScript(args_info.matlab_setup_arg, "%", "sl_", 1);
  }

  if (args_info.python_setup_given) {
    if (DebugLevel > 1) 
      fprintf(logfp, "Generating python setup script '%s'\n",
	      args_info.python_setup_arg);
    generateSetupScript(args_info.python_setup_arg, "#", "", 0);
  }

  if (args_info.header_flag) {
# endif
    fprintf(ofp,"%% Time\tNR\tRPQ\t");
    for(i=0;i<NSpecies;i++)
      fprintf(ofp,"%6s\t",SpeciesName[i]);
    fprintf(ofp,"%6s\t","Volume");
    for(i=0;i<NOperators; i++)
      fprintf(ofp,"%6s\t",&Operator[i].Name[8]);
#ifdef RMM_MODS
    /* RNAP counts */
    if (args_info.pops_given) {
      for (i = 0; i < NPromotors; ++i)
	fprintf(ofp, "%6s-RNAP\t", Promotor[i]->Name);
    }
    fprintf(ofp,"\n");
  }
# else
  fprintf(ofp,"\n");
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
    fprintf(logfp,"NR= %d\t",NReactions);    
    fprintf(logfp,"react= %d\trnap= %d\tribo= %d\tmrnap= %d\tmribo= %d\n",
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

  /* Print some information for the user */
  if (DebugLevel) {
    fprintf(stderr, "%g: \tCNT = %d, RPQ = %e\n", t, cnt, rpq);
  }

  fprintf(ofp,"%e\t%d\t%e\t",t,cnt,rpq);
  for(i=0; i<NSpecies; i++)
    fprintf(ofp,"%7d\t",Concentration[i]);

#ifdef RMM_MODS
  /* Use the reference cell size for normalization */
  fprintf(ofp,"%e\t",EColi->V/EColi->V0);
#else
  fprintf(ofp,"%e\t",EColi->V/EColi->VI);
#endif

  for(i=0; i<NOperators; i++)
    fprintf(ofp,"%7d\t",Operator[i].CurrentState);

#ifdef RMM_MODS
  if (args_info.pops_given) {
    /* Print out the number of RNApolymerases for each promoter */
    for (i = 0; i < NPromotors; ++i)
      fprintf(ofp, "%7d\t", Promotor[i]->RNAPCount);
  }
#endif

  fprintf(ofp,"\n");
  fflush(ofp);
}

void generateSetupScript(char *filename, char *comment, char *prefix, 
			 int offset)
{
  FILE *setup_fp;
  int i;

  /* Open up the file for writing (overwrite mode) */
  if ((setup_fp = fopen(filename, "w")) == NULL) { 
    perror(filename); 
    exit(-1); 
  }

  /* Print out some header information */
  time_t curtime = time(NULL);
  fprintf(setup_fp, "%s Simulac setup file\n", comment);
  fprintf(setup_fp, "%s Run generated: %s", comment, ctime(&curtime));
   
  /* Print out parameters governing the simulation */
  fprintf(setup_fp, "\n%s Simulation parameters\n", comment);
  fprintf(setup_fp, "%sconfig_file = '%s';\n", prefix, SystemFile);
  fprintf(setup_fp, "%smaximum_time = %g;\n", prefix, MaximumTime);
  fprintf(setup_fp, "%sprint_time = %g;\n", prefix, PrintTime);
  fprintf(setup_fp, "%sseed = %ld;\n", prefix, SEED);

  /* Print out information about the number of objects of each type */
  fprintf(setup_fp, "\n%s System size\n", comment);
  fprintf(setup_fp, "%sn_species = %d;\n", prefix, NSpecies);
  fprintf(setup_fp, "%sn_operators = %d;\n", prefix, NOperators);
  fprintf(setup_fp, "%sn_promoters = %d;\n", prefix, NPromotors);

  /* Print the column indices for the output file */
  fprintf(setup_fp, "\n%s Data indices\n", comment);
  int col = offset;
  fprintf(setup_fp, "%stime_index = %d;\n", prefix, col++);
  fprintf(setup_fp, "%sNR_index = %d;\n", prefix, col++);
  fprintf(setup_fp, "%sRPQ_index = %d;\n", prefix, col++);
  for (i = 0; i < NSpecies; ++i) 
    fprintf(setup_fp, "%sspecies_%s_index = %d;\n", prefix, 
	    SpeciesName[i], col++);
  fprintf(setup_fp, "%svolume_index = %d;\n", prefix, col++);
  for (i=0 ; i < NOperators; ++i)
    fprintf(setup_fp, "%soperator_%s_index = %d;\n", prefix, 
	    Operator[i].Name, col++);
  if (args_info.pops_given) {
    for (i = 0; i < NPromotors; ++i)
      fprintf(setup_fp, "%spromoter_%s_index = %d;\n", prefix,
	      Promotor[i]->Name, col++);
  }
    
  /* Close up the file */
  fclose(setup_fp);
}
