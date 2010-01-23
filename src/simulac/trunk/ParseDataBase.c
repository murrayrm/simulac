/**************
 *
 * These routines parse a mechanistic outline file
 *
 * An outline file contains a list of mech files each of
 * which contains either a kinetic reaction system, a genetic
 * reaction system. Two other special types of mech files are
 * RNAP and RIBSOME files which describe the dynamics of 
 * these particular complexes. These last files are necessary only
 * if a genetic reaction system is included in the outline file
 *
 **************/

/****************************/
/******* Includes ***********/
/****************************/

#ifndef _H_STDIO
 #include <stdio.h>
#endif

#ifndef _H_STLIB
 #include <stdlib.h>
#endif

#ifndef _H_STRINGS
 #include <strings.h>
#endif

#ifndef __H_CTYPE
 #include <ctype.h>
#endif

#ifndef DataStructures
 #include "DataStructures.h"
#endif

#ifndef UTILS
 #include "Util.h"
#endif

#define Success  0
#define Failure  1

#ifdef RMM_MODS
/* Create a variable for specifying the global MOI */
int global_MOI = 0;
#endif

/************************************/
/******* Simple Utilities ***********/
/************************************/

FILE *OpenFile(file,status)
char *file,*status;
{
  FILE *fp;

  if((fp=fopen(file,status))==NULL){
    fprintf(stderr,"%s: Unable to open file %s with status %s.\n",progid,file,status);
    exit(-1);
  }

  return(fp);
}

void GrowStomatReactions(n)
int n;
{
  if(StoMat1==NULL) StoMat1= (int **) rcalloc(n,sizeof(int *),"GrowStomatR.1");
  else              StoMat1= (int **) rrealloc(StoMat1,n,sizeof(int *),"GrowStomatR.1");

  StoMat1[n-1]= (int *) rcalloc(NSpecies,sizeof(int),"GrowStomatR.1.1");

  if(StoMat2==NULL) StoMat2= (int **) rcalloc(n,sizeof(int *),"GrowStomatR.2");
  else              StoMat2= (int **) rrealloc(StoMat2,n,sizeof(int *),"GrowStomatR.2");

  StoMat2[n-1]= (int *) rcalloc(NSpecies,sizeof(int),"GrowStomatR.2.1");
}

void GrowStomatSpecies(n)
int n;
{
  int i;

  fprintf(stderr,"NMassAction= %d NSpecies = %d\n",NMassAction, n);
  for(i=0; i<NMassAction; i++){

    if(StoMat1[i]==NULL) StoMat1[i]= (int *) rcalloc(n,sizeof(int),"GrowStomatS.1");
    else                 StoMat1[i]= (int *) rrealloc(StoMat1[i],n,sizeof(int),"GrowStomatS.1");

    if(StoMat2[i]==NULL) StoMat2[i]= (int *) rcalloc(n,sizeof(int),"GrowStomatS.2");
    else                 StoMat2[i]= (int *) rrealloc(StoMat2[i],n,sizeof(int),"GrowStomatS.2");
  }

}

int FindSpecies(name)
char *name;
{
  int i;

  if(strcmp(name,"___")==0) return(-1);

  for(i=0; i<NSpecies; i++)
    if(strcmp(name,SpeciesName[i])==0) break;

return(i);
}

void AddSpecies(name)
char *name;
{
  void GrowStomatSpecies();

  NSpecies++;
  GrowStomatSpecies(NSpecies);

  if(SpeciesName==NULL) SpeciesName= (char **) rcalloc(NSpecies, sizeof(char *),"AddSpecies.1");
  else                  SpeciesName= (char **) rrealloc(SpeciesName,NSpecies,sizeof(char *),"AddSpecies.1");

  SpeciesName[NSpecies-1]= (char *) rcalloc(strlen(name)+1,sizeof(char),"AddSpecies.2");
  strcpy(SpeciesName[NSpecies-1],name);

  if(Concentration==NULL) Concentration= (int *) rcalloc(NSpecies, sizeof(int),"AddSpecies.3");
  else                    Concentration= (int *) rrealloc(Concentration,NSpecies,sizeof(int),"AddSpecies.3");

  Concentration[NSpecies-1]=0;

  fprintf(stderr,"@@@ %s: Added Species: %s\n",progid,name);
}

/*** Base Unit is set to seconds ****/

char  *Units[20]= {"ms","msec","millisecond","s","sec","seconds","min","minute","hr","hour","NULL"};
double FindTimeUnit(unit)
char *unit;
{
  int i;

  i=0;
  while(strcmp(Units[i],"NULL")!=0 && strncasecmp(Units[i],unit,strlen(Units[i]))!=0) i++;

  if(strcmp(Units[i],"NULL")==0) return(-1.);

  switch(i){

  case 0:
  case 1:
  case 2:
    return(0.001);
  case 3:
  case 4:
  case 5:
    return(1.0);
  case 6:
  case 7:
    return(60.0);
  case 8:
  case 9:
    return(3600.0);
  }
return(-1.);
}
  
/************************************/
/******* Outline Parser *************/
/************************************/
  
short ParseOutline(file)
char *file;
{
  int i,len;
  short dflag=0,rnapflag=0,riboflag=0,cflag=0;
  char file1[80],buffer[200], token[80];
  FILE *fp,*fp1;

  FILE *OpenFile();
  void ReadCell();
  void ReadKinetics();
  void ReadRibosome();
  void ReadDNA();

  fp=OpenFile(file,"r");

  while(fgets(buffer,81,fp)!=NULL){
    
    i=0;
    len=strlen(buffer);

    while(i<len && isspace((int) buffer[i])) i++;
    if(i==len) continue;

    if(buffer[i]=='%') continue;

    sscanf(&buffer[i],"%s",file1);

    for(i=0; i<NMechanisms; i++)
      if(strcmp(Mechanism[i],file1)==0){
	fprintf(stderr,"%s: You cannot include the same mechanism twice (%s) in an outline file!\n",progid,file1);
	exit(-1);
      }
    
    if(Mechanism==NULL) Mechanism = (char **) rcalloc(1,sizeof(char *),"ParseOutLine.1");
    else                Mechanism = (char **) rrealloc(Mechanism,NMechanisms+1,sizeof(char *),"ParseOutLine.1");

    Mechanism[NMechanisms]= (char *) rcalloc(strlen(file1)+1,sizeof(char),"ParseOutLine.2");

    strcpy(Mechanism[NMechanisms],file1);
    NMechanisms++;

    fp1=OpenFile(file1,"r");
    while(fgets(buffer,81,fp1)!=NULL){
      
      sscanf(buffer,"%s",token);

      if(token[0]=='%' || isspace(token[0])) continue;

      if(strcasecmp(token,"Type=")!=0){
	fprintf(stderr,"%s: Incorrect token %s found in mechanism %s. Expected to find \"Type=\" token.\n",progid,token,file1);
	exit(-1);
      }

      sscanf(buffer,"%*s %s",token);

      break;
    }

    if(strcasecmp(token,"Cell")==0){
      ReadCell(fp1,file1);
      cflag++;
    } else if(strcasecmp(token,"Kinetic")==0){
      if(cflag==0)
	fprintf(stderr,"%s: At least one Cell object must be included in an outline file. All Cells must be included at the top of the outline file. (outline= %s)\n", progid,file1);
      ReadKinetics(fp1,file1);
    } else if(strcasecmp(token,"DNA")==0){
      if(cflag==0)
	fprintf(stderr,"%s: At least one Cell object must be included in an outline file. All Cells must be included at the top of the outline file. (outline= %s)\n", progid,file1);
      ReadDNA(fp1,file1);
      dflag=1;
    } else if(strcasecmp(token,"RNAP")==0){
      if(cflag==0)
	fprintf(stderr,"%s: At least one Cell object must be included in an outline file. All Cells must be included at the top of the outline file. (outline= %s)\n", progid,file1);
      if(rnapflag!=0){
	fprintf(stderr,"%s: No more than one RNAP may be included in an outline file. Extra RNAP file %s.\n",progid,file1);
	exit(-1);
      }
      ReadKinetics(fp1,file1);
      rnapflag++;
    }else if(strcasecmp(token,"Ribosome")==0){
      if(cflag==0)
	fprintf(stderr,"%s: At least one Cell object must be included in an outline file. All Cells must be included at the top of the outline file. (outline= %s)\n", progid,file1);
      if(riboflag!=0){
	fprintf(stderr,"%s: No more than one RIBOSOME may be included in an outline file. Extra RIBOSOME file %s.\n",progid,file1);
	exit(-1);
      }
      ReadKinetics(fp1,file1);
      riboflag++;
    }
  }

  if(dflag>0 && (rnapflag==0 || riboflag==0)){
    fprintf(stderr,"%s: An RNAP and Ribosome mechanism file must be included if a DNA file is present. A %s file is missing.\n",progid,(rnapflag==0 ? (riboflag==0 ? "RNAP and Ribosome":"RNAP") : "Ribsome"));
    exit(-1);
  }


  if(rnapflag>1 || riboflag>1){
    fprintf(stderr,"%s: No more than one RNAP file or Ribosome file may be included in an outline file.\n",progid);
    exit(-1);
  }

fclose(fp);
return(Success);
}

/**************************************/
/******* Specific Parsers *************/
/**************************************/

void ReadCell(fp,cell)
FILE *fp;
char *cell;
{
  double vol;
  char buffer[90];

  /**** Allocate EColi *****/

  EColi= (CELL *) calloc(1,sizeof(CELL));

  /* Read to Separator */
  while(fgets(buffer,81,fp)!=NULL){
    strcat(buffer,"\n");
    if(buffer[0]=='%') continue; /* '%' is a comment */
    if(buffer[0]=='\n') continue; /* Blank Line */
    if(buffer[0]=='-') break;
    fprintf(stderr,"%s: Syntax error before TYPE/Parameter Separator in Cell %s\n",progid,cell);
  }

  fscanf(fp,"%*s %lf",&vol);
  EColi->VI=vol;
  fscanf(fp,"%*s %lf",&vol);
  EColi->V0=vol;
  fscanf(fp,"%*s %lf",&vol);
  EColi->GrowthRate=vol;
  EColi->V=EColi->VI;
}


void ReadKinetics(fp,mech)
FILE *fp;
char *mech;
{
  int    i,j,spec,nreact,nspec;
  char   buffer[90],token[20];
  char  *tptr,*rightside,*leftside;
  int    molec;
  void   GrowStomatReactions();
  void   GrowStomatSpecies();
  int    FindSpecies();
  void   AddSpecies();
  double FindTimeUnit();

  while(fgets(buffer,81,fp)!=NULL){
    strcat(buffer,"\n");
    if(buffer[0]=='%') continue; /* '%' is a comment */
    if(buffer[0]=='\n') continue; /* Blank Line */
    if(buffer[0]=='-') break;
    fprintf(stderr,"%s: Syntax error before TYPE/Mechanism Separator in mechanism %s\n",progid,mech);
  }

  /* Read File */

  nreact=nspec=0;
  while(fgets(buffer,81,fp)!=NULL){
    if((tptr=strstr(buffer,"-->"))==NULL){
      if(buffer[0]=='-') break;
      if(buffer[0]=='%') continue;
      for(i=0; i<strlen(buffer); i++)
	if(isspace(buffer[i])==0) break;
      
      if(i!=strlen(buffer)){
	fprintf(stderr,"%s: Syntax Error before reaction %d in mechanism %s.\n",
		progid,nreact+1,mech);
	exit(-1);
      }
      
      continue;
    }
    
    tptr += 3;        /* Step past arrow */
    *(tptr-3)= '\0';   /* Put Null before arrow */
    
    rightside=buffer;
    i=0;
    while(isspace(tptr[i])) i++; /* Rid me of leading white space */
    tptr = &tptr[i];
    leftside=tptr;

    NMassAction++;
    nreact++;
    GrowStomatReactions(NMassAction);

    /* Analyze right size */

    do{

      sscanf(rightside,"%s",token);
      
      molec=1;
      if(isdigit(token[0]))
	sscanf(rightside,"%d %s",&molec,token);
      
      if(strcmp(token,"()")!=0){
	spec=FindSpecies(token);
	if(spec==NSpecies) AddSpecies(token);
	StoMat1[NMassAction-1][spec] += molec;
      }
      
      tptr=strstr(rightside,"+");
      if(tptr==NULL) break;
      tptr = &tptr[1]; /* Step past + sign */
      rightside=tptr;

    } while(rightside[0]!='\0' && rightside[0]!='\n');

    /* Analyze left side */

    do{

      sscanf(leftside,"%s",token);
      
      molec=1;
      if(isdigit(token[0]))
	sscanf(leftside,"%d %s",&molec,token);
      
      if(strcmp(token,"()")!=0){
	spec=FindSpecies(token);
	if(spec==NSpecies) AddSpecies(token);
	StoMat2[NMassAction-1][spec] += molec;
      }
      
      tptr=strstr(leftside,"+");
      if(tptr==NULL) break;
      tptr += 1; /* Step past + sign */
      leftside=tptr;
    } while(leftside[0]!='\0' && leftside[0]!='\n');
  }
  
  /* Read Reaction Probabilities */
  
  if(ReactionProbability==NULL) ReactionProbability= (double *) rcalloc(NMassAction,sizeof(double),"ReadKinetics.1");
  else                          ReactionProbability= (double *) rrealloc(ReactionProbability,NMassAction,sizeof(double),"ReadKinetics.1");

  for(i=0; i<nreact; i++){
    double mult;
    fscanf(fp,"%*s %*s %lf %s",&ReactionProbability[NMassAction-nreact+i],token);
    mult=FindTimeUnit(token);
    if(mult<0.0){
      fprintf(stderr,"%s: Illegal time scale found for reaction probability %d in mechanism %s.\n",
	      progid,i,mech);
      exit(-1);
    }
    ReactionProbability[NMassAction-nreact+i] /= mult;    
  }

  fscanf(fp,"%*s"); /* Read Spacer */

  /* Read Intial Concentrations */

  while(fscanf(fp,"%s %*s %d",token,&molec)!=EOF){
    
    j= FindSpecies(token);
    
    if(j==NSpecies){
      fprintf(stderr,"Found unknown species %s\n",token);
      exit(-1);
    }
    fprintf(stderr,"@@@ [%s]= %d\n",SpeciesName[j],molec);
    Concentration[j]=molec;
  }
  
fclose(fp);	 	 
}

void ReadDNA(fp, mech)
FILE *fp;
char *mech;
{
  int i,j;
  int rflag,nflag,moi,length,mult,firstcopy;
  DNA *sequence,*nextsequence;
  char token[80],buffer[90];
  char unit[30],type[80],direction[10];
  char **name,**parameters;
  SEGMENT  *seg;
  PROMOTOR *prom;
  void ReadPromotorData();
  void ReadTerminatorData();
  void ReadAntiTerminatorData();
  void ReadCodingData();
  void ReadNonCodingData();
  void *FindSegData();

  void SubmitTermination();
  void SubmitAntiTermination();
  void SubmitProduceTranscript();
  void SubmitSimpleJumpSegment();

  /**** Allocate memory for sequence *******/

  NSequences++;

  if(Sequence==NULL) Sequence= (DNA *) rcalloc(NSequences,sizeof(DNA),"ReadDNA.1");
  else               Sequence= (DNA *) rrealloc(Sequence,NSequences,sizeof(DNA),"ReadDNA.1");

  /****
   *
   * We've realloc'ed....we gotta reset 
   * the pointers to the first segment of 
   * the first sequence 
   ****/

  for(i=0; i<NSequences-1; i++){
    sequence= &Sequence[i];
  if(sequence->RightSegment != NULL)
    sequence->RightSegment->LeftSegment= sequence;
  }

  sequence=(DNA *) &Sequence[NSequences-1];
  
  /**** Read down to first separator ******/

  while(fgets(buffer,81,fp)!=NULL){
    strcat(buffer,"\n");
    if(buffer[0]=='%') continue; /* '%' is a comment */
    if(buffer[0]=='\n') continue; /* Blank Line */
    if(buffer[0]=='-') break;
    fprintf(stderr,"%s: Syntax error before TYPE/Mechanism Separator in mechanism %s\n",progid,mech);
  }

  /* Read DNA Segment */

  rflag=0;
  nflag=1;
  fscanf(fp,"%s",token);  
  do{
  
    if(strcmp(token,"-->")==0){
      if(rflag==0){
	fprintf(stderr,"%s: '-->' cannot be the first token in a dna sequence (%s).\n",progid,mech);
	exit(-1);
      }
      nflag=1;
    } else if(strncmp(token,"---",3)==0){
      if(rflag==0){
	fprintf(stderr,"%s: Null DNA Sequence encountered in mechanism %s.\n",progid,mech);
	fclose(fp);
	NSequences--;
	free(sequence);
	if(NSequences==0) free(Sequence);
	else              Sequence= (DNA *) rrealloc(Sequence,NSequences,sizeof(DNA),"ReadDNA.2");
	return;
      }
      if(nflag==1){
	fprintf(stderr,"%s: Last token in DNA sequence (%s) cannot be '-->'.\n", progid,mech);
	exit(-1);
      }
      nflag= -1;
      break;
    } else{ /* Next Piece of DNA */
      if(nflag!=1){
	fprintf(stderr,"%s: Expected to find '-->' token %s but found %s instead in mechanism %s.\n",progid,(rflag>0 ? "or a Separator" : ""), token,mech);
	exit(-1);
      }
      if(rflag!=0){
	nextsequence= (DNA *) rcalloc(1,sizeof(DNA),"ReadDNA.3");
	sequence->RightSegment=nextsequence;
	nextsequence->LeftSegment=sequence;
	nextsequence->RightSegment=NULL;
	sequence= nextsequence;
	sequence->Name= (char *) rcalloc(strlen(token)+1+4,sizeof(char),"ReadDNA.4");
	/* Extra 4 in length because we will append a copy-number */
	strcpy(sequence->Name,token);
      } else {
	sequence->RightSegment=NULL;
	sequence->LeftSegment=NULL;
	sequence->Name= (char *) rcalloc(strlen(token)+1+4,sizeof(char),"ReadDNA.5"); 
	/* Extra 4 in length because we will append a copy-number */
	strcpy(sequence->Name,token);
      }
      rflag++;
      nflag=0;
    }

  } while(fscanf(fp,"%s",token)!=EOF);
     
  if(nflag!= -1) {
    fprintf(stderr,"%s: Premature end-of-file in mechanism %s.\n",progid,mech);
    exit(-1);
  }
  
  /**** Now read the DNA parameters ******/

  /* First Get MOI */
  fscanf(fp,"%s %*s %d",token,&moi);
  if(strcmp(token,"MOI")!=0){
    fprintf(stderr,"%s: Expected an MOI token instead of %s in mechanism %s.\n"
	    ,progid,token,mech);
    exit(-1);
  }
#ifdef RMM_MODS
  /* Reset the MOI using the command line option */
  if (global_MOI != 0 && moi != 1) moi = global_MOI;
#endif


  name=        (char **) rcalloc(rflag,sizeof(char *),"ReadDNA.6");
  parameters=  (char **) rcalloc(rflag,sizeof(char *),"ReadDNA.7");

  for(i=0; i<rflag; i++){
    name[i]=       (char *) rcalloc(80,sizeof(char),"ReadDNA.8");
    parameters[i]= (char *) rcalloc(80,sizeof(char),"ReadDNA.9");

    if(fscanf(fp,"%s %d %s %s %s %s"
	      ,name[i],&length,unit,direction,type,parameters[i])==EOF){
      fprintf(stderr,"%s: Premature end-of-file in mechanism %s.\n",
	      progid,mech);
      exit(-1);
    }

    /**** Find name ****/

    sequence= &Sequence[NSequences-1];
    while(sequence!=NULL && strcmp(sequence->Name,name[i])!=0) sequence= sequence->RightSegment;

    if(sequence==NULL){
      fprintf(stderr,"%s: Parameters for unknown segment %s found in mechanism %s.\n",progid,name[i],mech);
      exit(-1);
    }

    if(strcasecmp(direction,"LEFT")==0)           sequence->Direction=LEFT;
    else if(strcasecmp(direction,"RIGHT")==0)     sequence->Direction= RIGHT;
    else {
      fprintf(stderr,"%s: Unknown direction type %s found for segment %s of mechanism %s.\n",
	      progid,direction,sequence->Name,mech);
      exit(-1);
    }

    if(strncasecmp(unit,"nt",2)==0 || strncasecmp(unit,"nucleotide",10)==0) mult=1;
    else
      if(strncasecmp(unit,"bp",2)==0 || strncasecmp(unit,"basepair",8)==0) mult=1;
      else
	if(strncasecmp(unit,"kb",2)==0 || strncasecmp(unit,"kilobase",8)==0) mult=1000;
	else {
	  fprintf(stderr,"%s: Unknown unit type %s found for segment %s of mechanism %s.\n",
		  progid,unit,sequence->Name,mech);
	  exit(-1);
	}

    sequence->Length= mult*length;      


    /**** When copying for MOI:
     ****
     **** New RNAP Queues for all!
     **** Full Copies of Shea data must be made
     ****/

    if(strcasecmp(type,"Promotor")==0)
      sequence->Type= DNA_Type_Promotor;
    else
      if(strcasecmp(type,"Terminator")==0)
	sequence->Type= DNA_Type_Terminator;
      else
	if(strcasecmp(type,"AntiTerminator")==0)
	  sequence->Type= DNA_Type_AntiTerminator;
	else
	  if(strcasecmp(type,"Gene")==0)
	    sequence->Type= DNA_Type_Coding;
          else
	    if(strcasecmp(type,"NonCoding")==0)
	      sequence->Type= DNA_Type_NonCoding;
	    else {
	      fprintf(stderr,"%s: Unknown segment type %s found for segment %s of mechanism %s.\n",
		      progid,type,sequence->Name,mech);
	      exit(-1);
	    }

  } /** End for rflag **/


  /**** Now we make multiple infections ******/

  firstcopy= NSequences-1;

  for(i=1; i<moi; i++){
    DNA *fseq;

    NSequences++;
    Sequence= (DNA *) rrealloc(Sequence,NSequences,sizeof(DNA),"ReadDNA.10");
    sequence=(DNA *) &Sequence[NSequences-1];

  /****
   *
   * We've realloc'ed....we gotta reset 
   * the pointers to the first segment of 
   * the first sequence 
   ****/
    
    for(j=0; j<NSequences-1; j++){
      fseq= &Sequence[j];
      if(fseq->RightSegment != NULL)
	fseq->RightSegment->LeftSegment=fseq;
    }

    fseq= &Sequence[firstcopy];

    sequence->Name= (char *) rcalloc(strlen(fseq->Name)+4+1,sizeof(char),"ReadDNA.11");
    sprintf(sequence->Name,"%s.%d",fseq->Name,i);
    sequence->Length       = fseq->Length;
    sequence->Direction    = fseq->Direction;
    sequence->Type         = fseq->Type;
    sequence->LeftSegment  = NULL;
    
    fseq=fseq->RightSegment;
    while(fseq!=NULL){

      sequence->RightSegment= (DNA *) rcalloc(1,sizeof(DNA),"ReadDNA.12");
      sequence->RightSegment->Name= (char *) rcalloc(strlen(fseq->Name)+4+1,sizeof(char),"ReadDNA.13");
      sprintf(sequence->RightSegment->Name,"%s.%d",fseq->Name,i);
      sequence->RightSegment->Length       = fseq->Length;
      sequence->RightSegment->Direction    = fseq->Direction;
      sequence->RightSegment->Type         = fseq->Type;

      sequence->RightSegment->LeftSegment=sequence;
      sequence=sequence->RightSegment;
      fseq=fseq->RightSegment;
    }
    sequence->RightSegment=NULL;
  }
      
  /**** Now we read in parameter files ******/

  /******
   ****** Currently, in all cases except promotors,
   ****** all copies of a given sequence can 
   ****** share common segment data-structures
   ****** in order to conserve memory.
   ****** 
   ****** Promotors on the other hand MUST have
   ****** separate ACKERSSHEA data structures since
   ****** the state of the promotor is time-dependent.
   ******
   ******/

  for(i=0; i<moi; i++){
   
   for(j=0; j<rflag; j++){

     if(i>0)
       sprintf(buffer,"%s.%d",name[j],i);
     else 
       sprintf(buffer,"%s",name[j]);

     sequence= &Sequence[firstcopy+i];
     while(sequence!=NULL && strcmp(sequence->Name,buffer)!=0)
       sequence= sequence->RightSegment;
     
     if(sequence==NULL){
       fprintf(stderr,"@@@ NULL sequence found when searching %s for %s.\n",
	       Sequence[firstcopy+i].Name,name[j]);
       exit(-1);
     }
     
     switch(sequence->Type){
       
     case DNA_Type_Promotor:
       
       prom=((PROMOTOR *) rcalloc(1,sizeof(PROMOTOR),"ReadDNA.14"));
       sequence->DNAStruct= (void *)  prom;
       prom->RNAPQueue=NULL;
       ReadPromotorData(i,prom,parameters[j]);
       break;
       
     case DNA_Type_Terminator:
       seg= (SEGMENT *) rcalloc(1,sizeof(SEGMENT),"ReadDNA.15");
       sequence->DNAStruct= (void *) seg; 
       seg->RNAPQueue=NULL;
       seg->SegmentFunc=SubmitTermination;
       if(i==0)
	 ReadTerminatorData(seg,parameters[j]);
       else
	 seg->SegmentData= (void *) FindSegData(&Sequence[firstcopy],name[j]);
       break;

     case DNA_Type_AntiTerminator:
       seg= (SEGMENT *) rcalloc(1,sizeof(SEGMENT),"ReadDNA.16");
       sequence->DNAStruct= (void *) seg; 
       seg->RNAPQueue=NULL;
       seg->SegmentFunc=SubmitAntiTermination;
       if(i==0)  
	 ReadAntiTerminatorData(seg,parameters[j]);
       else
	 seg->SegmentData= (void *) FindSegData(&Sequence[firstcopy],name[j]);
       break;

     case DNA_Type_Coding:
       seg= (SEGMENT *) rcalloc(1,sizeof(SEGMENT),"ReadDNA.17");
       sequence->DNAStruct= (void *) seg; 
       seg->RNAPQueue=NULL;
       seg->SegmentFunc=SubmitProduceTranscript;
       ReadCodingData(seg,parameters[j]);
/*
       if(i==0)  

       else
	 seg->SegmentData= (void *) FindSegData(&Sequence[firstcopy],name[j]);
*/
       break;

     case DNA_Type_NonCoding:
       seg= (SEGMENT *) rcalloc(1,sizeof(SEGMENT),"ReadDNA.19");
       sequence->DNAStruct= (void *) seg; 
       seg->RNAPQueue=NULL;
       seg->SegmentFunc=SubmitSimpleJumpSegment;
       if(i==0)  
	 ReadNonCodingData(seg,parameters[j]);
       else
	 seg->SegmentData= (void *) FindSegData(&Sequence[firstcopy],name[j]);
       break;
       
     default:
       fprintf(stderr,"%s: Unknown DNA_Type %d found in mechanism %s\n",
	       progid,sequence->Type,mech);
       exit(-1);
     }    
   } /* for j */
   
 }/* for i */

  for(i=0; i<rflag; i++){
    free(parameters[i]);
    free(name[i]);
  }
  free(parameters);
  free(name);
    

  fclose(fp);
}
  
void *FindSegData(Seq,name)
DNA *Seq;
char *name;
{
  DNA     *seq;
  SEGMENT *seg;

  seq=Seq;
  while(seq!=NULL){
    if(strcmp(seq->Name,name)==0){
      if(seq->Type==DNA_Type_Promotor){
	fprintf(stderr,"%s: FindSegData() passed an incorrect data type named %s\n",
		progid,name);
	exit(-1);
      }

      seg= (SEGMENT *) seq->DNAStruct;

      return((void *) seg->SegmentData);
    }
    seq=seq->RightSegment;
  }

  fprintf(stderr,"%s: FindSegData() unable to find segment named %s in sequence.\n",
	  progid,name);

  exit(-1);  
}

/***********************
 * 
 * Data Readers
 *
 ***********************/

void ReadPromotorData(copy,prom,params)
int        copy;
PROMOTOR  *prom;
char      *params;
{
  int i;
  char token1[85],token2[85],token3[85];
  
  FILE *fp;
  
  void ReadSheaAckers();
  void ReadIsoData();

  /***** Promotor Data files contain the following information:
   *****
   ***** 1) Transcription Direction
   ***** 2) Name of SHEA/ACKERS Data file
   ***** 3) Name of attendant Isomerzation rate file
   *****
   ***** There are three ways in which a S/A Data File might be referenced 
   ***** more than once:
   *****
   ***** 1) Two promotors might feed off the same set of operators (witness PR and PRM)
   ***** 2) The operator site might be genetically duplicated on the same piece of DNA
   ***** 3) The operator site might exist on multiple copies of the gene.
   *****
   ***** The user must help distiguish between the first 2 cases. If the operator is being
   ***** used by two or more promotion sites then it is sufficient for each promotor to
   ***** reference the same S/A file. On the other hand if the operator is duplicated on the 
   ***** same piece of DNA then the S/A files MUST BE GIVE DIFFERENT NAMES.
   *****
   ***** In the third case, duplicates will be made automatically be the program.
   *****/

  fp= OpenFile(params,"r");

  fgets(token3,81,fp);
  while(token3[0]== '%' || isspace(token3[0]))  fgets(token3,81,fp);
  
  sscanf(token3,"%s %*s %s",token1,token2);

  if(strcasecmp(token1,"TranscriptionDirection")!=0){
    fprintf(stderr,"%s: First Token in a promotor data file MUST be TranscriptionDirection.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }

  prom->TranscriptionDirection= ( strcasecmp(token2,"LEFT")==0  ? LEFT  :
				  strcasecmp(token2,"RIGHT")==0 ? RIGHT : -1);

  if(prom->TranscriptionDirection== -1){
    fprintf(stderr,"%s: Invalid Transcription Direction %s found in parameter file %s.\n",
	    progid,token2,params);
    exit(-1);
  }
 
  fscanf(fp,"%s %*s %s",token1,token2);

  if(strcasecmp(token1,"SheaAckers")!=0){
    fprintf(stderr,"%s: Second Token in a promotor data file MUST be SheaAckers.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }

  sprintf(token3,"%s.%d",token2,copy);

  for(i=0; i<NOperators; i++)
    if(strcmp(token3,Operator[i].Name)==0)
      break;
  

  if(i==NOperators){
    fprintf(stderr,"@@@ Reading SheaAckers for copy %d\n",copy);
    ReadSheaAckers(copy,token2);
  }

  prom->Data= i;

  fscanf(fp,"%s %*s %s",token1,token2);
  
  if(strcasecmp(token1,"IsoData")!=0){
    fprintf(stderr,"%s: Third Token in a promotor data file MUST be IsoData.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }

  ReadIsoData(prom,token2);

  fclose(fp);
}

#define kcal_per_joule (0.001/4.184)
#define RT (8.314*310.15*kcal_per_joule)            /* (J/(mol K))*K = J/mol */   

void ReadSheaAckers(copy,file)
int copy;
char *file;
{
  
  int       i,j;
  int       NSites,NConfigs;
  double   *deltaG;
  int     **configs;
  char      speciesname[80];
  FILE     *fp;
  SHEADATA *oper;

  int  FindSpecies();
  void AddSpecies();
  void MakeConfigList();

  if(NOperators==0) Operator= (SHEADATA *) rcalloc(1,sizeof(SHEADATA),"ReadSheaAckers.1");
  else              Operator= (SHEADATA *) rrealloc(Operator,NOperators+1,sizeof(SHEADATA),"ReadSheaAckers.1");

  oper= &Operator[NOperators];
  NOperators++;
    
  oper->Name= (char *) rcalloc(strlen(file)+1+4,sizeof(char),"ReadSheaAckers.2");
  sprintf(oper->Name,"%s.%d",file,copy);

  fp=OpenFile(file,"r");

  fscanf(fp,"%d %d",&NSites,&NConfigs);
  
  oper->NSites=NSites;
  oper->NConfigs=NConfigs;
  
  if(NSites<=0 || NConfigs <=0){
    fprintf(stderr,"%s: The number of binding sites and configurations must be more that 0 for this object.\n",progid);
    fprintf(stderr,"%s: Sites= %d, Configs= %d in file %s\n",progid,NSites,NConfigs,file);
    exit(-1);
  }
  
  deltaG= (double *)  rcalloc(NConfigs,sizeof(double),"ReadSheaAckers.3");
  configs= (int **)   rcalloc(NConfigs,sizeof(int *),"ReadSheaAckers.4");
  for(i=0; i<NConfigs; i++){

    configs[i]= (int *) rcalloc(NSites, sizeof(int),"ReadSheaAckers.5");
    
    for(j=0; j<NSites; j++){
      if(fscanf(fp,"%s",speciesname)==EOF){
        fprintf(stderr,"%s: Short Ackers/Shea Specification in %s. Last Full Config= %d\n",progid,file,i);
        exit(-1);
      }
      configs[i][j]=FindSpecies(speciesname);
      if(configs[i][j]== NSpecies)
	AddSpecies(speciesname);
    }
    if(fscanf(fp,"%lf",&deltaG[i])==EOF){
      fprintf(stderr,"%s: Short Ackers/Shea Specification in %s. Last Full Config= %d\n",progid,file,i);
      exit(-1);
    }
    /* Preconversion */
    deltaG[i]= exp(-deltaG[i]/RT);
  }
  
  oper->DeltaG=deltaG;
  oper->CurrentState=0;

  MakeConfigList(oper,configs);

  for(i=0; i<NConfigs; i++) free(configs[i]);
  free(configs);

  fclose(fp);
}

void MakeConfigList(data,configs)
SHEADATA *data;
int **configs;
{
  int i,j,pos;
  int nspec;
  int    *scnt=NULL,*sind=NULL;

  scnt= (int *) rcalloc(NSpecies,sizeof(int),"MakeConfigList.1");
  sind= (int *) rcalloc(NSpecies,sizeof(int),"MakeConfigList.2");

  data->CList= (int **) rcalloc(data->NConfigs,sizeof(int *),"MakeConfigList.3");

  for(i=0; i<data->NConfigs; i++){

    nspec=0;
    for(j=0; j<data->NSites; j++)
      if(configs[i][j]!= -1){
	scnt[configs[i][j]] += 1;
	if(sind[configs[i][j]]==0)
	  nspec++;
	sind[configs[i][j]]=1;
      }


    data->CList[i]= (int *) rcalloc(1+nspec*2,sizeof(int),"MakeConfigList.4");
    data->CList[i][0]=nspec;
    pos=1;
    for(j=0; j<NSpecies; j++)
      if(sind[j]==1){
	data->CList[i][pos]=j;
	data->CList[i][pos+1]=scnt[j];
	sind[j]=scnt[j]=0;
	pos+=2;
      }
  }

  free(scnt);
  free(sind);

}

void ReadIsoData(prom,file)
PROMOTOR *prom;
char     *file;
{
  int   i;
  int   NConfigs;
  FILE *fp;
  FILE *OpenFile();
  
  fp= OpenFile(file,"r");
  
  NConfigs= Operator[prom->Data].NConfigs;
  prom->IsoRate= (double *) rcalloc(NConfigs,sizeof(double),"ReadIsoData.1");

  for(i=0; i<NConfigs; i++)
    if(fscanf(fp,"%lf",&prom->IsoRate[i])==EOF){
      fprintf(stderr,"%s: Short IsoFile %s found linked to operator %s.\n",
	     progid,file,Operator[prom->Data].Name);
      exit(-1);
    }
  
  fclose(fp);  
}

void ReadTerminatorData(seg,params)
SEGMENT  *seg;
char     *params;
{
  char      token1[81],token2[81],token3[90];
  double    mult;
  TERMDATA *tdata;
  FILE     *fp,*OpenFile();
  double    FindTimeUnit();
  int       FindSpecies();
  void      AddSpecies();

  tdata= (TERMDATA *) rcalloc(1,sizeof(TERMDATA),"ReadTerminatorData.1");
  seg->SegmentData= (void *) tdata;

  fp= OpenFile(params,"r");
  
  fgets(token3,81,fp);
  while(token3[0]== '%' || isspace(token3[0]))  fgets(token3,81,fp);
  
  sscanf(token3,"%s %*s %s",token1,token2);

  if(strcasecmp(token1,"TermModifier")!=0){
    fprintf(stderr,"%s: First Token in a terminator data file MUST be TermModifier.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }


  tdata->SpeciesIndex= FindSpecies(token2);
  if(tdata->SpeciesIndex==NSpecies) AddSpecies(token2);


  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"BaseFallOffRate")!=0){
    fprintf(stderr,"%s: Second Token in a terminator data file MUST be BaseFallOffRate.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of BaseFallOffRate is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->BaseFallOffRate= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for BaseFallOffRate in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->BaseFallOffRate *= mult;


  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"BaseRNAPMotion")!=0){
    fprintf(stderr,"%s: Third Token in a terminator data file MUST be BaseRNAPMotion.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of BaseRNAPMotion is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->BaseRNAPMotion= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for BaseRNAPMotion in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->BaseRNAPMotion *= mult;


  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"AntiTerminatedFallOffRate")!=0){
    fprintf(stderr,"%s: Fourth Token in a terminator data file MUST be AntiTerminatedFallOffRate.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of AntiTerminatedFallOffRate is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->AntiTerminatedFallOffRate= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for AntiTerminatedFallOffRate in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->AntiTerminatedFallOffRate *= mult;

  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"AntiTerminatedRNAPMotion")!=0){
    fprintf(stderr,"%s: Fifth Token in a terminator data file MUST be AntiTerminatedRNAPMotion.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of AntiTerminatedRNAPMotion is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->AntiTerminatedRNAPMotion= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for AntiTerminatedRNAPMotion in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->AntiTerminatedRNAPMotion *= mult;
  
  fclose(fp);
}


void ReadAntiTerminatorData(seg,params)
SEGMENT  *seg;
char     *params;
{
  char          token1[85],token2[85],token3[90];
  double        mult;
  ANTITERMDATA *tdata;
  FILE         *fp,*OpenFile();
  double        FindTimeUnit();
  int           FindSpecies();
  void          AddSpecies();

  tdata= (ANTITERMDATA *) rcalloc(1,sizeof(ANTITERMDATA),"ReadAntiTerminatorDate");
  seg->SegmentData= (void *) tdata;

  fp= OpenFile(params,"r");

  fgets(token3,81,fp);
  while(token3[0]== '%' || isspace(token3[0]))  fgets(token3,81,fp);
  
  sscanf(token3,"%s %*s %s",token1,token2);

  if(strcasecmp(token1,"TermModifier")!=0){
    fprintf(stderr,"%s: First Token in an antiterminator data file MUST be TermModifier.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }

  tdata->SpeciesIndex= FindSpecies(token2);
  if(tdata->SpeciesIndex==NSpecies) AddSpecies(token2);

  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"UnBoundRNAPMotion")!=0){
    fprintf(stderr,"%s: Second Token in an antiterminator data file MUST be UnBoundRNAPMotion.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of UnBoundRNAPMotion is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->UnBoundRNAPMotion= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for UnBoundRNAPMotion in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->UnBoundRNAPMotion *= mult;


  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"BindingRate")!=0){
    fprintf(stderr,"%s: Third Token in an antiterminator data file MUST be BindingRate.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of BindingRate is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->BindingRate= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for BindingRate in antitermination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->BindingRate *= mult;


  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"BoundRNAPMotion")!=0){
    fprintf(stderr,"%s: Fourth Token in an antiterminator data file MUST be BoundRNAPMotion.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of BoundRNAPMotion is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->BoundRNAPMotion= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for BoundRNAPMotion in antitermination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->BoundRNAPMotion *= mult;

  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"UnBindingRate")!=0){
    fprintf(stderr,"%s: Fifth Token in an antiterminator data file MUST be UnBindingRate.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of UnBindingRate is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->UnBindingRate= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for UnBindingRate in antitermination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->UnBindingRate *= mult;
  
  fclose(fp);
}

void ReadCodingData(seg,params)
SEGMENT  *seg;
char     *params;
{

  char          token1[85],token2[85],token3[90];
  double        mult;
  CODINGDATA   *tdata;
  FILE         *fp,*OpenFile();
  double        FindTimeUnit();
  int           FindSpecies();
  void          AddSpecies();

  tdata= (CODINGDATA *) rcalloc(1,sizeof(CODINGDATA),"ReadCodingData.1");
  seg->SegmentData= (void *) tdata; 

  fp= OpenFile(params,"r");

  fgets(token3,81,fp);
  while(token3[0]== '%' || isspace(token3[0]))  fgets(token3,81,fp);
  
  sscanf(token3,"%s %*s %s",token1,token2);

  if(strcasecmp(token1,"Product")!=0){
    fprintf(stderr,"%s: First Token in a coding data file MUST be Product.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }

  tdata->SpeciesIndex= FindSpecies(token2);
  if(tdata->SpeciesIndex==NSpecies) AddSpecies(token2);

  fscanf(fp,"%s %*s %s %s",token1,token2,token3);

  if(strcasecmp(token1,"mRNADegradationRate")!=0){
    fprintf(stderr,"%s: Second Token in a coding data file MUST be mRNADegradationRate.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token1,params);
    exit(-1);
  }
  
  if(isdigit(token2[0])==0){
    fprintf(stderr,"%s: First argument of mRNADegradationRate is not a number.\n",progid);
    fprintf(stderr,"%s: You have a %s token in parameter file %s.\n",progid,token2,params);
    exit(-1);
  }

  tdata->mRNADegradationRate= atof(token2);

  mult= FindTimeUnit(token3);
  if(mult<0.0){
    fprintf(stderr,"%s: Illegal time unit %s found for mRNADegradationRate in termination data file %s.\n",
	    progid,token3,params);
    exit(-1);
  }

  tdata->mRNADegradationRate *= mult;
fclose(fp);
}

void ReadNonCodingData(seg,params)
SEGMENT  *seg;
char     *params;
{

  seg->SegmentData= (void *) NULL;

}








