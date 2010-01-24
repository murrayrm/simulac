# Makefile for simulac
# APA, circa 1996

# Compiler options
CC     = gcc
CFLAGS = -g -D_H_MALLOC -DRMM_MODS
LIBS   = -lm -lc

# Files needed to create the simulator
OBJS   = Util.o Memory.o Kinetics.o PromotorDynamics.o SegmentDynamics.o \
  ReactionManager.o ParseDataBase.o CellManager.o cmdline.o

# Rules for building the simulation
all : Simulac

# Simulation program
Simulac : $(OBJS) Main.o
	  $(CC) $(CFLAGS) -o Simulac $(OBJS) Main.o $(LIBS)

# Tar file for the source code
simulac.tgz:
	tar zcf simulac.tgz Makefile README *.c *.h 

# Files for command line options
cmdline.o: cmdline.c cmdline.h
cmdline.h cmdline.c: simulac.ggo
	gengetopt -u -i simulac.ggo 

# Source code dependencies
Util.o:			Util.c	DataStructures.h
Memory.o:		Memory.c	DataStructures.h
Kinetics.o :            Kinetics.c         DataStructures.h
CellManager.o :         CellManager.c      DataStructures.h
PromotorDynamics.o :    PromotorDynamics.c DataStructures.h
SegmentDynamics.o :     SegmentDynamics.c  DataStructures.h
ReactionManager.o :     ReactionManager.c  DataStructures.h
ParseDataBase.o :       ParseDataBase.c    DataStructures.h
Main.o :                Main.c             DataStructures.h cmdline.h

# Rules to clean things up
tidy:;		rm -f *.o *~ core
clean: tidy;	rm cmdline.[ch] Simulac 
