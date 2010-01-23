#
# 

CC     = gcc
#CFLAGS = -g -DNDEBUG -p -O3 
CFLAGS = -O 
LIBS   = -lm -lc
OBJS   = Util.o Memory.o Kinetics.o PromotorDynamics.o SegmentDynamics.o ReactionManager.o ParseDataBase.o CellManager.o Main.o
OOBJS  = Main.o
#

all : Simulac

Simulac : $(OBJS)
	  $(CC) $(CFLAGS) -o Simulac $(OBJS) $(LIBS)

Util.o:			Util.c	DataStructures.h
Memory.o:		Memory.c	DataStructures.h

Kinetics.o :            Kinetics.c         DataStructures.h

CellManager.o :         CellManager.c      DataStructures.h

PromotorDynamics.o :    PromotorDynamics.c DataStructures.h

SegmentDynamics.o :     SegmentDynamics.c  DataStructures.h

ReactionManager.o :     ReactionManager.c  DataStructures.h

ParseDataBase.o :       ParseDataBase.c    DataStructures.h

Main.o :                Main.c             DataStructures.h

clean :
	rm -f *.o *~ core
