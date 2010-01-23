#
# 

CC     = gcc
# CFLAGS = -g -D_H_MALLOC -DRMM_MODS -I/usr/local/include
CFLAGS = -O -D_H_MALLOC -I/usr/local/include
LIBS   = -lm -lc
OBJS   = Util.o Memory.o Kinetics.o PromotorDynamics.o SegmentDynamics.o ReactionManager.o ParseDataBase.o CellManager.o cmdline.o
OOBJS  = Main.o
#

all : Simulac

Simulac : $(OBJS) Main.o
	  $(CC) $(CFLAGS) -o Simulac $(OBJS) Main.o $(LIBS)

Util.o:			Util.c	DataStructures.h
Memory.o:		Memory.c	DataStructures.h

Kinetics.o :            Kinetics.c         DataStructures.h

CellManager.o :         CellManager.c      DataStructures.h

PromotorDynamics.o :    PromotorDynamics.c DataStructures.h

SegmentDynamics.o :     SegmentDynamics.c  DataStructures.h

ReactionManager.o :     ReactionManager.c  DataStructures.h

ParseDataBase.o :       ParseDataBase.c    DataStructures.h

Main.o :                Main.c             DataStructures.h cmdline.h

cmdline.o: cmdline.c cmdline.h
cmdline.h: simulac.ggo
	gengetopt -u -i simulac.ggo 

clean :
	rm -f *.o *~ core
