# Makefile.am - automake configuration file for simulac 
# RMM, 14 Mar 2010

# Flags to enable RMM modifications
AM_CPPFLAGS = -g -D_H_MALLOC -DRMM_MODS

# Programs and libraries built in this directory
bin_PROGRAMS = Simulac

# Rules for building simulator
Simulac_SOURCES = Main.c Util.c Memory.c Kinetics.c PromotorDynamics.c \
  SegmentDynamics.c ReactionManager.c ParseDataBase.c CellManager.c \
  DataStructures.h Memory.h Util.h param.c param.h \
  simulac.ggo cmdline.c cmdline.h
BUILT_SOURCES = cmdline.c cmdline.h

# Rule for creating gengetopt files
cmdline.h cmdline.c: simulac.ggo
	gengetopt -u -i simulac.ggo 
