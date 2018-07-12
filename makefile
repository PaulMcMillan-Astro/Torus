################################################################################
#
# Makefile for the torus project
#
# written by Walter Dehnen at Oxford University, 1995-97
#            Paul McMillan thereafter
# address: Theoretical Physics, 1 Keble Road, Oxford OX1 3NP, UK
# e-mail:  w.dehnen1.physics.oxford.ac.uk
#
################################################################################

all: mains


#
# MACRO DEFINITIONS that depend on the archetcture
#
# directories
INC		= src/
SRC		= src/
SRCPOT		= src/pot/
SRCUTILS	= src/utils/
EXE		= ${SRC}mains/
OBJ		= obj/
BIN		= bin/
INCL		= src/utils/
LIB		= obj/
EBFINC		= ../libebf_c_cpp-0.0.3/include/
EBFLIB		= ../libebf_c_cpp-0.0.3/lib/

# Compiler
CPP	= g++

# flags for compiler for optimized & debug code
CFLAGS 	= -c -o $@ -O3 -ffast-math -I$(INC) -I$(SRCPOT) -I$(SRCUTILS) -I$(EBFINC)
CFLAGS_WD = -c -o $@ -O3 -ffast-math -I$(SRCUTILS)
MFLAGS	= -O3 -ffast-math -I$(INC) -I$(SRCPOT) -I$(INCL) -I$(EBFINC)

#This block was added for portability across Mac OSX and Linux
# With thanks to Cecilia Mateu.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    MFLAGS += -Wl,-rpath=$(EBFLIB)
endif

# flags for linker for optimized & debug code
#LOTH        = -o $(BIN)$@ -L$(LIB)              -lOther -lWD -lm
#LPOT       = -o $(BIN)$@ -L$(LIB)        -lPot -lOther -lWD -lm
LFULL	    = -o $(BIN)$@ -L$(LIB) -L$(EBFLIB) -lTorus -lebf_cpp #-lPot -lOther -lWD -lm -lebf_cpp
LDFLAGS     = $(LFULL)


# commands to put file into library
AR            = ar r
RL            = ranlib

#
# include the macro definitions, dependencies, and rules that are
# unspecific to the machine used
#

include	make.torus

clean:
	- rm -rf obj/* lib/* bin/* *~ src/*~ src/mains/*~ \
	src/pot/*~ src/utils/*~

halfclean:
	- rm -rf obj/* lib/* bin/*

################################################################################
