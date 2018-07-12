EBF = ../../libebf_c_cpp-0.0.3/src/cpp/
EBFLIB = ../../libebf_c_cpp-0.0.3/lib/

MFLAGS = -I . -I pot/ -I utils/ -I$(EBF)/
#This block was added for portability across Mac OSX and Linux
# With thanks to Cecilia Mateu.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    MFLAGS += -Wl,-rpath=$(EBFLIB)
endif


.cc.o:
	g++ -c -fopenmp -I . -I utils/ -I pot/ -I$(EBF) $*.cc

#Torus.o: Torus.h Torus.cc Point_ClosedOrbitCheby.h

eTorus.o: eTorus.cc eTorus.h

iTorus.o: iTorus.cc iTorus.h

resTorus.o: resTorus.cc resTorus.h

resTorus_c.o: resTorus_c.cc resTorus_c.h

resTorus_L.o: resTorus_L.cc resTorus_L.h

#toy_isochrone.o: toy_isochrone.cc toy_isochrone.h

libtorus2.a : eTorus.o iTorus.o resTorus.o resTorus_c.o resTorus_L.o
	ar -qs libtorus2.a eTorus.o iTorus.o resTorus.o resTorus_c.o resTorus_L.o

lindblad.exe : lindblad.cc libtorus2.a bar_pot.o
	g++ lindblad.cc  bar_pot.o $(MFLAGS)  -L$(EBFLIB) -L../obj/ -L . -L pot/ -L utils/ -ltorus2 -lTorus -lebf_cpp  -o Lindblad.exe

corot.exe : corot.cc libtorus2.a bar_pot.o
	g++ corot.cc  bar_pot.o -I . $(MFLAGS) -L$(EBFLIB) -L../obj/ -L .  -ltorus2 -lTorus -lebf_cpp  -o corot.exe

nores.exe : nores.cc libtorus2.a bar_pot.o
	g++ nores.cc  bar_pot.o $(MFLAGS) -L$(EBFLIB) -L../obj/ -L .  -ltorus2 -lTorus -lebf_cpp  -o nores.exe
