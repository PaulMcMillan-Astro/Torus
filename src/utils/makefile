.cc.o:
	g++ -c -I../ -I . -I ../../../libebf_c_cpp-0.0.3/src/cpp/ $*.cc

Numerics.o : Numerics.cc Numerics.h Numerics.templates

libutil.a : CHB.o Compress.o Err.o Numerics.o PJMCoords.o PJMebf.o PJMNum.o Random.o jjb_utils.o WDMath.o
	ar -qs libutil.a CHB.o Compress.o Err.o Numerics.o PJMCoords.o PJMebf.o PJMNum.o Random.o jjb_utils.o WDMath.o
