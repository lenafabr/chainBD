OBJS = chainutil.o  keys.o readkey.o genutil.o mt19937.o inputparams.o browndyn.o

FC = gfortran
FLAGS = -O3 
DEBUGFLAGS =
LIBS =  -lblas -llapack

# executable name
EXE = ../chainBD.exe
TESTEXE = ../testrun.exe

.SUFFIXES: .o .f90 

# general suffix rules
%.o: %.f90 
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

test: $(OBJS) testmain.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) testmain.o -o $(TESTEXE) $(LIBS)

clean: 
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
chainutil.o: keys.o mt19937.o genutil.o
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o 
browndyn.o: mt19937.o chainutil.o
