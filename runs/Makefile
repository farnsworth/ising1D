
NRUN=1

MF=run$(NRUN)/run.cc
OUT=run$(NRUN)/out

MFO=$(MF:.cc=.o)

COMPILER=g++

DIRLIB=../external ../lib
DIRLIBFLAG=-L../external -L../lib
INCLUDEDIR=-I../external -I../lib

LIB=-lrandom.o -lmatrices.o -lerror.o -lio.o -lcommon.o -lising1D.o

MACFLAG=-framework veclib
FLAGS=-DDEBUG

default: $(MFO)
	@echo $(MF)
	@for dir in $(DIRLIB);do \
	$(MAKE) -C$${dir}; \
	done;
	$(COMPILER) $(MACFLAG) $(DIRLIBFLAG) $(LIB) $(MFO) -o $(OUT)


fast: FLAGS=-O3
fast: $(MFO)
	@echo $(MF)
	@for dir in $(DIRLIB);do \
	$(MAKE) -C$${dir} FLAGS=$(FLAGS); \
	done;
	$(COMPILER) $(MACFLAG) $(DIRLIBFLAG) $(LIB) $(MFO) -o $(OUT)

clean:
	@rm -f *.o out *~
	@rm -f run*/*.o run*/*out run*/*~

%.o: %.cc
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $< -o $*.o
