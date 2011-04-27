
NRUN=1

MF=run$(NRUN)/run.cc
MFO=$(MF:.cc=.o)
OUT=run$(NRUN)/out


COMPILER=g++

DIRLIB=../external ../lib
DIRLIBFLAG=-L../external -L../lib
INCLUDEDIR=-I../external -I../lib

LIB=-lrandom.o -lmatrices.o -lerror.o -lio.o -lcommon.o -lising1D.o

UNAME := $(shell uname)

FLAGS=-DDEBUG

ifeq ($(UNAME), Linux)

endif

ifeq ($(UNAME), Darwin)
   FLAGS += -framework veclib
endif



default: $(MFO)
	@echo $(MF)
	@for dir in $(DIRLIB);do \
	$(MAKE) -C$${dir}; \
	done;
	$(COMPILER) $(FLAGS) $(DIRLIBFLAG) $(LIB) $(MFO) -o $(OUT)


fast: FLAGS=-O3
fast: $(MFO)
	@echo $(MF)
	@for dir in $(DIRLIB);do \
	$(MAKE) -C$${dir} FLAGS=$(FLAGS); \
	done;
	$(COMPILER) $(FLAGS) $(DIRLIBFLAG) $(LIB) $(MFO) -o $(OUT)

clean:
	@rm -f *.o out *~
	@rm -f run*/*.o run*/*out run*/*~

%.o: %.cc
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $< -o $*.o
