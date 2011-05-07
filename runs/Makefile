
NRUN=1

MF=run$(NRUN)/run.cc
MFO=$(MF:.cc=.o)
OUT=run$(NRUN)/out

DIR=$(shell pwd)/..

ifneq ($(strip $(shell which g++)),)
	COMPILER=g++
endif

ifneq ($(strip $(shell which icpc)),)
	COMPILER=icpc
endif

INCLUDEDIR=-I../external -I../lib

ISING1DLIB=$(DIR)/lib/matrices.o $(DIR)/lib/error.o $(DIR)/lib/io.o $(DIR)/lib/common.o $(DIR)/lib/ising1D.o $(DIR)/external/random.o

UNAME := $(shell uname)
MACHINE := $(shell uname -m)

ifeq ($(UNAME), Linux)
   ifeq ($(MACHINE),x86_64)
	LIB=-llapack
#      LIB=-lmkl_solver_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
   endif
   ifeq ($(MACHINE),i686)
	LIB=-llapack
#      LIB=-lmkl_solver -Wl,--start-group -lmkl_intel -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
#   LIB=-lmkl_core -lmkl_lapack -lmkl_solver -lmkl_def -lmkl_intel -lmkl_intel_thread -liomp5 -lpthread
   endif
endif

ifeq ($(UNAME), Darwin)
   LIB += -framework veclib
endif

COMMAND:= $(COMPILER) $(FLAGS) $(LIB) $(ISING1DLIB) $(MFO) -o $(OUT)

default: FLAGS += -DDEBUG
default: $(MFO) $(ISING1DLIB)
	@echo $(MF)
	$(COMMAND)

fast: FLAGS += -O3
fast: $(MFO) $(ISING1DLIB)
	$(COMMAND)

clean:
	@rm -f *.o out *~
	@rm -f run*/*.o run*/*out run*/*~

%.o: %.cc
	echo $(DIR)
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $< -o $*.o


#@for dir in $(DIRLIB);do \
#$(MAKE) -C$${dir}; \
#done;