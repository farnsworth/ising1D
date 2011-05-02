
FILES=$(shell ls check*.cc)
FILESOUT=$(FILES:.cc=.o)

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
      LIB=-lmkl_solver_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
   endif
   ifeq ($(MACHINE),i686)
      LIB=-lmkl_solver -Wl,--start-group -lmkl_intel -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread
#   LIB=-lmkl_core -lmkl_lapack -lmkl_solver -lmkl_def -lmkl_intel -lmkl_intel_thread -liomp5 -lpthread
   endif
endif

ifeq ($(UNAME), Darwin)
   LIB += -framework veclib
endif

default: $(FILES) $(FILESOUT)
	@for file in $(FILESOUT); do \
	$(COMPILER) $(FLAGS) $(LIB) $(ISING1DLIB) $${file} -o $${file/o/out}; \
	done;

clean:
	@rm -f *.o *.out *~ *.in

%.o: %.cc
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $< -o $*.o