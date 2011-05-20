
ifneq ($(strip $(shell which g++)),)
        COMPILER=g++
	FLAGS = -fopenmp
endif

ifneq ($(strip $(shell which icpc)),)
        COMPILER=icpc
	FLAGS = -openmp
endif


SOURCE = $(shell ls *.cc)

OBJECTS = $(SOURCE:.cc=.o)

HEADERS = $(shell ls *.hh)
HEADERS_CODE = $(shell ls *.hpp)

INCLUDEDIR=-I../external/
FLAGS += -O3 -DUSUAL_ALG
#-DDK_ALG
#-Wall -pedantic

default: $(HEADERS) $(HEADERS_CODE) $(OBJECTS)

blas: FLAGS += -DBLAS
blas: default

clean:
	@rm -f *.o *~

%.o: %.cc
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $<