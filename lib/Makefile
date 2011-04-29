
ifneq ($(strip $(shell which g++)),)
        COMPILER=g++
endif

ifneq ($(strip $(shell which icpc)),)
        COMPILER=icpc
endif

SOURCE = $(shell ls *.cc)

OBJECTS = $(SOURCE:.cc=.o)

HEADERS = $(shell ls *.hh)
HEADERS_CODE = $(shell ls *.hpp)

INCLUDEDIR=-I../external/
#FLAGS=-framework veclib
FLAGS =-O3

default: $(HEADERS) $(HEADERS_CODE) $(OBJECTS)

clean:
	@rm -f *.o *~

%.o: %.cc
	$(COMPILER) $(FLAGS) $(INCLUDEDIR) -c $<