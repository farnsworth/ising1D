
ifneq ($(strip $(shell which g++)),)
        COMPILER=gcc
endif

ifneq ($(strip $(shell which icpc)),)
        COMPILER=icc
endif

default: random.o random.hh

%.o:%.c
	$(COMPILER) -DRISC -c random.c

clean:
	@rm -f *.o *~