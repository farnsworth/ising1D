
DIRLIB=external lib

RUNS=runs

ALLDIR=$(DIRLIB) $(RUNS)

#ALLDIR=$(shell ls -d */)

lib: force_look
	@for dir in $(DIRLIB);do \
	$(MAKE) -C$${dir}; \
	done;

runs: force_look
	$(MAKE) -C$(RUNS) all

all: force_look
	$(MAKE) lib
	$(MAKE) runs

clean:
	@rm -f *.o out *~
	@rm utilities/*~
	@for dir in $(ALLDIR);do \
	$(MAKE) -C$${dir} clean; \
	done;

force_look:
	@true