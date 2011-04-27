
DIRS=$(shell ls -d */)


#all:
#	@for dir in $(DIRLIB);do \
#	$(MAKE) -C$${dir}; \
#	done;

#all-fast:
#	@for dir in $(DIRLIB);do \
#	$(MAKE) -C$${dir}; \
#	done;

clean-all:
	@rm -f *.o out *~
	@for dir in $(DIRS);do \
	$(MAKE) -C$${dir} clean; \
	done;