
N=4

MF=run.cc

default: $(MF)
	make -C../ NRUN=$(N)

fast: $(MF)
	make -C../ NRUN=$(N) fast

clean:
	@rm -f *.o out *~
