
N=4

MF=run.cc

default: $(MF)
	make -C../ NRUN=$(N)

fast: $(MF)
	make -C../ NRUN=$(N) fast

omp:
	make -C../ NRUN=$(N) omp

clean:
	@rm -f *.o out *~
