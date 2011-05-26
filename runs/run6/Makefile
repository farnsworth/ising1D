
N=6

MF=run.cc

default: $(MF)
	make -C../ NRUN=$(N)

omp:
	make -C../ NRUN=$(N) omp

clean:
	@rm -f *.o out *~
