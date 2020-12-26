# production
USE_OMP_OPTIONS := $(shell lib/test_openmp.sh)
# no OpenMP debug
#USE_OMP_OPTIONS =

# production
GSL_OPTIONS := $(shell lib/test_gsl.sh)
# no GSL debug
#GSL_OPTIONS = -DPATPC_NOGSL

#production
CFITSIO_OPTIONS := $(shell lib/test_cfitsio.sh)
# no CFITSIO debug
#CFITSIO_OPTIONS = -DPATPC_NOCFITSIO

all: clean simulator patpc

simulator: simulator.c
	cc -Wall -o simulator simulator.c -lm $(GSL_OPTIONS)

patpc: patpc.c count_lines_in_ASCII_file.h
	# debug
	#cc -Wall -o patpc patpc.c -g -fsanitize=address,leak,undefined -lm $(GSL_OPTIONS) $(CFITSIO_OPTIONS) $(USE_OMP_OPTIONS)
	# production
	cc -Wall -o patpc patpc.c -lm $(GSL_OPTIONS) $(CFITSIO_OPTIONS) $(USE_OMP_OPTIONS)

clean:
	rm -f simulator patpc *~ lib/*~

test: patpc simulator test_patpc.sh
	./test_patpc.sh	
	
	