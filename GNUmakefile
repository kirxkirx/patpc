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

test: patpc
	@echo "Running tests"
	./patpc test_events.txt 2>&1 | grep 'The peak [[:alpha:]]* is at period' | awk '{if ( sqrt( ($$7-5711.5)*($$7-5711.5) ) < 15.0 ) print "PASS" ;else print "FAIL" }' | { grep --quiet "FAIL" && echo "TEST FAIED" && exit 1 || true; }
	./patpc test_events.txt 10000 100 0.2 2>&1 | grep 'The peak [[:alpha:]]* is at period' | awk '{if ( sqrt( ($$5711.5-0)*($$7-5711.5) ) < 15.0 ) print "PASS" ;else print "FAIL" }' | { grep --quiet "FAIL" && echo "TEST FAIED" && exit 1 || true; }
	./patpc test_events.txt 501.486 2>&1 | grep 'single-trial probability' | awk -F' is ' '{print $$2}' | awk '{if ( $$1 > 0.05 ) print "PASS" ;else print "FAIL" }' | { grep --quiet "FAIL" && echo "TEST FAIED" && exit 1 || true; }
	# end with this test as it will leave behind good looking artifacts
	./patpc test_events.txt 10000 100 2>&1 | grep 'The peak [[:alpha:]]* is at period' | awk '{if ( sqrt( ($$7-5711.5)*($$7-5711.5) ) < 15.0 ) print "PASS" ;else print "FAIL" }' | { grep --quiet "FAIL" && echo "TEST FAIED" && exit 1 || true; }

test_gsl: patpc simulator test_patpc.sh
	./test_patpc.sh	
	
	