CC   = g++
OPT  = -funroll-loops -fopenmp -O
INCL = -I$(MOCKINGDESI_BASE)/util -L$(MOCKINGDESI_BASE)/util

box2sky: box2sky.cpp
	$(CC) $(OPT) $(INCL) -o box2sky box2sky.cpp -lm

.PHONY: clean
clean:
	rm -f box2sky
