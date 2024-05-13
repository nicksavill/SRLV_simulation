DEBUG =
CC = gcc
WARNINGS = -Wall -Wno-unused-but-set-variable
ifdef DEBUG
	OPT = -g
else
	OPT = -O3
endif
FLAGS = $(OPT) $(WARNINGS) -fopenmp
LIBS = -lm -lgsl -lgslcblas

sim: sim.o ngm_biased_sale.o det.o ibm.o func.o
	$(CC) $(FLAGS) $+ -o $@ $(LIBS) 
sim.o: sim.c R0.h
	$(CC) -c $(FLAGS) $<
ngm_biased_sale.o: ngm_biased_sale.c R0.h
	$(CC) -c $(FLAGS) $<
det.o: det.c R0.h
	$(CC) -c $(FLAGS) $<
ibm.o: ibm.c R0.h
	$(CC) -c $(FLAGS) $<
func.o: func.c R0.h
	$(CC) -c $(FLAGS) $<


clean:
	rm -f *.o
