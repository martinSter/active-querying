SRC = .
CFLAGS = -W -Wall -DTIME -Ofast -march=native
LDFLAGS = 
CC = gcc

OBJ1 = o/run.o o/sim.o o/detsir.o o/misc.o o/inf.o o/query.o o/heap.o o/quick.o o/pcg_rnd.o

all : run

run: $(OBJ1)
	$(CC) $(LDFLAGS) -o $@ $^ -lm

o/run.o : $(SRC)/run.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/run.c -o $@
    
o/sim.o : $(SRC)/sim.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/sim.c -o $@

o/detsir.o : $(SRC)/detsir.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/detsir.c -o $@

o/misc.o : $(SRC)/misc.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/misc.c -o $@
    
o/inf.o : $(SRC)/inf.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/inf.c -o $@
    
o/query.o : $(SRC)/query.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/query.c -o $@

o/heap.o : $(SRC)/heap.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/heap.c -o $@

o/quick.o : $(SRC)/quick.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/quick.c -o $@

o/pcg_rnd.o : $(SRC)/pcg_rnd.c $(SRC)/run.h $(SRC)/Makefile
	$(CC) $(CFLAGS) -c $(SRC)/pcg_rnd.c -o $@
