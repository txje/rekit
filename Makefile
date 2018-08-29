CC     = gcc
LIBS   = -lz -lm -lhts

OBJECTS = rekit

all: $(OBJECTS)

rekit:
	$(CC) $(CFLAGS) -o rekit src/rekit.c src/io_utils.c src/bnx.c src/rmap.c src/lsh.c src/dtw.c src/hash.c src/sim.c src/digest.c src/cmap.c src/bam.c $(LIBS)

.PHONY: clean
clean:
	-rm $(OBJECTS)
