CC     = gcc
CFLAGS = -std=c99
LIBS   = -lz

OBJECTS = rekit

all: $(OBJECTS)

rekit:
	$(CC) $(CFLAGS) -o rekit src/rekit.c src/io_utils.c src/bnx.c src/rmap.c src/lsh.c src/dtw.c $(LIBS)

.PHONY: clean
clean:
	-rm $(OBJECTS)
