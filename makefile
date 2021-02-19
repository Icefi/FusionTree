SD = src/
ID = include/
BD = build/

OF = bin/FusionTree

OBJS = $(BD)main.o $(BD)FusionTree.o

.PHONY: all clean re

all: $(OF)

re: clean all

$(BD)main.o: $(SD)main.cpp
	g++ $^ -I include -c -o $@

$(BD)FusionTree.o: $(SD)FusionTree.cpp
	g++ $^ -I include -c -o $@

$(OF): $(OBJS)
	g++ $^ -o $@

clean:
	rm -f $(BD)*.o $(OF)
