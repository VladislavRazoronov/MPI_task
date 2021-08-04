
TempGrid.exe: main.o Functions.o
	mpic++	Functions.o	main.o	-o	TempGrid

Functions.o:
	mpic++	-c	../Functions.cpp

main.o:
	mpic++	-c	../main.cpp

.PHONY:	clean

clean:
	rm -f Functions.o
	rm -f main.o
