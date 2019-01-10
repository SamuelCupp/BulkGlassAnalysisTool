phony: all
all: glass.x

compiler = mpifort
options = -C

glass.x: main.o initMopac.o randomSeed.o structures.o vandcheck.o interfaces.o tally.o
	$(compiler) $(options) -o glass.x main.o structures.o initMopac.o vandcheck.o randomSeed.o interfaces.o tally.o
	
main.o: main.f90 initMopac.o randomSeed.o structures.o vandcheck.o interfaces.o tally.o
	$(compiler) -c $(options) main.f90
	
initMopac.o: initMopac.f90 structures.o interfaces.o vandcheck.o
	$(compiler) -c $(options) initMopac.f90

vandcheck.o: vandcheck.f90 structures.o
	$(compiler) -c $(options) vandcheck.f90

tally.o: tally.f90 structures.o
	$(compiler) -c $(options) tally.f90

interfaces.o: interfaces.f90 structures.o
	$(compiler) -c $(options) interfaces.f90

randomSeed.o: randomSeed.f90
	$(compiler) -c $(options) randomSeed.f90

structures.o: structures.f90
	$(compiler) -c $(options) structures.f90

phony: clean
clean:
	rm -f *.x
	rm -f *.o
	rm -f *.mod
	rm -f *.out
	rm -f *.aux
	rm -f *.arc

phony: rmdata
rmdata:
	rm *.dat
	rm init-*
