###########################################################################################################
#			Make file for artemide + harpy
#				modify the first section with your values
###########################################################################################################
# location of artemide
HOME       = $(PWD)

#PUT YOUR FORTRAN COMPILER
FCompilator=f95 
#PUT HERE extra flags for compilator (put "space" if not flags requared)
Fflags= -fopenmp -cpp
#Fflags=  

#options for COMILATOR to compile QCDinput. e.g. link to LHA

FOPT=-O3 -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -Wall
#### for optimization -O3 -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -Wall
#### for debuging -g -fbacktrace -fcheck=all -Wall -pedantic


################################################################### LIST OF FILES ####################################
SOURCEDIR       = $(HOME)/src
BIN		= $(HOME)/bin
OBJ		= $(HOME)/obj
MOD		= $(HOME)/mod

SOURCEFILES = \
$(SOURCEDIR)/IO_snowflake.f90 \
$(SOURCEDIR)/HexGrid.f90 \
$(SOURCEDIR)/EvolutionKernels.f90 \
$(SOURCEDIR)/SnowFlake.f90

CommonFiles=\
$(SOURCEDIR)/commonVariables.f90

ExtraFiles=\
$(SOURCEDIR)/ExpressionsForKernels.f90

OBJFILES = \
$(OBJ)/IO_snowflake.o \
$(OBJ)/HexGrid.o \
$(OBJ)/EvolutionKernels.o \
$(OBJ)/SnowFlake.o

################################################################### COMPILATION OF ARTEMIDE ####################################
FC=$(FCompilator) $(Fflags)

.PHONY: clean default obj program test

default: obj

obj: $(OBJFILES) $(SOURCEFILES) $(CommonFiles) $(ExtraFiles)

$(OBJ)/IO_snowflake.o: $(SOURCEDIR)/IO_snowflake.f90  $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/IO_snowflake.f90
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/HexGrid.o: $(SOURCEDIR)/HexGrid.f90 $(SOURCEDIR)/IO_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/HexGrid.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/EvolutionKernels.o: $(SOURCEDIR)/EvolutionKernels.f90 $(SOURCEDIR)/ExpressionsForKernels.f90 $(SOURCEDIR)/HexGrid.f90 $(SOURCEDIR)/IO_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/EvolutionKernels.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SnowFlake.o: $(SOURCEDIR)/SnowFlake.f90 $(SOURCEDIR)/EvolutionKernels.f90 $(SOURCEDIR)/HexGrid.f90 $(SOURCEDIR)/IO_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/SnowFlake.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

clean: 
	$(RM) a.out
	$(RM) count *.o *.mod
	$(RM) count $(OBJ)/*.o
	$(RM) count $(MOD)/*.mod

test:
	$(FC) $(HOME)/prog/TEST.f90 $(OBJFILES) $(FOPT) -I$(MOD)
	./a.out

program: 
	echo $(TARGET)
	$(FC) $(TARGET) $(OBJFILES) $(FOPT) -I$(MOD)
