###########################################################################################################
#			Make file for artemide + harpy
#				modify the first section with your values
###########################################################################################################
# location of artemide
HOME       = $(PWD)

#PUT YOUR FORTRAN COMPILER
FCompilator=f95 
#PUT HERE extra flags for compilator (put "space" if not flags requared)
#Fflags= -fopenmp -cpp
#Fflags=-O3 -fPIC -cpp -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -fopenmp
#FflagsPY= '-O3 -fPIC -cpp -march=native  -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -fopenmp'

Fpath=/usr/bin/f95
F77path=/usr/bin/f77

Fflags=-O3 -fPIC -cpp -fopenmp
FflagsPY= '-O3 -fPIC -cpp -fopenmp'

#options for COMILATOR to compile QCDinput. e.g. link to LHA
#### for optimization -O3 -fforce-addr -fstrength-reduce -fcaller-saves -funroll-loops -Wall
#### for debuging -g -fbacktrace -fcheck=all -Wall -pedantic


################################################################### LIST OF FILES ####################################
SOURCEDIR       = $(HOME)/src
BIN		= $(HOME)/bin
OBJ		= $(HOME)/obj
MOD		= $(HOME)/mod
PYDIR		= $(HOME)/PySnowflake

SOURCEFILES = \
$(SOURCEDIR)/IO_snowflake.f90 \
$(SOURCEDIR)/HexGrid.f90 \
$(SOURCEDIR)/EvolutionKernels.f90 \
$(SOURCEDIR)/LHA_alpha_snowflake.f90 \
$(SOURCEDIR)/SnowFlake_Model.f90 \
$(SOURCEDIR)/SnowFlake.f90

CommonFiles=\
$(SOURCEDIR)/commonVariables.f90

ExtraFiles=\
$(SOURCEDIR)/ExpressionsForKernels.f90\
$(SOURCEDIR)/ExpressionsForG2.f90\
$(SOURCEDIR)/ExpressionsForD2.f90

OBJFILES = \
$(OBJ)/IO_snowflake.o \
$(OBJ)/HexGrid.o \
$(OBJ)/EvolutionKernels.o \
$(OBJ)/LHA_alpha_snowflake.o \
$(OBJ)/SnowFlake_Model.o\
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

$(OBJ)/EvolutionKernels.o: $(SOURCEDIR)/EvolutionKernels.f90 $(ExtraFiles) $(SOURCEDIR)/HexGrid.f90 $(SOURCEDIR)/IO_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/EvolutionKernels.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/LHA_alpha_snowflake.o: $(SOURCEDIR)/LHA_alpha_snowflake.f90 $(SOURCEDIR)/IO_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/LHA_alpha_snowflake.f90 -I$(MOD)
	mv *.o $(OBJ)
	mv *.mod $(MOD)

$(OBJ)/SnowFlake_Model.o: $(SOURCEDIR)/SnowFlake_Model.f90 $(SOURCEDIR)/LHA_alpha_snowflake.f90 $(CommonFiles)
	$(FC) -c $(SOURCEDIR)/SnowFlake_Model.f90 -I$(MOD)
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
	$(FC) $(HOME)/prog/TEST.f90 $(OBJFILES) -I$(MOD)
	./a.out

program: 
	echo $(TARGET)
	$(FC) $(TARGET) $(OBJFILES) -I$(MOD)

################################################  Python PART  #######################################

pysnow-signature:
	f2py -h $(PYDIR)/pysnowflake.pyf --overwrite-signature $(PYDIR)/PySnowflake.f90
	sed -i '3i\\' $(PYDIR)/pysnowflake.pyf
	sed -i '3i interface' $(PYDIR)/pysnowflake.pyf
	sed -i '3i python module pysnowflake' $(PYDIR)/pysnowflake.pyf
	sed -i '3i\\' $(PYDIR)/pysnowflake.pyf
	echo 'end interface' >> $(PYDIR)/pysnowflake.pyf
	echo 'end python module pysnowflake' >> $(PYDIR)/pysnowflake.pyf

pysnow:
	f2py -c --f90exec=$(Fpath) --f77exec=$(F77path) --f90flags=$(FflagsPY) -lgomp -I$(MOD) $(SOURCEFILES) $(PYDIR)/PySnowflake.f90 $(PYDIR)/pysnowflake.pyf
	mv pysnowflake*.so $(PYDIR)
