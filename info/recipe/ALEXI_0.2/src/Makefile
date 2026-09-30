#____________________________________________________________________
#
#  Makefile for MSG ALEXI flux project. MA 03/23/2012
#____________________________________________________________________

#
# Define path macros
SRC     = .
SRCD    = .
ADV     = .
OBJ     = ./o
SRCU    = .
BIN     = ../bin


# Uses fixed leaf abs and soil refl fields
ALBD_OBJ  = $(OBJ)/USflux.o $(OBJ)/USflux_utl.o $(OBJ)/USflux_clear.o \
            $(OBJ)/USflux_cloud.o $(OBJ)/USflux_rad.o \
            $(OBJ)/USflux_cover.o $(OBJ)/USflux_run.o $(OBJ)/ALEXI.o \
            $(OBJ)/ALEXI_atmos.o $(OBJ)/ALEXI_rad.o \
            $(OBJ)/ALEXI_utl.o $(OBJ)/landcover.o \
            $(OBJ)/sfc_read.o $(OBJ)/USflux_main.o \
            $(OBJ)/pbl_read.o \
            $(OBJ)/ALEXI_water.o

# CDFLIB = -L/usr/local/lib
# CDFINC = /usr/local/include

INCLUS  = $(SRC)/date.inc $(SRCU)/ALEXIdata.inc
#FC     =  ifort  # LINUX
FC = gfortran # MAC
#FLAGS  = -c -g -qextname -qmaxmem=32768 -lX11 -lm -I$(SRC) -I$(SRCU)
#FLAGS  = -c -O2 -L/usr/X11R6/lib -lX11 -lm -I$(SRCU) -I$(SRC)
#FFLAGS  = -c -L/usr/X11R6/lib -lX11 -lm -I$(SRCU) -I$(SRC) -I$(ADV)
#FLAGS  = -c -g -I$(SRCU) -I$(SRC) -traceback -extend-source 132 -heap-arrays # LINUX
FLAGS  = -c -g -I$(SRCU) -I$(SRC) -fPIC -ffixed-line-length-132 # MAC
# FFLAGS  = -c -I$(SRCU) -I$(SRC) -I$(ADV)
#
# The programs:
# -------------------------------------------------------------
alexi_proc : $(ALBD_OBJ)
#       $(FC) -o $(BIN)/alexi_proc $(FLUX_OBJ)
#	$(FC) -g -O3 -traceback -assume byterecl -traceback -mcmodel=large -shared-intel -heap-arrays -o $(BIN)/alexi_proc $(ALBD_OBJ) # linux
	$(FC) -g -O3 -mcmodel=medium -o $(BIN)/alexi_proc $(ALBD_OBJ) # MAC

all: alexi_proc landcover interp_LAI
#
# The code:
# -------------------------------------------------------------
$(OBJ)/landcover.o: $(SRCU)/landcover.f
	$(FC) -o $(OBJ)/landcover.o $(FLAGS) $(SRCU)/landcover.f

$(OBJ)/USflux.o: $(SRCU)/USflux.f90
	$(FC) -o $(OBJ)/USflux.o $(FLAGS) $(SRCU)/USflux.f90
	
$(OBJ)/USflux_utl.o: $(SRCU)/USflux_utl.f
	$(FC) -o $(OBJ)/USflux_utl.o $(FLAGS) $(SRCU)/USflux_utl.f

$(OBJ)/USflux_rad.o: $(SRCU)/USflux_rad.f
	$(FC) -o $(OBJ)/USflux_rad.o $(FLAGS) $(SRCU)/USflux_rad.f

$(OBJ)/USflux_cover.o: $(SRCU)/USflux_cover.f
	$(FC) -o $(OBJ)/USflux_cover.o $(FLAGS) $(SRCU)/USflux_cover.f

$(OBJ)/USflux_run.o: $(SRCU)/USflux_run.f
	$(FC) -o $(OBJ)/USflux_run.o $(FLAGS) $(SRCU)/USflux_run.f

$(OBJ)/USflux_clear.o: $(SRCU)/USflux_clear.f
	$(FC) -o $(OBJ)/USflux_clear.o $(FLAGS) $(SRCU)/USflux_clear.f

$(OBJ)/USflux_cloud.o: $(SRCU)/USflux_cloud.f
	$(FC) -o $(OBJ)/USflux_cloud.o $(FLAGS) $(SRCU)/USflux_cloud.f
	
$(OBJ)/ALEXI.o: $(SRCU)/ALEXI.f
	$(FC) -o $(OBJ)/ALEXI.o $(FLAGS) $(SRCU)/ALEXI.f
	
$(OBJ)/ALEXI_atmos.o: $(SRCU)/ALEXI_atmos.f
	$(FC) -o $(OBJ)/ALEXI_atmos.o $(FLAGS) $(SRCU)/ALEXI_atmos.f
	
$(OBJ)/ALEXI_rad.o: $(SRCU)/ALEXI_rad.f
	$(FC) -o $(OBJ)/ALEXI_rad.o $(FLAGS) $(SRCU)/ALEXI_rad.f
	
$(OBJ)/ALEXI_utl.o: $(SRCU)/ALEXI_utl.f
	$(FC) -o $(OBJ)/ALEXI_utl.o $(FLAGS) $(SRCU)/ALEXI_utl.f

$(OBJ)/ALEXI_water.o: $(SRCU)/ALEXI_water.f
	$(FC) -o $(OBJ)/ALEXI_water.o $(FLAGS) $(SRCU)/ALEXI_water.f

$(OBJ)/pbl_read.o: $(SRCU)/pbl_read.f
	$(FC) -o $(OBJ)/pbl_read.o $(FLAGS) $(SRCU)/pbl_read.f

$(OBJ)/sfc_read.o: $(SRCU)/sfc_read.f
	$(FC) -o $(OBJ)/sfc_read.o $(FLAGS) $(SRCU)/sfc_read.f

$(OBJ)/USflux_main.o: $(SRCU)/USflux_main.f
	$(FC) -o $(OBJ)/USflux_main.o $(FLAGS) $(SRCU)/USflux_main.f

#
# Utilities:
# -------------------------------------------------------------
interp_LAI : $(OBJ)/interp_LAI.o
	$(FC) -o $(BIN)/interp_LAI $(OBJ)/interp_LAI.o

$(OBJ)/interp_LAI.o: $(SRCU)/interp_LAI.f
	$(FC) -o $(OBJ)/interp_LAI.o $(FLAGS) $(SRCU)/interp_LAI.f

#landcover:	$(OBJ)/landcover.o
#	$(FC) -o $(BIN)/landcover $(OBJ)/landcover.o

#$(OBJ)/landcover.o:	$(SRCU)/landcover.f
#	$(FC) $(FFLAGS) -g $(SRCU)/landcover.f -o $(OBJ)/landcover.o
		
#
#  Remove all object files to force remake.
clean:
	rm -f $(OBJ)/*.o


