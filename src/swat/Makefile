# Purpose: 
#	Compiles SWAT with gfortran
# 
# Author: Andre Dozier and Olaf David
#   Date: May 10, 2013 
# 
# 
.SUFFIXES: 
.SUFFIXES: .f .f90 .f03 .o .mod

SYS := $(shell gcc -dumpmachine) 
ifneq (, $(findstring mingw, $(SYS)))
OS := win
endif

# if not release, set the debug flag
ifneq ($(rel),true)
deb := -g
endif

ifeq ($(OS),win)

# windows 
RM = del /F
SWAT.OUT = swat.exe
SWAT.DLL = libswat.dll
EVENT.DIR = ..\event\\
EVENT.DLL = libevent.dll
LIB.DIR = ..\..\lib\\
BIN.DIR = ..\..\bin\\
FFLAGS = -fno-underscoring $(deb)
CP = copy /y

else

# linux 
RM = rm -f
SWAT.OUT = swat
SWAT.DLL = libswat.so
EVENT.DIR := ../event/
LIB.DIR = ../../lib/
BIN.DIR = ../../bin/
EVENT.DLL = libevent.so
FFLAGS = -fno-underscoring -fPIC $(deb)
CP = cp -f 

endif

EVENT = $(LIB.DIR)$(EVENT.DLL)
FC = gfortran
MODFS := modevent.f modparm.f
MODOBJ := $(patsubst %.f,%.o, $(MODFS))
OBJ := $(patsubst %.f,%.o, $(wildcard *.f)) $(patsubst %.f90,%.o, $(wildcard *.f90)) $(patsubst %.f03,%.o, $(wildcard *.f03))

all: swat swatdll

swat: $(EVENT) $(MODOBJ) $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) $(EVENT) -o $(SWAT.OUT)
	$(CP) $(SWAT.OUT) $(BIN.DIR)
	$(CP) *.mod $(LIB.DIR)

swatdll: $(EVENT) $(MODOBJ) $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) $(EVENT) -shared -o $(SWAT.DLL)
	$(CP) $(SWAT.DLL) $(LIB.DIR)
	$(CP) *.mod $(LIB.DIR)

clean: 
	$(MAKE) clean -C $(EVENT.DIR) 
	$(RM) *.o *.mod *.exe *.dll *.so $(SWAT.OUT)

.f.o:
	$(FC) $(FFLAGS) -c $<
.f90.o:
	$(FC) $(FFLAGS) -c $<
.f03.o:
	$(FC) $(FFLAGS) -c $<
$(LIB.DIR)$(EVENT.DIR): 
	$(MAKE) -C $(EVENT.DIR) 

	