# Makefile for GreenK
#
# Author: Kyungjin Ahn

# For compilation of GreenK, which precalculates and outputs k-space
# green function for LW intensity. Files will be used later in C2Ray
# when obtaining LW intensity by convolving with k-space source luminosity.

# Compiler
#FC = ifort # Intel compiler
FC = gfortran 
MPIFC = mpif90 # MPI compiler

# FFTW Chosun Teragon
#FFTW3 =/usr/local/fftw-3.2.2

# FFTW Ranger local
#FFTW3 =/share/home/01113/tg803913/fftw_3.2.2

# FFTW Chosun Teragon
#FFTW3 =/home/kjahn/fftw331
# This is for Marvin (ITIs laptop)
FFTW3=/p/software/default/stages/2025/software/FFTW/3.3.10-GCCcore-13.3.0

#FFTW3_INC =$(FFTW3)/include
FFTW3_INC =$(FFTW3)/include
FFTW3_LIB =$(FFTW3)/lib

# place FFTW3_FLAG after .o files when linking fails.
FFTW3_FLAG = -I$(FFTW3_INC) -L$(FFTW3_LIB) -lfftw3

# F90 options (ifort)
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT -shared-intel #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpp -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -DIFORT -check all -traceback
#IFORTFLAGS = -O2 -fpp  # Chosun Teragon

# F90 options (gfortran)
GFORTFLAGS = -O2  -ffree-line-length-none
# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS) 
#F90FLAGS1 = -xW $(IFORTFLAGS)
F90FLAGS1 = $(GFORTFLAGS)

#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#MPI_FLAGS = -DMPI -I/engrid/enhpc/mpich/intel/include

#OPENMP_FLAGS = -openmp -D_OPENMP# For Intel compiler


#-------------------------------------------------------

# Other F90 compilers
#FC = f90
F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)                 # MPI only

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS = 

#-------------------------------------------------------

UTILS=

CONSTANTS = mathconstants.o cgsconstants.o  cgsastroconstants.o cosmoparms.o


greenK: F90=$(MPIFC)
greenK: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS) $(FFTW3_FLAG)
greenK: precision.o $(CONSTANTS) $(UTILS) sizes.o file_admin.o mpi.o para_range.o readinput.o jLWgreen.o GreenK.o
	$(F90) precision.o $(UTILS) sizes.o file_admin.o mpi.o para_range.o readinput.o jLWgreen.o GreenK.o $(OPTIONS) -o $@ 


clean : 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o


