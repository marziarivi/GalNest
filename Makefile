###############################################################
#                                                             #
###############################################################

#------------------------------------------------------------------------- 
OMP  = -fopenmp  # OpenMP compiler switch

OPT += -DUSE_MPI

#------------------------------------------------------------------------- 

CASACORE_INC =  -I/share/splinter/cosmos/modules/sep_2019/install_dir/casacore-5.4.0/include 
CASACORE_LIB =  -L/share/splinter/cosmos/modules/sep_2019/install_dir/casacore-5.4.0/lib -lcasa_casa -lcasa_measures -lcasa_tables -lcasa_ms -std=c++11

MULTINEST_INC = -I/share/splinter/cosmos/modules/oct_2018/install_dir/MultiNest_v3.11_CMake/multinest/include
MULTINEST_LIB = -L/share/splinter/cosmos/modules/oct_2018/install_dir/MultiNest_v3.11_CMake/multinest/lib -lmultinest

//GSL_INC = -I/iranet/arcsoft/gsl/gsl-2.5/include
//GSL_LIB = -L/iranet/arcsoft/gsl/gsl-2.5/lib

SUP_INCL = -I. $(CASACORE_INC) $(MULTINEST_INC) $(GSL_INC)
OPTIMIZE = -O3 -g 

ifeq (USE_MPI,$(findstring USE_MPI,$(OPT)))
  CC  =  mpiCC -g
  EXEC = GalNest-mpi 
else
  CC  = g++
  EXEC = GalNest
endif

OPTIONS = $(OPTIMIZE) $(OPT)
OBJS  = Galnest.o utils.o distributions.o generate_random_values.o likelihood.o measurement_set.o 
 
 
INCL   = *.h Makefile
LIB_OPT =  -lgsl -lm $(CASACORE_LIB) $(MULTINEST_LIB)

CPPFLAGS = $(OPTIONS) $(SUP_INCL)  $(OMP)

LIBS   = $(LIB_OPT) $(OMP)

.SUFFIXES: .o .cc .cxx .cpp .cu

.cc.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cxx.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

.cpp.o:
	$(CC) -c $(CPPFLAGS) -o "$@" "$<"

GalNest: $(OBJS)
	$(CC)  $(OBJS)  $(OPTIONS) $(LIBS) -o $(EXEC)

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) 

realclean: clean
	rm -f $(EXEC)  

