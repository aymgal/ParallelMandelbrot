# Makefile for Mandelbrot set generator program

########## OPTIONS ##########

OPTS += -DOUTPUT_IMAGE
OPTS += -DOUTPUT_TIMINGS
# OPTS += -DUSE_STDCOMPLEX
# OPTS += -DUSE_DISTANCE_ESTIM
# OPTS += -DCOLOR_ITERATIONS
# OPTS += -DDOUBLE_PRECISION # do not use with MPI (for now)
# OPTS += -DSQUARE_TRICK


#-- only one of these two options at a time !
# OPTS += -DMPI_SIMPLE
OPTS += -DMPI_MASTER_WORKERS


#-- for debugging --#
OPTS += -DCOLOR_PRANK
# OPTS += -DVERBOSE

# DEBUG = -g

# PROFIL = -pg

#############################

# define any specific target machine
SYSTYPE = "mac_aym"
# SYSTYPE = "phpc_cluster"

ifeq ($(SYSTYPE), "mac_aym")
CXX = g++-7
endif
ifeq ($(SYSTYPE), "phpc_cluster")
CXX = g++
endif

LD = $(CXX)
OPTIMIZE = -Wall -Wextra -Werror -pedantic -std=c++11 -O3

CXXFLAGS = $(OPTIMIZE) $(DEBUG) $(PROFIL) $(OPTS)
LDFLAGS  = -lm $(CXXFLAGS)

OBJS     = main.o mandelbrot.o buffer.o grid.o dumpers.o \
		   utils.o timers.o
OBJS_OMP = main_omp.o mandelbrot_omp.o buffer.o grid.o dumpers.o \
		   utils.o timers.o
OBJS_MPI = main_mpi.o mandelbrot_mpi.o buffer_mpi.o grid_mpi.o dumpers_mpi.o \
		   utils_mpi.o timers_mpi.o

EXEC_SERIAL = mandel_srl
EXEC_OMP    = mandel_omp
EXEC_MPI    = mandel_mpi
EXEC_HYBRID = mandel_hyb

all: clean $(EXEC_SERIAL) $(EXEC_OMP) $(EXEC_MPI) $(EXEC_HYBRID)

srl: clean $(EXEC_SERIAL)
omp: clean $(EXEC_OMP)
mpi: clean $(EXEC_MPI)
hyb: clean $(EXEC_HYBRID)

#*** compilation of serial exec ***#
$(EXEC_SERIAL) : CXXFLAGS += -ftree-vectorize
$(EXEC_SERIAL) : $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)
#**********************************#

#*** compilation of OpenMP exec ***#
$(EXEC_OMP) : CXXFLAGS += -DPARALLEL_OPENMP -fopenmp
$(EXEC_OMP) : $(OBJS_OMP)
	$(LD) -o $@ $(OBJS_OMP) $(LDFLAGS)

main_omp.o: main.cc
	$(LD) $< $(CXXFLAGS) -c -o $@

mandelbrot_omp.o: mandelbrot.cc mandelbrot.hh
	$(LD) $< $(CXXFLAGS) -c -o $@
#**********************************#

#**** compilation of MPI exec *****#
$(EXEC_MPI) : CXXFLAGS += -DPARALLEL_MPI
$(EXEC_MPI) : CXX = mpic++
$(EXEC_MPI) : $(OBJS_MPI)
	$(LD) -o $@ $(OBJS_MPI) $(LDFLAGS)

main_mpi.o: main.cc
	$(LD) $< $(CXXFLAGS) -c -o $@

mandelbrot_mpi.o: mandelbrot.cc mandelbrot.hh
	$(LD) $< $(CXXFLAGS) -c -o $@

buffer_mpi.o: buffer.cc buffer.hh
	$(LD) $< $(CXXFLAGS) -c -o $@

grid_mpi.o: grid.cc grid.hh
	$(LD) $< $(CXXFLAGS) -c -o $@ 

dumpers_mpi.o: dumpers.cc dumpers.hh
	$(LD) $< $(CXXFLAGS) -c -o $@

utils_mpi.o: utils.cc utils.hh
	$(LD) $< $(CXXFLAGS) -c -o $@

timers_mpi.o: timers.cc timers.hh
	$(LD) $< $(CXXFLAGS) -c -o $@
#**********************************#

#**** compilation of hydrid OpenMP+MPI exec *****#
$(EXEC_HYBRID) : CXXFLAGS += -DPARALLEL_MPI -DPARALLEL_OPENMP -fopenmp
$(EXEC_HYBRID) : CXX = mpic++
$(EXEC_HYBRID) : $(OBJS_MPI)
	$(LD) -o $@ $(OBJS_MPI) $(LDFLAGS)
#**********************************#

# utility targets
clean:
	rm -rf $(EXECS) *.o *.dSYM *~

cleano:
	rm -rf *.o *~

cleane:
	rm -rf $(EXEC_SERIAL) $(EXEC_OMP) $(EXEC_MPI)

cleani:
	rm -rf *.pgm *.bmp

cleanall: clean cleano cleane cleani

