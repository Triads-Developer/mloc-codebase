OBJFILES=	mloc_allocate.o \
			mloc_calibration.o \
			mloc_commands.o \
			mloc_declare.o \
			mloc_hyposat.o \
			mloc_init.o \
			mloc_inv.o \
			mloc_lib.o \
			mloc_math.o \
			mloc_mnf.o \
			mloc_out.o \
			mloc_phases.o \
			mloc_plots.o \
			mloc_set.o \
			mloc_stations.o \
			mloc_taup.o \
			mloc_tt.o \
			mloc.o

# Link LAPACK and BLAS libraries during final compilation
mloc_g1103p: $(OBJFILES)
	gfortran $(OBJFILES) -o mloc_g1103p -llapack -lblas

# Flags for development and debugging
#FFLAGS=	-c -W -Wall -Wno-tabs -Wno-unused-parameter -fno-automatic -fstack-check -g \
#		-fbacktrace -ffpe-trap=zero,overflow,underflow -fno-align-commons \
#		-fcheck=all
# -c (compile but don't link)
# -W (warn more verbosely than usual)
# -Wall (enable commonly used warnings)
# -g (include debug data)
# -fbacktrace (only for main program)
# -ffpe-trap=list (Specify a list of IEEE exceptions when a Floating Point Exception (FPE) should be raised.)
#FFLAGS= -c -W -Wall -Wextra -Wpedantic

# Flags for development and debugging (including profiling with gprof)
FFLAGS=	-c -O2 -pg

# Location of source files
MLOCSRCDIR=/Users/jjwlove/Downloads/mloc_src_final_prototype-09-09-2024
SRCDIR=$(MLOCSRCDIR)

mloc_declare.o: $(SRCDIR)/mloc_declare.f90
	gfortran $(FFLAGS) $(SRCDIR)/mloc_declare.f90

mloc_init.o: $(SRCDIR)/mloc_init.f90
	gfortran $(FFLAGS) $(SRCDIR)/mloc_init.f90

mloc_taup.o: $(SRCDIR)/mloc_taup.f90 mloc_declare.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_taup.f90

mloc_allocate.o: $(SRCDIR)/mloc_allocate.f90 mloc_declare.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_allocate.f90

mloc_math.o: $(SRCDIR)/mloc_math2.f90 mloc_declare.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_math2.f90 -o mloc_math.o

mloc_lib.o: $(SRCDIR)/mloc_lib.f90 mloc_declare.o mloc_math.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_lib.f90

mloc_inv.o: $(SRCDIR)/mloc_inv2.f90  mloc_declare.o mloc_allocate.o mloc_math.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_inv2.f90 -o mloc_inv.o

mloc_commands.o: $(SRCDIR)/mloc_commands.f90 mloc_declare.o mloc_math.o mloc_lib.o \
		mloc_hyposat.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_commands.f90

mloc_hyposat.o: $(SRCDIR)/mloc_hyposat.f90 mloc_declare.o mloc_math.o mloc_lib.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_hyposat.f90

mloc_calibration.o: $(SRCDIR)/mloc_calibration.f90 mloc_declare.o mloc_math.o mloc_lib.o \
		mloc_allocate.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_calibration.f90

mloc_tt.o: $(SRCDIR)/mloc_tt.f90 mloc_declare.o mloc_lib.o mloc_taup.o \
		mloc_hyposat.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_tt.f90

mloc_phases.o: $(SRCDIR)/mloc_phases.f90 mloc_declare.o mloc_taup.o mloc_lib.o \
		mloc_math.o mloc_tt.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_phases.f90

mloc_plots.o: $(SRCDIR)/mloc_plots.f90 mloc_declare.o mloc_lib.o mloc_taup.o mloc_tt.o \
		mloc_hyposat.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_plots.f90

mloc_stations.o: $(SRCDIR)/mloc_stations.f90 mloc_declare.o mloc_math.o mloc_lib.o \
		mloc_lib.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_stations.f90

mloc_mnf.o: $(SRCDIR)/mloc_mnf.f90 mloc_declare.o mloc_lib.o mloc_math.o mloc_taup.o \
		mloc_phases.o mloc_stations.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_mnf.f90

mloc_out.o: $(SRCDIR)/mloc_out.f90 mloc_declare.o mloc_lib.o mloc_math.o mloc_taup.o \
		mloc_plots.o mloc_tt.o mloc_hyposat.o mloc_mnf.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_out.f90

mloc_set.o: $(SRCDIR)/mloc_set.f90 mloc_declare.o mloc_math.o mloc_plots.o \
		mloc_calibration.o mloc_taup.o mloc_tt.o mloc_out.o mloc_inv.o mloc_lib.o \
		mloc_phases.o mloc_mnf.o mloc_stations.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc_set.f90

mloc.o: $(SRCDIR)/mloc.f90 mloc_declare.o mloc_allocate.o mloc_calibration.o mloc_lib.o \
		mloc_set.o mloc_commands.o mloc_stations.o
	gfortran $(FFLAGS) $(SRCDIR)/mloc.f90

clean:
	rm -f mloc_g1103p *.o *.mod

