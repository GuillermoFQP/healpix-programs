# ======================
# Paths
# ======================
HEALPIX   := /Users/guillermo/Documents/Healpix_3.83
FITSDIR   := /opt/local/lib

INCDIR    := $(HEALPIX)/include_gfortran
LIBDIR    := $(HEALPIX)/lib_gfortran
SHARPLDIR := $(HEALPIX)/lib

LIBFITS   := cfitsio

# ======================
# Compiler
# ======================
FC := gfortran

# ======================
# Flags
# ======================
OPTFLAGS  := -O3
PRECFLAGS := -fdefault-real-8 -fdefault-double-8
OMPFLAGS  := -fopenmp
PICFLAGS  := -fPIC
WARNFLAGS := -w

FFLAGS := $(OPTFLAGS) $(PRECFLAGS) $(OMPFLAGS) $(PICFLAGS) $(WARNFLAGS) -I$(INCDIR)

# ======================
# Linker flags
# ======================
LIBS := -lhealpix -lhpxgif -lsharp -l$(LIBFITS) -lcurl -llapack -lblas

LDFLAGS := -L$(LIBDIR) -L$(FITSDIR) -L$(SHARPLDIR) $(LIBS) -Wl,-rpath,$(LIBDIR) -Wl,-rpath,$(FITSDIR) -Wl,-rpath,$(SHARPLDIR)

# ======================
# Modules and programs
# ======================
MODULES := hp_globals.f90 hp_utils.f90
PROG1   := emd
PROG2   := memd
PROG3   := minmax
PROG4   := imilc
PROG5   := mapstacker
PROG6   := mapops
PROG7   := ilc

# ======================
# Build rules
# ======================
all: $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7)

# Executables
$(PROG1): $(MODULES) $(PROG1).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

$(PROG2): $(MODULES) $(PROG2).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
	
$(PROG3): $(MODULES) $(PROG3).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
	
$(PROG4): $(MODULES) $(PROG4).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

$(PROG5): $(MODULES) $(PROG5).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
	
$(PROG6): $(MODULES) $(PROG6).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)
	
$(PROG7): $(MODULES) $(PROG7).f90
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

# ======================
# Utilities
# ======================
clean:
	@echo "Removing .o and .mod files..."
	@rm -f *.o *.mod

cleanout:
	@echo "Removing executables..."
	@rm -f $(PROG1) $(PROG2) $(PROG3) $(PROG4) $(PROG5) $(PROG6) $(PROG7)

cleanall: clean cleanout

.PHONY: all clean cleanout cleanall
