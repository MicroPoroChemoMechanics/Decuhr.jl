# Makefile — build the reference DECUHR (Espelid & Genz, 1994) Fortran code
# together with the comparison driver, producing `decuhr_cmp.exe`.
#
# Requires gfortran (e.g. from MinGW).  Usage:
#     make                       # build decuhr_cmp.exe
#     make run                   # build then run (prints the comparison table)
#     make clean                 # remove the executable
#     make FC=/path/to/gfortran  # override the compiler
#
# The reference routines are fixed-form Fortran 77 (*.f); the driver is
# free-form (*.f90).  gfortran selects the source form from the extension,
# so a single compile/link command handles both.

FC     ?= gfortran
FFLAGS ?= -std=legacy -w -O2 -static
SRCDIR := src_fortran

# 14 reference routines (the full DECUHR Fortran package).
FSRC := \
	$(SRCDIR)/d07hre.f $(SRCDIR)/d09hre.f $(SRCDIR)/d113re.f $(SRCDIR)/d132re.f \
	$(SRCDIR)/deadhr.f $(SRCDIR)/decalp.f $(SRCDIR)/dechhr.f $(SRCDIR)/decuhr.f \
	$(SRCDIR)/defshr.f $(SRCDIR)/deinhr.f $(SRCDIR)/derlhr.f $(SRCDIR)/desbhr.f \
	$(SRCDIR)/detrhr.f $(SRCDIR)/dexthr.f

DRIVER := driver.f90
EXE    := decuhr_cmp.exe

.PHONY: all run clean

all: $(EXE)

$(EXE): $(FSRC) $(DRIVER)
	$(FC) $(FFLAGS) $(FSRC) $(DRIVER) -o $(EXE)

run: $(EXE)
	./$(EXE)

clean:
	-rm -f $(EXE)
