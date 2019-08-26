CC  = icc -openmp
CXX = icpc -openmp
FF  = ifort
AR  = ar rc 
RANLIB = ranlib
RM  = rm -f
CP  = cp
MV  = mv
LIBEXT  = .a
LIBPREF = lib

OUTC    = -o
OUTF    = -o 

CFLAGS   = -O3 -mkl
OPTF = $(CFLAGS)
OPTC = $(CFLAGS)
LD_FLAGS = 
