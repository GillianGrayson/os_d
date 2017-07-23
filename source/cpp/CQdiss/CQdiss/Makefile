topdir = ..
NAME_APP = cdiss
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = CalcEigs.o      CalcRho.o       f_d_sym.o  \
CalcGs.o        genMatrix.o   Model.o \
CalcKs.o        Init_a1_a2.o  sortTensor.o \
calcODE.o       CalcRs.o        Init_f_d.o    stdToCrs.o \
calcODE_real.o  CalcTraseRO2.o  InitFs.o      testODE.o  \
CalcODE_rk.o    coef_coord.o    InitH.o       transpFs.o \
CalcQs.o        f_d_ch.o        linSolv.o


default: build

build: $(NAME_APP)

$(NAME_APP): $(OBJS)
		$(AR) $(prefix)lib$@.a $^

.cpp.o:
	$(CC) -I../SparceM_op/ -I../read_config/ $(OPTC) -c $*.cpp $(OUTC)$*.o

clean:
	$(RM) *.o 
	$(RM) *.a
	$(RM) $(prefix)$(NAME_APP)