topdir = ..
NAME_APP = cdiss
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = CalcEigs.o \
CalcGs.o \
CalcKs.o \
calcODE.o \
CalcQs.o \
CalcRho.o \
CalcRs.o \
coef_coord.o \
f_d_ch.o \
f_d_sym.o \
Init_a1_a2.o \
InitH.o \
linSolv.o \
Model.o \
sortTensor.o \
stdToCrs.o 



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