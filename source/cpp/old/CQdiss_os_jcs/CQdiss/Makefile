topdir = ..
NAME_APP = cdiss
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = CalcGs.o CalcKs.o calcODE.o CalcQs.o CalcRho.o \
CalcRs.o characteristics.o coef_coord.o data.o f_basis_init.o \
Init_a1_a2.o Init_f_d.o InitFs.o InitFs.o InitH.o \
Model.o stdToCrs.o transpFs.o utils.o main.o

default: build

build: $(NAME_APP)

$(NAME_APP): $(OBJS)
		$(CC) -I../SparceM_op/ -I../read_config/ \
		$(OUTC)$(prefix)$@ $(OPTC) $^ $(LIB_LINK) \
		-L$(prefix) -lread_config -lmatrix_op

.cpp.o:
	$(CC) -I../SparceM_op/ -I../read_config/ $(OPTC) -c $*.cpp $(OUTC) $*.o

clean:
	$(RM) *.o 
	$(RM) *.a
	$(RM) $(prefix)$(NAME_APP)