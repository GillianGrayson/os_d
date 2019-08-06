topdir = ..
NAME_APP = matrix_op
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = Matrix_add.o Matrix_mult.o Matrix_op.o Matrix.o

default: build

build: $(NAME_APP)

$(NAME_APP): $(OBJS)
		$(AR) $(prefix)lib$@.a $^

.cpp.o:
	$(CC) $(OPTC) -c $*.cpp $(OUTC)$*.o

clean:
	$(RM) *.o 
	$(RM) *.a
	$(RM) $(prefix)lib$@.a