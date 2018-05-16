topdir = ..
NAME_APP = cdiss_os_jcs
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = main.o

default: build

build: $(NAME_APP)

$(NAME_APP): $(OBJS)
		$(CC) -I../SparceM_op/ -I../CQdiss/ -I../read_config/ \
		$(OUTC)$(prefix)$@ $(OPTC) $^ $(LIB_LINK) \
		-L$(prefix) -lcdiss -lread_config -lmatrix_op

.cpp.o:
	$(CC) -I../SparceM_op/ -I../CQdiss/ -I../read_config/ $(OPTC) -c $*.cpp $(OUTC) $*.o

clean:
	$(RM) *.o 
	$(RM) *.a
	$(RM) $(prefix)$(NAME_APP)