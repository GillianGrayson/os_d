topdir = ..
NAME_APP = read_config
prefix   = ../bin/

include $(topdir)/Makefile.inc

OBJS = read_config.o


default: build

build: $(NAME_APP)

$(NAME_APP): $(OBJS)
		$(AR) $(prefix)lib$@.a $^

.cpp.o:
	$(CC) $(OPTC) -c $*.cpp $(OUTC) $*.o

clean:
	$(RM) *.o 
	$(RM) *.a
	$(RM) $(prefix)lib$@.a