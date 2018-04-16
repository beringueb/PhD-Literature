CC = gcc
AR = ar rvs 
CCFLAG = -O2 -DCAMB -fPIC -DHYRECPATH=\"$(PWD)/\" 
LDFLAG = -O2

%.o: %.c 
	$(CC) $(CCFLAG) -c $*.c -o $*.o

HYREC_SRC = hyrectools.o helium.o hydrogen.o history.o hyrec_21cm.o
HYREC_EXE = hyrec.o 

HYREC_SRC2 = hyrectools.c helium.c hydrogen.c history.c hyrec.c hyrec_21cm.c

default: libhyrec.a

clean: 
	rm *.o

hyrec: $(HYREC_SRC2) 
	$(CC) $(LDFLAG) -o hyrec  $(HYREC_SRC2) -lm


libhyrec.a: $(HYREC_SRC)
	$(AR) $@ $(HYREC_SRC)
