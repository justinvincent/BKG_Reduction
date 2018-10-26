PROG	= NSepoch2power
CSRCS	= ${PROG}.c
FSRCS	= fftsun.f
LOBJS	= 
OBJS	= ${CSRCS:.c=.o} ${FSRCS:.f=.o}

CFLAGS	= -O
CC	= cc ${CFLAGS}
FC	= f77 -e -I4 -O

.c.o:
	${CC} -c $<
.f.o:
	${FC} -c $<

${PROG}: ${OBJS} 
	${FC} -o $@ ${OBJS} -lm

clean:
	/bin/rm ${OBJS} ${PROG}
