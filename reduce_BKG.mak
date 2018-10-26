#$Header: /data/petsun4/data1/src_solaris/eeg/RCS/reduce_BKG.mak,v 1.4 2006/12/17 01:15:48 avi Exp $
#$Log: reduce_BKG.mak,v $
# Revision 1.4  2006/12/17  01:15:48  avi
# Revision 1.3  2006/12/17  01:13:39  avi
# Solaris 10
#
# Revision 1.2  2006/03/17  04:43:12  avi
# double precision matrix operations
#
# Revision 1.1  2006/02/01  04:07:07  avi
# Initial revision
#

PROG	= reduce_BKG
CSRCS	= ${PROG}.c
FSRCS	= basis_functions.f dmatinv.f deigen.f fftsun.f matopr.f npad.f eigen.f
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
