#$Header: /data/petsun4/data1/src_solaris/eeg/RCS/EKG_phase.mak,v 1.2 2006/12/16 23:06:38 avi Exp $
#$Log: EKG_phase.mak,v $
# Revision 1.2  2006/12/16  23:06:38  avi
# Solaris 10
#
# Revision 1.1  2006/04/16  22:36:01  avi
# Initial revision
#

PROG	= EKG_phase
CSRCS	= ${PROG}.c
FSRCS	= butt1d.f eigen.f fftsun.f matopr.f
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
