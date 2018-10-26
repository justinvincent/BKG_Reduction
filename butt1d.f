C$Header: /data/petsun4/data1/src_solaris/interp_4dfp/RCS/butt1d.f,v 1.3 2004/05/26 20:30:53 avi Exp avi $
C$Log: butt1d.f,v $
c Revision 1.3  2004/05/26  20:30:53  avi
c subroutine butt1db
c
c Revision 1.2  2002/06/27  05:12:21  avi
c correct code to compute factor
c
c Revision 1.1  2002/06/25  05:16:29  avi
c Initial revision
c
      subroutine butt1d_rcs
      write (*,"('$Id: butt1d.f,v 1.3 2004/05/26 20:30:53 avi Exp avi $')")
      return
      end

      subroutine butt1dh(data,n,delta,fhalf,iorder)
      real*4 data(n)
      parameter (nmax=1024)
      real*4 a(nmax/2+1),b(nmax/2+1)

      if(n.gt.nmax)then
        write(*,"('butt1dh: input array length',i,' exceeds',i)")n,nmax
        call exit(-1)
      endif
      if(mod(n,2))then
        write(*,"('butt1dh: illegal odd input array length',i)")n
        call exit(-1)
      endif

      i=1
      do 21 k=1,n,2
      a(i)=data(k)
      b(i)=data(k+1)
   21 i=i+1

      call FFT  (a,b,1,n/2,1,-1)
      call REALT(a,b,1,n/2,1,-1)
      do 31 i=1,n/2+1
      f=float(i-1)/(float(n)*delta)
      r=(f/fhalf)**(2*iorder)
      factor=r/(1.0+r)
      a(i)=factor*a(i)
   31 b(i)=factor*b(i)
      call REALT(a,b,1,n/2,1,+1)
      call FFT  (a,b,1,n/2,1,+1)

      i=1
      do 41 k=1,n,2
      data(k)  =a(i)
      data(k+1)=b(i)
   41 i=i+1

      return
      end

      subroutine butt1db(data,n,delta,fhalf_lo,iorder_lo,fhalf_hi,iorder_hi)
      real*4 data(n)
      parameter (nmax=1024)
      real*4 a(nmax/2+1),b(nmax/2+1)

      if(n.gt.nmax)then
        write(*,"('butt1db: input array length',i,' exceeds',i)")n,nmax
        call exit(-1)
      endif
      if(mod(n,2))then
        write(*,"('butt1db: illegal odd input array length',i)")n
        call exit(-1)
      endif
      if(iorder_lo.lt.0.or.iorder_hi.lt.0)then
        write(*,"('butt1db: negative Butterworth filter orders not allowed')")
        call exit(-1)
      endif

      i=1
      do 21 k=1,n,2
      a(i)=data(k)
      b(i)=data(k+1)
   21 i=i+1

      call FFT  (a,b,1,n/2,1,-1)
      call REALT(a,b,1,n/2,1,-1)
      do 31 i=1,n/2+1
      f=float(i-1)/(float(n)*delta)
      if(iorder_lo.gt.0)then
        r_lo=(f/fhalf_lo)**(2*iorder_lo)
        factor_lo=r_lo/(1.0+r_lo)
      else
        factor_lo=1.0
      endif
      if(iorder_hi.gt.0)then
        r_hi=(f/fhalf_hi)**(2*iorder_hi)
        factor_hi=1.0/(1.0+r_hi)
      else
        factor_hi=1.0
      endif
      a(i)=factor_lo*factor_hi*a(i)
   31 b(i)=factor_lo*factor_hi*b(i)
      call REALT(a,b,1,n/2,1,+1)
      call FFT  (a,b,1,n/2,1,+1)

      i=1
      do 41 k=1,n,2
      data(k)  =a(i)
      data(k+1)=b(i)
   41 i=i+1

      return
      end

      subroutine butt1dba(data,n,delta,fhalf_lo,iorder_lo,fhalf_hi,iorder_hi)
c     version of butt1db that allocates buffers on each call
      real*4 data(n)
      real*4 a(1),b(1)
      pointer (pa,a),(pb,b)

      if(mod(n,2))then
        write(*,"('butt1dba: illegal odd input array length',i)")n
        call exit(-1)
      endif
      if(iorder_lo.lt.0.or.iorder_hi.lt.0)then
        write(*,"('butt1dba: negative Butterworth filter orders not allowed')")
        call exit(-1)
      endif

      pa=malloc(4*(n/2+1))
      pb=malloc(4*(n/2+1))
      if(pa.eq.0.or.pb.eq.0)then
        write(*,"('butt1dba: memory allocation error')")
        call exit(-1)
      endif

      i=1
      do 21 k=1,n,2
      a(i)=data(k)
      b(i)=data(k+1)
   21 i=i+1

      call FFT  (a,b,1,n/2,1,-1)
      call REALT(a,b,1,n/2,1,-1)
      do 31 i=1,n/2+1
      f=float(i-1)/(float(n)*delta)
      if(iorder_lo.gt.0)then
        r_lo=(f/fhalf_lo)**(2*iorder_lo)
        factor_lo=r_lo/(1.0+r_lo)
      else
        factor_lo=1.0
      endif
      if(iorder_hi.gt.0)then
        r_hi=(f/fhalf_hi)**(2*iorder_hi)
        factor_hi=1.0/(1.0+r_hi)
      else
        factor_hi=1.0
      endif
      a(i)=factor_lo*factor_hi*a(i)
   31 b(i)=factor_lo*factor_hi*b(i)
      call REALT(a,b,1,n/2,1,+1)
      call FFT  (a,b,1,n/2,1,+1)

      i=1
      do 41 k=1,n,2
      data(k)  =a(i)
      data(k+1)=b(i)
   41 i=i+1

      call free(pa)
      call free(pb)
      return
      end
