cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005
c washington university, mallinckrodt institute of radiology.
c all rights reserved.
c this software may not be reproduced, copied, or distributed without writte
c permission of washington university. for further information contact a. z. snyder.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $Header: /data/petsun4/data1/src_solaris/librms/RCS/deigen.f,v 1.1 2005/09/04 05:41:40 avi Exp $
c $Log: deigen.f,v $
c Revision 1.1  2005/09/04  05:41:40  avi
c Initial revision
c
c     program deigentst
      subroutine deigentst
      implicit real*8 (a-h,o-z)
      parameter (n=256)
      real*8 q(n,n),w(n,n),t(n,n),r(n,n)
      do 1 i=1,n
      do 1 j=1,n
    1 t(i,j)=drand(0)
      do 2 i=1,n
      do 2 j=1,n
      q(i,j)=0.0
      do 2 k=1,n
      q(i,j)=q(i,j)+t(i,k)*t(j,k)
    2 r(i,j)=q(i,j)
      call deigen(q,w,n)
      if (n.eq.3) then
        do i=1,n
          write(*,101)(r(i,k),k=1,n),(q(i,k),k=1,n),(w(i,k),k=1,n)
        enddo
      endif
  101 format(3(3f16.12,4x))
      do 12 i=1,n
      do 12 j=1,n
      t(i,j)=0.0
      do 12 k=1,n
   12 t(i,j)=t(i,j)+w(i,k)*q(k,j)
      do 13 i=1,n
      do 13 j=1,n
      q(i,j)=0.0
      do 13 k=1,n
   13 q(i,j)=q(i,j)+t(i,k)*w(j,k)
      write(*,"()")
      if (n.eq.3)then
        do i=1,n
          write(*,101)(q(i,k),k=1,n),(q(i,k)-r(i,k),k=1,n)
         enddo
      endif
      derr=0.0
      ierr=0
      jerr=0
      do 6 i=1,n
      do 6 j=1,n
      if (dabs(q(i,j)-r(i,j)).gt.derr)then
        ierr=i
        jerr=j
        derr=dabs(q(i,j)-r(i,j))
      endif
    6 continue
      write(*,"('maximum error = ',e12.6)")derr
      end
      subroutine deigen(a,p,n)
      implicit real*8 (a-h,o-z)
      real*8 a(n,n),p(n,n)
      logical*4 lind
      data range/1.e-21/
      do 40 l=1,n
   40 p(l,l)=1.0
      do 41 l=1,n-1
      do 41 m=l+1,n
      p(l,m)=0.
   41 p(m,l)=0.
      thr=0.
      do 30 l=1,n-1
      do 30 m=l+1,n
   30 thr=thr+a(l,m)**2
      thr=dsqrt(2.*thr)
   45 thr=thr/dfloat(n)
   46 lind=.false.
      do 50 l=1,n-1
      do 55 m=l+1,n
      if(dabs(a(l,m)).le.thr)goto 55
      lind=.true.
      alamda=-a(l,m)
      amu=.5*(a(l,l)-a(m,m))
      omega=dsign(1.,amu)*alamda/dsqrt(alamda**2+amu**2)
      st=omega/dsqrt(2.*(1.+dsqrt(1.-omega**2)))
      st2=st**2
      ct2=1.-st2
      ct=dsqrt(ct2)
      stct=st*ct
      do 125 i=1,n
      if(i.eq.l.or.i.eq.m)goto 115
      x=ct*a(i,l)-st*a(i,m)
      a(i,m)=st*a(i,l)+ct*a(i,m)
      a(m,i)=a(i,m)
      a(i,l)=x
      a(l,i)=x
  115 x=ct*p(i,l)-st*p(i,m)
      p(i,m)=st*p(i,l)+ct*p(i,m)
      p(i,l)=x
  125 continue
      x=2.*a(l,m)*stct
      y=a(l,l)*ct2+a(m,m)*st2-x
      x=a(l,l)*st2+a(m,m)*ct2+x
      a(l,m)=(a(l,l)-a(m,m))*stct+a(l,m)*(ct2-st2)
      a(m,l)=a(l,m)
      a(l,l)=y
      a(m,m)=x
   55 continue
   50 continue
      if(lind)goto 46
      diag=dabs(a(n,n))
      offd=0.
      do 48 l=1,n-1
      diag=diag+dabs(a(l,l))
      do 48 m=l+1,n
   48 offd=offd+dabs(a(l,m))
      if(offd/(.5*dfloat(n-1)).gt.range*diag)goto 45
      do 185 l=1,n-1
      do 185 m=l+1,n
      if(a(l,l).ge.a(m,m))goto 185
      x=a(l,l)
      a(l,l)=a(m,m)
      a(m,m)=x
      do 180 i=1,n
      x=p(i,l)
      p(i,l)=p(i,m)
  180 p(i,m)=x
  185 continue
      return
      end
