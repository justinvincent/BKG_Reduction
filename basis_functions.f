      subroutine basis_function_test
      parameter (npts=500)
      parameter (nbasis=30)
      real*4 b(0:npts-1,0:nbasis),a(0:nbasis,0:nbasis)

      call sincbas(b,npts,nbasis)

      do 11 i=0,nbasis
      do 11 j=0,nbasis
      a(i,j)=0.
      do 12 k=0,npts-1
   12 a(i,j)=a(i,j)+b(k,i)*b(k,j)
   11 a(i,j)=a(i,j)/float(npts)

      if(0)then
      do 31 i=0,nbasis
   31 write(*,"(31f4.1)")(a(i,j),j=0,nbasis)
      else
      do 52 i=0,npts-1
   52 write(*,"(i5,10f10.4)")i,(b(i,j),j=15,24)
      endif
      end

      subroutine sinbas(b,npts,nbasis)
      real*4 b(0:npts-1,0:nbasis)

      if(mod(nbasis,2).ne.0)then
        write(*,"('sinbas: nbasis argument must be even')")
        call exit(-1)
      endif

      twopi=2.*atan2(0.,-1.)
      root2=sqrt(2.)
      do 2 i=0,npts-1
      t=twopi*float(i)/float(npts)
      ibasis=0
      b(i,ibasis)=1.
      ibasis=ibasis+1
      do 1 j=1,nbasis/2
      b(i,ibasis)=root2*cos(t*float(j))
      ibasis=ibasis+1
      b(i,ibasis)=root2*sin(t*float(j))
    1 ibasis=ibasis+1
    2 continue
      return
      end

      subroutine sincbas(b,npts,nbasis)
      parameter (alpha=0.2)
      real*4 b(0:npts-1,0:nbasis)

      real*4 f(0:npts-1,0:nbasis)
      real*4 x(0:nbasis,0:nbasis),w(0:nbasis,0:nbasis),a(0:nbasis,0:nbasis)
      real*4 r(0:nbasis,0:nbasis)
      pointer (pf,f),(px,x),(pw,w),(pa,a),(pr,r)

      n1=nbasis+1
      pf=malloc(4*n1*npts)
      px=malloc(4*n1*n1)
      pw=malloc(4*n1*n1)
      pa=malloc(4*n1*n1)
      pr=malloc(4*n1*n1)

      call sinbas(f,npts,nbasis)
      twopi=8.*atan(1.)
      do 2 i=0,npts-1
      t=twopi*float(i)/float(npts)
      e=exp(-alpha*t)
      do 2 j=0,nbasis
    2 f(i,j)=f(i,j)*e

      do 11 i=0,nbasis
      do 11 j=0,nbasis
      a(i,j)=0.
      do 12 k=0,npts-1
   12 a(i,j)=a(i,j)+f(k,i)*f(k,j)
   11 a(i,j)=a(i,j)/float(npts)

      call eigen(a,w,nbasis+1)

      do 42 i=0,nbasis
      do 43 j=0,nbasis
   43 r(i,j)=0.0
   42 r(i,i)=1./sqrt(a(i,i))
      call matmul(w,r,x,nbasis+1)

      do 51 i=0,npts-1
      do 51 j=0,nbasis
      b(i,j)=0.
      do 51 k=0,nbasis
   51 b(i,j)=b(i,j)+f(i,k)*x(k,j)

      do 53 j=0,nbasis
      imax=0
      do 54 i=0,npts-1
      if(abs(b(i,j)).gt.abs(b(imax,j)))imax=i
   54 continue
      if(b(imax,j).lt.0.)then
        do i=0,npts-1
          b(i,j)=-b(i,j)
        enddo
      endif
   53 continue

      call free(pf)
      call free(px)
      call free(pw)
      call free(pa)
      call free(pr)
      return
      end
