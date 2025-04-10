      subroutine rampfit(n,y,n1,n2,i1,i2,y0,smin,swap)
      implicit none
      integer*4 n,n1,n2,i1,i2
      real*8 y0,smin
      real*8 y(n),swap(n,0:1)
c
      integer*4 i,j1,j2,id,j1min,j1max,j2min,j2max
      real*8 sigma,delta
c
      swap(n,0)=y(n)
      do i=n-1,n1,-1
        swap(i,0)=swap(i+1,0)+y(i)
      enddo
c
      swap(1,1)=y(1)
      do i=2,n2
        swap(i,1)=swap(i-1,1)+dble(i)*y(i)
      enddo
c
      i2=n1
      i1=i2-1
      y0=swap(i2,0)/dble(1+n-i2)
      smin=-y0**2*dble(1+n-i2)
c
      j1min=n1
      j1max=n2
      j2min=n1
      j2max=n2
      id=1+(n2-n1)/50
c
10    continue
      do j2=j2min+1,j2max,id
        y0=swap(j2,0)/dble(1+n-j2)
        delta=-y0**2*dble(1+n-j2)
        do j1=j1min,min0(j1max,j2-1),id
          sigma=delta
     &         +(y0*dble(j2-j1-1)*dble(2*(j2-j1)-1)/6.d0
     &         +2.d0*(dble(j1)*(swap(j1,0)-swap(j2-1,0))
     &             -swap(j2-1,1)+swap(j1,1)))*y0/dble(j2-j1)
          if(smin.gt.sigma)then
            i1=j1
            i2=j2
            smin=sigma
          endif
        enddo
      enddo
      if(id.gt.1)then
        j1min=max0(n1,i1-5*id/2)
        j1max=min0(n2,i1+5*id/2)
        j2min=max0(n1,i2-5*id/2)
        j2max=min0(n2,i2+5*id/2)
        id=1+id/5
        goto 10
      endif
c
      do i=1,n
        smin=smin+y(i)**2
      enddo
c
      return
      end