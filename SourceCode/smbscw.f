      subroutine smbscw(ist,nwin)
      use smalloc
      implicit none
      integer*4 ist,nwin
c
      integer*4 i,j,k,l,ipre,ipga,ipgaj,isdw,iddw,ipst
      real*8 al,bl,preoff,pstoff,pga,pgaj,sigma
c
c     select pre-seismic time window
c
      ipre=1+idint((ponset(ist)-start(ist)-DTP)/dt)
      k=ipre-1-idint(6.d0*PREWIN/dt)
      if(k.gt.0)then
        ipre=ipre-k
        nwin=nwin-k
        start(ist)=start(ist)+dble(k)*dt
        length(ist)=length(ist)-dble(k)*dt
        do j=1,3
          do i=1,nwin
            acc(i,j)=acc(i+k,j)
          enddo
        enddo
      endif
c
c     remove pre-event sm offset
c
      do j=1,3
        preoff=0.d0
        do i=1,ipre
          preoff=preoff+acc(i,j)
        enddo
        preoff=preoff/dble(ipre)
        do i=1,nwin
          acc(i,j)=acc(i,j)-preoff
        enddo
      enddo
c
c     select signal and pst-seismic time windows
c
      ipga=1+idint((tpga(ist)-start(ist))/dt)
      isdw=1+idint((tsdw(ist)-start(ist))/dt)
      iddw=1+idint((tddw(ist)-start(ist))/dt)
c
      nwin=min0(nwin,iddw+isdw-ipre)
c
c     integrate accelerograms with event induced baseline errors to uncorrected velocigrams
c
      do j=1,3
        vel(1,j)=0.d0
        do i=2,nwin
          vel(i,j)=vel(i-1,j)+acc(i,j)*dt
        enddo
c
        call linefit(ipre,vel(1,j),al,bl)
c
c       update pre-seismic baseline correction
c
        do i=1,nwin
          vel(i,j)=vel(i,j)-al-(bl-al)*dble(i-1)/dble(ipre-1)
        enddo
        do i=1,nwin
          acc(i,j)=acc(i,j)-(bl-al)/(dble(ipre-1)*dt)
        enddo
c
        ipgaj=ipre
        pgaj=0.d0
        do i=ipre+1,isdw
          if(dabs(acc(i,j)).gt.pgaj)then
            pgaj=dabs(acc(i,j))
            ipgaj=i
          endif
        enddo
c
        call bscmono(nwin,ipre,min0(ipga,ipgaj),isdw,vel(1,j),
     &               err(1,j),dt,offset(j,ist),rbserr(j,ist))
c
        dis(1,j)=0.d0
        do i=2,nwin
          dis(i,j)=dis(i-1,j)+vel(i,j)*dt
        enddo
      enddo
c
      return
      end