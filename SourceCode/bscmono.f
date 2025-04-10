      subroutine bscmono(nw,ipp,ipga,isdw,vel,
     &                   bse,dt,offset,rbserr)
      implicit none
c
c     nw: length of velocity time series
c     ipp: time of P arrival
c     ipga: time of pga
c     isdw: start of post-seismic period
c     dt: sampling interval
c     vel(nw): input and returned velocity time series
c     bse(nw): errors estimated based on empirical baseline correction
c     offset: returned estimate of permanent displacement offset
c     rbserr: ratio of the baseline effect to co-seismic average acceleration
c
      integer*4 nw,ipp,ipga,isdw
      real*8 dt,offset,rbserr
      real*8 vel(nw),bse(nw)
c
      integer*4 i,j,i0,i1,i2,i3,i4,i5,i6,iw,ierr
      integer*4 k,k1,k2,k3,k4,k5,k6,id,nd,ip,npeak
      real*8 a0,a1,delta,smin,smax,sigma,gamma,alpha
      real*8 d1,d2,d3,d4,d5,d6,d1opt,d2opt,d3opt,d4opt,d5opt
      logical*2 postsmooth
c
      integer*4 ndeg
      parameter(ndeg=2)
c
      real*8 pi
      data pi/3.14159265358979d+00/
c
      real*8, allocatable:: poly(:,:),bat(:),dis(:),swp(:)
c
      allocate (poly(nw,0:ndeg),stat=ierr)
      if(ierr.ne.0)stop ' Error in bscmono: poly not allocated!'
      allocate (bat(0:ndeg),stat=ierr)
      if(ierr.ne.0)stop ' Error in bscmono: bat not allocated!'
      allocate (dis(nw),stat=ierr)
      if(ierr.ne.0)stop ' Error in bscmono: dis not allocated!'
      allocate (swp(nw),stat=ierr)
      if(ierr.ne.0)stop ' Error in bscmono: swp not allocated!'
c
c     remove pre-event trend
c     ip = time of p arrival
c
      ip=ipp
      call d2dfit(ip,vel,poly,bat,1,bse)
c
      a0=bse(1)
      a1=(bse(ip)-bse(1))/dble(ip-1)
c
      do i=1,nw
        vel(i)=vel(i)-(a0+a1*dble(i-1))
      enddo
c
      iw=isdw
c
c     estimate the earliest start time of baseline shift that is
c     defined by the latest zero crossing of uncorrected
c     displacogram, swp, before i = ip+(ipga-ip)/3
c
      dis(1)=0.d0
      do i=2,nw
        dis(i)=dis(i-1)+vel(i)*dt
      enddo
c
      i0=ip
      do i=1+(2*ip+ipga)/3,ip+1,-1
        if(dis(i)*dis(i-1).le.0.d0)then
          i0=i
          goto 200
        endif
      enddo
200   ip=i0
c
c     get postseismic linear trend of the uncorrected velocitogram
c     through 2nd order polynomial regression of displacogram within
c     the time window i0 <= i <= nw, where i0 = iw+2*(nw-iw)/3
c
      i0=1+(nw+2*iw)/3
      call d2dfit(2+nw-i0,dis(i0-1),poly,bat,2,swp(i0-1))
      do i=i0,nw
        bse(i)=(swp(i)-swp(i-1))/dt
      enddo
c
c     get baseline effect in velocitogram by smoothing the time series bse(i)
c     which is equal to 0 for 1 <= i <= ip, vel(i) for ip+1 <= i <= i0 and
c     the sum of right tappered vel(i) and left tappered postseismic linear
c     trend of vel(i) for i0+1 <= i <= nw
c
      i0=1+(nw+2*iw)/3
      do i=1,ip
        bse(i)=0.d0
      enddo
      do i=ip+1,i0
        bse(i)=vel(i)
      enddo
      do i=i0+1,nw
        bse(i)=vel(i)*dcos(0.5d0*pi*dble(i-i0)/dble(nw-i0))**2
     &        +bse(i)*dsin(0.5d0*pi*dble(i-i0)/dble(nw-i0))**2
      enddo
c
c     smooth bse(i) iteratively by fixing bse(ip) at i = ip and i = nw
c
      postsmooth=.true.
      k=0
400   k=k+1
      do i=ip,nw
        swp(i)=bse(i)
      enddo
      do i=ip+1,nw-1
        bse(i)=(swp(i-1)+swp(i)+swp(i+1))/3.d0
      enddo
c
c     account peaks of current bse(i) (after subtracting the linear trend)
c     within the post-seismic window and repeat the smoothing procedure if
c     there exists more than one peak
c
      if(postsmooth)then
        do i=iw+1,nw-1
          if(bse(i).lt.bse(i-1).and.bse(i).lt.bse(i+1).or.
     &       bse(i).gt.bse(i-1).and.bse(i).gt.bse(i+1))then
            postsmooth=.true.
            goto 400
          endif
        enddo
        postsmooth=.false.
      endif
      npeak=0
      do i=ip+1,iw
        if(bse(i).lt.bse(i-1).and.bse(i).lt.bse(i+1).or.
     &     bse(i).gt.bse(i-1).and.bse(i).gt.bse(i+1))npeak=npeak+1
      enddo
      if(npeak.gt.1)goto 400
c
c     shift ip rightward if the average gradient of bse(i) within
c     ip <= i <= iw and that after iw have the same sign and the former
c     is smaller than the latter
c
      call d2dfit(1+nw-iw,bse(iw),poly,bat,1,swp(iw))
c
      d1=(swp(nw)-swp(iw))/dble(nw-iw)
      d2=swp(iw)/dble(iw-ip)
      if(d2*d1.ge.0.d0.and.dabs(d2).lt.dabs(d1))then
c
c       transient baseline shift has the same sign as the permanent one
c
        ip=(ip+2*max0(ip,iw-idint(swp(iw)/d1)))/3
      endif
c
c     shift iw leftward
c
      i0=1+(2*iw+nw)/3
      call d2dfit(1+nw-i0,bse(i0),poly,bat,1,swp(i0))
      do i=ip,i0-1
        swp(i)=swp(i0)+(swp(nw)-swp(i0))*dble(i-i0)/dble(nw-i0)
      enddo
      sigma=0.d0
      gamma=0.d0
      do i=iw,nw
        sigma=dmax1(sigma,bse(i)-swp(i))
        gamma=dmin1(gamma,bse(i)-swp(i))
      enddo
      delta=sigma-gamma
c
      i0=iw
      sigma=0.d0
      gamma=0.d0
      do i=i0,(ip+iw)/2,-1
        sigma=dmax1(sigma,bse(i)-swp(i))
        gamma=dmin1(gamma,bse(i)-swp(i))
        if(sigma-gamma.gt.delta)then
          iw=i
          goto 300
        endif
      enddo
300   continue
c
c     fix bse(i) as the baseline effect in the uncorrected
c     velocitogram for iw < i <= nw
c
c     get a monotonic polyline function (6 segmants) bse(i) best fitting the
c     uncorrected velocitogram vel (i) within ip < i <= iw by grid search under
c     the condition that bse(i) is a monotone function varies between
c     bse(iw)*[(i-ip)/(iw-ip)]^2 and bse(iw)*(i-ip)/(iw-ip) for ip < i < iw
c
      id=1+(iw-ip)/6
      if(ip+6*id.ge.nw)id=id-1
      i1=ip+id
      i2=i1+id
      i3=i2+id
      i4=i3+id
      i5=i4+id
      i6=i5+id
      smin=0.d0
c
      nd=max0(1,id/5)
c
      k6=36
      d6=bse(i6)
      delta=d6/36.d0
c
      do k5=25,30
        d5=dble(k5)*delta
        do k4=16,min0(24,k5)
          d4=dble(k4)*delta
          do k3=9,min0(18,k4)
            d3=dble(k3)*delta
            do k2=4,min0(12,k3)
              d2=dble(k2)*delta
              do k1=1,min0(6,k2)
                d1=dble(k1)*delta
                sigma=0.d0
                do i=ip+1,i1,nd
                  sigma=sigma+(bse(i)
     &                 -d1*(dble(i-ip)/dble(id))**2)**2
                enddo
                do i=i1+1,i2,nd
                  sigma=sigma+(bse(i)
     &                 -(d1+(d2-d1)*dble(i-i1)/dble(id)))**2
                enddo
                do i=i2+1,i3,nd
                  sigma=sigma+(bse(i)
     &                 -(d2+(d3-d2)*dble(i-i2)/dble(id)))**2
                enddo
                do i=i3+1,i4,nd
                  sigma=sigma+(bse(i)
     &                 -(d3+(d4-d3)*dble(i-i3)/dble(id)))**2
                enddo
                do i=i4+1,i5,nd
                  sigma=sigma+(bse(i)
     &                 -(d4+(d5-d4)*dble(i-i4)/dble(id)))**2
                enddo
                do i=i5+1,i6,nd
                  sigma=sigma+(bse(i)
     &                 -(d5+(d6-d5)*dble(i-i5)/dble(id)))**2
                enddo
                if(smin.le.0.d0.or.sigma.le.smin)then
                  smin=sigma
                  d1opt=d1
                  d2opt=d2
                  d3opt=d3
                  d4opt=d4
                  d5opt=d5
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
c
      do i=1,ip
        bse(i)=0.d0
      enddo
      do i=ip+1,i1
        bse(i)=d1opt*(dble(i-ip)/dble(id))**2
      enddo
      do i=i1+1,i2
        bse(i)=d1opt+(d2opt-d1opt)*dble(i-i1)/dble(id)
      enddo
      do i=i2+1,i3
        bse(i)=d2opt+(d3opt-d2opt)*dble(i-i2)/dble(id)
      enddo
      do i=i3+1,i4
        bse(i)=d3opt+(d4opt-d3opt)*dble(i-i3)/dble(id)
      enddo
      do i=i4+1,i5
        bse(i)=d4opt+(d5opt-d4opt)*dble(i-i4)/dble(id)
      enddo
      do i=i5+1,i6-1
        bse(i)=d5opt+(d6-d5opt)*dble(i-i5)/dble(id)
      enddo
c
600   continue
c
c     re-smoothing
c
      k=0
800   k=k+1
      do i=ipp,nw
        swp(i)=bse(i)
      enddo
      do i=ipp+1,nw-1
        bse(i)=(swp(i-1)+swp(i)+swp(i+1))/3.d0
      enddo
      npeak=0
      do i=ipp+1,nw-1
        if(bse(i).lt.bse(i-1).and.bse(i).lt.bse(i+1).or.
     &     bse(i).gt.bse(i-1).and.bse(i).gt.bse(i+1))npeak=npeak+1
      enddo
      if(k.lt.(iw-ipp)/2.or.npeak.gt.1)goto 800
c
      do i=1,nw
        vel(i)=vel(i)-bse(i)
      enddo
c
      dis(1)=0.d0
      do i=2,nw
        dis(i)=dis(i-1)+vel(i)*dt
      enddo
c
      call rampfit(nw,dis,ipp,iw,i1,i2,offset,sigma,poly)
c
c     final optimization
c
      i0=1+(2*iw+nw)/3
      call d2dfit(1+nw-i0,dis(i0),poly,bat,1,swp(i0))
      delta=(swp(nw)-swp(i0))/dble(nw-i0)/dt
c
      do i=ipga+1,iw
        sigma=dsin(0.5d0*pi*dble(i-ipga)/dble(iw-ipga))**2
        vel(i)=vel(i)-delta*sigma
        bse(i)=bse(i)+delta*sigma
      enddo
      do i=iw+1,nw
        vel(i)=vel(i)-delta
        bse(i)=bse(i)+delta
      enddo
c
      i0=(2*iw+nw)/3
      do i=i0,nw
        swp(i)=vel(i)*dsin(0.5d0*pi*dble(i-i0)/dble(nw-i0))**2
        vel(i)=vel(i)-swp(i)
        bse(i)=bse(i)+swp(i)
      enddo
c
      dis(1)=0.d0
      do i=2,nw
        dis(i)=dis(i-1)+vel(i)*dt
      enddo
c
      call rampfit(nw,dis,ipp,iw,i1,i2,offset,sigma,poly)
c
      if(offset.lt.-99.999d0)then
        offset=-99.999d0
      else if(offset.gt.99.999d0)then
        offset=99.999d0
      endif
c
      sigma=vel(ipp)
      gamma=vel(ipp)
      do i=ipp+1,isdw
        sigma=dmax1(sigma,vel(i))
        gamma=dmin1(gamma,vel(i))
      enddo
      delta=(sigma-gamma)/(dble(isdw-ipp)*dt)
c
      sigma=0.d0
      do i=ipp+1,isdw
        sigma=sigma+dabs(bse(i))*dt
      enddo
      sigma=2.d0*sigma/(dble(isdw-ipp)*dt)**2
c
      gamma=0.d0
      do i=isdw+1,nw
        gamma=gamma+dabs(bse(i))*dt
      enddo
      gamma=2.d0*gamma/(dble(nw-isdw)*dt)**2
c
      alpha=0.d0
      do i=1,ip
        alpha=dmax1(alpha,dabs(vel(i)))
      enddo
      alpha=alpha/(dble(ip)*dt)
c
      rbserr=(sigma+gamma+alpha)/delta
      if(rbserr.gt.99.9999d0)rbserr=99.9999d0
c
c     add the initial pre-seismic correction
c
      do i=1,nw
        bse(i)=bse(i)+a0+a1*dble(i-1)
      enddo
c
      deallocate(poly,bat,swp)
      return
      end