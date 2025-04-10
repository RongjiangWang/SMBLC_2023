      subroutine smgetout(ierr)
      use smalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,k,l,m,ist,ip,ipga,isdw,ipre,ipst,iddw,nwin,nsam
      real*8 tplus,pga,delta,sigma
      real*8 accoff(3)
c
      print *,' Read strong-motion data ...'
      print *,' Make baseline correction ...'
c
      write(*,'(a)')'   Station  Lat[deg]  Lon[deg] Epdis[km]'
     &            //'   East[m]  North[m]     Up[m]'
     &            //'   RbserrE   RbserrN   RbserrU'
c
      do ist=1,nst
c
c       read strong-motion data
c
        open(21,file=datadir(1:datadirlen)//'/'
     &      //stcode(ist)(1:stclen(ist))//'.dat',status='old')
        k=0
        do i=1,nwinmax
          read(21,*,end=101)(dat(i,icmp(j)),j=1,3)
          k=k+1
        enddo
101     nwin=k
        length(ist)=dble(nwin-1)*sample(ist)
        close(21)
c
c       initial pre-seismic baseline correction
c
        ipre=1+idint((ponset(ist)-start(ist)-DTP)/sample(ist))
        do j=1,3
          delta=0.d0
          do i=1,ipre
            delta=delta+dat(i,j)
          enddo
          delta=delta/dble(ipre)
          do i=1,nwin
            dat(i,j)=dat(i,j)-delta
          enddo
          accoff(j)=0.d0
          do i=ipre,nwin
            accoff(j)=accoff(j)+dat(i,j)
          enddo
          accoff(j)=accoff(j)/dble(1+nwin-ipre)
        enddo
c
c       determine time of pga and start of post-seismic period
c
        ipga=ipre
        pga=0.d0
        do i=1,ipre
          ene(i)=0.d0
        enddo
        do i=ipre,nwin
          sigma=dsqrt((dat(i,1)-accoff(1))**2
     &               +(dat(i,2)-accoff(2))**2
     &               +(dat(i,3)-accoff(3))**2)
          if(pga.lt.sigma)then
            pga=sigma
            ipga=i
          endif
          ene(i)=ene(i-1)+sigma
        enddo
        tpga(ist)=start(ist)+dble(ipga)*sample(ist)
        nwin=min0(nwin,
     &            ipre+20*idnint((tpga(ist)-ponset(ist))/sample(ist)))
c
        isdw=ipre
        iddw=ipre
        do i=1+ipre,nwin
          if(ene(i).le.SDW*ene(nwin))isdw=i
          if(ene(i).le.DDW*ene(nwin))iddw=i
        enddo
        tsdw(ist)=start(ist)+dble(isdw)*sample(ist)
        tddw(ist)=start(ist)+dble(iddw)*sample(ist)
c
        if(tsdw(ist).lt.tpga(ist))then
          tpga(ist)=tsdw(ist)
        endif
c
        length(ist)=dmin1(length(ist),
     &        tddw(ist)-start(ist)+dmin1(tpga(ist)-start(ist),PSTWIN))
        nwin=idint(length(ist)/sample(ist))
c
        okay(ist)=tsdw(ist).le.
     &    start(ist)+length(ist)-dmin1(tpga(ist)-start(ist),PSTWIN)
        if(.not.okay(ist))then
          write(*,'(a10,a)')stcode(ist)(1:stclen(ist)),
     &                     '   ... data length not enough ...'
          goto 200
        endif
c
c       down sampling
c
        if(sample(ist).lt.dt)then
          nsam=idnint(dt/sample(ist))
        else
          nsam=1
        endif
        ipre=1+idint((ponset(ist)-start(ist)-DTP)/dt)
        nwin=nwin/nsam
c
        do i=1,nwin
          l=max0(0,idint((dble(i-1)-0.5d0)*dt/sample(ist)))
          do j=1,3
            acc(i,j)=0.d0
            do m=l+1,l+nsam
              acc(i,j)=acc(i,j)+dat(m,j)
            enddo
            acc(i,j)=accunit*acc(i,j)/dble(nsam)
          enddo
        enddo
        sample(ist)=dt
c
c       make baseline correction
c
        call smbscw(ist,nwin)
c
        write(*,1001)stcode(ist)(1:stclen(ist)),
     &         lat(ist),lon(ist),epidis(ist)/KM2M,
     &         (offset(j,ist),j=1,3),(rbserr(j,ist),j=1,3)
c
        open(31,file=outdir(1:outdirlen)//'/'
     &     //stcode(ist)(1:stclen(ist))//'_blc.dat',status='unknown')
        write(31,'(a)')'        Time'
     &               //'          VdatE          VdatN          VdatZ'
     &               //'         BlerrE         BlerrN         BlerrZ'
     &               //'      VelocityE      VelocityN      VelocityZ'
     &               //'  DisplacementE  DisplacementN  DisplacementZ'
        do i=1,nwin
          write(31,'(f12.3,$)')start(ist)+dble(i-1)*dt
          write(31,'(3E15.7,$)')(vel(i,j)+err(i,j),j=1,3)
          write(31,'(3E15.7,$)')(err(i,j),j=1,3)
          write(31,'(3E15.7,$)')(vel(i,j),j=1,3)
          write(31,'(3E15.7)')(dis(i,j),j=1,3)
        enddo
        close(31)
200     continue
      enddo
c
      open(32,file=coseis,status='unknown')
      write(32,'(a)')'   Station  Lat[deg]  Lon[deg] Epdis[km]'
     &            //'   East[m]  North[m]     Up[m]'
     &            //'   RbserrE   RbserrN   RbserrU'
      k=0
      do ist=1,nst
        if(okay(ist))then
          k=k+1
          write(32,1001)stcode(ist)(1:stclen(ist)),
     &         lat(ist),lon(ist),epidis(ist)/KM2M,
     &         (offset(j,ist),j=1,3),(rbserr(j,ist),j=1,3)
        endif
      enddo
      nst=k
      close(32)
1001  format(a10,2f10.4,4f10.3,4f10.4)
c
      write(*,'(a,i4,a)')' ====== Baseline correction for ',nst,
     &                   ' stations completed ======='
      return
      end