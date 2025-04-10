      subroutine smgetinp(ierr)
      use smalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,k,ind,ist
      real*8 dnorth,deast,datswp
c
      print *,' Read input file ...'
c
      open(10,file=inputfile,status='old')
c
c     Earthquake parameters
c
      call skipdoc(10)
      read(10,*)year,month,day,hour,minute,hyptime
      call skipdoc(10)
      read(10,*)hyplat,hyplon,hypdep
      hypdep=hypdep*KM2M
c
c     Strong-motion data folder
c
      call skipdoc(10)
      read(10,*)datadir
c
      call skipdoc(10)
      read(10,*)stdismin,stdismax
      stdismin=stdismin*KM2M
      stdismax=stdismax*KM2M
c
c     Output folder
c
      call skipdoc(10)
      read(10,*)outdir
c
      do i=1,80
        if(outdir(i:i).eq.' ')goto 101
      enddo
101   ind=i-1
      if(outdir(ind:ind).eq.'/'.or.outdir(ind:ind).eq.'\')then
        ind=ind-1
        outdir=outdir(1:ind)
      endif
      outdirlen=ind
c
      call skipdoc(10)
      read(10,*)coseis
c
      call skipdoc(10)
      read(10,*)dt
c
      coseis=outdir(1:outdirlen)//'/'//coseis
c
      close(10)
c
      do i=1,80
        if(datadir(i:i).eq.' ')goto 102
      enddo
102   ind=i-1
      if(datadir(ind:ind).eq.'/'.or.datadir(ind:ind).eq.'\')then
        ind=ind-1
        datadir=datadir(1:ind)
      endif
      datadirlen=ind
c
      open(20,file=datadir(1:datadirlen)//'/'//'SMDataInfo.dat',
     &     status='old')
      call skipdoc(20)
      read(20,*)year0,month0,day0,hour0,minute0,hyptime0
      if(year0.ne.year.or.month0.ne.month.or.day0.ne.day.or.
     &   hour0.ne.hour.or.minute0.ne.minute.or.hyptime0.ne.hyptime)then
        print *,' EQ origin times are not identical!'
        stop
      endif
      call skipdoc(20)
      read(20,*)hyplat0,hyplon0,hypdep0
      hypdep0=hypdep0*KM2M
      if(hyplat0.ne.hyplat.or.hyplon0.ne.hyplon.or.
     &   hypdep0.ne.hypdep)then
        print *,' EQ hypocenters are not identical!'
        stop
      endif
c
      call skipdoc(20)
      read(20,*)nst,accunit
c
      call skipdoc(20)
      read(20,*)(icmp(i),i=1,3)
c
      allocate(stcode(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: stcode not allocated!'
      allocate(stclen(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: stclen not allocated!'
      allocate(lat(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: lat not allocated!'
      allocate(lon(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: lon not allocated!'
      allocate(start(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: start not allocated!'
      allocate(ponset(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: ponset not allocated!'
      allocate(length(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: length not allocated!'
      allocate(epidis(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: epidis not allocated!'
      allocate(sample(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: sample not allocated!'
      allocate(offset(3,nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: offset not allocated!'
      allocate(rbserr(3,nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: bspst not allocated!'
      allocate(tpga(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: tpga not allocated!'
      allocate(tsdw(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: tsdw not allocated!'
      allocate(tddw(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: tddw not allocated!'
      allocate(okay(nst),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: okay not allocated!'
c
      k=0
      do ist=1,nst
        call skipdoc(20)
        read(20,*)stcode(ist),lat(ist),lon(ist),start(ist),ponset(ist),
     &            length(ist),sample(ist)
        if(sample(ist).le.0.d0)then
          print *,' SMDataInfo: Bad sampling interval!'
          stop
        endif
        if(ponset(ist).lt.start(ist)+PREWIN)then
          write(*,'(a)')stcode(ist)
     &                //' ... pre-seismic window not enough ...'
        endif
        call disazi(REARTH,hyplat,hyplon,lat(ist),lon(ist),dnorth,deast)
        epidis(ist)=dsqrt(dnorth**2+deast**2)
        if(epidis(ist).ge.stdismin.and.epidis(ist).le.stdismax.and.
     &     ponset(ist).ge.start(ist)+PREWIN)then
          k=k+1
          stcode(k)=stcode(ist)
          lat(k)=lat(ist)
          lon(k)=lon(ist)
          start(k)=start(ist)
          ponset(k)=ponset(ist)
          length(k)=length(ist)
          epidis(k)=epidis(ist)
          sample(k)=sample(ist)
        endif
      enddo
c
      nst=k
      if(nst.le.0)then
        print *,' No data usable!'
        stop
      endif
c
      close(20)
c
c     sort stations by epicentral distance
c
      do i=1,nst
        do j=i+1,nst
          if(epidis(j).lt.epidis(i))then
            stswp=stcode(i)
            stcode(i)=stcode(j)
            stcode(j)=stswp
c
            datswp=lat(i)
            lat(i)=lat(j)
            lat(j)=datswp
c
            datswp=lon(i)
            lon(i)=lon(j)
            lon(j)=datswp
c
            datswp=start(i)
            start(i)=start(j)
            start(j)=datswp
c
            datswp=ponset(i)
            ponset(i)=ponset(j)
            ponset(j)=datswp
c
            datswp=length(i)
            length(i)=length(j)
            length(j)=datswp
c
            datswp=epidis(i)
            epidis(i)=epidis(j)
            epidis(j)=datswp
c
            datswp=sample(i)
            sample(i)=sample(j)
            sample(j)=datswp
          endif
        enddo
      enddo
c
      nwinmax=0
      do ist=1,nst
        nwinmax=max0(nwinmax,1+2*idint(length(ist)/sample(ist)))
        do i=1,80
          if(stcode(ist)(i:i).eq.' ')goto 103
        enddo
103     stclen(ist)=i-1
      enddo
c
      allocate(acc(nwinmax,3),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: acc not allocated!'
      allocate(vel(nwinmax,3),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: vel not allocated!'
      allocate(dis(nwinmax,3),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: dis not allocated!'
      allocate(err(nwinmax,3),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: err not allocated!'
      allocate(swp(nwinmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: swp not allocated!'
c
      allocate(ene(nwinmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: ene not allocated!'
      allocate(dat(nwinmax,3),stat=ierr)
      if(ierr.ne.0)stop ' Error in smgetinp: dat not allocated!'
c
      return
      end