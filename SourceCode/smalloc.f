      module smalloc
c===================================================================
c     global constants
c===================================================================
      real*8 REARTH,KM2M
      parameter(REARTH=6.371d+06,KM2M=1.0d+03)
      real*8 DEG2RAD
      parameter(DEG2RAD=1.745329251994328d-02)
      real*8 PREWIN,PSTWIN,DTP
      parameter(PREWIN=5.0d+00,PSTWIN=2.5d+01,DTP=1.5d+00)
      real*8 SDW,DDW
      parameter(SDW=0.85d+00,DDW=0.95d+00)
c===================================================================
c     global variables
c===================================================================
      integer*4 nwinmax,nst
      integer*4 year,month,day,hour,minute
      integer*4 year0,month0,day0,hour0,minute0
      integer*4 datadirlen,outdirlen
      integer*4 icmp(3)
      real*8 hyptime,hyptime0
      real*8 hyplat,hyplon,hypdep,hyplat0,hyplon0,hypdep0
      real*8 stdismin,stdismax
      real*8 dt,accunit
      character*10 stswp
      character*80 datadir,outdir,inputfile,coseis
c===================================================================
c     global allocatable variables
c===================================================================
      integer*4, allocatable:: stclen(:)
      real*8, allocatable:: lat(:),lon(:),start(:),ponset(:)
      real*8, allocatable:: epidis(:),tpga(:),tsdw(:),tddw(:)
      real*8, allocatable:: length(:),sample(:),swp(:)
      real*8, allocatable:: acc(:,:),vel(:,:),dis(:,:),err(:,:)
      real*8, allocatable:: dat(:,:),ene(:),offset(:,:),rbserr(:,:)
      logical*2, allocatable:: okay(:)
      character*10, allocatable:: stcode(:)
c
      end module
