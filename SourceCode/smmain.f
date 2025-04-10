      program smmain
      use smalloc
      implicit none
c
      integer*4 ierr
c
c     read data
c
      write(*,'(a,$)')'  Input file name: '
      read(*,'(a)')inputfile
c
      call smgetinp(ierr)
c
c     make baseline correction
c
      call smgetout(ierr)
c
      stop
      end