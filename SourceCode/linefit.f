      subroutine linefit(n,vel,al,bl)
      implicit none
      integer*4 n
      real*8 vel(n)
c
      real*8 al,bl
c
      integer*4 i
      real*8 det
      real*8 mat(2,2),bat(2)
c
      bat(1)=0.d0
      bat(2)=0.d0
      mat(1,1)=dble(n)
      mat(1,2)=0.d0
      mat(2,2)=0.d0
      do i=1,n
        bat(1)=bat(1)+vel(i)
        bat(2)=bat(2)+vel(i)*dble(i-1)
        mat(1,2)=mat(1,2)+dble(i-1)
        mat(2,2)=mat(2,2)+dble(i-1)**2
      enddo
      mat(2,1)=mat(1,2)
c
      det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      if(det.eq.0.d0)then
        print *,' Error in linefit (singularity problem)!'
        stop
      endif
c
      al=(mat(2,2)*bat(1)-mat(1,2)*bat(2))/det
      bl=(mat(1,1)*bat(2)-mat(2,1)*bat(1))/det
c
      bl=al+bl*dble(n-1)
c
      return
      end