program main
  implicit none

  integer, parameter :: xblocks = 4
  integer, parameter :: yblocks = 1
  integer, parameter :: zblocks = 1
  integer, parameter :: levelx = 2
  integer, parameter :: levely = 2
  integer :: imx, jmx, kmx
  integer, dimension(xblocks):: imin, imax
  integer, dimension(yblocks):: jmin, jmax
  integer, dimension(zblocks):: kmin, kmax
  integer :: xb,yb, zb
  integer :: xnp,ynp, znp
  !--------------------------------------------------------
  ! nx is multiple of 17, min 34 
  ! ny is multiple of 4, min  24 
  ! 35 x 25 -> coarest mesh
  ! read mesh obatin from nasa archieve
  !---------------------------------------------------------

  real*8, dimension(:, :, :), allocatable :: x, y, z
  character(len=*), parameter :: infilename="grid.txt"
  character(len=32) :: outfilename
  character(len=*), parameter :: outfolder="grid"
  integer :: count=0
  integer :: i, j, k, n

  !-----------------------------------------------------------
  !blocking in x according to x=0.0
  !-----------------------------------------------------------

  open(1, file= infilename, status='old', action='read')

  read(1,*) imx, jmx, kmx
  print*, imx, jmx, kmx
  allocate(x(1:imx,1:jmx,1:kmx))
  allocate(y(1:imx,1:jmx,1:kmx))
  allocate(z(1:imx,1:jmx,1:kmx))
  do k =1, kmx
    do j = 1,jmx
      do i = 1,imx
        read(1,*) x(i,j,k),&
                  y(i,j,k),&
                  z(i,j,k)
        end do
      end do
    end do

  xnp = (imx-1)/xblocks + 1
  ynp = (jmx-1)/yblocks + 1
  znp = (kmx-1)/zblocks + 1

  call find_indices()
  call make_blocks()



  contains

    subroutine find_indices()
      implicit none
      imin(1) = 1
      jmin(1) = 1
      kmin(1) = 1
      imax(xblocks) = imx
      jmax(yblocks) = jmx
      kmax(zblocks) = kmx
      j = 1
      k = 1
      do i= 1,imx
        if (x(i,j,k)<=0.0 .and. xblocks > 1) then
          imin(1)       = 1
          imax(1)       = i
          imin(2)       = i
          imax(xblocks) = imx
          imax(xblocks) = imx
          do n = 2,xblocks-1
            imax(n)     = imin(n)+((imx-imax(1))/(xblocks-1))
            imin(n+1)     = imax(n)
          end do
        end if
      end do

      do n = 1,yblocks-1
        jmax(n)     = jmin(n)+((jmx-jmin(1))/(yblocks))
        jmin(n+1)   = jmax(n)
      end do

      do n = 1,zblocks-1
        kmax(n)     = kmin(n)+((kmx-kmin(1))/(zblocks))
        kmin(n+1)   = kmax(n)
      end do
      print*, imin,imax
      print*, jmin,jmax
      print*, kmin,kmax

    end subroutine find_indices


    subroutine make_blocks()
      implicit none
      
        call system('mkdir -p '//trim(outfolder))
        do zb = 1,zblocks
          do yb = 1,yblocks
            do xb = 1,xblocks
!              print*, imin, imax
!              print*, jmin, jmax
!              print*, kmin, kmax
              write(outfilename,'(a,i2.2,a)') trim(outfolder)//"/grid_", count, ".txt"
              !write(outfilename,'(a,i4.4)') trim(outfolder)//"/fort.", count+1
              open(2, file=outfilename)
              write(2, '(3I4.2)') (imax(xb)-imin(xb))/levelx+ 1, &
                                  (jmax(yb)-jmin(yb))/levely+ 1, &
                                  (kmax(zb)-kmin(zb))+ 1
              do k = kmin(zb),kmax(zb)
                do j = jmin(yb),jmax(yb), levely
                  do i = imin(xb),imax(xb), levelx
                    write(2,*) x(i,j,k), y(i,j,k), z(i,j,k)
                  end do
                end do
              end do
              close(2)
              count= count + 1
            end do
          end do
        end do
      
        close(1)
      end subroutine make_blocks
end program main
