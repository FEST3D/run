program main
  implicit none

  integer, parameter :: xblocks = 6
  !< Number of block in I-directions
  integer, parameter :: yblocks = 1
  !< Number of block in J-directions
  integer, parameter :: zblocks = 1
  !< Number of block in Z-directions
  integer, parameter :: levelx = 1
  !< Coarse the mesh in I-direction by levelx amount
  integer, parameter :: levely = 1
  !< Coarse the mesh in J-direction by levely amount
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
  character(len=*), parameter :: infilename="raetaf.x"
  character(len=32) :: outfilename
  character(len=*), parameter :: outfolder="grid"
  integer :: count=0
  integer :: i, j, k, n

  !-----------------------------------------------------------
  !blocking in x according to x=0.0
  !-----------------------------------------------------------

  open(1, file= infilename, status='old', action='read')

  read(1,*) imx, jmx
  kmx = 2
 ! print*, imx, jmx, kmx
  allocate(x(1:imx,1:jmx,1:kmx))
  allocate(y(1:imx,1:jmx,1:kmx))
  allocate(z(1:imx,1:jmx,1:kmx))
  read(1,*) ((x(i,j,1), i=1,imx), j=1,jmx), ((y(i,j,1), i=1,imx), j=1,jmx)
  z(:,:,1) = 0.0
  z(:,:,2) = 1.0
  x(:,:,2) = x(:,:,1)
  y(:,:,2) = y(:,:,1)

  xnp = (imx-1)/xblocks + 1
  ynp = (jmx-1)/yblocks + 1
  znp = (kmx-1)/zblocks + 1

  call find_indices()
  call make_blocks()
  print*, "---> Created blocked mesh in 'grid/' folder"



  contains

    subroutine find_indices()
      implicit none
      imin(1) = 1
      jmin(1) = 1
      kmin(1) = 1
      imax(xblocks) = imx
      jmax(yblocks) = jmx
      kmax(zblocks) = kmx

      !i direction
      Imax(1) = 33
      Imin(2) = 33
      
      do n = 2,xblocks-1
        imax(n)     = imin(n)+((304)/(xblocks-1))+1
        imin(n+1)   = imax(n)
      end do
      Imax(xblocks) = imx
      Imin(xblocks) = imx-32
      Imax(xblocks-1) = imx-32


      do n = 1,yblocks-1
        jmax(n)     = jmin(n)+((jmx-jmin(1))/(yblocks))
        jmin(n+1)   = jmax(n)
      end do

      do n = 1,zblocks-1
        kmax(n)     = kmin(n)+((kmx-kmin(1))/(zblocks))
        kmin(n+1)   = kmax(n)
      end do
      !print*, imin,imax
      !print*, jmin,jmax
      !print*, kmin,kmax

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
