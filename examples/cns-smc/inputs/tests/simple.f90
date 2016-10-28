module advance_module

!  use bl_error_module
!  use multifab_module

  implicit none

  private
      
  ! These index constants are shared with the initial data routine.                                                                    
  integer, parameter, public :: irho = 1
  integer, parameter, public :: imx  = 2
  integer, parameter, public :: imy  = 3
  integer, parameter, public :: imz  = 4
  integer, parameter, public :: iene = 5

  integer, parameter :: qu    = 2
  integer, parameter :: qv    = 3
  integer, parameter :: qw    = 4
  integer, parameter :: qpres = 5

  double precision, parameter :: ALP =  0.8d0
  double precision, parameter :: BET = -0.2d0
  double precision, parameter :: GAM =  4.d0/105.d0
  double precision, parameter :: DEL = -1.d0/280.d0

  public :: advance

contains

  subroutine advance (lo,hi,ng,dx,dt,a,b,Unew,U,Q,D,F,eta,alam,courno)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3),dt,a,b
    integer                       :: i,j,k

    double precision, intent(out) :: Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in ) :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(inout) :: Q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)
    double precision, intent(inout) :: D(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)
    double precision, intent(inout) :: F(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

    double precision, intent(inout), optional :: courno

    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    double precision :: divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz,txx,tyy,tzz
    double precision :: mechwork, uxy,uxz,vyz,wzx,wzy,vyx

    do k = lo(3),hi(3)
       do j = lo(2),hi(3)

          Unew(i,j,k,irho) = a * U(i,j,k,irho) + &
               dt * (D(i,j-1,k,irho) + F(-1,j,k,irho))

          Unew(i,j,k,irho) = a * Unew(i-1,j,k,irho) 
          
          do i = lo(1),hi(1)
             Unew(i,j,k,imx)  = a * U(i-1,j,k,imx) + b * (Unew(i+1,j,k,imx) + &
                                dt * (D(i+1,j,k,imx) + F(i,j,k,imx)))

             Unew(i,j,k,imx)  = a * U(i-1,j,k,imx) + b * (Unew(i+1,j,k,imx) + &
                                dt * (D(i+1,j,k,imx) + F(i,j,k,imx)))
          end do

          do i = lo(1),hi(1)
             Unew(i,j,k,imx)  = a * U(i-1,j,k,imx) + b * (Unew(i+1,j,k,imx) + &
                                dt * (D(i+1,j,k,imx) + F(i,j,k,imx)))

          
          end do
       end do
       do j = lo(2),hi(3)

          Unew(i,j,k,irho) = a * U(i,j,k,irho) + &
               dt * (D(i,j-1,k,irho) + F(-1,j,k,irho))

          Unew(i,j,k,irho) = a * Unew(i-1,j,k,irho) 
          
          do i = lo(1),hi(1)
             Unew(i,j,k,imx)  = a * U(i-1,j,k,imx) + b * (Unew(i+1,j,k,imx) + &
                                dt * (D(i+1,j,k,imx) + F(i,j,k,imx)))

             Unew(i,j,k,imx)  = a * U(i-1,j,k,imx) + b * (Unew(i+1,j,k,imx) + &
                                dt * (D(i+1,j,k,imx) + F(i,j,k,imx)))
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine advance 

end module advance_module

