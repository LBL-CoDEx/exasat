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

  subroutine advance ()

  end subroutine advance 


  subroutine diffterm (lo,hi,ng,dx,q,difflux,eta,alam, ux, wx, vx)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
!    double precision, intent(out ) :: ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
 !   double precision, intent(out ) :: vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
  !  double precision, intent(out ) :: wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(out) :: difflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)
    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: dxinv(3)
    double precision :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    double precision :: divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz,txx,tyy,tzz
    double precision :: mechwork, uxy,uxz,vyz,wzx,wzy,vyx
    integer          :: i,j,k

    double precision, parameter :: OneThird   = 1.0d0/3.0d0
    double precision, parameter :: TwoThirds  = 2.0d0/3.0d0
    double precision, parameter :: FourThirds = 4.0d0/3.0d0

    double precision, parameter :: CENTER = -205.d0/72.d0
    double precision, parameter :: OFF1   =    8.d0/5.d0
    double precision, parameter :: OFF2   =   -0.2d0
    double precision, parameter :: OFF3   =    8.d0/315.d0
    double precision, parameter :: OFF4   =   -1.d0/560.d0

    allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    difflux(:,:,:,irho) = 0.0d0

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             ux(i,j,k)= &
                   (ALP*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                  + BET*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                  + GAM*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                  + DEL*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))*dxinv(1)

             vx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                  + BET*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                  + GAM*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                  + DEL*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))*dxinv(1)

             wx(i,j,k)= &
                   (ALP*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                  + BET*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                  + GAM*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                  + DEL*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))*dxinv(1)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

  end subroutine diffterm

end module advance_module

