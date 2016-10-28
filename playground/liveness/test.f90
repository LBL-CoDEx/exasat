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

  subroutine hypterm (lo,hi,ng,dx,cons,q,flux)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(inout ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(inout ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(out) :: flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

    integer          :: i,j,k
    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
     do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)

             flux(i,j,k,irho)= flux(i,j,k,irho)- unm3

             flux(i,j,k, irho) = unp1

             cons(i,j,k,imx) = flux(i,j,k,irho)

          enddo
     enddo
    enddo

    q(1,1,1,imx) = cons(1,1,1,irho)

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             flux(i,j,k,irho)= flux(i,j,k,irho)- unm3

             q(i,j,k,imx) = cons(i,j,k,irho)

          enddo
       enddo
    enddo

    !$OMP END PARALLEL DO
end subroutine hypterm

end module advance_module

