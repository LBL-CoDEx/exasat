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

  subroutine ctoprim (lo,hi,u,q,dx,ng,courno)

    integer,          intent(in ) :: lo(3), hi(3), ng
    double precision, intent(in ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in ) :: dx(3)
    double precision, intent(out) :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)

    double precision, intent(inout), optional :: courno

    integer          :: i, j, k
    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    !$OMP PARALLEL DO PRIVATE(i,j,k,eint,rhoinv)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             rhoinv     = 1.0d0/u(i,j,k,1)
             q(i,j,k,1) = u(i,j,k,1)
             q(i,j,k,2) = u(i,j,k,2)*rhoinv
             q(i,j,k,3) = u(i,j,k,3)*rhoinv
             q(i,j,k,4) = u(i,j,k,4)*rhoinv

             eint = u(i,j,k,5)*rhoinv - 0.5d0*(q(i,j,k,2)**2 + q(i,j,k,3)**2 + q(i,j,k,4)**2)

             q(i,j,k,5) = (GAMMA-1.d0)*eint*u(i,j,k,1)
             q(i,j,k,6) = eint/CV

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if ( present(courno) ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,c,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                c     = sqrt(GAMMA*q(i,j,k,5)/q(i,j,k,1))
                courx = ( c+abs(q(i,j,k,2)) ) / dx(1)
                coury = ( c+abs(q(i,j,k,3)) ) / dx(2) 
                courz = ( c+abs(q(i,j,k,4)) ) / dx(3)

                courmx = max( courmx, courx )
                courmy = max( courmy, coury )
                courmz = max( courmz, courz )

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
       !
       ! Compute running max of Courant number over grids.
       !
       courno = max( courmx, courmy, courmz , courno )

    end if

  end subroutine ctoprim
                  
end module advance_module

