module kernels_module
!  use bc_module
!  use chemistry_module, only : nspecies, molecular_weight
!  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
!       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb, &
!       M8, M6, M4, M2, BRB, BLB
!  use variables_module
  implicit none

  private

  public :: hypterm_3d

  ! for 8th-order first derivatives
  double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

  
contains
  function first_deriv_8(u) result(du)
    double precision :: du
    double precision, intent(in) :: u(-4:4)
    du =   D8(1)*(u(1)-u(-1)) &
         + D8(2)*(u(2)-u(-2)) &
         + D8(3)*(u(3)-u(-3)) &
         + D8(4)*(u(4)-u(-4))
  end function first_deriv_8

  subroutine hypterm_3d (lo,hi,ng,dx,cons,q,rhs,dlo,dhi,bclo,bchi, nspecies, &
    stencil_ng, qu, imx, qpres, iry1, irho, imy, imz, ncons, nprim, iene)

    integer,          intent(in ) :: lo(3),hi(3),ng, stencil_ng, nspecies, iene,qu, imx, qpres, iry1, irho, imy, imz, ncons,nprim
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) ::  rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
    integer          , intent(in) :: dlo(3),dhi(3),bclo(3),bchi(3)

    integer          :: i,j,k,n
    double precision :: un(-4:4)
    double precision :: dxinv(3)
    integer :: slo(3), shi(3) 

    double precision :: rose_temp__4
    double precision :: u__2




    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    slo = dlo + stencil_ng
    shi = dhi - stencil_ng

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    rhs = 0.d0

    !$omp parallel private(i,j,k,n,un)
    !$omp do 
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)

             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) = dxinv(n) 
             end do
             
             un = q(i-4:i+4,j,k,qu)
             !temp = cons(i-4:i+4 ,j ,k ,imx) * un             
             !du = 
             !

           !  DOUBLE PRECISION, PUBLIC :: rose_temp__4
           !  DOUBLE PRECISION, DIMENSION((-4):4), PUBLIC :: u__2 = cons((i - 4):(i + 4),j,k,imx)
           !  DOUBLE PRECISION :: du
           !  DOUBLE PRECISION, DIMENSION((-4):4), INTENT(IN) :: u
           !  du = D8(1) *(u__2(1) - u__2(-1)) + D8(2) *(u__2(2) - u__2(-2)) + D8(3) *(u__2(3) - u__2(-3)) + D8(4) *(u__2(4) - u__2(-4))
           ! rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * rose_temp__4
      
!             u__2= cons((i - 4):(i + 4),j,k,imx)
             !  DOUBLE PRECISION :: du
             !  DOUBLE PRECISION, DIMENSION((-4):4), INTENT(IN) :: u
!            
!             rose_temp__4 = D8(1) *(u__2(1) - u__2(-1)) + D8(2) *(u__2(2) - u__2(-2)) + D8(3) *(u__2(3) - u__2(-3)) + D8(4) *(u__2(4) - u__2(-4))
!             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * rose_temp__4
           


             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx) ) 

!             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
 !                 first_deriv_8( cons(i-4:i+4,j,k,imx)*un+q(i-4:i+4,j,k,qpres) )
!
          !   rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
 !                 first_deriv_8( cons(i-4:i+4,j,k,imy)*un ) 

  !           rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
   !               first_deriv_8( cons(i-4:i+4,j,k,imz)*un ) 

    !         rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
     !             first_deriv_8( (cons(i-4:i+4,j,k,iene)+q(i-4:i+4,j,k,qpres))*un )

      !       do n = iry1, iry1+nspecies-1
       !         rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
        !             first_deriv_8( cons(i-4:i+4,j,k,n)*un )
         !    end do
             
          enddo
       enddo
    enddo
    !$omp end do

  end subroutine hypterm_3d

  
end module kernels_module
