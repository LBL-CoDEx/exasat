module kernels_module
  use bc_module
  use chemistry_module, only : nspecies, molecular_weight
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb, &
       M8, M6, M4, M2, BRB, BLB
  use variables_module
  implicit none

  private

  public :: hypterm_3d, compact_diffterm_3d, chemterm_3d, comp_courno_3d, &
       S3D_diffterm_1, S3D_diffterm_2

contains

  subroutine hypterm_3d (lo,hi,ng,dx,cons,q,rhs,dlo,dhi,bclo,bchi)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(out) ::  rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
    integer          , intent(in) :: dlo(3),dhi(3),bclo(3),bchi(3)

    integer          :: i,j,k,n
    double precision :: un(-4:4)
    double precision :: dxinv(3)
    integer :: slo(3), shi(3) 

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

             un = q(i-4:i+4,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx)*un+q(i-4:i+4,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imy)*un ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imz)*un ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_8( (cons(i-4:i+4,j,k,iene)+q(i-4:i+4,j,k,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_8( cons(i-4:i+4,j,k,n)*un )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)

             un = q(i,j-4:j+4,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imx)*un )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy)*un+q(i,j-4:j+4,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imz)*un )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_8( (cons(i,j-4:j+4,k,iene)+q(i,j-4:j+4,k,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_8( cons(i,j-4:j+4,k,n)*un )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un = q(i,j,k-4:k+4,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imx)*un )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imy)*un )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz)*un+q(i,j,k-4:k+4,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_8( (cons(i,j,k-4:k+4,iene)+q(i,j,k-4:k+4,qpres))*un )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_8(cons(i,j,k-4:k+4,n)*un)
             end do

          enddo
       enddo
    enddo
    !$omp end do

    ! ----------------- boundary -----------------------

    ! ----- lo-x boundary -----
    if (dlo(1) .eq. lo(1)) then 
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)

             ! if (bclo(1) .eq. WALL???) then
             !    i = lo(1)
             !    ! use completely right-biased stencil
                
             !    un(0:3) = q(i:i+3,j,k,qu)
                
             !    rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
             !         first_deriv_rb( cons(i:i+3,j,k,imx) ) 
                
             !    rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
             !         first_deriv_rb( cons(i:i+3,j,k,imx)*un(0:3)+q(i:i+3,j,k,qpres) )
                
             !    rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
             !         first_deriv_rb( cons(i:i+3,j,k,imy)*un(0:3) ) 
                
             !    rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
             !         first_deriv_rb( cons(i:i+3,j,k,imz)*un(0:3) ) 
                
             !    rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
             !         first_deriv_rb( (cons(i:i+3,j,k,iene)+q(i:i+3,j,k,qpres))*un(0:3) )
                
             !    do n = iry1, iry1+nspecies-1
             !       rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
             !            first_deriv_rb( cons(i:i+3,j,k,n)*un(0:3) )
             !    end do
             ! end if

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil

             un(-1:2) = q(i-1:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imx)*un(-1:2)+q(i-1:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imy)*un(-1:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imz)*un(-1:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_r3( (cons(i-1:i+2,j,k,iene)+q(i-1:i+2,j,k,qpres))*un(-1:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_r3( cons(i-1:i+2,j,k,n)*un(-1:2) )
             end do

             i = lo(1)+2
             ! use 4th-order stencil

             un(-2:2) = q(i-2:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx)*un(-2:2)+q(i-2:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imy)*un(-2:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imz)*un(-2:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_4( (cons(i-2:i+2,j,k,iene)+q(i-2:i+2,j,k,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_4( cons(i-2:i+2,j,k,n)*un(-2:2) )
             end do

             i = lo(1)+3
             ! use 6th-order stencil

             un(-3:3) = q(i-3:i+3,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx)*un(-3:3)+q(i-3:i+3,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imy)*un(-3:3) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imz)*un(-3:3) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_6( (cons(i-3:i+3,j,k,iene)+q(i-3:i+3,j,k,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_6( cons(i-3:i+3,j,k,n)*un(-3:3) )
             end do

          enddo
       enddo
       !$omp end do
    end if

    ! ----- hi-x boundary -----
    if (dhi(1) .eq. hi(1)) then 
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)

             i = hi(1)-3
             ! use 6th-order stencil

             un(-3:3) = q(i-3:i+3,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx)*un(-3:3)+q(i-3:i+3,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imy)*un(-3:3) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imz)*un(-3:3) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_6( (cons(i-3:i+3,j,k,iene)+q(i-3:i+3,j,k,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_6( cons(i-3:i+3,j,k,n)*un(-3:3) )
             end do

             i = hi(1)-2
             ! use 4th-order stencil

             un(-2:2) = q(i-2:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx)*un(-2:2)+q(i-2:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imy)*un(-2:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imz)*un(-2:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_4( (cons(i-2:i+2,j,k,iene)+q(i-2:i+2,j,k,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_4( cons(i-2:i+2,j,k,n)*un(-2:2) )
             end do

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil

             un(-2:1) = q(i-2:i+1,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imx)*un(-2:1)+q(i-2:i+1,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imy)*un(-2:1) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imz)*un(-2:1) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_l3( (cons(i-2:i+1,j,k,iene)+q(i-2:i+1,j,k,qpres))*un(-2:1) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_l3( cons(i-2:i+1,j,k,n)*un(-2:1) )
             end do

             ! if (bchi(1) .eq. WALL???) then
             !    i = hi(1)
             !    ! use completely left-biased stencil
                
             !    un(-3:0) = q(i-3:i,j,k,qu)
                
             !    rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
             !         first_deriv_lb( cons(i-3:i,j,k,imx) ) 
                
             !    rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
             !         first_deriv_lb( cons(i-3:i,j,k,imx)*un(-3:0)+q(i-3:i,j,k,qpres) )
                
             !    rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
             !         first_deriv_lb( cons(i-3:i,j,k,imy)*un(-3:0) ) 
                
             !    rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
             !         first_deriv_lb( cons(i-3:i,j,k,imz)*un(-3:0) ) 
                
             !    rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
             !         first_deriv_lb( (cons(i-3:i,j,k,iene)+q(i-3:i,j,k,qpres))*un(-3:0) )
                
             !    do n = iry1, iry1+nspecies-1
             !       rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
             !            first_deriv_lb( cons(i-3:i,j,k,n)*un(-3:0) )
             !    end do
             ! end if

          enddo
       enddo
       !$omp end do
    end if

    ! ----- lo-y boundary -----
    if (dlo(2) .eq. lo(2)) then 
       !$omp do
       do k=lo(3),hi(3)

          ! if (bclo(2) .eq. WALL???) then
          !    j = lo(2)
          !    ! use completely right-biased stencil

          !    do i=lo(1),hi(1)
                
          !       un(0:3) = q(i,j:j+3,k,qv)
                
          !       rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
          !            first_deriv_rb( cons(i,j:j+3,k,imy) )
                
          !       rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
          !            first_deriv_rb( cons(i,j:j+3,k,imx)*un(0:3) )
                
          !       rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
          !            first_deriv_rb( cons(i,j:j+3,k,imy)*un(0:3)+q(i,j:j+3,k,qpres) )
                
          !       rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
          !            first_deriv_rb( cons(i,j:j+3,k,imz)*un(0:3) )
                
          !       rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
          !            first_deriv_rb( (cons(i,j:j+3,k,iene)+q(i,j:j+3,k,qpres))*un(0:3) )
                
          !       do n = iry1, iry1+nspecies-1
          !          rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
          !               first_deriv_rb( cons(i,j:j+3,k,n)*un(0:3) )
          !       end do

          !    enddo
          ! end if

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil

          do i=lo(1),hi(1)

             un(-1:2) = q(i,j-1:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imx)*un(-1:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imy)*un(-1:2)+q(i,j-1:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imz)*un(-1:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_r3( (cons(i,j-1:j+2,k,iene)+q(i,j-1:j+2,k,qpres))*un(-1:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_r3( cons(i,j-1:j+2,k,n)*un(-1:2) )
             end do

          enddo

          j = lo(2)+2
          ! use 4th-order stencil

          do i=lo(1),hi(1)

             un(-2:2) = q(i,j-2:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy)*un(-2:2)+q(i,j-2:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imz)*un(-2:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_4( (cons(i,j-2:j+2,k,iene)+q(i,j-2:j+2,k,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_4( cons(i,j-2:j+2,k,n)*un(-2:2) )
             end do

          enddo

          j = lo(2)+3
          ! use 6th-order stencil

          do i=lo(1),hi(1)

             un(-3:3) = q(i,j-3:j+3,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy)*un(-3:3)+q(i,j-3:j+3,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imz)*un(-3:3) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_6( (cons(i,j-3:j+3,k,iene)+q(i,j-3:j+3,k,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_6( cons(i,j-3:j+3,k,n)*un(-3:3) )
             end do

          enddo

       enddo
       !$omp end do
    end if

    ! ----- hi-y boundary -----
    if (dhi(2) .eq. hi(2)) then 
       !$omp do
       do k=lo(3),hi(3)

          j = hi(2)-3
          ! use 6th-order stencil

          do i=lo(1),hi(1)

             un(-3:3) = q(i,j-3:j+3,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy)*un(-3:3)+q(i,j-3:j+3,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imz)*un(-3:3) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_6( (cons(i,j-3:j+3,k,iene)+q(i,j-3:j+3,k,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_6( cons(i,j-3:j+3,k,n)*un(-3:3) )
             end do

          enddo

          j = hi(2)-2
          ! use 4th-order stencil

          do i=lo(1),hi(1)

             un(-2:2) = q(i,j-2:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy)*un(-2:2)+q(i,j-2:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imz)*un(-2:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_4( (cons(i,j-2:j+2,k,iene)+q(i,j-2:j+2,k,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_4( cons(i,j-2:j+2,k,n)*un(-2:2) )
             end do

          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil

          do i=lo(1),hi(1)

             un(-2:1) = q(i,j-2:j+1,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imx)*un(-2:1) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imy)*un(-2:1)+q(i,j-2:j+1,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imz)*un(-2:1) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_l3( (cons(i,j-2:j+1,k,iene)+q(i,j-2:j+1,k,qpres))*un(-2:1) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_l3( cons(i,j-2:j+1,k,n)*un(-2:1) )
             end do

          enddo

          ! if (bchi(2) .eq. WALL???) then
          !    j = hi(2)
          !    ! use completely left-biased stencil
             
          !    do i=lo(1),hi(1)
                
          !       un(-3:0) = q(i,j-3:j,k,qv)
                
          !       rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
          !            first_deriv_lb( cons(i,j-3:j,k,imy) )
                
          !       rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
          !            first_deriv_lb( cons(i,j-3:j,k,imx)*un(-3:0) )
                
          !       rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
          !            first_deriv_lb( cons(i,j-3:j,k,imy)*un(-3:0)+q(i,j-3:j,k,qpres) )
                
          !       rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
          !            first_deriv_lb( cons(i,j-3:j,k,imz)*un(-3:0) )
                
          !       rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
          !            first_deriv_lb( (cons(i,j-3:j,k,iene)+q(i,j-3:j,k,qpres))*un(-3:0) )
                
          !       do n = iry1, iry1+nspecies-1
          !          rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
          !               first_deriv_lb( cons(i,j-3:j,k,n)*un(-3:0) )
          !       end do

          !    enddo
          ! end if

       enddo
       !$omp end do
    end if


    ! ----- lo-z boundary -----
    if (dlo(3) .eq. lo(3)) then

       ! if (bclo(3) .eq. WALL???) then
       !    k = lo(3)
       !    ! use completely right-biased stencil
       !    !$omp do
       !    do j=lo(2),hi(2)
       !       do i=lo(1),hi(1)
                
       !          un(0:3) = q(i,j,k:k+3,qw)
                
       !          rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
       !               first_deriv_rb( cons(i,j,k:k+3,imz) )
                
       !          rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
       !               first_deriv_rb( cons(i,j,k:k+3,imx)*un(0:3) )
                
       !          rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
       !               first_deriv_rb( cons(i,j,k:k+3,imy)*un(0:3) )
                
       !          rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
       !               first_deriv_rb( cons(i,j,k:k+3,imz)*un(0:3)+q(i,j,k:k+3,qpres) )
                
       !          rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
       !               first_deriv_rb( (cons(i,j,k:k+3,iene)+q(i,j,k:k+3,qpres))*un(0:3) )
                
       !          do n = iry1, iry1+nspecies-1
       !             rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
       !                  first_deriv_rb(cons(i,j,k:k+3,n)*un(0:3))
       !          end do
                
       !       enddo
       !    enddo
       !    !$omp end do nowait
       ! end if

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-1:2) = q(i,j,k-1:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imx)*un(-1:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imy)*un(-1:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imz)*un(-1:2)+q(i,j,k-1:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_r3( (cons(i,j,k-1:k+2,iene)+q(i,j,k-1:k+2,qpres))*un(-1:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_r3( cons(i,j,k-1:k+2,n)*un(-1:2) )
             end do

          enddo
       enddo
       !$omp end do nowait

       k = lo(3)+2
       ! use 4th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:2) = q(i,j,k-2:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imy)*un(-2:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz)*un(-2:2)+q(i,j,k-2:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_4( (cons(i,j,k-2:k+2,iene)+q(i,j,k-2:k+2,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_4( cons(i,j,k-2:k+2,n)*un(-2:2) )
             end do

          enddo
       enddo
       !$omp end do nowait

       k = lo(3)+3
       ! use 6th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-3:3) = q(i,j,k-3:k+3,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imy)*un(-3:3) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz)*un(-3:3)+q(i,j,k-3:k+3,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_6( (cons(i,j,k-3:k+3,iene)+q(i,j,k-3:k+3,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_6( cons(i,j,k-3:k+3,n)*un(-3:3) )
             end do

          enddo
       enddo
       !$omp end do

    end if

    ! ----- hi-z boundary -----
    if (dhi(3) .eq. hi(3)) then

       k = hi(3)-3
       ! use 6th-order stencil
       
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-3:3) = q(i,j,k-3:k+3,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imy)*un(-3:3) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz)*un(-3:3)+q(i,j,k-3:k+3,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_6( (cons(i,j,k-3:k+3,iene)+q(i,j,k-3:k+3,qpres))*un(-3:3) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_6(cons(i,j,k-3:k+3,n)*un(-3:3))
             end do

          enddo
       enddo
       !$omp end do nowait

       k = hi(3)-2
       ! use 4th-order stencil
       
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:2) = q(i,j,k-2:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imy)*un(-2:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz)*un(-2:2)+q(i,j,k-2:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_4( (cons(i,j,k-2:k+2,iene)+q(i,j,k-2:k+2,qpres))*un(-2:2) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_4(cons(i,j,k-2:k+2,n)*un(-2:2))
             end do

          enddo
       enddo
       !$omp end do nowait

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:1) = q(i,j,k-2:k+1,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imx)*un(-2:1) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imy)*un(-2:1) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imz)*un(-2:1)+q(i,j,k-2:k+1,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_l3( (cons(i,j,k-2:k+1,iene)+q(i,j,k-2:k+1,qpres))*un(-2:1) )

             do n = iry1, iry1+nspecies-1
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_l3(cons(i,j,k-2:k+1,n)*un(-2:1))
             end do

          enddo
       enddo
       !$omp end do nowait

       ! if (bchi(3) .eq. WALL???) then
       !    k = hi(3)
       !    ! use completely left-biased stencil
          
       !    !$omp do
       !    do j=lo(2),hi(2)
       !       do i=lo(1),hi(1)
                
       !          un(-3:0) = q(i,j,k-3:k,qw)
                
       !          rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
       !               first_deriv_lb( cons(i,j,k-3:k,imz) )
                
       !          rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
       !               first_deriv_lb( cons(i,j,k-3:k,imx)*un(-3:0) )
                
       !          rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
       !               first_deriv_lb( cons(i,j,k-3:k,imy)*un(-3:0) )
                
       !          rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
       !               first_deriv_lb( cons(i,j,k-3:k,imz)*un(-3:0)+q(i,j,k-3:k,qpres) )
                
       !          rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
       !               first_deriv_lb( (cons(i,j,k-3:k,iene)+q(i,j,k-3:k,qpres))*un(-3:0) )
                
       !          do n = iry1, iry1+nspecies-1
       !             rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
       !                  first_deriv_lb(cons(i,j,k-3:k,n)*un(-3:0))
       !          end do
                
       !       enddo
       !    enddo
       !    !$omp end do
       ! end if

    end if

    !$omp end parallel

  end subroutine hypterm_3d


  subroutine compact_diffterm_3d (lo,hi,ng,dx,q,rhs,mu,xi,lam,dxy,dlo,dhi,bclo,bchi)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(out) :: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
    integer          , intent(in) :: dlo(3),dhi(3),bclo(3),bchi(3)

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable, dimension(:,:,:) :: vsp,vsm, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hg, dpy, dxe
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3), dx2inv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies), hhalf
    integer          :: i,j,k,n, qxn, qyn, qhn
    integer :: slo(3), shi(3)

    double precision :: muM8(8), M8p(8), M8X(8)
    double precision :: muM6(6), M6p(6), M6X(6)
    double precision :: muM4(4), M4p(4), M4X(4)
    double precision :: muM2(2), M2p(2), M2X(2)
    double precision :: muBB(4), BBp(4), BBX(4)
    double precision :: rhstmp(nspecies), rhstot, rhsene
    double precision :: Hcell(0:1,2:ncons)
    integer :: iface

    ! used to turn off some terms
    double precision :: finlo(3), finhi(3)
    double precision :: foulo(3), fouhi(3)

    finlo = 1.d0 
    finhi = 1.d0
    foulo = 1.d0 
    fouhi = 1.d0

    do i=1,3
       if (bclo(i) .eq. INLET) then
          finlo(i) = 0.d0
       else if (bclo(i) .eq. OUTLET) then
          foulo(i) = 0.d0
       end if

       if (bchi(i) .eq. INLET) then
          finhi(i) = 0.d0
       else if (bchi(i) .eq. OUTLET) then
          fouhi(i) = 0.d0
       end if
    end do

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    slo = dlo + stencil_ng
    shi = dhi - stencil_ng
    
    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i)**2
    end do

    allocate(ux( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(vx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(wx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(uy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(vy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(wy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))

    allocate(uz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(vz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(wz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))

    allocate(vsp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(vsm(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    !$omp parallel &
    !$omp private(i,j,k,tauxx,tauyy,tauzz,divu) &
    !$omp private(dmvxdy,dmwxdz,dmvywzdx,dmuydx,dmwydz,dmuxwzdy,dmuzdx,dmvzdy,dmuxvydz)

    !$omp workshare
    rhs = 0.d0
    !$omp end workshare

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp do
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=slo(1),shi(1)
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qw))
          enddo

          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qw))

             i = lo(1)+2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = lo(1)+3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))

             i = hi(1)-2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qw))

             i = hi(1)
             ! use completely left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
          end if
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=slo(2),shi(2)   
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          enddo
       enddo

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qw))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qw))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo
       end if

       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qw))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qw))
          enddo
       end if
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=slo(3),shi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+2
       ! use 4th-order stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = lo(3)+3
       ! use 6th-order stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)-2
       ! use 4th-order stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT

       k = hi(3)
       ! use completely left-biased stencil
       !$OMP DO
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          enddo
       enddo
       !$OMP END DO NOWAIT
    end if

    !$omp barrier

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             divu = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*ux(i,j,k) + divu
             tauyy = 2.d0*mu(i,j,k)*vy(i,j,k) + divu
             tauzz = 2.d0*mu(i,j,k)*wz(i,j,k) + divu
             
             ! change in internal energy
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  tauxx*ux(i,j,k) + tauyy*vy(i,j,k) + tauzz*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    ! d()/dx
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=slo(1),shi(1)
             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_8( vsm(i-4:i+4,j,k)*(vy(i-4:i+4,j,k)+wz(i-4:i+4,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uy(i-4:i+4,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uz(i-4:i+4,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx
          end do

          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             i = lo(1)
             ! use completely right-biased stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_rb( vsm(i:i+3,j,k)*(vy(i:i+3,j,k)+wz(i:i+3,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_rb( mu(i:i+3,j,k)*uy(i:i+3,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_rb( mu(i:i+3,j,k)*uz(i:i+3,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx*finlo(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx  *foulo(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx  *foulo(1)

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_r3( vsm(i-1:i+2,j,k)*(vy(i-1:i+2,j,k)+wz(i-1:i+2,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_r3( mu(i-1:i+2,j,k)*uy(i-1:i+2,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_r3( mu(i-1:i+2,j,k)*uz(i-1:i+2,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx

             i = lo(1)+2
             ! use 4th-order stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_4( vsm(i-2:i+2,j,k)*(vy(i-2:i+2,j,k)+wz(i-2:i+2,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_4( mu(i-2:i+2,j,k)*uy(i-2:i+2,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_4( mu(i-2:i+2,j,k)*uz(i-2:i+2,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx

             i = lo(1)+3
             ! use 6th-order stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_6( vsm(i-3:i+3,j,k)*(vy(i-3:i+3,j,k)+wz(i-3:i+3,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_6( mu(i-3:i+3,j,k)*uy(i-3:i+3,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_6( mu(i-3:i+3,j,k)*uz(i-3:i+3,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_6( vsm(i-3:i+3,j,k)*(vy(i-3:i+3,j,k)+wz(i-3:i+3,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_6( mu(i-3:i+3,j,k)*uy(i-3:i+3,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_6( mu(i-3:i+3,j,k)*uz(i-3:i+3,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx

             i = hi(1)-2
             ! use 4th-order stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_4( vsm(i-2:i+2,j,k)*(vy(i-2:i+2,j,k)+wz(i-2:i+2,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_4( mu(i-2:i+2,j,k)*uy(i-2:i+2,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_4( mu(i-2:i+2,j,k)*uz(i-2:i+2,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_l3( vsm(i-2:i+1,j,k)*(vy(i-2:i+1,j,k)+wz(i-2:i+1,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_l3( mu(i-2:i+1,j,k)*uy(i-2:i+1,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_l3( mu(i-2:i+1,j,k)*uz(i-2:i+1,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx

             i = hi(1)
             ! use completely left-biased stencil

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_lb( vsm(i-3:i,j,k)*(vy(i-3:i,j,k)+wz(i-3:i,j,k)) )
             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_lb( mu(i-3:i,j,k)*uy(i-3:i,j,k) )
             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_lb( mu(i-3:i,j,k)*uz(i-3:i,j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvywzdx*finhi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuydx  *fouhi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuzdx  *fouhi(1)
          end if

       end do
    end do
    !$omp end do

    ! d()/dy
    !$omp do
    do k=lo(3),hi(3)

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vx(i,j-4:j+4,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_8( vsm(i,j-4:j+4,k)*(ux(i,j-4:j+4,k)+wz(i,j-4:j+4,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vz(i,j-4:j+4,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do
       end do
       
       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_rb( mu(i,j:j+3,k)*vx(i,j:j+3,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_rb( vsm(i,j:j+3,k)*(ux(i,j:j+3,k)+wz(i,j:j+3,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_rb( mu(i,j:j+3,k)*vz(i,j:j+3,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy  *foulo(2)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy*finlo(2)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy  *foulo(2)
          end do

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_r3( mu(i,j-1:j+2,k)*vx(i,j-1:j+2,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_r3( vsm(i,j-1:j+2,k)*(ux(i,j-1:j+2,k)+wz(i,j-1:j+2,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_r3( mu(i,j-1:j+2,k)*vz(i,j-1:j+2,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_4( mu(i,j-2:j+2,k)*vx(i,j-2:j+2,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_4( vsm(i,j-2:j+2,k)*(ux(i,j-2:j+2,k)+wz(i,j-2:j+2,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_4( mu(i,j-2:j+2,k)*vz(i,j-2:j+2,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_6( mu(i,j-3:j+3,k)*vx(i,j-3:j+3,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_6( vsm(i,j-3:j+3,k)*(ux(i,j-3:j+3,k)+wz(i,j-3:j+3,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_6( mu(i,j-3:j+3,k)*vz(i,j-3:j+3,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do
       end if
       
       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_6( mu(i,j-3:j+3,k)*vx(i,j-3:j+3,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_6( vsm(i,j-3:j+3,k)*(ux(i,j-3:j+3,k)+wz(i,j-3:j+3,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_6( mu(i,j-3:j+3,k)*vz(i,j-3:j+3,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_4( mu(i,j-2:j+2,k)*vx(i,j-2:j+2,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_4( vsm(i,j-2:j+2,k)*(ux(i,j-2:j+2,k)+wz(i,j-2:j+2,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_4( mu(i,j-2:j+2,k)*vz(i,j-2:j+2,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_l3( mu(i,j-2:j+1,k)*vx(i,j-2:j+1,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_l3( vsm(i,j-2:j+1,k)*(ux(i,j-2:j+1,k)+wz(i,j-2:j+1,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_l3( mu(i,j-2:j+1,k)*vz(i,j-2:j+1,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy
          end do

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_lb( mu(i,j-3:j,k)*vx(i,j-3:j,k) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_lb( vsm(i,j-3:j,k)*(ux(i,j-3:j,k)+wz(i,j-3:j,k)) )
             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_lb( mu(i,j-3:j,k)*vz(i,j-3:j,k) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmvxdy  *fouhi(2)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmuxwzdy*finhi(2)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmvzdy  *fouhi(2)
          end do
       end if
    end do
    !$omp end do 
    
    ! d()/dz
    !$omp do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wx(i,j,k-4:k+4) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wy(i,j,k-4:k+4) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_8( vsm(i,j,k-4:k+4)*(ux(i,j,k-4:k+4)+vy(i,j,k-4:k+4)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
    end do
    !$omp end do nowait

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_rb( mu(i,j,k:k+3)*wx(i,j,k:k+3) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_rb( mu(i,j,k:k+3)*wy(i,j,k:k+3) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_rb( vsm(i,j,k:k+3)*(ux(i,j,k:k+3)+vy(i,j,k:k+3)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz  *foulo(3)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz  *foulo(3)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz*finlo(3)
          end do
       end do
       !$omp end do nowait

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_r3( mu(i,j,k-1:k+2)*wx(i,j,k-1:k+2) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_r3( mu(i,j,k-1:k+2)*wy(i,j,k-1:k+2) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_r3( vsm(i,j,k-1:k+2)*(ux(i,j,k-1:k+2)+vy(i,j,k-1:k+2)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait

       k = lo(3)+2
       ! use 4th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_4( mu(i,j,k-2:k+2)*wx(i,j,k-2:k+2) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_4( mu(i,j,k-2:k+2)*wy(i,j,k-2:k+2) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_4( vsm(i,j,k-2:k+2)*(ux(i,j,k-2:k+2)+vy(i,j,k-2:k+2)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait

       k = lo(3)+3
       ! use 6th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_6( mu(i,j,k-3:k+3)*wx(i,j,k-3:k+3) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_6( mu(i,j,k-3:k+3)*wy(i,j,k-3:k+3) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_6( vsm(i,j,k-3:k+3)*(ux(i,j,k-3:k+3)+vy(i,j,k-3:k+3)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_6( mu(i,j,k-3:k+3)*wx(i,j,k-3:k+3) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_6( mu(i,j,k-3:k+3)*wy(i,j,k-3:k+3) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_6( vsm(i,j,k-3:k+3)*(ux(i,j,k-3:k+3)+vy(i,j,k-3:k+3)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait

       k = hi(3)-2
       ! use 4th-order stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_4( mu(i,j,k-2:k+2)*wx(i,j,k-2:k+2) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_4( mu(i,j,k-2:k+2)*wy(i,j,k-2:k+2) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_4( vsm(i,j,k-2:k+2)*(ux(i,j,k-2:k+2)+vy(i,j,k-2:k+2)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_l3( mu(i,j,k-2:k+1)*wx(i,j,k-2:k+1) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_l3( mu(i,j,k-2:k+1)*wy(i,j,k-2:k+1) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_l3( vsm(i,j,k-2:k+1)*(ux(i,j,k-2:k+1)+vy(i,j,k-2:k+1)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz
          end do
       end do
       !$omp end do nowait

       k = hi(3)
       ! use completely left-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_lb( mu(i,j,k-3:k)*wx(i,j,k-3:k) )
             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_lb( mu(i,j,k-3:k)*wy(i,j,k-3:k) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_lb( vsm(i,j,k-3:k)*(ux(i,j,k-3:k)+vy(i,j,k-3:k)) )
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dmwxdz  *fouhi(3)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dmwydz  *fouhi(3)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dmuxvydz*finhi(3)
          end do
       end do
       !$omp end do nowait
    end if

    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz)

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    !$omp parallel &
    !$omp private(i,j,k,n,qxn,qyn,qhn,Htot,Htmp,Ytmp,hhalf,muM8,M8p,M8X) &
    !$omp private(muM6,M6p,M6X,muM4,M4p,M4X,muM2,M2p,M2X,muBB,BBp,BBX) &
    !$omp private(rhstmp,rhstot,rhsene,Hcell,iface)

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$OMP DO
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             do i=dlo(1),dhi(1)
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp barrier

    ! ------- BEGIN x-direction -------
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=slo(1),shi(1)+1
             muM8 = matmul(   mu(i-4:i+3,j,k)      , M8)
             M8p  = matmul(M8, q(i-4:i+3,j,k,qpres))

             Hg(i,j,k,imx) = dot_product(matmul(vsp(i-4:i+3,j,k   ), M8), &
                  &                               q(i-4:i+3,j,k,qu))

             Hg(i,j,k,imy) = dot_product(muM8, q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imz) = dot_product(muM8, q(i-4:i+3,j,k,qw))

             Hg(i,j,k,iene) = dot_product(matmul(lam(i-4:i+3,j,k      ), M8), &
                  &                                q(i-4:i+3,j,k,qtemp))      &
                  &                + dot_product(dpe(i-4:i+3,j,k), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i-4:i+3,j,k,qxn))
                
                Htmp(n) = dot_product(dpy(i-4:i+3,j,k,n), M8p) &
                     +    dot_product(dxy(i-4:i+3,j,k,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i-4:i+3,j,k,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do
          end do

          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             i = lo(1)
             ! use completely right-biased stencil
             muBB = matmul(    mu(i:i+3,j,k) , BRB)
             BBp  = matmul(BRB, q(i:i+3,j,k,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + finlo(1) * dx2inv(1) * &
                  dot_product(matmul(vsp(i:i+3,j,k), BRB), &
                  &                    q(i:i+3,j,k,qu) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) + foulo(1)*dx2inv(1)*dot_product(muBB,q(i:i+3,j,k,qv))
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + foulo(1)*dx2inv(1)*dot_product(muBB,q(i:i+3,j,k,qw))

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(1)*dx2inv(1) * &
                  ( dot_product(matmul(lam(i:i+3,j,k), BRB), &
                  &                      q(i:i+3,j,k,qtemp)) &
                  + dot_product(       dpe(i:i+3,j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BRB, q(i:i+3,j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i:i+3,j,k,n), BBp) &
                     +      dot_product(dxy(i:i+3,j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i:i+3,j,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(1)*dx2inv(1) * rhsene


             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 2nd-order stencil for cell lo(1)+1,j,k
             do iface=0,1 
                i = lo(1)+1 + iface
                muM2 = matmul(   mu(i-1:i,j,k)      , M2)
                M2p  = matmul(M2, q(i-1:i,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-1:i,j,k   ), M2), &
                     &                                  q(i-1:i,j,k,qu))

                Hcell(iface,imy) = dot_product(muM2, q(i-1:i,j,k,qv))
                Hcell(iface,imz) = dot_product(muM2, q(i-1:i,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-1:i,j,k      ), M2), &
                     &                                   q(i-1:i,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-1:i,j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i-1:i,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-1:i,j,k,n), M2p) &
                        +    dot_product(dxy(i-1:i,j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-1:i,j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = lo(1)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 4th-order stencil for cell lo(1)+2,j,k
             do iface=0,1 
                i = lo(1)+2 + iface
                muM4 = matmul(   mu(i-2:i+1,j,k)      , M4)
                M4p  = matmul(M4, q(i-2:i+1,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-2:i+1,j,k   ), M4), &
                     &                                  q(i-2:i+1,j,k,qu))

                Hcell(iface,imy) = dot_product(muM4, q(i-2:i+1,j,k,qv))
                Hcell(iface,imz) = dot_product(muM4, q(i-2:i+1,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-2:i+1,j,k      ), M4), &
                     &                                   q(i-2:i+1,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-2:i+1,j,k), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i-2:i+1,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-2:i+1,j,k,n), M4p) &
                        +    dot_product(dxy(i-2:i+1,j,k,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-2:i+1,j,k,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = lo(1)+2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 6th-order stencil for cell lo(1)+3,j,k
             do iface=0,1 
                i = lo(1)+3 + iface
                muM6 = matmul(   mu(i-3:i+2,j,k)      , M6)
                M6p  = matmul(M6, q(i-3:i+2,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-3:i+2,j,k   ), M6), &
                     &                                  q(i-3:i+2,j,k,qu))

                Hcell(iface,imy) = dot_product(muM6, q(i-3:i+2,j,k,qv))
                Hcell(iface,imz) = dot_product(muM6, q(i-3:i+2,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-3:i+2,j,k      ), M6), &
                     &                                   q(i-3:i+2,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-3:i+2,j,k), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i-3:i+2,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-3:i+2,j,k,n), M6p) &
                        +    dot_product(dxy(i-3:i+2,j,k,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-3:i+2,j,k,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = lo(1)+3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 6th-order stencil for cell hi(1)-3,j,k
             do iface=0,1  ! two faces of 
                i = hi(1)-3 + iface
                muM6 = matmul(   mu(i-3:i+2,j,k)      , M6)
                M6p  = matmul(M6, q(i-3:i+2,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-3:i+2,j,k   ), M6), &
                     &                                  q(i-3:i+2,j,k,qu))

                Hcell(iface,imy) = dot_product(muM6, q(i-3:i+2,j,k,qv))
                Hcell(iface,imz) = dot_product(muM6, q(i-3:i+2,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-3:i+2,j,k      ), M6), &
                     &                                   q(i-3:i+2,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-3:i+2,j,k), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i-3:i+2,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-3:i+2,j,k,n), M6p) &
                        +    dot_product(dxy(i-3:i+2,j,k,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-3:i+2,j,k,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = hi(1)-3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 4th-order stencil for cell hi(1)-2,j,k
             do iface=0,1 
                i = hi(1)-2 + iface
                muM4 = matmul(   mu(i-2:i+1,j,k)      , M4)
                M4p  = matmul(M4, q(i-2:i+1,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-2:i+1,j,k   ), M4), &
                     &                                  q(i-2:i+1,j,k,qu))

                Hcell(iface,imy) = dot_product(muM4, q(i-2:i+1,j,k,qv))
                Hcell(iface,imz) = dot_product(muM4, q(i-2:i+1,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-2:i+1,j,k      ), M4), &
                     &                                   q(i-2:i+1,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-2:i+1,j,k), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i-2:i+1,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-2:i+1,j,k,n), M4p) &
                        +    dot_product(dxy(i-2:i+1,j,k,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-2:i+1,j,k,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = hi(1)-2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 2nd-order stencil for cell hi(1)-1,j,k
             do iface=0,1 
                i = hi(1)-1 + iface
                muM2 = matmul(   mu(i-1:i,j,k)      , M2)
                M2p  = matmul(M2, q(i-1:i,j,k,qpres))

                Hcell(iface,imx) = dot_product(matmul(vsp(i-1:i,j,k   ), M2), &
                     &                                  q(i-1:i,j,k,qu))

                Hcell(iface,imy) = dot_product(muM2, q(i-1:i,j,k,qv))
                Hcell(iface,imz) = dot_product(muM2, q(i-1:i,j,k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i-1:i,j,k      ), M2), &
                     &                                   q(i-1:i,j,k,qtemp))      &
                     &            + dot_product(       dpe(i-1:i,j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i-1:i,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-1:i,j,k,n), M2p) &
                        +    dot_product(dxy(i-1:i,j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-1:i,j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i-1,j,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i-1,j,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             i = hi(1)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             i = hi(1)
             ! use completely left-biased stencil
             muBB = matmul(    mu(i-3:i,j,k) , BLB)
             BBp  = matmul(BLB, q(i-3:i,j,k,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + finhi(1) * dx2inv(1) * &
                  dot_product(matmul(vsp(i-3:i,j,k), BLB), &
                  &                    q(i-3:i,j,k,qu) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) + fouhi(1)*dx2inv(1)*dot_product(muBB,q(i-3:i,j,k,qv))
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + fouhi(1)*dx2inv(1)*dot_product(muBB,q(i-3:i,j,k,qw))

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(1)*dx2inv(1) * &
                  ( dot_product(matmul(lam(i-3:i,j,k), BLB), &
                  &                      q(i-3:i,j,k,qtemp)) &
                  + dot_product(       dpe(i-3:i,j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i-3:i,j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i-3:i,j,k,n), BBp) &
                     +      dot_product(dxy(i-3:i,j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i-3:i,j,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(1)*dx2inv(1) * rhsene
          end if

       end do
    end do
    !$omp end do

    ! add x-direction rhs
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END x-direction -------

    !$omp barrier

    ! ------- BEGIN y-direction -------
    !$omp do
    do k=lo(3),hi(3)

       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)
             
             muM8 = matmul(   mu(i,j-4:j+3,k)      , M8)
             M8p  = matmul(M8, q(i,j-4:j+3,k,qpres))

             Hg(i,j,k,imx) = dot_product(muM8, q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imz) = dot_product(muM8, q(i,j-4:j+3,k,qw))

             Hg(i,j,k,imy) = dot_product(matmul(vsp(i,j-4:j+3,k   ), M8), &
                  &                               q(i,j-4:j+3,k,qv))

             Hg(i,j,k,iene) = dot_product(matmul(lam(i,j-4:j+3,k      ), M8), &
                  &                                q(i,j-4:j+3,k,qtemp))      &
                  +                  dot_product(dpe(i,j-4:j+3,k), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i,j-4:j+3,k,qxn))

                Htmp(n) = dot_product(dpy(i,j-4:j+3,k,n), M8P) &
                     +    dot_product(dxy(i,j-4:j+3,k,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i,j-4:j+3,k,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             muBB = matmul(    mu(i,j:j+3,k) , BRB)
             BBp  = matmul(BRB, q(i,j:j+3,k,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + foulo(2)*dx2inv(2)*dot_product(muBB,q(i,j:j+3,k,qu))
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + foulo(2)*dx2inv(2)*dot_product(muBB,q(i,j:j+3,k,qw))

             rhs(i,j,k,imy) = rhs(i,j,k,imy) + finlo(2) * dx2inv(2) * &
                  dot_product(matmul(vsp(i,j:j+3,k), BRB), &
                  &                    q(i,j:j+3,k,qv) )

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(2)*dx2inv(2) * &
                  ( dot_product(matmul(lam(i,j:j+3,k), BRB), &
                  &                      q(i,j:j+3,k,qtemp)) &
                  + dot_product(       dpe(i,j:j+3,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BRB, q(i,j:j+3,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j:j+3,k,n), BBp) &
                     +      dot_product(dxy(i,j:j+3,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j:j+3,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(2)*dx2inv(2) * rhsene
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell i,lo(2)+1,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+1 + iface
                muM2 = matmul(   mu(i,j-1:j,k)      , M2)
                M2p  = matmul(M2, q(i,j-1:j,k,qpres))

                Hcell(iface,imx) = dot_product(muM2, q(i,j-1:j,k,qu))
                Hcell(iface,imz) = dot_product(muM2, q(i,j-1:j,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-1:j,k   ), M2), &
                     &                                  q(i,j-1:j,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-1:j,k      ), M2), &
                     &                                   q(i,j-1:j,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-1:j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j-1:j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-1:j,k,n), M2p) &
                        +    dot_product(dxy(i,j-1:j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-1:j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = lo(2)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell i,lo(2)+2,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+2 + iface
                muM4 = matmul(   mu(i,j-2:j+1,k)      , M4)
                M4p  = matmul(M4, q(i,j-2:j+1,k,qpres))

                Hcell(iface,imx) = dot_product(muM4, q(i,j-2:j+1,k,qu))
                Hcell(iface,imz) = dot_product(muM4, q(i,j-2:j+1,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-2:j+1,k   ), M4), &
                     &                                  q(i,j-2:j+1,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-2:j+1,k      ), M4), &
                     &                                   q(i,j-2:j+1,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-2:j+1,k), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i,j-2:j+1,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-2:j+1,k,n), M4p) &
                        +    dot_product(dxy(i,j-2:j+1,k,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-2:j+1,k,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = lo(2)+2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell i,lo(2)+3,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+3 + iface
                muM6 = matmul(   mu(i,j-3:j+2,k)      , M6)
                M6p  = matmul(M6, q(i,j-3:j+2,k,qpres))

                Hcell(iface,imx) = dot_product(muM6, q(i,j-3:j+2,k,qu))
                Hcell(iface,imz) = dot_product(muM6, q(i,j-3:j+2,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-3:j+2,k   ), M6), &
                     &                                  q(i,j-3:j+2,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-3:j+2,k      ), M6), &
                     &                                   q(i,j-3:j+2,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-3:j+2,k), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i,j-3:j+2,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-3:j+2,k,n), M6p) &
                        +    dot_product(dxy(i,j-3:j+2,k,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-3:j+2,k,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = lo(2)+3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

       end if

       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell i,hi(2)-3,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-3 + iface
                muM6 = matmul(   mu(i,j-3:j+2,k)      , M6)
                M6p  = matmul(M6, q(i,j-3:j+2,k,qpres))

                Hcell(iface,imx) = dot_product(muM6, q(i,j-3:j+2,k,qu))
                Hcell(iface,imz) = dot_product(muM6, q(i,j-3:j+2,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-3:j+2,k   ), M6), &
                     &                                  q(i,j-3:j+2,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-3:j+2,k      ), M6), &
                     &                                   q(i,j-3:j+2,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-3:j+2,k), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i,j-3:j+2,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-3:j+2,k,n), M6p) &
                        +    dot_product(dxy(i,j-3:j+2,k,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-3:j+2,k,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = hi(2)-3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell i,hi(2)-2,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-2 + iface
                muM4 = matmul(   mu(i,j-2:j+1,k)      , M4)
                M4p  = matmul(M4, q(i,j-2:j+1,k,qpres))

                Hcell(iface,imx) = dot_product(muM4, q(i,j-2:j+1,k,qu))
                Hcell(iface,imz) = dot_product(muM4, q(i,j-2:j+1,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-2:j+1,k   ), M4), &
                     &                                  q(i,j-2:j+1,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-2:j+1,k      ), M4), &
                     &                                   q(i,j-2:j+1,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-2:j+1,k), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i,j-2:j+1,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-2:j+1,k,n), M4p) &
                        +    dot_product(dxy(i,j-2:j+1,k,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-2:j+1,k,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = hi(2)-2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell i,hi(2)-1,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-1 + iface
                muM2 = matmul(   mu(i,j-1:j,k)      , M2)
                M2p  = matmul(M2, q(i,j-1:j,k,qpres))

                Hcell(iface,imx) = dot_product(muM2, q(i,j-1:j,k,qu))
                Hcell(iface,imz) = dot_product(muM2, q(i,j-1:j,k,qw))

                Hcell(iface,imy) = dot_product(matmul(vsp(i,j-1:j,k   ), M2), &
                     &                                  q(i,j-1:j,k,qv))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j-1:j,k      ), M2), &
                     &                                   q(i,j-1:j,k,qtemp))      &
                     &            + dot_product(       dpe(i,j-1:j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j-1:j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-1:j,k,n), M2p) &
                        +    dot_product(dxy(i,j-1:j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-1:j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j-1,k,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j-1,k,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             j = hi(2)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          j = hi(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             muBB = matmul(    mu(i,j-3:j,k) , BLB)
             BBp  = matmul(BLB, q(i,j-3:j,k,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + fouhi(2)*dx2inv(2)*dot_product(muBB,q(i,j-3:j,k,qu))
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + fouhi(2)*dx2inv(2)*dot_product(muBB,q(i,j-3:j,k,qw))

             rhs(i,j,k,imy) = rhs(i,j,k,imy) + finhi(2) * dx2inv(2) * &
                  dot_product(matmul(vsp(i,j-3:j,k), BLB), &
                  &                    q(i,j-3:j,k,qv) )

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(2)*dx2inv(2) * &
                  ( dot_product(matmul(lam(i,j-3:j,k), BLB), &
                  &                      q(i,j-3:j,k,qtemp)) &
                  + dot_product(       dpe(i,j-3:j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i,j-3:j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j-3:j,k,n), BBp) &
                     +      dot_product(dxy(i,j-3:j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j-3:j,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(2)*dx2inv(2) * rhsene
          end do

       end if

    end do
    !$omp end do

    ! add y-direction rhs
    do n=2,ncons
       !$omp do
       do k=lo(3),hi(3)
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END y-direction -------

    !$omp barrier

    ! ------- BEGIN z-direction -------
    !$omp do
    do k=slo(3),shi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             muM8 = matmul(   mu(i,j,k-4:k+3)      , M8)
             M8p  = matmul(M8, q(i,j,k-4:k+3,qpres))

             Hg(i,j,k,imx) = dot_product(muM8, q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imy) = dot_product(muM8, q(i,j,k-4:k+3,qv))

             Hg(i,j,k,imz) = dot_product(matmul(vsp(i,j,k-4:k+3   ), M8), &
                  &                               q(i,j,k-4:k+3,qw))

             Hg(i,j,k,iene) = dot_product(matmul(lam(i,j,k-4:k+3      ), M8), &
                  &                                q(i,j,k-4:k+3,qtemp))      &
                  +                  dot_product(dpe(i,j,k-4:k+3), M8p)

             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                M8X = matmul(M8, q(i,j,k-4:k+3,qxn))

                Htmp(n) = dot_product(dpy(i,j,k-4:k+3,n), M8P) &
                     +    dot_product(dxy(i,j,k-4:k+3,n), M8X)

                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
                     +    dot_product(dxe(i,j,k-4:k+3,n), M8X)

                Htot = Htot + Htmp(n)
                Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       k = lo(3)
       ! use completely right-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             muBB = matmul(    mu(i,j,k:k+3) , BRB)
             BBp  = matmul(BRB, q(i,j,k:k+3,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + foulo(3)*dx2inv(3)*dot_product(muBB,q(i,j,k:k+3,qu))
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + foulo(3)*dx2inv(3)*dot_product(muBB,q(i,j,k:k+3,qv))

             rhs(i,j,k,imz) = rhs(i,j,k,imz) + finlo(3) * dx2inv(3) * &
                  dot_product(matmul(vsp(i,j,k:k+3), BRB), &
                  &                    q(i,j,k:k+3,qw) )

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(3)*dx2inv(3) * &
                  ( dot_product(matmul(lam(i,j,k:k+3), BRB), &
                  &                      q(i,j,k:k+3,qtemp)) &
                  + dot_product(       dpe(i,j,k:k+3), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BRB, q(i,j,k:k+3,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j,k:k+3,n), BBp) &
                     +      dot_product(dxy(i,j,k:k+3,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j,k:k+3,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(3)*dx2inv(3) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(3)*dx2inv(3) * rhsene
             
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,j,lo(3)+1
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+1 + iface
                muM2 = matmul(   mu(i,j,k-1:k)      , M2)
                M2p  = matmul(M2, q(i,j,k-1:k,qpres))

                Hcell(iface,imx) = dot_product(muM2, q(i,j,k-1:k,qu))
                Hcell(iface,imy) = dot_product(muM2, q(i,j,k-1:k,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-1:k   ), M2), &
                     &                                  q(i,j,k-1:k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-1:k      ), M2), &
                     &                                   q(i,j,k-1:k,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-1:k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j,k-1:k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-1:k,n), M2p) &
                        +    dot_product(dxy(i,j,k-1:k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-1:k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = lo(3)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,j,lo(3)+2
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+2 + iface
                muM4 = matmul(   mu(i,j,k-2:k+1)      , M4)
                M4p  = matmul(M4, q(i,j,k-2:k+1,qpres))

                Hcell(iface,imx) = dot_product(muM4, q(i,j,k-2:k+1,qu))
                Hcell(iface,imy) = dot_product(muM4, q(i,j,k-2:k+1,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-2:k+1   ), M4), &
                     &                                  q(i,j,k-2:k+1,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-2:k+1      ), M4), &
                     &                                   q(i,j,k-2:k+1,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-2:k+1), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i,j,k-2:k+1,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-2:k+1,n), M4p) &
                        +    dot_product(dxy(i,j,k-2:k+1,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-2:k+1,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = lo(3)+2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,j,lo(3)+3
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+3 + iface
                muM6 = matmul(   mu(i,j,k-3:k+2)      , M6)
                M6p  = matmul(M6, q(i,j,k-3:k+2,qpres))

                Hcell(iface,imx) = dot_product(muM6, q(i,j,k-3:k+2,qu))
                Hcell(iface,imy) = dot_product(muM6, q(i,j,k-3:k+2,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-3:k+2   ), M6), &
                     &                                  q(i,j,k-3:k+2,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-3:k+2      ), M6), &
                     &                                   q(i,j,k-3:k+2,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-3:k+2), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i,j,k-3:k+2,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-3:k+2,n), M6p) &
                        +    dot_product(dxy(i,j,k-3:k+2,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-3:k+2,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = lo(3)+3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do
       
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,j,hi(3)-3
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-3 + iface
                muM6 = matmul(   mu(i,j,k-3:k+2)      , M6)
                M6p  = matmul(M6, q(i,j,k-3:k+2,qpres))

                Hcell(iface,imx) = dot_product(muM6, q(i,j,k-3:k+2,qu))
                Hcell(iface,imy) = dot_product(muM6, q(i,j,k-3:k+2,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-3:k+2   ), M6), &
                     &                                  q(i,j,k-3:k+2,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-3:k+2      ), M6), &
                     &                                   q(i,j,k-3:k+2,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-3:k+2), M6p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M6X = matmul(M6, q(i,j,k-3:k+2,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-3:k+2,n), M6p) &
                        +    dot_product(dxy(i,j,k-3:k+2,n), M6X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-3:k+2,n), M6X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = hi(3)-3
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,j,hi(3)-2
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-2 + iface
                muM4 = matmul(   mu(i,j,k-2:k+1)      , M4)
                M4p  = matmul(M4, q(i,j,k-2:k+1,qpres))

                Hcell(iface,imx) = dot_product(muM4, q(i,j,k-2:k+1,qu))
                Hcell(iface,imy) = dot_product(muM4, q(i,j,k-2:k+1,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-2:k+1   ), M4), &
                     &                                  q(i,j,k-2:k+1,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-2:k+1      ), M4), &
                     &                                   q(i,j,k-2:k+1,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-2:k+1), M4p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M4X = matmul(M4, q(i,j,k-2:k+1,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-2:k+1,n), M4p) &
                        +    dot_product(dxy(i,j,k-2:k+1,n), M4X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-2:k+1,n), M4X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = hi(3)-2
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,j,hi(3)-1
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-1 + iface
                muM2 = matmul(   mu(i,j,k-1:k)      , M2)
                M2p  = matmul(M2, q(i,j,k-1:k,qpres))

                Hcell(iface,imx) = dot_product(muM2, q(i,j,k-1:k,qu))
                Hcell(iface,imy) = dot_product(muM2, q(i,j,k-1:k,qv))

                Hcell(iface,imz) = dot_product(matmul(vsp(i,j,k-1:k   ), M2), &
                     &                                  q(i,j,k-1:k,qw))

                Hcell(iface,iene) = dot_product(matmul(lam(i,j,k-1:k      ), M2), &
                     &                                   q(i,j,k-1:k,qtemp))      &
                     &            + dot_product(       dpe(i,j,k-1:k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j,k-1:k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-1:k,n), M2p) &
                        +    dot_product(dxy(i,j,k-1:k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-1:k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (q(i,j,k-1,qyn) + q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (q(i,j,k-1,qhn) + q(i,j,k,qhn)) / 2.d0
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n) * hhalf * Htot
                end do
             end do

             k = hi(3)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       k = hi(3)
       ! use completely right-biased stencil
       !$omp do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             muBB = matmul(    mu(i,j,k-3:k) , BLB)
             BBp  = matmul(BLB, q(i,j,k-3:k,qpres))

             rhs(i,j,k,imx) = rhs(i,j,k,imx) + fouhi(3)*dx2inv(3)*dot_product(muBB,q(i,j,k-3:k,qu))
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + fouhi(3)*dx2inv(3)*dot_product(muBB,q(i,j,k-3:k,qv))

             rhs(i,j,k,imz) = rhs(i,j,k,imz) + finhi(3) * dx2inv(3) * &
                  dot_product(matmul(vsp(i,j,k-3:k), BLB), &
                  &                    q(i,j,k-3:k,qw) )

             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(3)*dx2inv(3) * &
                  ( dot_product(matmul(lam(i,j,k-3:k), BLB), &
                  &                      q(i,j,k-3:k,qtemp)) &
                  + dot_product(       dpe(i,j,k-3:k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i,j,k-3:k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j,k-3:k,n), BBp) &
                     +      dot_product(dxy(i,j,k-3:k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j,k-3:k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(3)*dx2inv(3) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(3)*dx2inv(3) * rhsene
             
          end do
       end do
       !$omp end do
    end if

    ! add z-direction rhs
    do n=2,ncons
       !$omp do
       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do
       !$omp end do nowait
    end do
    ! ------- END z-direction -------
    
    !$omp barrier

    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do 
    
    !$omp end parallel

    deallocate(Hg,dpy,dxe,dpe,vsp,vsm)

  end subroutine compact_diffterm_3d


  subroutine chemterm_3d(lo,hi,ng,q,up) ! up is UPrime that has no ghost cells
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in )   :: q (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: up(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer :: iwrk, i,j,k
    double precision :: Yt(nspecies), wdot(nspecies), rwrk

    !$omp parallel do private(i,j,k,iwrk,rwrk,Yt,wdot)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Yt = q(i,j,k,qy1:qy1+nspecies-1)
             call ckwyr(q(i,j,k,qrho), q(i,j,k,qtemp), Yt, iwrk, rwrk, wdot)
             up(i,j,k,iry1:) = up(i,j,k,iry1:) + wdot * molecular_weight
             
          end do
       end do
    end do
    !$omp end parallel do 

  end subroutine chemterm_3d

  subroutine comp_courno_3d(lo,hi,ng,dx,Q,courno)
    integer, intent(in) :: lo(3), hi(3), ng
    double precision, intent(in) :: dx(3)
    double precision, intent(in) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: courno

    integer :: i,j,k, iwrk
    double precision :: dxinv(3), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury, courz

    double precision, parameter :: Ru = 8.31451d7

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel do private(i,j,k,iwrk,rwrk,Tt,X,gamma,Cv,Cp,c) &
    !$omp private(courx,coury,courz) &
    !$omp reduction(max:courno)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Tt = q(i,j,k,qtemp)
             X  = q(i,j,k,qx1:qx1+nspecies-1)
             call ckcvbl(Tt, X, iwrk, rwrk, Cv)
             Cp = Cv + Ru
             gamma = Cp / Cv
             c = sqrt(gamma*q(i,j,k,qpres)/q(i,j,k,qrho))

             courx = (c+abs(q(i,j,k,qu))) * dxinv(1)
             coury = (c+abs(q(i,j,k,qv))) * dxinv(2)
             courz = (c+abs(q(i,j,k,qw))) * dxinv(3)
             
             courno = max( courx, coury, courz , courno )

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine comp_courno_3d

  subroutine S3D_diffterm_1(lo,hi,ng,ndq,dx,q,rhs,mu,xi,qx,qy,qz)
 
    integer,          intent(in ) :: lo(3),hi(3),ng,ndq
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(out) :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(out) :: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    double precision, allocatable, dimension(:,:,:) :: vsm

    double precision :: dxinv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    integer :: i,j,k,n, qxn, qdxn

    allocate(vsm(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qdxn,divu,tauxx,tauyy,tauzz) &
    !$omp   private(dmvxdy,dmwxdz,dmvywzdx,dmuydx,dmwydz,dmuxwzdy,dmuzdx,dmvzdy,dmuxvydz)

    !$omp workshare
    rhs(:,:,:,irho) = 0.d0
    !$omp end workshare

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)
             qx(i,j,k,idu) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qu) )
             qx(i,j,k,idv) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qv) )
             qx(i,j,k,idw) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
             qy(i,j,k,idu) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qu) )
             qy(i,j,k,idv) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qv) )
             qy(i,j,k,idw) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             qz(i,j,k,idu) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qu) )
             qz(i,j,k,idv) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qv) )
             qz(i,j,k,idw) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qw) )
          enddo
       enddo
    enddo
    !$OMP END DO


    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! d(mu*dv/dx)/dy
             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qx(i,j-4:j+4,k,idv) )

             ! d(mu*dw/dx)/dz
             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qx(i,j,k-4:k+4,idw) )

             ! d((xi-2/3*mu)*(vy+wz))/dx
             dmvywzdx = dxinv(1) * &
                  first_deriv_8( vsm(i-4:i+4,j,k)*(qy(i-4:i+4,j,k,idv)+qz(i-4:i+4,j,k,idw)) )

             ! d(mu*du/dy)/dx
             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qy(i-4:i+4,j,k,idu) )

             ! d(mu*dw/dy)/dz
             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qy(i,j,k-4:k+4,idw) )

             ! d((xi-2/3*mu)*(ux+wz))/dy
             dmuxwzdy = dxinv(2) * &
                  first_deriv_8( vsm(i,j-4:j+4,k)*(qx(i,j-4:j+4,k,idu)+qz(i,j-4:j+4,k,idw)) )

             ! d(mu*du/dz)/dx
             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qz(i-4:i+4,j,k,idu) )

             ! d(mu*dv/dz)/dy
             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qz(i,j-4:j+4,k,idv) )

             ! d((xi-2/3*mu)*(ux+vy))/dz
             dmuxvydz = dxinv(3) * &
                  first_deriv_8( vsm(i,j,k-4:k+4)*(qx(i,j,k-4:k+4,idu)+qy(i,j,k-4:k+4,idv)) )

             rhs(i,j,k,imx) = dmvxdy + dmwxdz + dmvywzdx
             rhs(i,j,k,imy) = dmuydx + dmwydz + dmuxwzdy
             rhs(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu
             tauyy = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu
             tauzz = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu
             
             ! change in internal energy
             rhs(i,j,k,iene) = tauxx*qx(i,j,k,idu) + tauyy*qy(i,j,k,idv) + tauzz*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    rhs(:,:,:,iry1:) = 0.d0
    !$omp end workshare

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             qx(i,j,k,idT) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qtemp) )
             qx(i,j,k,idp) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qpres) )

             qy(i,j,k,idT) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qtemp) )
             qy(i,j,k,idp) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qpres) )

             qz(i,j,k,idT) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qtemp) )
             qz(i,j,k,idp) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qpres) )

          enddo
       enddo
    enddo
    !$omp end do nowait

    do n=1,nspecies
       qxn = qx1 + n - 1
       qdxn = idX1 + n -1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qxn) )
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qxn) )
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qxn) )
             enddo
          enddo
       enddo
       !$omp end do nowait
    enddo

    !$omp end parallel

    deallocate(vsm)

  end subroutine S3D_diffterm_1


  subroutine S3D_diffterm_2(lo,hi,ng,ndq,dx,q,rhs,mu,xi,lam,dxy,qx,qy,qz)

    integer,          intent(in )  :: lo(3),hi(3),ng,ndq
    double precision, intent(in )  :: dx(3)
    double precision, intent(in )  :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in )  :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in )  :: dxy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
    double precision, intent(in)   :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(in)   :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(in)   :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ndq)
    double precision, intent(inout):: rhs(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
 
    double precision, allocatable, dimension(:,:,:) :: vp, dpe, FE
    double precision, allocatable, dimension(:,:,:,:) :: dpy, FY
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! NOT USING ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3), rhoVc
    integer          :: i,j,k,n, qxn, qyn, qhn, idXn, iryn

    allocate(vp(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(dpy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(dpe(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    allocate(FY(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate(FE(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qyn,qhn,idXn,iryn,rhoVc)

    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             vp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do nowait

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       !$OMP DO
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             do i=lo(1)-ng,hi(1)+ng
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
       !$omp end do nowait
    end do

    ! ===== mx =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) &
                  + dxinv(1) * first_deriv_8( vp(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idu)) &
                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idu)) &
                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idu))
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== my =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) &
                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idv) ) &
                  + dxinv(2) * first_deriv_8( vp(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idv) ) &
                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idv) ) 
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== mz =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) &
                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idw) ) &
                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idw) ) &
                  + dxinv(3) * first_deriv_8( vp(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idw) )
          end do
       end do
    end do
    !$omp end do
    
    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do

    ! thermal conduction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + dxinv(1) * first_deriv_8( lam(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idT) ) &
                  + dxinv(2) * first_deriv_8( lam(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idT) ) &
                  + dxinv(3) * first_deriv_8( lam(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idT) )
          end do
       end do
    end do
    !$omp end do nowait

    ! x-direction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)-ng,hi(1)+ng

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qx(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qx(i,j,k,idXn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qx(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(1) * first_deriv_8( FY(i-4:i+4,j,k,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(1) * first_deriv_8( FE(i-4:i+4,j,k) )
          end do
       end do
    end do
    !$omp end do

    ! y-direction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qy(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qy(i,j,k,idXn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qy(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(2) * first_deriv_8( FY(i,j-4:j+4,k,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(2) * first_deriv_8( FE(i,j-4:j+4,k) )
          end do
       end do
    end do
    !$omp end do

    ! z-direction
    !$omp do
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             rhoVc = 0.d0
             FE(i,j,k) = dpe(i,j,k) * qz(i,j,k,idp)

             do n=1,nspecies
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = dxy(i,j,k,n)*qz(i,j,k,idXn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qz(i,j,k,idXn)*q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do

    do n=1,nspecies    
       iryn = iry1+n-1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(3) * first_deriv_8( FY(i,j,k-4:k+4,n) )
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(3) * first_deriv_8( FE(i,j,k-4:k+4) )
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(vp,dpy,dpe,FY,FE)

  end subroutine S3D_diffterm_2

end module kernels_module
