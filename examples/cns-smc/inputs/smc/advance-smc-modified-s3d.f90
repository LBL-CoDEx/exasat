module advance_module

!  use bl_error_module
!  use derivative_stencil_module
!  use kernels_module
!  use multifab_module
!  use time_module
!  use transport_properties
!  use variables_module

!  use chemistry_module, only : nspecies

  implicit none

  private
  ! These index constants are shared with the initial data routine.                                                                    
  integer, parameter, public :: irho = 1
  integer, parameter, public :: imx  = 2
  integer, parameter, public :: imy  = 3
  integer, parameter, public :: imz  = 4
  integer, parameter, public :: iene = 5
  integer, parameter, public :: qtemp = 6
  integer, parameter, public :: qe = 7
  integer, parameter, public :: qy1 = 8

  integer, parameter :: qrho    = 1
  integer, parameter :: qu    = 2
  integer, parameter :: qv    = 3
  integer, parameter :: qw    = 4
  integer, parameter :: qpres = 5

  ! Arithmetic constants 
  double precision, parameter :: Zero          = 0.d0
  double precision, parameter :: One           = 1.d0
  double precision, parameter :: OneThird      = 1.d0/3.d0
  double precision, parameter :: TwoThirds     = 2.d0/3.d0
  double precision, parameter :: FourThirds    = 4.d0/3.d0
  double precision, parameter :: OneQuarter    = 1.d0/4.d0
  double precision, parameter :: ThreeQuarters = 3.d0/4.d0

  ! Indices for S3D style first-derivatives
  integer, parameter :: idu=1, idv=2, idw=3, idT=4, idp=5, idX1=6 

  public advance

contains
  subroutine multifab_fill_boundary(U, lo, hi, ncons, ng)
    
    integer, intent(in) :: lo(3),hi(3),ng
    integer, intent(in) :: ncons
    
    double precision, intent(in   ) :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)
    
  end subroutine multifab_fill_boundary
  
  subroutine chemterm_3d(lo,hi,ng, Q, Uprime, ncons, nprim, nspecies)
    integer, intent(in)::ncons, nprim, nspecies
    
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: Q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: Uprime(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)
    
  end subroutine chemterm_3d

  subroutine ctoprim(U, Q, ng, lo, hi, nprim, ncons, courno, dx, nspecies)

    integer,          intent(in ) :: lo(3),hi(3),ng, nprim, ncons, nspecies
    double precision, intent(inout),optional ::courno
    double precision, intent(in) :: dx(3)

    double precision, intent(in ) :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)
    double precision, intent(inout ) :: Q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)

    integer :: iwrk, i, j, k , n, qxn, qx1, qh1, iry1
    double precision :: dxinv(3), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), Y(nspecies), h(nspecies), gamma
    double precision :: courx, coury, courz

    double precision ::rho, rhoinv, ei, Pt
    double precision, parameter :: Ru = 8.31451d7

    !call ctoprim_3d(lo, hi, U, Q, ng)
    !$omp parallel do private(i, j, k, n, iwrk, rho, rhoinv, rwrk) &
    !$omp private(X, Y, h, ei, Tt, Pt)
    do k = lo(3)-ng, hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1) + ng

             rho = U(i,j,k,irho)
             rhoinv = 1.d0/rho
             Q(i,j,k,qrho) = rho
             Q(i,j,k,qu) = U(i,j,k,imx) * rhoinv
             Q(i,j,k,qv) = U(i,j,k,imy) * rhoinv
             Q(i,j,k,qw) = U(i,j,k,imz) * rhoinv

             do n=1,nspecies
                Y(n) = U(i,j,k,iry1+n-1) * rhoinv
                Q(i,j,k,qy1+n-1) = Y(n)
             end do

!             call ckytx(Y, iwrk, rwrk, X)

             do n=1,nspecies               
                Q(i,j,k,qx1+n-1) = X(n)
             end do

             ei = rhoinv*U(i,j,k,iene) - 0.5d0*(Q(i,j,k,qu)**2+Q(i,j,k,qv)**2+Q(i,j,k,qw)**2)
             Q(i,j,k,qe) = ei

!             call feeytt(ei, Y, iwrk, rwrk, Tt)
             Q(i,j,k,qtemp) = Tt

!             call CKPY(rho, Tt, Y, iwrk, rwrk, Pt)
             Q(i,j,k,qpres) = Pt

!             call ckhms(Tt, iwrk, rwrk, h)

             do n=1,nspecies
                Q(i,j,k,qh1+n-1) = h(n)
             end do
          enddo
       enddo
    enddo
    !$omp end parallel do

  
    !call build(bpt_courno, "courno")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    if (present(courno)) then
       do i=1,3
          dxinv(i) = 1.0d0 / dx(i)
       end do
       
       !$omp parallel do private(i,j,k,iwrk,rwrk,Tt,X,gamma,Cv,Cp,c) &
       !$omp private(courx,coury,courz) &
       !$omp reduction(max:courno)
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                
                Tt = Q(i,j,k,qtemp)
                
                do n=1,nspecies
                   qxn = qx1 + n -1
                   X(n)  = Q(i,j,k,qxn)
                end do

!                call ckcvbl(Tt, X, iwrk, rwrk, Cv)
                Cp = Cv + Ru
                gamma = Cp / Cv
                c = sqrt(gamma*Q(i,j,k,qpres)/Q(i,j,k,qrho))
                
                courx = (c+abs(Q(i,j,k,qu))) * dxinv(1)
                coury = (c+abs(Q(i,j,k,qv))) * dxinv(2)
                courz = (c+abs(Q(i,j,k,qw))) * dxinv(3)
                
                courno = max( courx, coury, courz , courno )
                
             end do
          end do
       end do
       !$omp end parallel do
       
       !   call compute_courno(Q, dx, courno)
    end if

  end subroutine ctoprim

 subroutine get_transport_properties(Q, mu, xi, lam, Ddiag, lo, hi, ng, nprim, nspecies)

   integer, intent(in) :: lo(3), hi(3), ng, nprim, nspecies
   double precision, intent(in ) :: Q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
   
   
   double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
   double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
   double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
   double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

  end subroutine get_transport_properties

  subroutine advance (lo,hi,ng,dx,dt,Uprime,U, Unew, istep, nspecies, ncons, nprim)
!  subroutine advance(U, dt, dx, istep)

    integer, intent(in) :: istep
    integer, intent(in) :: lo(3),hi(3),ng
    integer, intent(in) :: nspecies, ncons, nprim

    double precision, intent(in   ) :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)
    double precision, intent(inout   ) :: Uprime(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)
    double precision, intent(inout   ) :: Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)

    double precision,  intent(  out) :: dt
    double precision,  intent(in   ) :: dx(3)

    double precision, allocatable, dimension(:,:,:,:) :: Q, Ddiag, Fdif, Fhyp
    double precision, allocatable, dimension(:,:,:) :: mu, xi, lam

    double precision, allocatable, dimension(:,:,:) :: vsm, dpe, vp, FE
    double precision, allocatable, dimension(:,:,:,:) :: qx, qy, qz, dpy, FY

!    double precision :: un(-4:4)

    double precision :: dxinv(3), dx2inv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies), hhalf

    double precision :: rho, rhoVc
!    double precision, parameter :: Ru = 8.31451d7

    double precision :: M8p(8), M8X(8), mmtmp(8)
    double precision, save, dimension(8,8) :: M8
    ! for 8th-order first derivatives                                                                                                                                                                                                                                             
    double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

    double precision :: courno

    integer :: qx1, iry1, qh1
    integer :: i,j,k, m, n, dm
    integer :: iryn, qxn, qyn, qhn
    integer :: dlo(3), dhi(3)

    allocate( Q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim))
    allocate(Fhyp(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons))
    allocate(mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(Ddiag(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate( Fdif(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons))

    ! RK Step 1
    !call build(bpt_rkstep1, "rkstep1")   !! vvvvvvvvvvvvvvvvvvvvvvv timer

    courno = 1.0d-50

    !call build(bpt_rkstep1dUdt, "rkstep1dUdt")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call dUdt(U, Uprime, dx, courno=courno_proc)

    call multifab_fill_boundary(U, lo, hi, ncons, ng)

    ! Calculate primitive variables based on U
    !
    !call build(bpt_ctoprim, "ctoprim")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call ctoprim(U, Q, ng)
    call ctoprim(U, Q, ng, lo, hi, nprim, ncons, courno, dx, nspecies)
    !call destroy(bpt_ctoprim)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !call destroy(bpt_courno)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !call build(bpt_gettrans, "gettrans")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call get_transport_properties(Q, mu, xi, lam, Ddiag)
    call get_transport_properties(Q, mu, xi, lam, Ddiag, lo, hi, ng, nprim, nspecies)
    !call destroy(bpt_gettrans)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !
    ! Hyperbolic terms
    !
    !call build(bpt_hypterm, "hypterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
   
    !up => dataptr(U,n)
    !qp => dataptr(Q,n)
    !fhp=> dataptr(Fhyp,n)
    
    !call hypterm_3d(lo,hi,ng,dx,up,qp,fhp)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

!EXPAND    Fhyp = 0.d0
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fhyp(i,j,k,irho) = 0.d0
             Fhyp(i,j,k,imx ) = 0.d0
             Fhyp(i,j,k,imy ) = 0.d0
             Fhyp(i,j,k,imz ) = 0.d0
             Fhyp(i,j,k,iene) = 0.d0
             do n=1,nspecies
                iryn = iry1+n-1
                Fhyp(i,j,k,iryn) = 0.d0
             enddo
          enddo
       enddo
    enddo
    
    !$omp parallel private(i,j,k,n,un)
    !$omp do 
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

!             un = Q(i-4:i+4,j,k,qu)

!EXPAND             Fhyp(i,j,k,irho) = Fhyp(i,j,k,irho) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imx) ) 
             Fhyp(i,j,k,irho) = Fhyp(i,j,k,irho) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imx)-U(i-1,j,k,imx)) &
                + D8(2)*(U(i+2,j,k,imx)-U(i-2,j,k,imx)) &
                + D8(3)*(U(i+3,j,k,imx)-U(i-3,j,k,imx)) &
                + D8(4)*(U(i+4,j,k,imx)-U(i-4,j,k,imx)) )

!EXPAND             Fhyp(i,j,k,imx) = Fhyp(i,j,k,imx) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imx)*un+Q(i-4:i+4,j,k,qpres) )
             Fhyp(i,j,k,imx) = Fhyp(i,j,k,imx) - dxinv(1) * &
                ( D8(1)*((U(i+1,j,k,imx)*Q(i+1,j,k,qu)+Q(i+1,j,k,qpres))-(U(i-1,j,k,imx)*Q(i-1,j,k,qu)+Q(i-1,j,k,qpres))) &
                + D8(2)*((U(i+2,j,k,imx)*Q(i+2,j,k,qu)+Q(i+2,j,k,qpres))-(U(i-2,j,k,imx)*Q(i-2,j,k,qu)+Q(i-2,j,k,qpres))) &
                + D8(3)*((U(i+3,j,k,imx)*Q(i+3,j,k,qu)+Q(i+3,j,k,qpres))-(U(i-3,j,k,imx)*Q(i-3,j,k,qu)+Q(i-3,j,k,qpres))) &
                + D8(4)*((U(i+4,j,k,imx)*Q(i+4,j,k,qu)+Q(i+4,j,k,qpres))-(U(i-4,j,k,imx)*Q(i-4,j,k,qu)+Q(i-4,j,k,qpres))) )

!EXPAND             Fhyp(i,j,k,imy) = Fhyp(i,j,k,imy) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imy)*un ) 
             Fhyp(i,j,k,imy) = Fhyp(i,j,k,imy) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imy)*Q(i+1,j,k,qu)-U(i-1,j,k,imy)*Q(i-1,j,k,qu)) &
                + D8(2)*(U(i+2,j,k,imy)*Q(i+2,j,k,qu)-U(i-2,j,k,imy)*Q(i-2,j,k,qu)) &
                + D8(3)*(U(i+3,j,k,imy)*Q(i+3,j,k,qu)-U(i-3,j,k,imy)*Q(i-3,j,k,qu)) &
                + D8(4)*(U(i+4,j,k,imy)*Q(i+4,j,k,qu)-U(i-4,j,k,imy)*Q(i-4,j,k,qu)) )

!EXPAND             Fhyp(i,j,k,imz) = Fhyp(i,j,k,imz) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imz)*un ) 
             Fhyp(i,j,k,imz) = Fhyp(i,j,k,imz) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imz)*Q(i+1,j,k,qu)-U(i-1,j,k,imz)*Q(i-1,j,k,qu)) &
                + D8(2)*(U(i+2,j,k,imz)*Q(i+2,j,k,qu)-U(i-2,j,k,imz)*Q(i-2,j,k,qu)) &
                + D8(3)*(U(i+3,j,k,imz)*Q(i+3,j,k,qu)-U(i-3,j,k,imz)*Q(i-3,j,k,qu)) &
                + D8(4)*(U(i+4,j,k,imz)*Q(i+4,j,k,qu)-U(i-4,j,k,imz)*Q(i-4,j,k,qu)) )

!EXPAND             Fhyp(i,j,k,iene) = Fhyp(i,j,k,iene) - dxinv(1) * &
!EXPAND                  first_deriv_8( (U(i-4:i+4,j,k,iene)+Q(i-4:i+4,j,k,qpres))*un )
             Fhyp(i,j,k,iene) = Fhyp(i,j,k,iene) - dxinv(1) * &
                ( D8(1)*((U(i+1,j,k,iene)+Q(i+1,j,k,qpres))*Q(i+1,j,k,qu)-(U(i-1,j,k,iene)+Q(i-1,j,k,qpres))*Q(i-1,j,k,qu)) &
                + D8(2)*((U(i+2,j,k,iene)+Q(i+2,j,k,qpres))*Q(i+2,j,k,qu)-(U(i-2,j,k,iene)+Q(i-2,j,k,qpres))*Q(i-2,j,k,qu)) &
                + D8(3)*((U(i+3,j,k,iene)+Q(i+3,j,k,qpres))*Q(i+3,j,k,qu)-(U(i-3,j,k,iene)+Q(i-3,j,k,qpres))*Q(i-3,j,k,qu)) &
                + D8(4)*((U(i+4,j,k,iene)+Q(i+4,j,k,qpres))*Q(i+4,j,k,qu)-(U(i-4,j,k,iene)+Q(i-4,j,k,qpres))*Q(i-4,j,k,qu)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(1) * &
!EXPAND                     first_deriv_8( U(i-4:i+4,j,k,n)*un )
                Fhyp(i,j,k,iryn) = Fhyp(i,j,k,iryn) - dxinv(1) * &
                   ( D8(1)*(U(i+1,j,k,iryn)*Q(i+1,j,k,qu)-U(i-1,j,k,iryn)*Q(i-1,j,k,qu)) &
                   + D8(2)*(U(i+2,j,k,iryn)*Q(i+2,j,k,qu)-U(i-2,j,k,iryn)*Q(i-2,j,k,qu)) &
                   + D8(3)*(U(i+3,j,k,iryn)*Q(i+3,j,k,qu)-U(i-3,j,k,iryn)*Q(i-3,j,k,qu)) &
                   + D8(4)*(U(i+4,j,k,iryn)*Q(i+4,j,k,qu)-U(i-4,j,k,iryn)*Q(i-4,j,k,qu)) )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

!             un = Q(i,j-4:j+4,k,qv)

!EXPAND             Fhyp(i,j,k,irho)=Fhyp(i,j,k,irho) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imy) )
             Fhyp(i,j,k,irho) = Fhyp(i,j,k,irho) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imy)-U(i,j-1,k,imy)) &
                + D8(2)*(U(i,j+2,k,imy)-U(i,j-2,k,imy)) &
                + D8(3)*(U(i,j+3,k,imy)-U(i,j-3,k,imy)) &
                + D8(4)*(U(i,j+4,k,imy)-U(i,j-4,k,imy)) )

!EXPAND             Fhyp(i,j,k,imx)=Fhyp(i,j,k,imx) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imx)*un )
             Fhyp(i,j,k,imx) = Fhyp(i,j,k,imx) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imx)*Q(i,j+1,k,qv)-U(i,j-1,k,imx)*Q(i,j-1,k,qv)) &
                + D8(2)*(U(i,j+2,k,imx)*Q(i,j+2,k,qv)-U(i,j-2,k,imx)*Q(i,j-2,k,qv)) &
                + D8(3)*(U(i,j+3,k,imx)*Q(i,j+3,k,qv)-U(i,j-3,k,imx)*Q(i,j-3,k,qv)) &
                + D8(4)*(U(i,j+4,k,imx)*Q(i,j+4,k,qv)-U(i,j-4,k,imx)*Q(i,j-4,k,qv)) )

!EXPAND             Fhyp(i,j,k,imy)=Fhyp(i,j,k,imy) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imy)*un+Q(i,j-4:j+4,k,qpres) )
             Fhyp(i,j,k,imy) = Fhyp(i,j,k,imy) - dxinv(2) * &
                ( D8(1)*((U(i,j+1,k,imy)*Q(i,j+1,k,qv)+Q(i,j+1,k,qpres))-(U(i,j-1,k,imy)*Q(i,j-1,k,qv)+Q(i,j-1,k,qpres))) &
                + D8(2)*((U(i,j+2,k,imy)*Q(i,j+2,k,qv)+Q(i,j+2,k,qpres))-(U(i,j-2,k,imy)*Q(i,j-2,k,qv)+Q(i,j-2,k,qpres))) &
                + D8(3)*((U(i,j+3,k,imy)*Q(i,j+3,k,qv)+Q(i,j+3,k,qpres))-(U(i,j-3,k,imy)*Q(i,j-3,k,qv)+Q(i,j-3,k,qpres))) &
                + D8(4)*((U(i,j+4,k,imy)*Q(i,j+4,k,qv)+Q(i,j+4,k,qpres))-(U(i,j-4,k,imy)*Q(i,j-4,k,qv)+Q(i,j-4,k,qpres))) )

!EXPAND             Fhyp(i,j,k,imz)=Fhyp(i,j,k,imz) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imz)*un )
             Fhyp(i,j,k,imz) = Fhyp(i,j,k,imz) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imz)*Q(i,j+1,k,qv)-U(i,j-1,k,imz)*Q(i,j-1,k,qv)) &
                + D8(2)*(U(i,j+2,k,imz)*Q(i,j+2,k,qv)-U(i,j-2,k,imz)*Q(i,j-2,k,qv)) &
                + D8(3)*(U(i,j+3,k,imz)*Q(i,j+3,k,qv)-U(i,j-3,k,imz)*Q(i,j-3,k,qv)) &
                + D8(4)*(U(i,j+4,k,imz)*Q(i,j+4,k,qv)-U(i,j-4,k,imz)*Q(i,j-4,k,qv)) )

!EXPAND             Fhyp(i,j,k,iene)=Fhyp(i,j,k,iene) - dxinv(2) * &
!EXPAND                  first_deriv_8( (U(i,j-4:j+4,k,iene)+Q(i,j-4:j+4,k,qpres))*un )
             Fhyp(i,j,k,iene) = Fhyp(i,j,k,iene) - dxinv(2) * &
                ( D8(1)*((U(i,j+1,k,iene)+Q(i,j+1,k,qpres))*Q(i,j+1,k,qv)-(U(i,j-1,k,iene)+Q(i,j-1,k,qpres))*Q(i,j-1,k,qv)) &
                + D8(2)*((U(i,j+2,k,iene)+Q(i,j+2,k,qpres))*Q(i,j+2,k,qv)-(U(i,j-2,k,iene)+Q(i,j-2,k,qpres))*Q(i,j-2,k,qv)) &
                + D8(3)*((U(i,j+3,k,iene)+Q(i,j+3,k,qpres))*Q(i,j+3,k,qv)-(U(i,j-3,k,iene)+Q(i,j-3,k,qpres))*Q(i,j-3,k,qv)) &
                + D8(4)*((U(i,j+4,k,iene)+Q(i,j+4,k,qpres))*Q(i,j+4,k,qv)-(U(i,j-4,k,iene)+Q(i,j-4,k,qpres))*Q(i,j-4,k,qv)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(2) * &
!EXPAND                     first_deriv_8( U(i,j-4:j+4,k,n)*un )
                Fhyp(i,j,k,iryn) = Fhyp(i,j,k,iryn) - dxinv(2) * &
                   ( D8(1)*(U(i,j+1,k,iryn)*Q(i,j+1,k,qv)-U(i,j-1,k,iryn)*Q(i,j-1,k,qv)) &
                   + D8(2)*(U(i,j+2,k,iryn)*Q(i,j+2,k,qv)-U(i,j-2,k,iryn)*Q(i,j-2,k,qv)) &
                   + D8(3)*(U(i,j+3,k,iryn)*Q(i,j+3,k,qv)-U(i,j-3,k,iryn)*Q(i,j-3,k,qv)) &
                   + D8(4)*(U(i,j+4,k,iryn)*Q(i,j+4,k,qv)-U(i,j-4,k,iryn)*Q(i,j-4,k,qv)) )
             end do

          enddo
       enddo
    enddo
    !$omp end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

!             qwUN = Q(i,j,k-4:k+4,qw)

!EXPAND             Fhyp(i,j,k,irho)=Fhyp(i,j,k,irho) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imz) )
             Fhyp(i,j,k,irho) = Fhyp(i,j,k,irho) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imz)-U(i,j,k-1,imz)) &
                + D8(2)*(U(i,j,k+2,imz)-U(i,j,k-2,imz)) &
                + D8(3)*(U(i,j,k+3,imz)-U(i,j,k-3,imz)) &
                + D8(4)*(U(i,j,k+4,imz)-U(i,j,k-4,imz)) )

!EXPAND             Fhyp(i,j,k,imx)=Fhyp(i,j,k,imx) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imx)*qwUN )
             Fhyp(i,j,k,imx) = Fhyp(i,j,k,imx) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imx)*Q(i,j,k+1,qw)-U(i,j,k-1,imx)*Q(i,j,k-1,qw)) &
                + D8(2)*(U(i,j,k+2,imx)*Q(i,j,k+2,qw)-U(i,j,k-2,imx)*Q(i,j,k-2,qw)) &
                + D8(3)*(U(i,j,k+3,imx)*Q(i,j,k+3,qw)-U(i,j,k-3,imx)*Q(i,j,k-3,qw)) &
                + D8(4)*(U(i,j,k+4,imx)*Q(i,j,k+4,qw)-U(i,j,k-4,imx)*Q(i,j,k-4,qw)) )

!EXPAND             Fhyp(i,j,k,imy)=Fhyp(i,j,k,imy) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imy)*qwUN )
             Fhyp(i,j,k,imy) = Fhyp(i,j,k,imy) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imy)*Q(i,j,k+1,qw)-U(i,j,k-1,imy)*Q(i,j,k-1,qw)) &
                + D8(2)*(U(i,j,k+2,imy)*Q(i,j,k+2,qw)-U(i,j,k-2,imy)*Q(i,j,k-2,qw)) &
                + D8(3)*(U(i,j,k+3,imy)*Q(i,j,k+3,qw)-U(i,j,k-3,imy)*Q(i,j,k-3,qw)) &
                + D8(4)*(U(i,j,k+4,imy)*Q(i,j,k+4,qw)-U(i,j,k-4,imy)*Q(i,j,k-4,qw)) )

!EXPAND             Fhyp(i,j,k,imz)=Fhyp(i,j,k,imz) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imz)*qwUN+Q(i,j,k-4:k+4,qpres) )
             Fhyp(i,j,k,imz) = Fhyp(i,j,k,imz) - dxinv(3) * &
                ( D8(1)*((U(i,j,k+1,imz)*Q(i,j,k+1,qw)+Q(i,j,k+1,qpres))-(U(i,j,k-1,imz)*Q(i,j,k-1,qw)+Q(i,j,k-1,qpres))) &
                + D8(2)*((U(i,j,k+2,imz)*Q(i,j,k+2,qw)+Q(i,j,k+2,qpres))-(U(i,j,k-2,imz)*Q(i,j,k-2,qw)+Q(i,j,k-2,qpres))) &
                + D8(3)*((U(i,j,k+3,imz)*Q(i,j,k+3,qw)+Q(i,j,k+3,qpres))-(U(i,j,k-3,imz)*Q(i,j,k-3,qw)+Q(i,j,k-3,qpres))) &
                + D8(4)*((U(i,j,k+4,imz)*Q(i,j,k+4,qw)+Q(i,j,k+4,qpres))-(U(i,j,k-4,imz)*Q(i,j,k-4,qw)+Q(i,j,k-4,qpres))) )

!EXPAND             Fhyp(i,j,k,iene)=Fhyp(i,j,k,iene) - dxinv(3) * &
!EXPAND                  first_deriv_8( (U(i,j,k-4:k+4,iene)+Q(i,j,k-4:k+4,qpres))*qwUN )
             Fhyp(i,j,k,iene) = Fhyp(i,j,k,iene) - dxinv(3) * &
                ( D8(1)*((U(i,j,k+1,iene)+Q(i,j,k+1,qpres))*Q(i,j,k+1,qw)-(U(i,j,k-1,iene)+Q(i,j,k-1,qpres))*Q(i,j,k-1,qw)) &
                + D8(2)*((U(i,j,k+2,iene)+Q(i,j,k+2,qpres))*Q(i,j,k+2,qw)-(U(i,j,k-2,iene)+Q(i,j,k-2,qpres))*Q(i,j,k-2,qw)) &
                + D8(3)*((U(i,j,k+3,iene)+Q(i,j,k+3,qpres))*Q(i,j,k+3,qw)-(U(i,j,k-3,iene)+Q(i,j,k-3,qpres))*Q(i,j,k-3,qw)) &
                + D8(4)*((U(i,j,k+4,iene)+Q(i,j,k+4,qpres))*Q(i,j,k+4,qw)-(U(i,j,k-4,iene)+Q(i,j,k-4,qpres))*Q(i,j,k-4,qw)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(3) * &
!EXPAND                     first_deriv_8(U(i,j,k-4:k+4,n)*qwUN)
                Fhyp(i,j,k,iryn) = Fhyp(i,j,k,iryn) - dxinv(3) * &
                   ( D8(1)*(U(i,j,k+1,iryn)*Q(i,j,k+1,qw)-U(i,j,k-1,iryn)*Q(i,j,k-1,qw)) &
                   + D8(2)*(U(i,j,k+2,iryn)*Q(i,j,k+2,qw)-U(i,j,k-2,iryn)*Q(i,j,k-2,qw)) &
                   + D8(3)*(U(i,j,k+3,iryn)*Q(i,j,k+3,qw)-U(i,j,k-3,iryn)*Q(i,j,k-3,qw)) &
                   + D8(4)*(U(i,j,k+4,iryn)*Q(i,j,k+4,qw)-U(i,j,k-4,iryn)*Q(i,j,k-4,qw)) )
             end do

          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !
    ! Transport terms
    !
    !call build(bpt_diffterm, "diffterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call compact_diffterm_3d(lo,hi,ng,dx,qp,fdp,mup,xip,lamp,Ddp)


    allocate(vsm(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$omp parallel private(i,j,k,n,qxn,qryn,divu,tauxx,tauyy,tauzz) &
    !$omp   private(dmvxdy,dmwxdz,dmvywzdx,dmuydx,dmwydz,dmuxwzdy,dmuzdx,dmvzdy,dmuxvydz)

    !$omp workshare
    Fdif(:,:,:,irho) = 0.d0
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
!EXPAND             qx(i,j,k,idu) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qu) )
             qx(i,j,k,idu) = dxinv(1) * &
                ( D8(1)*(Q(i+1,j,k,qu)-Q(i-1,j,k,qu)) &
                + D8(2)*(Q(i+2,j,k,qu)-Q(i-2,j,k,qu)) &
                + D8(3)*(Q(i+3,j,k,qu)-Q(i-3,j,k,qu)) &
                + D8(4)*(Q(i+4,j,k,qu)-Q(i-4,j,k,qu)) )
!EXPAND             qx(i,j,k,idv) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qv) )
             qx(i,j,k,idv) = dxinv(1) * &
                ( D8(1)*(Q(i+1,j,k,qv)-Q(i-1,j,k,qv)) &
                + D8(2)*(Q(i+2,j,k,qv)-Q(i-2,j,k,qv)) &
                + D8(3)*(Q(i+3,j,k,qv)-Q(i-3,j,k,qv)) &
                + D8(4)*(Q(i+4,j,k,qv)-Q(i-4,j,k,qv)) )
!EXPAND             qx(i,j,k,idw) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qw) )
             qx(i,j,k,idw) = dxinv(1) * &
                ( D8(1)*(Q(i+1,j,k,qw)-Q(i-1,j,k,qw)) &
                + D8(2)*(Q(i+2,j,k,qw)-Q(i-2,j,k,qw)) &
                + D8(3)*(Q(i+3,j,k,qw)-Q(i-3,j,k,qw)) &
                + D8(4)*(Q(i+4,j,k,qw)-Q(i-4,j,k,qw)) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng
!EXPAND             qy(i,j,k,idu) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qu) )
             qy(i,j,k,idu) = dxinv(2) * &
                ( D8(1)*(Q(i,j+1,k,qu)-Q(i,j-1,k,qu)) &
                + D8(2)*(Q(i,j+2,k,qu)-Q(i,j-2,k,qu)) &
                + D8(3)*(Q(i,j+3,k,qu)-Q(i,j-3,k,qu)) &
                + D8(4)*(Q(i,j+4,k,qu)-Q(i,j-4,k,qu)) )
!EXPAND             qy(i,j,k,idv) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qv) )
             qy(i,j,k,idv) = dxinv(2) * &
                ( D8(1)*(Q(i,j+1,k,qv)-Q(i,j-1,k,qv)) &
                + D8(2)*(Q(i,j+2,k,qv)-Q(i,j-2,k,qv)) &
                + D8(3)*(Q(i,j+3,k,qv)-Q(i,j-3,k,qv)) &
                + D8(4)*(Q(i,j+4,k,qv)-Q(i,j-4,k,qv)) )
!EXPAND             qy(i,j,k,idw) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qw) )
             qy(i,j,k,idw) = dxinv(2) * &
                ( D8(1)*(Q(i,j+1,k,qw)-Q(i,j-1,k,qw)) &
                + D8(2)*(Q(i,j+2,k,qw)-Q(i,j-2,k,qw)) &
                + D8(3)*(Q(i,j+3,k,qw)-Q(i,j-3,k,qw)) &
                + D8(4)*(Q(i,j+4,k,qw)-Q(i,j-4,k,qw)) )
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
!EXPAND             qz(i,j,k,idu) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qu) )
             qz(i,j,k,idu) = dxinv(3) * &
                ( D8(1)*(Q(i,j,k+1,qu)-Q(i,j,k-1,qu)) &
                + D8(2)*(Q(i,j,k+2,qu)-Q(i,j,k-2,qu)) &
                + D8(3)*(Q(i,j,k+3,qu)-Q(i,j,k-3,qu)) &
                + D8(4)*(Q(i,j,k+4,qu)-Q(i,j,k-4,qu)) )
!EXPAND             qz(i,j,k,idv) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qv) )
             qz(i,j,k,idv) = dxinv(3) * &
                ( D8(1)*(Q(i,j,k+1,qv)-Q(i,j,k-1,qv)) &
                + D8(2)*(Q(i,j,k+2,qv)-Q(i,j,k-2,qv)) &
                + D8(3)*(Q(i,j,k+3,qv)-Q(i,j,k-3,qv)) &
                + D8(4)*(Q(i,j,k+4,qv)-Q(i,j,k-4,qv)) )
!EXPAND             qz(i,j,k,idw) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qw) )
             qz(i,j,k,idw) = dxinv(3) * &
                ( D8(1)*(Q(i,j,k+1,qw)-Q(i,j,k-1,qw)) &
                + D8(2)*(Q(i,j,k+2,qw)-Q(i,j,k-2,qw)) &
                + D8(3)*(Q(i,j,k+3,qw)-Q(i,j,k-3,qw)) &
                + D8(4)*(Q(i,j,k+4,qw)-Q(i,j,k-4,qw)) )
          enddo
       enddo
    enddo
    !$OMP END DO


    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! d(mu*dv/dx)/dy
!EXPAND             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qx(i,j-4:j+4,k,idv) )
             dmvxdy = dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*qx(i,j+1,k,idv)-mu(i,j-1,k)*qx(i,j-1,k,idv)) &
                + D8(2)*(mu(i,j+2,k)*qx(i,j+2,k,idv)-mu(i,j-2,k)*qx(i,j-2,k,idv)) &
                + D8(3)*(mu(i,j+3,k)*qx(i,j+3,k,idv)-mu(i,j-3,k)*qx(i,j-3,k,idv)) &
                + D8(4)*(mu(i,j+4,k)*qx(i,j+4,k,idv)-mu(i,j-4,k)*qx(i,j-4,k,idv)) )

             ! d(mu*dw/dx)/dz
!EXPAND             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qx(i,j,k-4:k+4,idw) )
             dmwxdz = dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*qx(i,j,k+1,idw)-mu(i,j,k-1)*qx(i,j,k-1,idw)) &
                + D8(2)*(mu(i,j,k+2)*qx(i,j,k+2,idw)-mu(i,j,k-2)*qx(i,j,k-2,idw)) &
                + D8(3)*(mu(i,j,k+3)*qx(i,j,k+3,idw)-mu(i,j,k-3)*qx(i,j,k-3,idw)) &
                + D8(4)*(mu(i,j,k+4)*qx(i,j,k+4,idw)-mu(i,j,k-4)*qx(i,j,k-4,idw)) )

             ! d((xi-2/3*mu)*(vy+wz))/dx
!EXPAND             dmvywzdx = dxinv(1) * &
!EXPAND                  first_deriv_8( vsm(i-4:i+4,j,k)*(qy(i-4:i+4,j,k,idv)+qz(i-4:i+4,j,k,idw)) )
             dmvywzdx = dxinv(1) * &
                ( D8(1)*(vsm(i+1,j,k)*(qy(i+1,j,k,idv)+qz(i+1,j,k,idw))-vsm(i-1,j,k)*(qy(i-1,j,k,idv)+qz(i-1,j,k,idw))) &
                + D8(2)*(vsm(i+2,j,k)*(qy(i+2,j,k,idv)+qz(i+2,j,k,idw))-vsm(i-2,j,k)*(qy(i-2,j,k,idv)+qz(i-2,j,k,idw))) &
                + D8(3)*(vsm(i+3,j,k)*(qy(i+3,j,k,idv)+qz(i+3,j,k,idw))-vsm(i-3,j,k)*(qy(i-3,j,k,idv)+qz(i-3,j,k,idw))) &
                + D8(4)*(vsm(i+4,j,k)*(qy(i+4,j,k,idv)+qz(i+4,j,k,idw))-vsm(i-4,j,k)*(qy(i-4,j,k,idv)+qz(i-4,j,k,idw))) )

             ! d(mu*du/dy)/dx
!EXPAND             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qy(i-4:i+4,j,k,idu) )
             dmuydx = dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*qy(i+1,j,k,idu)-mu(i-1,j,k)*qy(i-1,j,k,idu)) &
                + D8(2)*(mu(i+2,j,k)*qy(i+2,j,k,idu)-mu(i-2,j,k)*qy(i-2,j,k,idu)) &
                + D8(3)*(mu(i+3,j,k)*qy(i+3,j,k,idu)-mu(i-3,j,k)*qy(i-3,j,k,idu)) &
                + D8(4)*(mu(i+4,j,k)*qy(i+4,j,k,idu)-mu(i-4,j,k)*qy(i-4,j,k,idu)) )

             ! d(mu*dw/dy)/dz
!EXPAND             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qy(i,j,k-4:k+4,idw) )
             dmwydz = dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*qy(i,j,k+1,idw)-mu(i,j,k-1)*qy(i,j,k-1,idw)) &
                + D8(2)*(mu(i,j,k+2)*qy(i,j,k+2,idw)-mu(i,j,k-2)*qy(i,j,k-2,idw)) &
                + D8(3)*(mu(i,j,k+3)*qy(i,j,k+3,idw)-mu(i,j,k-3)*qy(i,j,k-3,idw)) &
                + D8(4)*(mu(i,j,k+4)*qy(i,j,k+4,idw)-mu(i,j,k-4)*qy(i,j,k-4,idw)) )

             ! d((xi-2/3*mu)*(ux+wz))/dy
!EXPAND             dmuxwzdy = dxinv(2) * &
!EXPAND                  first_deriv_8( vsm(i,j-4:j+4,k)*(qx(i,j-4:j+4,k,idu)+qz(i,j-4:j+4,k,idw)) )
             dmuxwzdy = dxinv(2) * &
                ( D8(1)*(vsm(i,j+1,k)*(qx(i,j+1,k,idu)+qz(i,j+1,k,idw))-vsm(i,j-1,k)*(qx(i,j-1,k,idu)+qz(i,j-1,k,idw))) &
                + D8(2)*(vsm(i,j+2,k)*(qx(i,j+2,k,idu)+qz(i,j+2,k,idw))-vsm(i,j-2,k)*(qx(i,j-2,k,idu)+qz(i,j-2,k,idw))) &
                + D8(3)*(vsm(i,j+3,k)*(qx(i,j+3,k,idu)+qz(i,j+3,k,idw))-vsm(i,j-3,k)*(qx(i,j-3,k,idu)+qz(i,j-3,k,idw))) &
                + D8(4)*(vsm(i,j+4,k)*(qx(i,j+4,k,idu)+qz(i,j+4,k,idw))-vsm(i,j-4,k)*(qx(i,j-4,k,idu)+qz(i,j-4,k,idw))) )

             ! d(mu*du/dz)/dx
!EXPAND             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qz(i-4:i+4,j,k,idu) )
             dmuzdx = dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*qz(i+1,j,k,idu)-mu(i-1,j,k)*qz(i-1,j,k,idu)) &
                + D8(2)*(mu(i+2,j,k)*qz(i+2,j,k,idu)-mu(i-2,j,k)*qz(i-2,j,k,idu)) &
                + D8(3)*(mu(i+3,j,k)*qz(i+3,j,k,idu)-mu(i-3,j,k)*qz(i-3,j,k,idu)) &
                + D8(4)*(mu(i+4,j,k)*qz(i+4,j,k,idu)-mu(i-4,j,k)*qz(i-4,j,k,idu)) )

             ! d(mu*dv/dz)/dy
!EXPAND             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qz(i,j-4:j+4,k,idv) )
             dmvzdy = dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*qz(i,j+1,k,idv)-mu(i,j-1,k)*qz(i,j-1,k,idv)) &
                + D8(2)*(mu(i,j+2,k)*qz(i,j+2,k,idv)-mu(i,j-2,k)*qz(i,j-2,k,idv)) &
                + D8(3)*(mu(i,j+3,k)*qz(i,j+3,k,idv)-mu(i,j-3,k)*qz(i,j-3,k,idv)) &
                + D8(4)*(mu(i,j+4,k)*qz(i,j+4,k,idv)-mu(i,j-4,k)*qz(i,j-4,k,idv)) )

             ! d((xi-2/3*mu)*(ux+vy))/dz
!EXPAND             dmuxvydz = dxinv(3) * &
!EXPAND                  first_deriv_8( vsm(i,j,k-4:k+4)*(qx(i,j,k-4:k+4,idu)+qy(i,j,k-4:k+4,idv)) )
             dmuxvydz = dxinv(3) * &
                ( D8(1)*(vsm(i,j,k+1)*(qx(i,j,k+1,idu)+qy(i,j,k+1,idv))-vsm(i,j,k-1)*(qx(i,j,k-1,idu)+qy(i,j,k-1,idv))) &
                + D8(2)*(vsm(i,j,k+2)*(qx(i,j,k+2,idu)+qy(i,j,k+2,idv))-vsm(i,j,k-2)*(qx(i,j,k-2,idu)+qy(i,j,k-2,idv))) &
                + D8(3)*(vsm(i,j,k+3)*(qx(i,j,k+3,idu)+qy(i,j,k+3,idv))-vsm(i,j,k-3)*(qx(i,j,k-3,idu)+qy(i,j,k-3,idv))) &
                + D8(4)*(vsm(i,j,k+4)*(qx(i,j,k+4,idu)+qy(i,j,k+4,idv))-vsm(i,j,k-4)*(qx(i,j,k-4,idu)+qy(i,j,k-4,idv))) )

             Fdif(i,j,k,imx) = dmvxdy + dmwxdz + dmvywzdx
             Fdif(i,j,k,imy) = dmuydx + dmwydz + dmuxwzdy
             Fdif(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu
             tauyy = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu
             tauzz = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu
             
             ! change in internal energy
             Fdif(i,j,k,iene) = tauxx*qx(i,j,k,idu) + tauyy*qy(i,j,k,idv) + tauzz*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    !$omp workshare
    Fdif(:,:,:,iry1:) = 0.d0
    !$omp end workshare

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

!EXPAND             qx(i,j,k,idT) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qtemp) )
             qx(i,j,k,idT) = dxinv(1) * &
                ( D8(1)*(Q(i+1,j,k,qtemp)-Q(i-1,j,k,qtemp)) &
                + D8(2)*(Q(i+2,j,k,qtemp)-Q(i-2,j,k,qtemp)) &
                + D8(3)*(Q(i+3,j,k,qtemp)-Q(i-3,j,k,qtemp)) &
                + D8(4)*(Q(i+4,j,k,qtemp)-Q(i-4,j,k,qtemp)) )
!EXPAND             qx(i,j,k,idp) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qpres) )
             qx(i,j,k,idp) = dxinv(1) * &
                ( D8(1)*(Q(i+1,j,k,qpres)-Q(i-1,j,k,qpres)) &
                + D8(2)*(Q(i+2,j,k,qpres)-Q(i-2,j,k,qpres)) &
                + D8(3)*(Q(i+3,j,k,qpres)-Q(i-3,j,k,qpres)) &
                + D8(4)*(Q(i+4,j,k,qpres)-Q(i-4,j,k,qpres)) )

!EXPAND             qy(i,j,k,idT) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qtemp) )
             qy(i,j,k,idT) = dxinv(2) * &
                ( D8(1)*(Q(i,j+1,k,qtemp)-Q(i,j-1,k,qtemp)) &
                + D8(2)*(Q(i,j+2,k,qtemp)-Q(i,j-2,k,qtemp)) &
                + D8(3)*(Q(i,j+3,k,qtemp)-Q(i,j-3,k,qtemp)) &
                + D8(4)*(Q(i,j+4,k,qtemp)-Q(i,j-4,k,qtemp)) )
!EXPAND             qy(i,j,k,idp) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qpres) )
             qy(i,j,k,idp) = dxinv(2) * &
                ( D8(1)*(Q(i,j+1,k,qpres)-Q(i,j-1,k,qpres)) &
                + D8(2)*(Q(i,j+2,k,qpres)-Q(i,j-2,k,qpres)) &
                + D8(3)*(Q(i,j+3,k,qpres)-Q(i,j-3,k,qpres)) &
                + D8(4)*(Q(i,j+4,k,qpres)-Q(i,j-4,k,qpres)) )

!EXPAND             qz(i,j,k,idT) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qtemp) )
             qz(i,j,k,idT) = dxinv(3) * &
                ( D8(1)*(Q(i,j,k+1,qtemp)-Q(i,j,k-1,qtemp)) &
                + D8(2)*(Q(i,j,k+2,qtemp)-Q(i,j,k-2,qtemp)) &
                + D8(3)*(Q(i,j,k+3,qtemp)-Q(i,j,k-3,qtemp)) &
                + D8(4)*(Q(i,j,k+4,qtemp)-Q(i,j,k-4,qtemp)) )
!EXPAND             qz(i,j,k,idp) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qpres) )
             qz(i,j,k,idp) = dxinv(3) * &
                ( D8(1)*(Q(i,j,k+1,qpres)-Q(i,j,k-1,qpres)) &
                + D8(2)*(Q(i,j,k+2,qpres)-Q(i,j,k-2,qpres)) &
                + D8(3)*(Q(i,j,k+3,qpres)-Q(i,j,k-3,qpres)) &
                + D8(4)*(Q(i,j,k+4,qpres)-Q(i,j,k-4,qpres)) )

          enddo
       enddo
    enddo
    !$omp end do nowait

    do n=1,nspecies
       qxn = qx1 + n - 1
       iryn = iry1 + n -1
       !$omp do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                qx(i,j,k,iryn) = dxinv(1) * first_deriv_8( Q(i-4:i+4,j,k,qxn) )
                qx(i,j,k,iryn) = dxinv(1) * &
                   ( D8(1)*(Q(i+1,j,k,qxn)-Q(i-1,j,k,qxn)) &
                   + D8(2)*(Q(i+2,j,k,qxn)-Q(i-2,j,k,qxn)) &
                   + D8(3)*(Q(i+3,j,k,qxn)-Q(i-3,j,k,qxn)) &
                   + D8(4)*(Q(i+4,j,k,qxn)-Q(i-4,j,k,qxn)) )
!EXPAND                qy(i,j,k,iryn) = dxinv(2) * first_deriv_8( Q(i,j-4:j+4,k,qxn) )
                qy(i,j,k,iryn) = dxinv(2) * &
                   ( D8(1)*(Q(i,j+1,k,qxn)-Q(i,j-1,k,qxn)) &
                   + D8(2)*(Q(i,j+2,k,qxn)-Q(i,j-2,k,qxn)) &
                   + D8(3)*(Q(i,j+3,k,qxn)-Q(i,j-3,k,qxn)) &
                   + D8(4)*(Q(i,j+4,k,qxn)-Q(i,j-4,k,qxn)) )
!EXPAND                qz(i,j,k,iryn) = dxinv(3) * first_deriv_8( Q(i,j,k-4:k+4,qxn) )
                qz(i,j,k,iryn) = dxinv(3) * &
                   ( D8(1)*(Q(i,j,k+1,qxn)-Q(i,j,k-1,qxn)) &
                   + D8(2)*(Q(i,j,k+2,qxn)-Q(i,j,k-2,qxn)) &
                   + D8(3)*(Q(i,j,k+3,qxn)-Q(i,j,k-3,qxn)) &
                   + D8(4)*(Q(i,j,k+4,qxn)-Q(i,j,k-4,qxn)) )
             enddo
          enddo
       enddo
       !$omp end do nowait
    enddo

    !$omp end parallel

    deallocate(vsm)

!  end subroutine S3D_diffterm_1


!  subroutine S3D_diffterm_2(lo,hi,ng,cons,dx,q,Fdif,mu,xi,lam,Ddiag,qx,qy,qz)

!    integer,          intent(in )  :: lo(3),hi(3),ng,cons
!    double precision, intent(in )  :: dx(3)
!    double precision, intent(in )  :: q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
!    double precision, intent(in )  :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
!    double precision, intent(in )  :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
!    double precision, intent(in )  :: lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)

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
                dpy(i,j,k,n) = Ddiag(i,j,k,n)/Q(i,j,k,qpres)*(Q(i,j,k,qxn)-Q(i,j,k,qyn))
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*Q(i,j,k,qhn)
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
!EXPAND             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) &
!EXPAND                  + dxinv(1) * first_deriv_8( vp(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idu)) 
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dxinv(1) * &
                ( D8(1)*(vp(i+1,j,k)*qx(i+1,j,k,idu)-vp(i-1,j,k)*qx(i-1,j,k,idu)) &
                + D8(2)*(vp(i+2,j,k)*qx(i+2,j,k,idu)-vp(i-2,j,k)*qx(i-2,j,k,idu)) &
                + D8(3)*(vp(i+3,j,k)*qx(i+3,j,k,idu)-vp(i-3,j,k)*qx(i-3,j,k,idu)) &
                + D8(4)*(vp(i+4,j,k)*qx(i+4,j,k,idu)-vp(i-4,j,k)*qx(i-4,j,k,idu)) )
!EXPAND             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) &
!EXPAND                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idu)) 
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*qy(i,j+1,k,idu)-mu(i,j-1,k)*qy(i,j-1,k,idu)) &
                + D8(2)*(mu(i,j+2,k)*qy(i,j+2,k,idu)-mu(i,j-2,k)*qy(i,j-2,k,idu)) &
                + D8(3)*(mu(i,j+3,k)*qy(i,j+3,k,idu)-mu(i,j-3,k)*qy(i,j-3,k,idu)) &
                + D8(4)*(mu(i,j+4,k)*qy(i,j+4,k,idu)-mu(i,j-4,k)*qy(i,j-4,k,idu)) )
!EXPAND             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) &
!EXPAND                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idu))
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*qz(i,j,k+1,idu)-mu(i,j,k-1)*qz(i,j,k-1,idu)) &
                + D8(2)*(mu(i,j,k+2)*qz(i,j,k+2,idu)-mu(i,j,k-2)*qz(i,j,k-2,idu)) &
                + D8(3)*(mu(i,j,k+3)*qz(i,j,k+3,idu)-mu(i,j,k-3)*qz(i,j,k-3,idu)) &
                + D8(4)*(mu(i,j,k+4)*qz(i,j,k+4,idu)-mu(i,j,k-4)*qz(i,j,k-4,idu)) )
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== my =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) &
!EXPAND                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idv) ) 
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*qx(i+1,j,k,idv)-mu(i-1,j,k)*qx(i-1,j,k,idv)) &
                + D8(2)*(mu(i+2,j,k)*qx(i+2,j,k,idv)-mu(i-2,j,k)*qx(i-2,j,k,idv)) &
                + D8(3)*(mu(i+3,j,k)*qx(i+3,j,k,idv)-mu(i-3,j,k)*qx(i-3,j,k,idv)) &
                + D8(4)*(mu(i+4,j,k)*qx(i+4,j,k,idv)-mu(i-4,j,k)*qx(i-4,j,k,idv)) )
!EXPAND             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) &
!EXPAND                  + dxinv(2) * first_deriv_8( vp(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idv) ) 
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dxinv(2) * &
                ( D8(1)*(vp(i,j+1,k)*qy(i,j+1,k,idv)-vp(i,j-1,k)*qy(i,j-1,k,idv)) &
                + D8(2)*(vp(i,j+2,k)*qy(i,j+2,k,idv)-vp(i,j-2,k)*qy(i,j-2,k,idv)) &
                + D8(3)*(vp(i,j+3,k)*qy(i,j+3,k,idv)-vp(i,j-3,k)*qy(i,j-3,k,idv)) &
                + D8(4)*(vp(i,j+4,k)*qy(i,j+4,k,idv)-vp(i,j-4,k)*qy(i,j-4,k,idv)) )
!EXPAND             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) &
!EXPAND                  + dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idv) ) 
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*qz(i,j,k+1,idv)-mu(i,j,k-1)*qz(i,j,k-1,idv)) &
                + D8(2)*(mu(i,j,k+2)*qz(i,j,k+2,idv)-mu(i,j,k-2)*qz(i,j,k-2,idv)) &
                + D8(3)*(mu(i,j,k+3)*qz(i,j,k+3,idv)-mu(i,j,k-3)*qz(i,j,k-3,idv)) &
                + D8(4)*(mu(i,j,k+4)*qz(i,j,k+4,idv)-mu(i,j,k-4)*qz(i,j,k-4,idv)) )
          end do
       end do
    end do
    !$omp end do nowait
    
    ! ===== mz =====
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) &
!EXPAND                  + dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idw) ) 
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*qx(i+1,j,k,idw)-mu(i-1,j,k)*qx(i-1,j,k,idw)) &
                + D8(2)*(mu(i+2,j,k)*qx(i+2,j,k,idw)-mu(i-2,j,k)*qx(i-2,j,k,idw)) &
                + D8(3)*(mu(i+3,j,k)*qx(i+3,j,k,idw)-mu(i-3,j,k)*qx(i-3,j,k,idw)) &
                + D8(4)*(mu(i+4,j,k)*qx(i+4,j,k,idw)-mu(i-4,j,k)*qx(i-4,j,k,idw)) )
!EXPAND             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) &
!EXPAND                  + dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idw) ) 
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*qy(i,j+1,k,idw)-mu(i,j-1,k)*qy(i,j-1,k,idw)) &
                + D8(2)*(mu(i,j+2,k)*qy(i,j+2,k,idw)-mu(i,j-2,k)*qy(i,j-2,k,idw)) &
                + D8(3)*(mu(i,j+3,k)*qy(i,j+3,k,idw)-mu(i,j-3,k)*qy(i,j-3,k,idw)) &
                + D8(4)*(mu(i,j+4,k)*qy(i,j+4,k,idw)-mu(i,j-4,k)*qy(i,j-4,k,idw)) )
!EXPAND             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) &
!EXPAND                  + dxinv(3) * first_deriv_8( vp(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idw) )
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dxinv(3) * &
                ( D8(1)*(vp(i,j,k+1)*qz(i,j,k+1,idw)-vp(i,j,k-1)*qz(i,j,k-1,idw)) &
                + D8(2)*(vp(i,j,k+2)*qz(i,j,k+2,idw)-vp(i,j,k-2)*qz(i,j,k-2,idw)) &
                + D8(3)*(vp(i,j,k+3)*qz(i,j,k+3,idw)-vp(i,j,k-3)*qz(i,j,k-3,idw)) &
                + D8(4)*(vp(i,j,k+4)*qz(i,j,k+4,idw)-vp(i,j,k-4)*qz(i,j,k-4,idw)) )
          end do
       end do
    end do
    !$omp end do
    
    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) &
                  + Fdif(i,j,k,imx)*Q(i,j,k,qu) &
                  + Fdif(i,j,k,imy)*Q(i,j,k,qv) &
                  + Fdif(i,j,k,imz)*Q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do

    ! thermal conduction
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) &
!EXPAND                  + dxinv(1) * first_deriv_8( lam(i-4:i+4,j,k)*qx(i-4:i+4,j,k,idT) ) 
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(1) * &
                ( D8(1)*(lam(i+1,j,k)*qx(i+1,j,k,idT)-lam(i-1,j,k)*qx(i-1,j,k,idT)) &
                + D8(2)*(lam(i+2,j,k)*qx(i+2,j,k,idT)-lam(i-2,j,k)*qx(i-2,j,k,idT)) &
                + D8(3)*(lam(i+3,j,k)*qx(i+3,j,k,idT)-lam(i-3,j,k)*qx(i-3,j,k,idT)) &
                + D8(4)*(lam(i+4,j,k)*qx(i+4,j,k,idT)-lam(i-4,j,k)*qx(i-4,j,k,idT)) )
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) &
!EXPAND                  + dxinv(2) * first_deriv_8( lam(i,j-4:j+4,k)*qy(i,j-4:j+4,k,idT) ) 
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(2) * &
                ( D8(1)*(lam(i,j+1,k)*qy(i,j+1,k,idT)-lam(i,j-1,k)*qy(i,j-1,k,idT)) &
                + D8(2)*(lam(i,j+2,k)*qy(i,j+2,k,idT)-lam(i,j-2,k)*qy(i,j-2,k,idT)) &
                + D8(3)*(lam(i,j+3,k)*qy(i,j+3,k,idT)-lam(i,j-3,k)*qy(i,j-3,k,idT)) &
                + D8(4)*(lam(i,j+4,k)*qy(i,j+4,k,idT)-lam(i,j-4,k)*qy(i,j-4,k,idT)) )
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) &
!EXPAND                  + dxinv(3) * first_deriv_8( lam(i,j,k-4:k+4)*qz(i,j,k-4:k+4,idT) )
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(3) * &
                ( D8(1)*(lam(i,j,k+1)*qz(i,j,k+1,idT)-lam(i,j,k-1)*qz(i,j,k-1,idT)) &
                + D8(2)*(lam(i,j,k+2)*qz(i,j,k+2,idT)-lam(i,j,k-2)*qz(i,j,k-2,idT)) &
                + D8(3)*(lam(i,j,k+3)*qz(i,j,k+3,idT)-lam(i,j,k-3)*qz(i,j,k-3,idT)) &
                + D8(4)*(lam(i,j,k+4)*qz(i,j,k+4,idT)-lam(i,j,k-4)*qz(i,j,k-4,idT)) )
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
                iryn = iry1 + n - 1

                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qx(i,j,k,iryn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qx(i,j,k,iryn)*Q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*Q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*Q(i,j,k,qyn)*Q(i,j,k,qhn)
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
!EXPAND                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + &
!EXPAND                     dxinv(1) * first_deriv_8( FY(i-4:i+4,j,k,n) )
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + dxinv(1) * &
                   ( D8(1)*(FY(i+1,j,k,n)-FY(i-1,j,k,n)) &
                   + D8(2)*(FY(i+2,j,k,n)-FY(i-2,j,k,n)) &
                   + D8(3)*(FY(i+3,j,k,n)-FY(i-3,j,k,n)) &
                   + D8(4)*(FY(i+4,j,k,n)-FY(i-4,j,k,n)) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + &
!EXPAND                  dxinv(1) * first_deriv_8( FE(i-4:i+4,j,k) )
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(1) * &
                ( D8(1)*(FE(i+1,j,k)-FE(i-1,j,k)) &
                + D8(2)*(FE(i+2,j,k)-FE(i-2,j,k)) &
                + D8(3)*(FE(i+3,j,k)-FE(i-3,j,k)) &
                + D8(4)*(FE(i+4,j,k)-FE(i-4,j,k)) )
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
                iryn = iry1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qy(i,j,k,iryn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qy(i,j,k,iryn)*Q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*Q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*Q(i,j,k,qyn)*Q(i,j,k,qhn)
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
!EXPAND                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + &
!EXPAND                     dxinv(2) * first_deriv_8( FY(i,j-4:j+4,k,n) )
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + dxinv(2) * &
                   ( D8(1)*(FY(i,j+1,k,n)-FY(i,j-1,k,n)) &
                   + D8(2)*(FY(i,j+2,k,n)-FY(i,j-2,k,n)) &
                   + D8(3)*(FY(i,j+3,k,n)-FY(i,j-3,k,n)) &
                   + D8(4)*(FY(i,j+4,k,n)-FY(i,j-4,k,n)) )
             end do
          end do
       end do
       !$omp end do nowait
    end do
    
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + &
!EXPAND                  dxinv(2) * first_deriv_8( FE(i,j-4:j+4,k) )
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(2) * &
                ( D8(1)*(FE(i,j+1,k)-FE(i,j-1,k)) &
                + D8(2)*(FE(i,j+2,k)-FE(i,j-2,k)) &
                + D8(3)*(FE(i,j+3,k)-FE(i,j-3,k)) &
                + D8(4)*(FE(i,j+4,k)-FE(i,j-4,k)) )
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
                iryn = iry1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qz(i,j,k,iryn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qz(i,j,k,iryn)*Q(i,j,k,qhn)
                rhoVc = rhoVc + FY(i,j,k,n)
             end do

             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = FY(i,j,k,n) - rhoVc*Q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rhoVc*Q(i,j,k,qyn)*Q(i,j,k,qhn)
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
!EXPAND                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + &
!EXPAND                     dxinv(3) * first_deriv_8( FY(i,j,k-4:k+4,n) )
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + dxinv(3) * &
                   ( D8(1)*(FY(i,j,k+1,n)-FY(i,j,k-1,n)) &
                   + D8(2)*(FY(i,j,k+2,n)-FY(i,j,k-2,n)) &
                   + D8(3)*(FY(i,j,k+3,n)-FY(i,j,k-3,n)) &
                   + D8(4)*(FY(i,j,k+4,n)-FY(i,j,k-4,n)) )
             end do
          end do
       end do
       !$omp end do nowait
    end do

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + &
!EXPAND                  dxinv(3) * first_deriv_8( FE(i,j,k-4:k+4) )
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + dxinv(3) * &
                ( D8(1)*(FE(i,j,k+1)-FE(i,j,k-1)) &
                + D8(2)*(FE(i,j,k+2)-FE(i,j,k-2)) &
                + D8(3)*(FE(i,j,k+3)-FE(i,j,k-3)) &
                + D8(4)*(FE(i,j,k+4)-FE(i,j,k-4)) )
          end do
       end do
    end do
    !$omp end do

    !$omp end parallel

    deallocate(vp,dpy,dpe,FY,FE)


    !call destroy(bpt_diffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !
    ! Calculate U'
    !
    !call build(bpt_calcU, "calcU")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Uprime(i,j,k,irho) = Fhyp(i,j,k,irho) + Fdif(i,j,k,irho)
             Uprime(i,j,k,imx ) = Fhyp(i,j,k,imx ) + Fdif(i,j,k,imx )
             Uprime(i,j,k,imy ) = Fhyp(i,j,k,imy ) + Fdif(i,j,k,imy )
             Uprime(i,j,k,imz ) = Fhyp(i,j,k,imz ) + Fdif(i,j,k,imz )
             Uprime(i,j,k,iene) = Fhyp(i,j,k,iene) + Fdif(i,j,k,iene)
             do n = 1, nspecies
                iryn = iry1+n-1
                Uprime(i,j,k,iryn) = Fhyp(i,j,k,iryn) + Fdif(i,j,k,iryn)
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !call destroy(bpt_calcU)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! 
    ! Add chemistry
    !
    !call build(bpt_chemterm, "chemterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    call chemterm_3d(lo,hi,ng,Q,Uprime, ncons, nprim, nspecies)
    

    !call destroy(bpt_chemterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !call destroy(Q)
    !call destroy(Fhyp)
    !call destroy(Fdif)

    !call destroy(mu)
    !call destroy(xi)
    !call destroy(lam)
    !call destroy(Ddiag)

  !end subroutine dUdt_compact

    !call destroy(bpt_rkstep1dUdt)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

!compute new time step dt 
    !call set_dt(dt, courno_proc, istep)

!    call update_rk3(Zero,Unew, One,U, dt,Uprime)
       !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Unew(i,j,k,irho) = Zero * Unew(i,j,k,irho) + One * U(i,j,k,irho) + dt * Uprime(i,j,k,irho)
             Unew(i,j,k,imx ) = Zero * Unew(i,j,k,imx ) + One * U(i,j,k,imx ) + dt * Uprime(i,j,k,imx )
             Unew(i,j,k,imy ) = Zero * Unew(i,j,k,imy ) + One * U(i,j,k,imy ) + dt * Uprime(i,j,k,imy )
             Unew(i,j,k,imz ) = Zero * Unew(i,j,k,imz ) + One * U(i,j,k,imz ) + dt * Uprime(i,j,k,imz )
             Unew(i,j,k,iene) = Zero * Unew(i,j,k,iene) + One * U(i,j,k,iene) + dt * Uprime(i,j,k,iene)
             do n = 1, nspecies
                iryn = iry1+n-1
                Unew(i,j,k,iryn) = Zero * Unew(i,j,k,iryn) + One * U(i,j,k,iryn) + dt * Uprime(i,j,k,iryn)
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    

    !call reset_density(Unew)
    !$omp parallel do private(i,j,k,n,rho)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rho = 0.d0
             do n=1, nspecies
                rho = rho + Unew(i,j,k,iry1+n-1)
             end do
             Unew(i,j,k,irho) = rho
          end do
       end do
    end do
    !$omp end parallel do

    !call build(bpt_rkstep1dUdt, "rkstep1dUdt")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call dUdt(U, Uprime, dx, courno=courno_proc)
    !call set_dt(dt, courno_proc, istep)
    !call update_rk3(Zero,Unew, One,U, dt,Uprime)
    !call reset_density(Unew)

    !call destroy(bpt_rkstep1)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 2
    !call build(bpt_rkstep2, "rkstep2")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call dUdt(Unew, Uprime, dx)
    !call update_rk3(OneQuarter, Unew, ThreeQuarters, U, OneQuarter*dt, Uprime)
    !call reset_density(Unew)
    !call destroy(bpt_rkstep2)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    ! RK Step 3
    !call build(bpt_rkstep3, "rkstep3")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call dUdt(Unew, Uprime, dx)
    !call update_rk3(OneThird, U, TwoThirds, Unew, TwoThirds*dt, Uprime)
    !call reset_density(U)
    !call destroy(bpt_rkstep3)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

    !call destroy(Unew)
    !call destroy(Uprime)

    !if (contains_nan(U)) then
    !   call bl_error("U contains nan")
    !end if

  end subroutine advance



end module advance_module
