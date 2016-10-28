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

    double precision, allocatable, dimension(:,:,:,:) :: Q, Ddiag
    double precision, allocatable, dimension(:,:,:) :: mu, xi, lam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable, dimension(:,:,:) :: vsp,vsm, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hgx, Hgy, Hgz, dpy, dxe

!    double precision :: un(-4:4)

    double precision :: dxinv(3), dx2inv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies), hhalf

    double precision :: rho
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
    allocate(mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(Ddiag(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))

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
    ! Transport terms
    !
    !call build(bpt_diffterm, "diffterm")   !! vvvvvvvvvvvvvvvvvvvvvvv timer
    !call compact_diffterm_3d(lo,hi,ng,dx,qp,fdp,mup,xip,lamp,Ddp)
    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i) ** 2
    end do

    dlo = lo - ng
    dhi = hi + ng

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
          do i=dlo(1),dhi(1)
             if (i >= lo(1) .and. i <= hi(1)) then
   !EXPAND             ux(i,j,k) = dxinv(1)*first_deriv_8(Q(i-4:i+4,j,k,qu))
                ux(i,j,k) = dxinv(1) * &
                   ( D8(1)*(Q(i+1,j,k,qu)-Q(i-1,j,k,qu)) &
                   + D8(2)*(Q(i+2,j,k,qu)-Q(i-2,j,k,qu)) &
                   + D8(3)*(Q(i+3,j,k,qu)-Q(i-3,j,k,qu)) &
                   + D8(4)*(Q(i+4,j,k,qu)-Q(i-4,j,k,qu)) )
   !EXPAND             vx(i,j,k) = dxinv(1)*first_deriv_8(Q(i-4:i+4,j,k,qv))
                vx(i,j,k) = dxinv(1) * &
                   ( D8(1)*(Q(i+1,j,k,qv)-Q(i-1,j,k,qv)) &
                   + D8(2)*(Q(i+2,j,k,qv)-Q(i-2,j,k,qv)) &
                   + D8(3)*(Q(i+3,j,k,qv)-Q(i-3,j,k,qv)) &
                   + D8(4)*(Q(i+4,j,k,qv)-Q(i-4,j,k,qv)) )
   !EXPAND             wx(i,j,k) = dxinv(1)*first_deriv_8(Q(i-4:i+4,j,k,qw))
                wx(i,j,k) = dxinv(1) * &
                   ( D8(1)*(Q(i+1,j,k,qw)-Q(i-1,j,k,qw)) &
                   + D8(2)*(Q(i+2,j,k,qw)-Q(i-2,j,k,qw)) &
                   + D8(3)*(Q(i+3,j,k,qw)-Q(i-3,j,k,qw)) &
                   + D8(4)*(Q(i+4,j,k,qw)-Q(i-4,j,k,qw)) )
             end if
             if (j >= lo(2) .and. j <= hi(2)) then
   !EXPAND             uy(i,j,k) = dxinv(2)*first_deriv_8(Q(i,j-4:j+4,k,qu))
                uy(i,j,k) = dxinv(2) * &
                   ( D8(1)*(Q(i,j+1,k,qu)-Q(i,j-1,k,qu)) &
                   + D8(2)*(Q(i,j+2,k,qu)-Q(i,j-2,k,qu)) &
                   + D8(3)*(Q(i,j+3,k,qu)-Q(i,j-3,k,qu)) &
                   + D8(4)*(Q(i,j+4,k,qu)-Q(i,j-4,k,qu)) )
   !EXPAND             vy(i,j,k) = dxinv(2)*first_deriv_8(Q(i,j-4:j+4,k,qv))
                vy(i,j,k) = dxinv(2) * &
                   ( D8(1)*(Q(i,j+1,k,qv)-Q(i,j-1,k,qv)) &
                   + D8(2)*(Q(i,j+2,k,qv)-Q(i,j-2,k,qv)) &
                   + D8(3)*(Q(i,j+3,k,qv)-Q(i,j-3,k,qv)) &
                   + D8(4)*(Q(i,j+4,k,qv)-Q(i,j-4,k,qv)) )
   !EXPAND             wy(i,j,k) = dxinv(2)*first_deriv_8(Q(i,j-4:j+4,k,qw))
                wy(i,j,k) = dxinv(2) * &
                   ( D8(1)*(Q(i,j+1,k,qw)-Q(i,j-1,k,qw)) &
                   + D8(2)*(Q(i,j+2,k,qw)-Q(i,j-2,k,qw)) &
                   + D8(3)*(Q(i,j+3,k,qw)-Q(i,j-3,k,qw)) &
                   + D8(4)*(Q(i,j+4,k,qw)-Q(i,j-4,k,qw)) )
             end if
             if (k >= lo(3) .and. k <= hi(3)) then
   !EXPAND             uz(i,j,k) = dxinv(3)*first_deriv_8(Q(i,j,k-4:k+4,qu))
                uz(i,j,k) = dxinv(3) * &
                   ( D8(1)*(Q(i,j,k+1,qu)-Q(i,j,k-1,qu)) &
                   + D8(2)*(Q(i,j,k+2,qu)-Q(i,j,k-2,qu)) &
                   + D8(3)*(Q(i,j,k+3,qu)-Q(i,j,k-3,qu)) &
                   + D8(4)*(Q(i,j,k+4,qu)-Q(i,j,k-4,qu)) )
   !EXPAND             vz(i,j,k) = dxinv(3)*first_deriv_8(Q(i,j,k-4:k+4,qv))
                vz(i,j,k) = dxinv(3) * &
                   ( D8(1)*(Q(i,j,k+1,qv)-Q(i,j,k-1,qv)) &
                   + D8(2)*(Q(i,j,k+2,qv)-Q(i,j,k-2,qv)) &
                   + D8(3)*(Q(i,j,k+3,qv)-Q(i,j,k-3,qv)) &
                   + D8(4)*(Q(i,j,k+4,qv)-Q(i,j,k-4,qv)) )
   !EXPAND             wz(i,j,k) = dxinv(3)*first_deriv_8(Q(i,j,k-4:k+4,qw))
                wz(i,j,k) = dxinv(3) * &
                   ( D8(1)*(Q(i,j,k+1,qw)-Q(i,j,k-1,qw)) &
                   + D8(2)*(Q(i,j,k+2,qw)-Q(i,j,k-2,qw)) &
                   + D8(3)*(Q(i,j,k+3,qw)-Q(i,j,k-3,qw)) &
                   + D8(4)*(Q(i,j,k+4,qw)-Q(i,j,k-4,qw)) )
             end if
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$omp barrier

    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ! d((xi-2/3*mu)*(vy+wz))/dx
!EXPAND             dmvywzdx = dxinv(1) * &
!EXPAND                  first_deriv_8( vsm(i-4:i+4,j,k)*(vy(i-4:i+4,j,k)+wz(i-4:i+4,j,k)) )
             dmvywzdx = dxinv(1) * &
                ( D8(1)*(vsm(i+1,j,k)*(vy(i+1,j,k)+wz(i+1,j,k))-vsm(i-1,j,k)*(vy(i-1,j,k)+wz(i-1,j,k))) &
                + D8(2)*(vsm(i+2,j,k)*(vy(i+2,j,k)+wz(i+2,j,k))-vsm(i-2,j,k)*(vy(i-2,j,k)+wz(i-2,j,k))) &
                + D8(3)*(vsm(i+3,j,k)*(vy(i+3,j,k)+wz(i+3,j,k))-vsm(i-3,j,k)*(vy(i-3,j,k)+wz(i-3,j,k))) &
                + D8(4)*(vsm(i+4,j,k)*(vy(i+4,j,k)+wz(i+4,j,k))-vsm(i-4,j,k)*(vy(i-4,j,k)+wz(i-4,j,k))) )
             ! d(mu*du/dy)/dx
!EXPAND             dmuydx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uy(i-4:i+4,j,k) )
             dmuydx = dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*uy(i+1,j,k)-mu(i-1,j,k)*uy(i-1,j,k)) &
                + D8(2)*(mu(i+2,j,k)*uy(i+2,j,k)-mu(i-2,j,k)*uy(i-2,j,k)) &
                + D8(3)*(mu(i+3,j,k)*uy(i+3,j,k)-mu(i-3,j,k)*uy(i-3,j,k)) &
                + D8(4)*(mu(i+4,j,k)*uy(i+4,j,k)-mu(i-4,j,k)*uy(i-4,j,k)) )
             ! d(mu*du/dz)/dx
!EXPAND             dmuzdx = dxinv(1) * first_deriv_8( mu(i-4:i+4,j,k)*uz(i-4:i+4,j,k) )
             dmuzdx = dxinv(1) * &
                ( D8(1)*(mu(i+1,j,k)*uz(i+1,j,k)-mu(i-1,j,k)*uz(i-1,j,k)) &
                + D8(2)*(mu(i+2,j,k)*uz(i+2,j,k)-mu(i-2,j,k)*uz(i-2,j,k)) &
                + D8(3)*(mu(i+3,j,k)*uz(i+3,j,k)-mu(i-3,j,k)*uz(i-3,j,k)) &
                + D8(4)*(mu(i+4,j,k)*uz(i+4,j,k)-mu(i-4,j,k)*uz(i-4,j,k)) )

             ! d(mu*dv/dx)/dy
!EXPAND             dmvxdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vx(i,j-4:j+4,k) )
             dmvxdy = dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*vx(i,j+1,k)-mu(i,j-1,k)*vx(i,j-1,k)) &
                + D8(2)*(mu(i,j+2,k)*vx(i,j+2,k)-mu(i,j-2,k)*vx(i,j-2,k)) &
                + D8(3)*(mu(i,j+3,k)*vx(i,j+3,k)-mu(i,j-3,k)*vx(i,j-3,k)) &
                + D8(4)*(mu(i,j+4,k)*vx(i,j+4,k)-mu(i,j-4,k)*vx(i,j-4,k)) )
             ! d((xi-2/3*mu)*(ux+wz))/dy
!EXPAND             dmuxwzdy = dxinv(2) * &
!EXPAND                  first_deriv_8( vsm(i,j-4:j+4,k)*(ux(i,j-4:j+4,k)+wz(i,j-4:j+4,k)) )
             dmuxwzdy = dxinv(2) * &
                ( D8(1)*(vsm(i,j+1,k)*(ux(i,j+1,k)+wz(i,j+1,k))-vsm(i,j-1,k)*(ux(i,j-1,k)+wz(i,j-1,k))) &
                + D8(2)*(vsm(i,j+2,k)*(ux(i,j+2,k)+wz(i,j+2,k))-vsm(i,j-2,k)*(ux(i,j-2,k)+wz(i,j-2,k))) &
                + D8(3)*(vsm(i,j+3,k)*(ux(i,j+3,k)+wz(i,j+3,k))-vsm(i,j-3,k)*(ux(i,j-3,k)+wz(i,j-3,k))) &
                + D8(4)*(vsm(i,j+4,k)*(ux(i,j+4,k)+wz(i,j+4,k))-vsm(i,j-4,k)*(ux(i,j-4,k)+wz(i,j-4,k))) )
             ! d(mu*dv/dz)/dy
!EXPAND             dmvzdy = dxinv(2) * first_deriv_8( mu(i,j-4:j+4,k)*vz(i,j-4:j+4,k) )
             dmvzdy = dxinv(2) * &
                ( D8(1)*(mu(i,j+1,k)*vz(i,j+1,k)-mu(i,j-1,k)*vz(i,j-1,k)) &
                + D8(2)*(mu(i,j+2,k)*vz(i,j+2,k)-mu(i,j-2,k)*vz(i,j-2,k)) &
                + D8(3)*(mu(i,j+3,k)*vz(i,j+3,k)-mu(i,j-3,k)*vz(i,j-3,k)) &
                + D8(4)*(mu(i,j+4,k)*vz(i,j+4,k)-mu(i,j-4,k)*vz(i,j-4,k)) )

             ! d(mu*dw/dx)/dz
!EXPAND             dmwxdz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wx(i,j,k-4:k+4) )
             dmwxdz = dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*wx(i,j,k+1)-mu(i,j,k-1)*wx(i,j,k-1)) &
                + D8(2)*(mu(i,j,k+2)*wx(i,j,k+2)-mu(i,j,k-2)*wx(i,j,k-2)) &
                + D8(3)*(mu(i,j,k+3)*wx(i,j,k+3)-mu(i,j,k-3)*wx(i,j,k-3)) &
                + D8(4)*(mu(i,j,k+4)*wx(i,j,k+4)-mu(i,j,k-4)*wx(i,j,k-4)) )
             ! d(mu*dw/dy)/dz
!EXPAND             dmwydz = dxinv(3) * first_deriv_8( mu(i,j,k-4:k+4)*wy(i,j,k-4:k+4) )
             dmwydz = dxinv(3) * &
                ( D8(1)*(mu(i,j,k+1)*wy(i,j,k+1)-mu(i,j,k-1)*wy(i,j,k-1)) &
                + D8(2)*(mu(i,j,k+2)*wy(i,j,k+2)-mu(i,j,k-2)*wy(i,j,k-2)) &
                + D8(3)*(mu(i,j,k+3)*wy(i,j,k+3)-mu(i,j,k-3)*wy(i,j,k-3)) &
                + D8(4)*(mu(i,j,k+4)*wy(i,j,k+4)-mu(i,j,k-4)*wy(i,j,k-4)) )
             ! d((xi-2/3*mu)*(ux+vy))/dz
!EXPAND             dmuxvydz = dxinv(3) * &
!EXPAND                  first_deriv_8( vsm(i,j,k-4:k+4)*(ux(i,j,k-4:k+4)+vy(i,j,k-4:k+4)) )
             dmuxvydz = dxinv(3) * &
                ( D8(1)*(vsm(i,j,k+1)*(ux(i,j,k+1)+vy(i,j,k+1))-vsm(i,j,k-1)*(ux(i,j,k-1)+vy(i,j,k-1))) &
                + D8(2)*(vsm(i,j,k+2)*(ux(i,j,k+2)+vy(i,j,k+2))-vsm(i,j,k-2)*(ux(i,j,k-2)+vy(i,j,k-2))) &
                + D8(3)*(vsm(i,j,k+3)*(ux(i,j,k+3)+vy(i,j,k+3))-vsm(i,j,k-3)*(ux(i,j,k-3)+vy(i,j,k-3))) &
                + D8(4)*(vsm(i,j,k+4)*(ux(i,j,k+4)+vy(i,j,k+4))-vsm(i,j,k-4)*(ux(i,j,k-4)+vy(i,j,k-4))) )

             Uprime(i,j,k,imx) = dmvywzdx + dmvxdy + dmwxdz
             Uprime(i,j,k,imy) = dmuydx + dmuxwzdy + dmwydz
             Uprime(i,j,k,imz) = dmuzdx + dmvzdy + dmuxvydz

             divu = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx = 2.d0*mu(i,j,k)*ux(i,j,k) + divu
             tauyy = 2.d0*mu(i,j,k)*vy(i,j,k) + divu
             tauzz = 2.d0*mu(i,j,k)*wz(i,j,k) + divu
             
             ! change in internal energy
             Uprime(i,j,k,iene) = tauxx*ux(i,j,k) + tauyy*vy(i,j,k) + tauzz*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do
    !$omp end do nowait

    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz)

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hgx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))
    allocate(Hgy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))
    allocate(Hgz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    !$omp parallel &
    !$omp private(i,j,k,n,qxn,qyn,qhn,Htot,Htmp,Ytmp,hhalf,M8p,M8X,mmtmp)

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             dpe(i,j,k) = 0.d0
             do n=1,nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                qhn = qh1+n-1
                dpy(i,j,k,n) = Ddiag(i,j,k,n)/Q(i,j,k,qpres)*(Q(i,j,k,qxn)-Q(i,j,k,qyn))
                dxe(i,j,k,n) = Ddiag(i,j,k,n)*Q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*Q(i,j,k,qhn)
             end do
          end do
       end do
    end do
    !$omp end do nowait

    !$omp barrier

    ! ------- BEGIN x-direction -------
    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             if (j <= hi(2) .and. k <= hi(3)) then
   !EXPAND             mmtmp = matmul(vsp(i-4:i+3,j,k), M8)
                mmtmp(1) = vsp(i-4,j,k) * M8(1,1) &
                         + vsp(i-3,j,k) * M8(2,1) &
                         + vsp(i-2,j,k) * M8(3,1) &
                         + vsp(i-1,j,k) * M8(4,1) &
                         + vsp(i  ,j,k) * M8(5,1)
                mmtmp(2) = vsp(i-4,j,k) * M8(1,2) &
                         + vsp(i-3,j,k) * M8(2,2) &
                         + vsp(i-2,j,k) * M8(3,2) &
                         + vsp(i-1,j,k) * M8(4,2) &
                         + vsp(i  ,j,k) * M8(5,2) &
                         + vsp(i+1,j,k) * M8(6,2)
                mmtmp(3) = vsp(i-4,j,k) * M8(1,3) &
                         + vsp(i-3,j,k) * M8(2,3) &
                         + vsp(i-2,j,k) * M8(3,3) &
                         + vsp(i-1,j,k) * M8(4,3) &
                         + vsp(i  ,j,k) * M8(5,3) &
                         + vsp(i+1,j,k) * M8(6,3) &
                         + vsp(i+2,j,k) * M8(7,3)
                mmtmp(4) = vsp(i-4,j,k) * M8(1,4) &
                         + vsp(i-3,j,k) * M8(2,4) &
                         + vsp(i-2,j,k) * M8(3,4) &
                         + vsp(i-1,j,k) * M8(4,4) &
                         + vsp(i  ,j,k) * M8(5,4) &
                         + vsp(i+1,j,k) * M8(6,4) &
                         + vsp(i+2,j,k) * M8(7,4) &
                         + vsp(i+3,j,k) * M8(8,4)
                mmtmp(5) =-vsp(i-4,j,k) * M8(8,4) &
                         - vsp(i-3,j,k) * M8(7,4) &
                         - vsp(i-2,j,k) * M8(6,4) &
                         - vsp(i-1,j,k) * M8(5,4) &
                         - vsp(i  ,j,k) * M8(4,4) &
                         - vsp(i+1,j,k) * M8(3,4) &
                         - vsp(i+2,j,k) * M8(2,4) &
                         - vsp(i+3,j,k) * M8(1,4)
                mmtmp(6) =-vsp(i-3,j,k) * M8(7,3) &
                         - vsp(i-2,j,k) * M8(6,3) &
                         - vsp(i-1,j,k) * M8(5,3) &
                         - vsp(i  ,j,k) * M8(4,3) &
                         - vsp(i+1,j,k) * M8(3,3) &
                         - vsp(i+2,j,k) * M8(2,3) &
                         - vsp(i+3,j,k) * M8(1,3)
                mmtmp(7) =-vsp(i-2,j,k) * M8(6,2) &
                         - vsp(i-1,j,k) * M8(5,2) &
                         - vsp(i  ,j,k) * M8(4,2) &
                         - vsp(i+1,j,k) * M8(3,2) &
                         - vsp(i+2,j,k) * M8(2,2) &
                         - vsp(i+3,j,k) * M8(1,2)
                mmtmp(8) =-vsp(i-1,j,k) * M8(5,1) &
                         - vsp(i  ,j,k) * M8(4,1) &
                         - vsp(i+1,j,k) * M8(3,1) &
                         - vsp(i+2,j,k) * M8(2,1) &
                         - vsp(i+3,j,k) * M8(1,1)
   !EXPAND             Hgx(i,j,k,imx) = dot_product(mmtmp, Q(i-4:i+3,j,k,qu))
                Hgx(i,j,k,imx) =  &
                   ( Q(i-4,j,k,qu)*mmtmp(1) + Q(i-3,j,k,qu)*mmtmp(2) &
                   + Q(i-2,j,k,qu)*mmtmp(3) + Q(i-1,j,k,qu)*mmtmp(4) &
                   + Q(i  ,j,k,qu)*mmtmp(5) + Q(i+1,j,k,qu)*mmtmp(6) &
                   + Q(i+2,j,k,qu)*mmtmp(7) + Q(i+3,j,k,qu)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(mu(i-4:i+3,j,k), M8)
                mmtmp(1) = mu(i-4,j,k) * M8(1,1) &
                         + mu(i-3,j,k) * M8(2,1) &
                         + mu(i-2,j,k) * M8(3,1) &
                         + mu(i-1,j,k) * M8(4,1) &
                         + mu(i  ,j,k) * M8(5,1)
                mmtmp(2) = mu(i-4,j,k) * M8(1,2) &
                         + mu(i-3,j,k) * M8(2,2) &
                         + mu(i-2,j,k) * M8(3,2) &
                         + mu(i-1,j,k) * M8(4,2) &
                         + mu(i  ,j,k) * M8(5,2) &
                         + mu(i+1,j,k) * M8(6,2)
                mmtmp(3) = mu(i-4,j,k) * M8(1,3) &
                         + mu(i-3,j,k) * M8(2,3) &
                         + mu(i-2,j,k) * M8(3,3) &
                         + mu(i-1,j,k) * M8(4,3) &
                         + mu(i  ,j,k) * M8(5,3) &
                         + mu(i+1,j,k) * M8(6,3) &
                         + mu(i+2,j,k) * M8(7,3)
                mmtmp(4) = mu(i-4,j,k) * M8(1,4) &
                         + mu(i-3,j,k) * M8(2,4) &
                         + mu(i-2,j,k) * M8(3,4) &
                         + mu(i-1,j,k) * M8(4,4) &
                         + mu(i  ,j,k) * M8(5,4) &
                         + mu(i+1,j,k) * M8(6,4) &
                         + mu(i+2,j,k) * M8(7,4) &
                         + mu(i+3,j,k) * M8(8,4)
                mmtmp(5) =-mu(i-4,j,k) * M8(8,4) &
                         - mu(i-3,j,k) * M8(7,4) &
                         - mu(i-2,j,k) * M8(6,4) &
                         - mu(i-1,j,k) * M8(5,4) &
                         - mu(i  ,j,k) * M8(4,4) &
                         - mu(i+1,j,k) * M8(3,4) &
                         - mu(i+2,j,k) * M8(2,4) &
                         - mu(i+3,j,k) * M8(1,4)
                mmtmp(6) =-mu(i-3,j,k) * M8(7,3) &
                         - mu(i-2,j,k) * M8(6,3) &
                         - mu(i-1,j,k) * M8(5,3) &
                         - mu(i  ,j,k) * M8(4,3) &
                         - mu(i+1,j,k) * M8(3,3) &
                         - mu(i+2,j,k) * M8(2,3) &
                         - mu(i+3,j,k) * M8(1,3)
                mmtmp(7) =-mu(i-2,j,k) * M8(6,2) &
                         - mu(i-1,j,k) * M8(5,2) &
                         - mu(i  ,j,k) * M8(4,2) &
                         - mu(i+1,j,k) * M8(3,2) &
                         - mu(i+2,j,k) * M8(2,2) &
                         - mu(i+3,j,k) * M8(1,2)
                mmtmp(8) =-mu(i-1,j,k) * M8(5,1) &
                         - mu(i  ,j,k) * M8(4,1) &
                         - mu(i+1,j,k) * M8(3,1) &
                         - mu(i+2,j,k) * M8(2,1) &
                         - mu(i+3,j,k) * M8(1,1)
   !EXPAND             Hgx(i,j,k,imy) = dot_product(mmtmp, Q(i-4:i+3,j,k,qv))
                Hgx(i,j,k,imy) =  &
                   ( Q(i-4,j,k,qv)*mmtmp(1) + Q(i-3,j,k,qv)*mmtmp(2) &
                   + Q(i-2,j,k,qv)*mmtmp(3) + Q(i-1,j,k,qv)*mmtmp(4) &
                   + Q(i  ,j,k,qv)*mmtmp(5) + Q(i+1,j,k,qv)*mmtmp(6) &
                   + Q(i+2,j,k,qv)*mmtmp(7) + Q(i+3,j,k,qv)*mmtmp(8) )
   !EXPAND             Hgx(i,j,k,imz) = dot_product(mmtmp, Q(i-4:i+3,j,k,qw))
                Hgx(i,j,k,imz) =  &
                   ( Q(i-4,j,k,qw)*mmtmp(1) + Q(i-3,j,k,qw)*mmtmp(2) &
                   + Q(i-2,j,k,qw)*mmtmp(3) + Q(i-1,j,k,qw)*mmtmp(4) &
                   + Q(i  ,j,k,qw)*mmtmp(5) + Q(i+1,j,k,qw)*mmtmp(6) &
                   + Q(i+2,j,k,qw)*mmtmp(7) + Q(i+3,j,k,qw)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(lam(i-4:i+3,j,k), M8)
                mmtmp(1) = lam(i-4,j,k) * M8(1,1) &
                         + lam(i-3,j,k) * M8(2,1) &
                         + lam(i-2,j,k) * M8(3,1) &
                         + lam(i-1,j,k) * M8(4,1) &
                         + lam(i  ,j,k) * M8(5,1)
                mmtmp(2) = lam(i-4,j,k) * M8(1,2) &
                         + lam(i-3,j,k) * M8(2,2) &
                         + lam(i-2,j,k) * M8(3,2) &
                         + lam(i-1,j,k) * M8(4,2) &
                         + lam(i  ,j,k) * M8(5,2) &
                         + lam(i+1,j,k) * M8(6,2)
                mmtmp(3) = lam(i-4,j,k) * M8(1,3) &
                         + lam(i-3,j,k) * M8(2,3) &
                         + lam(i-2,j,k) * M8(3,3) &
                         + lam(i-1,j,k) * M8(4,3) &
                         + lam(i  ,j,k) * M8(5,3) &
                         + lam(i+1,j,k) * M8(6,3) &
                         + lam(i+2,j,k) * M8(7,3)
                mmtmp(4) = lam(i-4,j,k) * M8(1,4) &
                         + lam(i-3,j,k) * M8(2,4) &
                         + lam(i-2,j,k) * M8(3,4) &
                         + lam(i-1,j,k) * M8(4,4) &
                         + lam(i  ,j,k) * M8(5,4) &
                         + lam(i+1,j,k) * M8(6,4) &
                         + lam(i+2,j,k) * M8(7,4) &
                         + lam(i+3,j,k) * M8(8,4)
                mmtmp(5) =-lam(i-4,j,k) * M8(8,4) &
                         - lam(i-3,j,k) * M8(7,4) &
                         - lam(i-2,j,k) * M8(6,4) &
                         - lam(i-1,j,k) * M8(5,4) &
                         - lam(i  ,j,k) * M8(4,4) &
                         - lam(i+1,j,k) * M8(3,4) &
                         - lam(i+2,j,k) * M8(2,4) &
                         - lam(i+3,j,k) * M8(1,4)
                mmtmp(6) =-lam(i-3,j,k) * M8(7,3) &
                         - lam(i-2,j,k) * M8(6,3) &
                         - lam(i-1,j,k) * M8(5,3) &
                         - lam(i  ,j,k) * M8(4,3) &
                         - lam(i+1,j,k) * M8(3,3) &
                         - lam(i+2,j,k) * M8(2,3) &
                         - lam(i+3,j,k) * M8(1,3)
                mmtmp(7) =-lam(i-2,j,k) * M8(6,2) &
                         - lam(i-1,j,k) * M8(5,2) &
                         - lam(i  ,j,k) * M8(4,2) &
                         - lam(i+1,j,k) * M8(3,2) &
                         - lam(i+2,j,k) * M8(2,2) &
                         - lam(i+3,j,k) * M8(1,2)
                mmtmp(8) =-lam(i-1,j,k) * M8(5,1) &
                         - lam(i  ,j,k) * M8(4,1) &
                         - lam(i+1,j,k) * M8(3,1) &
                         - lam(i+2,j,k) * M8(2,1) &
                         - lam(i+3,j,k) * M8(1,1)
   !EXPAND             M8p = matmul(M8, Q(i-4:i+3,j,k,qpres))
                M8p(1) = M8(1,1) * Q(i-4,j,k,qpres) &
                       + M8(1,2) * Q(i-3,j,k,qpres) &
                       + M8(1,3) * Q(i-2,j,k,qpres) &
                       + M8(1,4) * Q(i-1,j,k,qpres) &
                       - M8(8,4) * Q(i  ,j,k,qpres)
                M8p(2) = M8(2,1) * Q(i-4,j,k,qpres) &
                       + M8(2,2) * Q(i-3,j,k,qpres) &
                       + M8(2,3) * Q(i-2,j,k,qpres) &
                       + M8(2,4) * Q(i-1,j,k,qpres) &
                       - M8(7,4) * Q(i  ,j,k,qpres) &
                       - M8(7,3) * Q(i+1,j,k,qpres)
                M8p(3) = M8(3,1) * Q(i-4,j,k,qpres) &
                       + M8(3,2) * Q(i-3,j,k,qpres) &
                       + M8(3,3) * Q(i-2,j,k,qpres) &
                       + M8(3,4) * Q(i-1,j,k,qpres) &
                       - M8(6,4) * Q(i  ,j,k,qpres) &
                       - M8(6,3) * Q(i+1,j,k,qpres) &
                       - M8(6,2) * Q(i+2,j,k,qpres)
                M8p(4) = M8(4,1) * Q(i-4,j,k,qpres) &
                       + M8(4,2) * Q(i-3,j,k,qpres) &
                       + M8(4,3) * Q(i-2,j,k,qpres) &
                       + M8(4,4) * Q(i-1,j,k,qpres) &
                       - M8(5,4) * Q(i  ,j,k,qpres) &
                       - M8(5,3) * Q(i+1,j,k,qpres) &
                       - M8(5,2) * Q(i+2,j,k,qpres) &
                       - M8(5,1) * Q(i+3,j,k,qpres)
                M8p(5) = M8(5,1) * Q(i-4,j,k,qpres) &
                       + M8(5,2) * Q(i-3,j,k,qpres) &
                       + M8(5,3) * Q(i-2,j,k,qpres) &
                       + M8(5,4) * Q(i-1,j,k,qpres) &
                       - M8(4,4) * Q(i  ,j,k,qpres) &
                       - M8(4,3) * Q(i+1,j,k,qpres) &
                       - M8(4,2) * Q(i+2,j,k,qpres) &
                       - M8(4,1) * Q(i+3,j,k,qpres)
                M8p(6) = M8(6,2) * Q(i-3,j,k,qpres) &
                       + M8(6,3) * Q(i-2,j,k,qpres) &
                       + M8(6,4) * Q(i-1,j,k,qpres) &
                       - M8(3,4) * Q(i  ,j,k,qpres) &
                       - M8(3,3) * Q(i+1,j,k,qpres) &
                       - M8(3,2) * Q(i+2,j,k,qpres) &
                       - M8(3,1) * Q(i+3,j,k,qpres)
                M8p(7) = M8(7,3) * Q(i-2,j,k,qpres) &
                       + M8(7,4) * Q(i-1,j,k,qpres) &
                       - M8(2,4) * Q(i  ,j,k,qpres) &
                       - M8(2,3) * Q(i+1,j,k,qpres) &
                       - M8(2,2) * Q(i+2,j,k,qpres) &
                       - M8(2,1) * Q(i+3,j,k,qpres)
                M8p(8) = M8(8,4) * Q(i-1,j,k,qpres) &
                       - M8(1,4) * Q(i  ,j,k,qpres) &
                       - M8(1,3) * Q(i+1,j,k,qpres) &
                       - M8(1,2) * Q(i+2,j,k,qpres) &
                       - M8(1,1) * Q(i+3,j,k,qpres)
   !EXPAND             Hgx(i,j,k,iene) = dot_product(mmtmp, Q(i-4:i+3,j,k,qtemp)) &
   !EXPAND                  &         + dot_product(     dpe(i-4:i+3,j,k), M8p)
                Hgx(i,j,k,iene) =  &
                   ( ( Q(i-4,j,k,qtemp)*mmtmp(1) + Q(i-3,j,k,qtemp)*mmtmp(2) &
                     + Q(i-2,j,k,qtemp)*mmtmp(3) + Q(i-1,j,k,qtemp)*mmtmp(4) &
                     + Q(i  ,j,k,qtemp)*mmtmp(5) + Q(i+1,j,k,qtemp)*mmtmp(6) &
                     + Q(i+2,j,k,qtemp)*mmtmp(7) + Q(i+3,j,k,qtemp)*mmtmp(8) ) &
                   + ( dpe(i-4,j,k)*M8p(1) + dpe(i-3,j,k)*M8p(2) &
                     + dpe(i-2,j,k)*M8p(3) + dpe(i-1,j,k)*M8p(4) &
                     + dpe(i  ,j,k)*M8p(5) + dpe(i+1,j,k)*M8p(6) &
                     + dpe(i+2,j,k)*M8p(7) + dpe(i+3,j,k)*M8p(8) ) )

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

   !EXPAND                M8X = matmul(M8, Q(i-4:i+3,j,k,qxn))
                   M8X(1) = M8(1,1) * Q(i-4,j,k,qxn) &
                          + M8(1,2) * Q(i-3,j,k,qxn) &
                          + M8(1,3) * Q(i-2,j,k,qxn) &
                          + M8(1,4) * Q(i-1,j,k,qxn) &
                          - M8(8,4) * Q(i  ,j,k,qxn)
                   M8X(2) = M8(2,1) * Q(i-4,j,k,qxn) &
                          + M8(2,2) * Q(i-3,j,k,qxn) &
                          + M8(2,3) * Q(i-2,j,k,qxn) &
                          + M8(2,4) * Q(i-1,j,k,qxn) &
                          - M8(7,4) * Q(i  ,j,k,qxn) &
                          - M8(7,3) * Q(i+1,j,k,qxn)
                   M8X(3) = M8(3,1) * Q(i-4,j,k,qxn) &
                          + M8(3,2) * Q(i-3,j,k,qxn) &
                          + M8(3,3) * Q(i-2,j,k,qxn) &
                          + M8(3,4) * Q(i-1,j,k,qxn) &
                          - M8(6,4) * Q(i  ,j,k,qxn) &
                          - M8(6,3) * Q(i+1,j,k,qxn) &
                          - M8(6,2) * Q(i+2,j,k,qxn)
                   M8X(4) = M8(4,1) * Q(i-4,j,k,qxn) &
                          + M8(4,2) * Q(i-3,j,k,qxn) &
                          + M8(4,3) * Q(i-2,j,k,qxn) &
                          + M8(4,4) * Q(i-1,j,k,qxn) &
                          - M8(5,4) * Q(i  ,j,k,qxn) &
                          - M8(5,3) * Q(i+1,j,k,qxn) &
                          - M8(5,2) * Q(i+2,j,k,qxn) &
                          - M8(5,1) * Q(i+3,j,k,qxn)
                   M8X(5) = M8(5,1) * Q(i-4,j,k,qxn) &
                          + M8(5,2) * Q(i-3,j,k,qxn) &
                          + M8(5,3) * Q(i-2,j,k,qxn) &
                          + M8(5,4) * Q(i-1,j,k,qxn) &
                          - M8(4,4) * Q(i  ,j,k,qxn) &
                          - M8(4,3) * Q(i+1,j,k,qxn) &
                          - M8(4,2) * Q(i+2,j,k,qxn) &
                          - M8(4,1) * Q(i+3,j,k,qxn)
                   M8X(6) = M8(6,2) * Q(i-3,j,k,qxn) &
                          + M8(6,3) * Q(i-2,j,k,qxn) &
                          + M8(6,4) * Q(i-1,j,k,qxn) &
                          - M8(3,4) * Q(i  ,j,k,qxn) &
                          - M8(3,3) * Q(i+1,j,k,qxn) &
                          - M8(3,2) * Q(i+2,j,k,qxn) &
                          - M8(3,1) * Q(i+3,j,k,qxn)
                   M8X(7) = M8(7,3) * Q(i-2,j,k,qxn) &
                          + M8(7,4) * Q(i-1,j,k,qxn) &
                          - M8(2,4) * Q(i  ,j,k,qxn) &
                          - M8(2,3) * Q(i+1,j,k,qxn) &
                          - M8(2,2) * Q(i+2,j,k,qxn) &
                          - M8(2,1) * Q(i+3,j,k,qxn)
                   M8X(8) = M8(8,4) * Q(i-1,j,k,qxn) &
                          - M8(1,4) * Q(i  ,j,k,qxn) &
                          - M8(1,3) * Q(i+1,j,k,qxn) &
                          - M8(1,2) * Q(i+2,j,k,qxn) &
                          - M8(1,1) * Q(i+3,j,k,qxn)
                   
   !EXPAND                Htmp(n) = dot_product(dpy(i-4:i+3,j,k,n), M8p) &
   !EXPAND                     +    dot_product(Ddiag(i-4:i+3,j,k,n), M8X)
                   Htmp(n) =  &
                      ( ( dpy(i-4,j,k,n)*M8p(1) + dpy(i-3,j,k,n)*M8p(2) &
                        + dpy(i-2,j,k,n)*M8p(3) + dpy(i-1,j,k,n)*M8p(4) &
                        + dpy(i  ,j,k,n)*M8p(5) + dpy(i+1,j,k,n)*M8p(6) &
                        + dpy(i+2,j,k,n)*M8p(7) + dpy(i+3,j,k,n)*M8p(8) ) &
                      + ( Ddiag(i-4,j,k,n)*M8X(1) + Ddiag(i-3,j,k,n)*M8X(2) &
                        + Ddiag(i-2,j,k,n)*M8X(3) + Ddiag(i-1,j,k,n)*M8X(4) &
                        + Ddiag(i  ,j,k,n)*M8X(5) + Ddiag(i+1,j,k,n)*M8X(6) &
                        + Ddiag(i+2,j,k,n)*M8X(7) + Ddiag(i+3,j,k,n)*M8X(8) ) )

   !EXPAND                Hgx(i,j,k,iene) = Hgx(i,j,k,iene) &
   !EXPAND                     +    dot_product(dxe(i-4:i+3,j,k,n), M8X)
                   Hgx(i,j,k,iene) = Hgx(i,j,k,iene)+ &
                      ( dxe(i-4,j,k,n)*M8X(1) + dxe(i-3,j,k,n)*M8X(2) &
                      + dxe(i-2,j,k,n)*M8X(3) + dxe(i-1,j,k,n)*M8X(4) &
                      + dxe(i  ,j,k,n)*M8X(5) + dxe(i+1,j,k,n)*M8X(6) &
                      + dxe(i+2,j,k,n)*M8X(7) + dxe(i+3,j,k,n)*M8X(8) )

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (Q(i-1,j,k,qyn) + Q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hgx(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (Q(i-1,j,k,qhn) + Q(i,j,k,qhn)) / 2.d0
                   Hgx(i,j,k,iene) =  Hgx(i,j,k,iene) - Ytmp(n) * hhalf * Htot
                end do

                ! add x-direction Fdif
                if (i >= lo(1)+1) then
                   Uprime(i-1,j,k,imx ) = Uprime(i-1,j,k,imx ) + (Hgx(i,j,k,imx ) - Hgx(i-1,j,k,imx )) * dx2inv(1)
                   Uprime(i-1,j,k,imy ) = Uprime(i-1,j,k,imy ) + (Hgx(i,j,k,imy ) - Hgx(i-1,j,k,imy )) * dx2inv(1)
                   Uprime(i-1,j,k,imz ) = Uprime(i-1,j,k,imz ) + (Hgx(i,j,k,imz ) - Hgx(i-1,j,k,imz )) * dx2inv(1)
                   Uprime(i-1,j,k,iene) = Uprime(i-1,j,k,iene) + (Hgx(i,j,k,iene) - Hgx(i-1,j,k,iene)) * dx2inv(1)
                   do n=1,nspecies
                      iryn = iry1+n-1
                      Uprime(i-1,j,k,iryn) = (Hgx(i,j,k,iryn) - Hgx(i-1,j,k,iryn)) * dx2inv(1)
                   end do
                end if
             end if
    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------
             if (i <= hi(1) .and. k <= hi(3)) then
   !EXPAND             mmtmp = matmul(mu(i,j-4:j+3,k), M8)
                mmtmp(1) = mu(i,j-4,k) * M8(1,1) &
                         + mu(i,j-3,k) * M8(2,1) &
                         + mu(i,j-2,k) * M8(3,1) &
                         + mu(i,j-1,k) * M8(4,1) &
                         + mu(i,j  ,k) * M8(5,1)
                mmtmp(2) = mu(i,j-4,k) * M8(1,2) &
                         + mu(i,j-3,k) * M8(2,2) &
                         + mu(i,j-2,k) * M8(3,2) &
                         + mu(i,j-1,k) * M8(4,2) &
                         + mu(i,j  ,k) * M8(5,2) &
                         + mu(i,j+1,k) * M8(6,2)
                mmtmp(3) = mu(i,j-4,k) * M8(1,3) &
                         + mu(i,j-3,k) * M8(2,3) &
                         + mu(i,j-2,k) * M8(3,3) &
                         + mu(i,j-1,k) * M8(4,3) &
                         + mu(i,j  ,k) * M8(5,3) &
                         + mu(i,j+1,k) * M8(6,3) &
                         + mu(i,j+2,k) * M8(7,3)
                mmtmp(4) = mu(i,j-4,k) * M8(1,4) &
                         + mu(i,j-3,k) * M8(2,4) &
                         + mu(i,j-2,k) * M8(3,4) &
                         + mu(i,j-1,k) * M8(4,4) &
                         + mu(i,j  ,k) * M8(5,4) &
                         + mu(i,j+1,k) * M8(6,4) &
                         + mu(i,j+2,k) * M8(7,4) &
                         + mu(i,j+3,k) * M8(8,4)
                mmtmp(5) =-mu(i,j-4,k) * M8(8,4) &
                         - mu(i,j-3,k) * M8(7,4) &
                         - mu(i,j-2,k) * M8(6,4) &
                         - mu(i,j-1,k) * M8(5,4) &
                         - mu(i,j  ,k) * M8(4,4) &
                         - mu(i,j+1,k) * M8(3,4) &
                         - mu(i,j+2,k) * M8(2,4) &
                         - mu(i,j+3,k) * M8(1,4)
                mmtmp(6) =-mu(i,j-3,k) * M8(7,3) &
                         - mu(i,j-2,k) * M8(6,3) &
                         - mu(i,j-1,k) * M8(5,3) &
                         - mu(i,j  ,k) * M8(4,3) &
                         - mu(i,j+1,k) * M8(3,3) &
                         - mu(i,j+2,k) * M8(2,3) &
                         - mu(i,j+3,k) * M8(1,3)
                mmtmp(7) =-mu(i,j-2,k) * M8(6,2) &
                         - mu(i,j-1,k) * M8(5,2) &
                         - mu(i,j  ,k) * M8(4,2) &
                         - mu(i,j+1,k) * M8(3,2) &
                         - mu(i,j+2,k) * M8(2,2) &
                         - mu(i,j+3,k) * M8(1,2)
                mmtmp(8) =-mu(i,j-1,k) * M8(5,1) &
                         - mu(i,j  ,k) * M8(4,1) &
                         - mu(i,j+1,k) * M8(3,1) &
                         - mu(i,j+2,k) * M8(2,1) &
                         - mu(i,j+3,k) * M8(1,1)
   !EXPAND             Hgy(i,j,k,imx) = dot_product(mmtmp, Q(i,j-4:j+3,k,qu))
                Hgy(i,j,k,imx) =  &
                   ( Q(i,j-4,k,qu)*mmtmp(1) + Q(i,j-3,k,qu)*mmtmp(2) &
                   + Q(i,j-2,k,qu)*mmtmp(3) + Q(i,j-1,k,qu)*mmtmp(4) &
                   + Q(i,j  ,k,qu)*mmtmp(5) + Q(i,j+1,k,qu)*mmtmp(6) &
                   + Q(i,j+2,k,qu)*mmtmp(7) + Q(i,j+3,k,qu)*mmtmp(8) )
   !EXPAND             Hgy(i,j,k,imz) = dot_product(mmtmp, Q(i,j-4:j+3,k,qw))
                Hgy(i,j,k,imz) =  &
                   ( Q(i,j-4,k,qw)*mmtmp(1) + Q(i,j-3,k,qw)*mmtmp(2) &
                   + Q(i,j-2,k,qw)*mmtmp(3) + Q(i,j-1,k,qw)*mmtmp(4) &
                   + Q(i,j  ,k,qw)*mmtmp(5) + Q(i,j+1,k,qw)*mmtmp(6) &
                   + Q(i,j+2,k,qw)*mmtmp(7) + Q(i,j+3,k,qw)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(vsp(i,j-4:j+3,k), M8)
                mmtmp(1) = vsp(i,j-4,k) * M8(1,1) &
                         + vsp(i,j-3,k) * M8(2,1) &
                         + vsp(i,j-2,k) * M8(3,1) &
                         + vsp(i,j-1,k) * M8(4,1) &
                         + vsp(i,j  ,k) * M8(5,1)
                mmtmp(2) = vsp(i,j-4,k) * M8(1,2) &
                         + vsp(i,j-3,k) * M8(2,2) &
                         + vsp(i,j-2,k) * M8(3,2) &
                         + vsp(i,j-1,k) * M8(4,2) &
                         + vsp(i,j  ,k) * M8(5,2) &
                         + vsp(i,j+1,k) * M8(6,2)
                mmtmp(3) = vsp(i,j-4,k) * M8(1,3) &
                         + vsp(i,j-3,k) * M8(2,3) &
                         + vsp(i,j-2,k) * M8(3,3) &
                         + vsp(i,j-1,k) * M8(4,3) &
                         + vsp(i,j  ,k) * M8(5,3) &
                         + vsp(i,j+1,k) * M8(6,3) &
                         + vsp(i,j+2,k) * M8(7,3)
                mmtmp(4) = vsp(i,j-4,k) * M8(1,4) &
                         + vsp(i,j-3,k) * M8(2,4) &
                         + vsp(i,j-2,k) * M8(3,4) &
                         + vsp(i,j-1,k) * M8(4,4) &
                         + vsp(i,j  ,k) * M8(5,4) &
                         + vsp(i,j+1,k) * M8(6,4) &
                         + vsp(i,j+2,k) * M8(7,4) &
                         + vsp(i,j+3,k) * M8(8,4)
                mmtmp(5) =-vsp(i,j-4,k) * M8(8,4) &
                         - vsp(i,j-3,k) * M8(7,4) &
                         - vsp(i,j-2,k) * M8(6,4) &
                         - vsp(i,j-1,k) * M8(5,4) &
                         - vsp(i,j  ,k) * M8(4,4) &
                         - vsp(i,j+1,k) * M8(3,4) &
                         - vsp(i,j+2,k) * M8(2,4) &
                         - vsp(i,j+3,k) * M8(1,4)
                mmtmp(6) =-vsp(i,j-3,k) * M8(7,3) &
                         - vsp(i,j-2,k) * M8(6,3) &
                         - vsp(i,j-1,k) * M8(5,3) &
                         - vsp(i,j  ,k) * M8(4,3) &
                         - vsp(i,j+1,k) * M8(3,3) &
                         - vsp(i,j+2,k) * M8(2,3) &
                         - vsp(i,j+3,k) * M8(1,3)
                mmtmp(7) =-vsp(i,j-2,k) * M8(6,2) &
                         - vsp(i,j-1,k) * M8(5,2) &
                         - vsp(i,j  ,k) * M8(4,2) &
                         - vsp(i,j+1,k) * M8(3,2) &
                         - vsp(i,j+2,k) * M8(2,2) &
                         - vsp(i,j+3,k) * M8(1,2)
                mmtmp(8) =-vsp(i,j-1,k) * M8(5,1) &
                         - vsp(i,j  ,k) * M8(4,1) &
                         - vsp(i,j+1,k) * M8(3,1) &
                         - vsp(i,j+2,k) * M8(2,1) &
                         - vsp(i,j+3,k) * M8(1,1)
   !EXPAND             Hgy(i,j,k,imy) = dot_product(mmtmp, Q(i,j-4:j+3,k,qv))
                Hgy(i,j,k,imy) =  &
                   ( Q(i,j-4,k,qv)*mmtmp(1) + Q(i,j-3,k,qv)*mmtmp(2) &
                   + Q(i,j-2,k,qv)*mmtmp(3) + Q(i,j-1,k,qv)*mmtmp(4) &
                   + Q(i,j  ,k,qv)*mmtmp(5) + Q(i,j+1,k,qv)*mmtmp(6) &
                   + Q(i,j+2,k,qv)*mmtmp(7) + Q(i,j+3,k,qv)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(lam(i,j-4:j+3,k), M8)
                mmtmp(1) = lam(i,j-4,k) * M8(1,1) &
                         + lam(i,j-3,k) * M8(2,1) &
                         + lam(i,j-2,k) * M8(3,1) &
                         + lam(i,j-1,k) * M8(4,1) &
                         + lam(i,j  ,k) * M8(5,1)
                mmtmp(2) = lam(i,j-4,k) * M8(1,2) &
                         + lam(i,j-3,k) * M8(2,2) &
                         + lam(i,j-2,k) * M8(3,2) &
                         + lam(i,j-1,k) * M8(4,2) &
                         + lam(i,j  ,k) * M8(5,2) &
                         + lam(i,j+1,k) * M8(6,2)
                mmtmp(3) = lam(i,j-4,k) * M8(1,3) &
                         + lam(i,j-3,k) * M8(2,3) &
                         + lam(i,j-2,k) * M8(3,3) &
                         + lam(i,j-1,k) * M8(4,3) &
                         + lam(i,j  ,k) * M8(5,3) &
                         + lam(i,j+1,k) * M8(6,3) &
                         + lam(i,j+2,k) * M8(7,3)
                mmtmp(4) = lam(i,j-4,k) * M8(1,4) &
                         + lam(i,j-3,k) * M8(2,4) &
                         + lam(i,j-2,k) * M8(3,4) &
                         + lam(i,j-1,k) * M8(4,4) &
                         + lam(i,j  ,k) * M8(5,4) &
                         + lam(i,j+1,k) * M8(6,4) &
                         + lam(i,j+2,k) * M8(7,4) &
                         + lam(i,j+3,k) * M8(8,4)
                mmtmp(5) =-lam(i,j-4,k) * M8(8,4) &
                         - lam(i,j-3,k) * M8(7,4) &
                         - lam(i,j-2,k) * M8(6,4) &
                         - lam(i,j-1,k) * M8(5,4) &
                         - lam(i,j  ,k) * M8(4,4) &
                         - lam(i,j+1,k) * M8(3,4) &
                         - lam(i,j+2,k) * M8(2,4) &
                         - lam(i,j+3,k) * M8(1,4)
                mmtmp(6) =-lam(i,j-3,k) * M8(7,3) &
                         - lam(i,j-2,k) * M8(6,3) &
                         - lam(i,j-1,k) * M8(5,3) &
                         - lam(i,j  ,k) * M8(4,3) &
                         - lam(i,j+1,k) * M8(3,3) &
                         - lam(i,j+2,k) * M8(2,3) &
                         - lam(i,j+3,k) * M8(1,3)
                mmtmp(7) =-lam(i,j-2,k) * M8(6,2) &
                         - lam(i,j-1,k) * M8(5,2) &
                         - lam(i,j  ,k) * M8(4,2) &
                         - lam(i,j+1,k) * M8(3,2) &
                         - lam(i,j+2,k) * M8(2,2) &
                         - lam(i,j+3,k) * M8(1,2)
                mmtmp(8) =-lam(i,j-1,k) * M8(5,1) &
                         - lam(i,j  ,k) * M8(4,1) &
                         - lam(i,j+1,k) * M8(3,1) &
                         - lam(i,j+2,k) * M8(2,1) &
                         - lam(i,j+3,k) * M8(1,1)
   !EXPAND             M8p = matmul(M8, Q(i,j-4:j+3,k,qpres))
                M8p(1) = M8(1,1) * Q(i,j-4,k,qpres) &
                       + M8(1,2) * Q(i,j-3,k,qpres) &
                       + M8(1,3) * Q(i,j-2,k,qpres) &
                       + M8(1,4) * Q(i,j-1,k,qpres) &
                       - M8(8,4) * Q(i,j  ,k,qpres)
                M8p(2) = M8(2,1) * Q(i,j-4,k,qpres) &
                       + M8(2,2) * Q(i,j-3,k,qpres) &
                       + M8(2,3) * Q(i,j-2,k,qpres) &
                       + M8(2,4) * Q(i,j-1,k,qpres) &
                       - M8(7,4) * Q(i,j  ,k,qpres) &
                       - M8(7,3) * Q(i,j+1,k,qpres)
                M8p(3) = M8(3,1) * Q(i,j-4,k,qpres) &
                       + M8(3,2) * Q(i,j-3,k,qpres) &
                       + M8(3,3) * Q(i,j-2,k,qpres) &
                       + M8(3,4) * Q(i,j-1,k,qpres) &
                       - M8(6,4) * Q(i,j  ,k,qpres) &
                       - M8(6,3) * Q(i,j+1,k,qpres) &
                       - M8(6,2) * Q(i,j+2,k,qpres)
                M8p(4) = M8(4,1) * Q(i,j-4,k,qpres) &
                       + M8(4,2) * Q(i,j-3,k,qpres) &
                       + M8(4,3) * Q(i,j-2,k,qpres) &
                       + M8(4,4) * Q(i,j-1,k,qpres) &
                       - M8(5,4) * Q(i,j  ,k,qpres) &
                       - M8(5,3) * Q(i,j+1,k,qpres) &
                       - M8(5,2) * Q(i,j+2,k,qpres) &
                       - M8(5,1) * Q(i,j+3,k,qpres)
                M8p(5) = M8(5,1) * Q(i,j-4,k,qpres) &
                       + M8(5,2) * Q(i,j-3,k,qpres) &
                       + M8(5,3) * Q(i,j-2,k,qpres) &
                       + M8(5,4) * Q(i,j-1,k,qpres) &
                       - M8(4,4) * Q(i,j  ,k,qpres) &
                       - M8(4,3) * Q(i,j+1,k,qpres) &
                       - M8(4,2) * Q(i,j+2,k,qpres) &
                       - M8(4,1) * Q(i,j+3,k,qpres)
                M8p(6) = M8(6,2) * Q(i,j-3,k,qpres) &
                       + M8(6,3) * Q(i,j-2,k,qpres) &
                       + M8(6,4) * Q(i,j-1,k,qpres) &
                       - M8(3,4) * Q(i,j  ,k,qpres) &
                       - M8(3,3) * Q(i,j+1,k,qpres) &
                       - M8(3,2) * Q(i,j+2,k,qpres) &
                       - M8(3,1) * Q(i,j+3,k,qpres)
                M8p(7) = M8(7,3) * Q(i,j-2,k,qpres) &
                       + M8(7,4) * Q(i,j-1,k,qpres) &
                       - M8(2,4) * Q(i,j  ,k,qpres) &
                       - M8(2,3) * Q(i,j+1,k,qpres) &
                       - M8(2,2) * Q(i,j+2,k,qpres) &
                       - M8(2,1) * Q(i,j+3,k,qpres)
                M8p(8) = M8(8,4) * Q(i,j-1,k,qpres) &
                       - M8(1,4) * Q(i,j  ,k,qpres) &
                       - M8(1,3) * Q(i,j+1,k,qpres) &
                       - M8(1,2) * Q(i,j+2,k,qpres) &
                       - M8(1,1) * Q(i,j+3,k,qpres)
   !EXPAND             Hgy(i,j,k,iene) = dot_product(mmtmp, Q(i,j-4:j+3,k,qtemp)) &
   !EXPAND                  +           dot_product(     dpe(i,j-4:j+3,k), M8p)
                Hgy(i,j,k,iene) =  &
                   ( ( Q(i,j-4,k,qtemp)*mmtmp(1) + Q(i,j-3,k,qtemp)*mmtmp(2) &
                     + Q(i,j-2,k,qtemp)*mmtmp(3) + Q(i,j-1,k,qtemp)*mmtmp(4) &
                     + Q(i,j  ,k,qtemp)*mmtmp(5) + Q(i,j+1,k,qtemp)*mmtmp(6) &
                     + Q(i,j+2,k,qtemp)*mmtmp(7) + Q(i,j+3,k,qtemp)*mmtmp(8) ) &
                   + ( dpe(i,j-4,k)*M8p(1) + dpe(i,j-3,k)*M8p(2) &
                     + dpe(i,j-2,k)*M8p(3) + dpe(i,j-1,k)*M8p(4) &
                     + dpe(i,j  ,k)*M8p(5) + dpe(i,j+1,k)*M8p(6) &
                     + dpe(i,j+2,k)*M8p(7) + dpe(i,j+3,k)*M8p(8) ) )

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

   !EXPAND                M8X = matmul(M8, Q(i,j-4:j+3,k,qxn))
                   M8X(1) = M8(1,1) * Q(i,j-4,k,qxn) &
                          + M8(1,2) * Q(i,j-3,k,qxn) &
                          + M8(1,3) * Q(i,j-2,k,qxn) &
                          + M8(1,4) * Q(i,j-1,k,qxn) &
                          - M8(8,4) * Q(i,j  ,k,qxn)
                   M8X(2) = M8(2,1) * Q(i,j-4,k,qxn) &
                          + M8(2,2) * Q(i,j-3,k,qxn) &
                          + M8(2,3) * Q(i,j-2,k,qxn) &
                          + M8(2,4) * Q(i,j-1,k,qxn) &
                          - M8(7,4) * Q(i,j  ,k,qxn) &
                          - M8(7,3) * Q(i,j+1,k,qxn)
                   M8X(3) = M8(3,1) * Q(i,j-4,k,qxn) &
                          + M8(3,2) * Q(i,j-3,k,qxn) &
                          + M8(3,3) * Q(i,j-2,k,qxn) &
                          + M8(3,4) * Q(i,j-1,k,qxn) &
                          - M8(6,4) * Q(i,j  ,k,qxn) &
                          - M8(6,3) * Q(i,j+1,k,qxn) &
                          - M8(6,2) * Q(i,j+2,k,qxn)
                   M8X(4) = M8(4,1) * Q(i,j-4,k,qxn) &
                          + M8(4,2) * Q(i,j-3,k,qxn) &
                          + M8(4,3) * Q(i,j-2,k,qxn) &
                          + M8(4,4) * Q(i,j-1,k,qxn) &
                          - M8(5,4) * Q(i,j  ,k,qxn) &
                          - M8(5,3) * Q(i,j+1,k,qxn) &
                          - M8(5,2) * Q(i,j+2,k,qxn) &
                          - M8(5,1) * Q(i,j+3,k,qxn)
                   M8X(5) = M8(5,1) * Q(i,j-4,k,qxn) &
                          + M8(5,2) * Q(i,j-3,k,qxn) &
                          + M8(5,3) * Q(i,j-2,k,qxn) &
                          + M8(5,4) * Q(i,j-1,k,qxn) &
                          - M8(4,4) * Q(i,j  ,k,qxn) &
                          - M8(4,3) * Q(i,j+1,k,qxn) &
                          - M8(4,2) * Q(i,j+2,k,qxn) &
                          - M8(4,1) * Q(i,j+3,k,qxn)
                   M8X(6) = M8(6,2) * Q(i,j-3,k,qxn) &
                          + M8(6,3) * Q(i,j-2,k,qxn) &
                          + M8(6,4) * Q(i,j-1,k,qxn) &
                          - M8(3,4) * Q(i,j  ,k,qxn) &
                          - M8(3,3) * Q(i,j+1,k,qxn) &
                          - M8(3,2) * Q(i,j+2,k,qxn) &
                          - M8(3,1) * Q(i,j+3,k,qxn)
                   M8X(7) = M8(7,3) * Q(i,j-2,k,qxn) &
                          + M8(7,4) * Q(i,j-1,k,qxn) &
                          - M8(2,4) * Q(i,j  ,k,qxn) &
                          - M8(2,3) * Q(i,j+1,k,qxn) &
                          - M8(2,2) * Q(i,j+2,k,qxn) &
                          - M8(2,1) * Q(i,j+3,k,qxn)
                   M8X(8) = M8(8,4) * Q(i,j-1,k,qxn) &
                          - M8(1,4) * Q(i,j  ,k,qxn) &
                          - M8(1,3) * Q(i,j+1,k,qxn) &
                          - M8(1,2) * Q(i,j+2,k,qxn) &
                          - M8(1,1) * Q(i,j+3,k,qxn)

   !EXPAND                Htmp(n) = dot_product(dpy(i,j-4:j+3,k,n), M8P) &
   !EXPAND                     +    dot_product(Ddiag(i,j-4:j+3,k,n), M8X)
                   Htmp(n) =  &
                      ( ( dpy(i,j-4,k,n)*M8P(1) + dpy(i,j-3,k,n)*M8P(2) &
                        + dpy(i,j-2,k,n)*M8P(3) + dpy(i,j-1,k,n)*M8P(4) &
                        + dpy(i,j  ,k,n)*M8P(5) + dpy(i,j+1,k,n)*M8P(6) &
                        + dpy(i,j+2,k,n)*M8P(7) + dpy(i,j+3,k,n)*M8P(8) ) &
                      + ( Ddiag(i,j-4,k,n)*M8X(1) + Ddiag(i,j-3,k,n)*M8X(2) &
                        + Ddiag(i,j-2,k,n)*M8X(3) + Ddiag(i,j-1,k,n)*M8X(4) &
                        + Ddiag(i,j  ,k,n)*M8X(5) + Ddiag(i,j+1,k,n)*M8X(6) &
                        + Ddiag(i,j+2,k,n)*M8X(7) + Ddiag(i,j+3,k,n)*M8X(8) ) )

   !EXPAND                Hgy(i,j,k,iene) = Hgy(i,j,k,iene) &
   !EXPAND                     +    dot_product(dxe(i,j-4:j+3,k,n), M8X)
                   Hgy(i,j,k,iene) = Hgy(i,j,k,iene)+ &
                      ( dxe(i,j-4,k,n)*M8X(1) + dxe(i,j-3,k,n)*M8X(2) &
                      + dxe(i,j-2,k,n)*M8X(3) + dxe(i,j-1,k,n)*M8X(4) &
                      + dxe(i,j  ,k,n)*M8X(5) + dxe(i,j+1,k,n)*M8X(6) &
                      + dxe(i,j+2,k,n)*M8X(7) + dxe(i,j+3,k,n)*M8X(8) )

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (Q(i,j-1,k,qyn) + Q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hgy(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (Q(i,j-1,k,qhn) + Q(i,j,k,qhn)) / 2.d0
                   Hgy(i,j,k,iene) =  Hgy(i,j,k,iene) - Ytmp(n) * hhalf * Htot
                end do

                ! add y-direction Fdif
                if (j >= lo(2)+1) then
                   Uprime(i,j-1,k,imx ) = Uprime(i,j-1,k,imx ) + (Hgy(i,j,k,imx ) - Hgy(i,j-1,k,imx )) * dx2inv(2)
                   Uprime(i,j-1,k,imy ) = Uprime(i,j-1,k,imy ) + (Hgy(i,j,k,imy ) - Hgy(i,j-1,k,imy )) * dx2inv(2)
                   Uprime(i,j-1,k,imz ) = Uprime(i,j-1,k,imz ) + (Hgy(i,j,k,imz ) - Hgy(i,j-1,k,imz )) * dx2inv(2)
                   Uprime(i,j-1,k,iene) = Uprime(i,j-1,k,iene) + (Hgy(i,j,k,iene) - Hgy(i,j-1,k,iene)) * dx2inv(2)
                   do n=1,nspecies
                      iryn = iry1+n-1
                      Uprime(i,j-1,k,iryn) = Uprime(i,j-1,k,iryn) + (Hgy(i,j,k,iryn) - Hgy(i,j-1,k,iryn)) * dx2inv(2)
                   end do
                end if
             end if
    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------
             if (i <= hi(1) .and. j <= hi(2)) then
   !EXPAND             mmtmp = matmul(mu(i,j,k-4:k+3), M8)
                mmtmp(1) = mu(i,j,k-4) * M8(1,1) &
                         + mu(i,j,k-3) * M8(2,1) &
                         + mu(i,j,k-2) * M8(3,1) &
                         + mu(i,j,k-1) * M8(4,1) &
                         + mu(i,j,k  ) * M8(5,1)
                mmtmp(2) = mu(i,j,k-4) * M8(1,2) &
                         + mu(i,j,k-3) * M8(2,2) &
                         + mu(i,j,k-2) * M8(3,2) &
                         + mu(i,j,k-1) * M8(4,2) &
                         + mu(i,j,k  ) * M8(5,2) &
                         + mu(i,j,k+1) * M8(6,2)
                mmtmp(3) = mu(i,j,k-4) * M8(1,3) &
                         + mu(i,j,k-3) * M8(2,3) &
                         + mu(i,j,k-2) * M8(3,3) &
                         + mu(i,j,k-1) * M8(4,3) &
                         + mu(i,j,k  ) * M8(5,3) &
                         + mu(i,j,k+1) * M8(6,3) &
                         + mu(i,j,k+2) * M8(7,3)
                mmtmp(4) = mu(i,j,k-4) * M8(1,4) &
                         + mu(i,j,k-3) * M8(2,4) &
                         + mu(i,j,k-2) * M8(3,4) &
                         + mu(i,j,k-1) * M8(4,4) &
                         + mu(i,j,k  ) * M8(5,4) &
                         + mu(i,j,k+1) * M8(6,4) &
                         + mu(i,j,k+2) * M8(7,4) &
                         + mu(i,j,k+3) * M8(8,4)
                mmtmp(5) =-mu(i,j,k-4) * M8(8,4) &
                         - mu(i,j,k-3) * M8(7,4) &
                         - mu(i,j,k-2) * M8(6,4) &
                         - mu(i,j,k-1) * M8(5,4) &
                         - mu(i,j,k  ) * M8(4,4) &
                         - mu(i,j,k+1) * M8(3,4) &
                         - mu(i,j,k+2) * M8(2,4) &
                         - mu(i,j,k+3) * M8(1,4)
                mmtmp(6) =-mu(i,j,k-3) * M8(7,3) &
                         - mu(i,j,k-2) * M8(6,3) &
                         - mu(i,j,k-1) * M8(5,3) &
                         - mu(i,j,k  ) * M8(4,3) &
                         - mu(i,j,k+1) * M8(3,3) &
                         - mu(i,j,k+2) * M8(2,3) &
                         - mu(i,j,k+3) * M8(1,3)
                mmtmp(7) =-mu(i,j,k-2) * M8(6,2) &
                         - mu(i,j,k-1) * M8(5,2) &
                         - mu(i,j,k  ) * M8(4,2) &
                         - mu(i,j,k+1) * M8(3,2) &
                         - mu(i,j,k+2) * M8(2,2) &
                         - mu(i,j,k+3) * M8(1,2)
                mmtmp(8) =-mu(i,j,k-1) * M8(5,1) &
                         - mu(i,j,k  ) * M8(4,1) &
                         - mu(i,j,k+1) * M8(3,1) &
                         - mu(i,j,k+2) * M8(2,1) &
                         - mu(i,j,k+3) * M8(1,1)
   !EXPAND             Hgz(i,j,k,imx) = dot_product(mmtmp, Q(i,j,k-4:k+3,qu))
                Hgz(i,j,k,imx) =  &
                   ( Q(i,j,k-4,qu)*mmtmp(1) + Q(i,j,k-3,qu)*mmtmp(2) &
                   + Q(i,j,k-2,qu)*mmtmp(3) + Q(i,j,k-1,qu)*mmtmp(4) &
                   + Q(i,j,k  ,qu)*mmtmp(5) + Q(i,j,k+1,qu)*mmtmp(6) &
                   + Q(i,j,k+2,qu)*mmtmp(7) + Q(i,j,k+3,qu)*mmtmp(8) )
   !EXPAND             Hgz(i,j,k,imy) = dot_product(mmtmp, Q(i,j,k-4:k+3,qv))
                Hgz(i,j,k,imy) =  &
                   ( Q(i,j,k-4,qv)*mmtmp(1) + Q(i,j,k-3,qv)*mmtmp(2) &
                   + Q(i,j,k-2,qv)*mmtmp(3) + Q(i,j,k-1,qv)*mmtmp(4) &
                   + Q(i,j,k  ,qv)*mmtmp(5) + Q(i,j,k+1,qv)*mmtmp(6) &
                   + Q(i,j,k+2,qv)*mmtmp(7) + Q(i,j,k+3,qv)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(vsp(i,j,k-4:k+3), M8)
                mmtmp(1) = vsp(i,j,k-4) * M8(1,1) &
                         + vsp(i,j,k-3) * M8(2,1) &
                         + vsp(i,j,k-2) * M8(3,1) &
                         + vsp(i,j,k-1) * M8(4,1) &
                         + vsp(i,j,k  ) * M8(5,1)
                mmtmp(2) = vsp(i,j,k-4) * M8(1,2) &
                         + vsp(i,j,k-3) * M8(2,2) &
                         + vsp(i,j,k-2) * M8(3,2) &
                         + vsp(i,j,k-1) * M8(4,2) &
                         + vsp(i,j,k  ) * M8(5,2) &
                         + vsp(i,j,k+1) * M8(6,2)
                mmtmp(3) = vsp(i,j,k-4) * M8(1,3) &
                         + vsp(i,j,k-3) * M8(2,3) &
                         + vsp(i,j,k-2) * M8(3,3) &
                         + vsp(i,j,k-1) * M8(4,3) &
                         + vsp(i,j,k  ) * M8(5,3) &
                         + vsp(i,j,k+1) * M8(6,3) &
                         + vsp(i,j,k+2) * M8(7,3)
                mmtmp(4) = vsp(i,j,k-4) * M8(1,4) &
                         + vsp(i,j,k-3) * M8(2,4) &
                         + vsp(i,j,k-2) * M8(3,4) &
                         + vsp(i,j,k-1) * M8(4,4) &
                         + vsp(i,j,k  ) * M8(5,4) &
                         + vsp(i,j,k+1) * M8(6,4) &
                         + vsp(i,j,k+2) * M8(7,4) &
                         + vsp(i,j,k+3) * M8(8,4)
                mmtmp(5) =-vsp(i,j,k-4) * M8(8,4) &
                         - vsp(i,j,k-3) * M8(7,4) &
                         - vsp(i,j,k-2) * M8(6,4) &
                         - vsp(i,j,k-1) * M8(5,4) &
                         - vsp(i,j,k  ) * M8(4,4) &
                         - vsp(i,j,k+1) * M8(3,4) &
                         - vsp(i,j,k+2) * M8(2,4) &
                         - vsp(i,j,k+3) * M8(1,4)
                mmtmp(6) =-vsp(i,j,k-3) * M8(7,3) &
                         - vsp(i,j,k-2) * M8(6,3) &
                         - vsp(i,j,k-1) * M8(5,3) &
                         - vsp(i,j,k  ) * M8(4,3) &
                         - vsp(i,j,k+1) * M8(3,3) &
                         - vsp(i,j,k+2) * M8(2,3) &
                         - vsp(i,j,k+3) * M8(1,3)
                mmtmp(7) =-vsp(i,j,k-2) * M8(6,2) &
                         - vsp(i,j,k-1) * M8(5,2) &
                         - vsp(i,j,k  ) * M8(4,2) &
                         - vsp(i,j,k+1) * M8(3,2) &
                         - vsp(i,j,k+2) * M8(2,2) &
                         - vsp(i,j,k+3) * M8(1,2)
                mmtmp(8) =-vsp(i,j,k-1) * M8(5,1) &
                         - vsp(i,j,k  ) * M8(4,1) &
                         - vsp(i,j,k+1) * M8(3,1) &
                         - vsp(i,j,k+2) * M8(2,1) &
                         - vsp(i,j,k+3) * M8(1,1)
   !EXPAND             Hgz(i,j,k,imz) = dot_product(mmtmp, Q(i,j,k-4:k+3,qw))
                Hgz(i,j,k,imz) =  &
                   ( Q(i,j,k-4,qw)*mmtmp(1) + Q(i,j,k-3,qw)*mmtmp(2) &
                   + Q(i,j,k-2,qw)*mmtmp(3) + Q(i,j,k-1,qw)*mmtmp(4) &
                   + Q(i,j,k  ,qw)*mmtmp(5) + Q(i,j,k+1,qw)*mmtmp(6) &
                   + Q(i,j,k+2,qw)*mmtmp(7) + Q(i,j,k+3,qw)*mmtmp(8) )

   !EXPAND             mmtmp = matmul(lam(i,j,k-4:k+3), M8)
                mmtmp(1) = lam(i,j,k-4) * M8(1,1) &
                         + lam(i,j,k-3) * M8(2,1) &
                         + lam(i,j,k-2) * M8(3,1) &
                         + lam(i,j,k-1) * M8(4,1) &
                         + lam(i,j,k  ) * M8(5,1)
                mmtmp(2) = lam(i,j,k-4) * M8(1,2) &
                         + lam(i,j,k-3) * M8(2,2) &
                         + lam(i,j,k-2) * M8(3,2) &
                         + lam(i,j,k-1) * M8(4,2) &
                         + lam(i,j,k  ) * M8(5,2) &
                         + lam(i,j,k+1) * M8(6,2)
                mmtmp(3) = lam(i,j,k-4) * M8(1,3) &
                         + lam(i,j,k-3) * M8(2,3) &
                         + lam(i,j,k-2) * M8(3,3) &
                         + lam(i,j,k-1) * M8(4,3) &
                         + lam(i,j,k  ) * M8(5,3) &
                         + lam(i,j,k+1) * M8(6,3) &
                         + lam(i,j,k+2) * M8(7,3)
                mmtmp(4) = lam(i,j,k-4) * M8(1,4) &
                         + lam(i,j,k-3) * M8(2,4) &
                         + lam(i,j,k-2) * M8(3,4) &
                         + lam(i,j,k-1) * M8(4,4) &
                         + lam(i,j,k  ) * M8(5,4) &
                         + lam(i,j,k+1) * M8(6,4) &
                         + lam(i,j,k+2) * M8(7,4) &
                         + lam(i,j,k+3) * M8(8,4)
                mmtmp(5) =-lam(i,j,k-4) * M8(8,4) &
                         - lam(i,j,k-3) * M8(7,4) &
                         - lam(i,j,k-2) * M8(6,4) &
                         - lam(i,j,k-1) * M8(5,4) &
                         - lam(i,j,k  ) * M8(4,4) &
                         - lam(i,j,k+1) * M8(3,4) &
                         - lam(i,j,k+2) * M8(2,4) &
                         - lam(i,j,k+3) * M8(1,4)
                mmtmp(6) =-lam(i,j,k-3) * M8(7,3) &
                         - lam(i,j,k-2) * M8(6,3) &
                         - lam(i,j,k-1) * M8(5,3) &
                         - lam(i,j,k  ) * M8(4,3) &
                         - lam(i,j,k+1) * M8(3,3) &
                         - lam(i,j,k+2) * M8(2,3) &
                         - lam(i,j,k+3) * M8(1,3)
                mmtmp(7) =-lam(i,j,k-2) * M8(6,2) &
                         - lam(i,j,k-1) * M8(5,2) &
                         - lam(i,j,k  ) * M8(4,2) &
                         - lam(i,j,k+1) * M8(3,2) &
                         - lam(i,j,k+2) * M8(2,2) &
                         - lam(i,j,k+3) * M8(1,2)
                mmtmp(8) =-lam(i,j,k-1) * M8(5,1) &
                         - lam(i,j,k  ) * M8(4,1) &
                         - lam(i,j,k+1) * M8(3,1) &
                         - lam(i,j,k+2) * M8(2,1) &
                         - lam(i,j,k+3) * M8(1,1)
   !EXPAND             M8p = matmul(M8, Q(i,j,k-4:k+3,qpres))
                M8p(1) = M8(1,1) * Q(i,j,k-4,qpres) &
                       + M8(1,2) * Q(i,j,k-3,qpres) &
                       + M8(1,3) * Q(i,j,k-2,qpres) &
                       + M8(1,4) * Q(i,j,k-1,qpres) &
                       - M8(8,4) * Q(i,j,k  ,qpres)
                M8p(2) = M8(2,1) * Q(i,j,k-4,qpres) &
                       + M8(2,2) * Q(i,j,k-3,qpres) &
                       + M8(2,3) * Q(i,j,k-2,qpres) &
                       + M8(2,4) * Q(i,j,k-1,qpres) &
                       - M8(7,4) * Q(i,j,k  ,qpres) &
                       - M8(7,3) * Q(i,j,k+1,qpres)
                M8p(3) = M8(3,1) * Q(i,j,k-4,qpres) &
                       + M8(3,2) * Q(i,j,k-3,qpres) &
                       + M8(3,3) * Q(i,j,k-2,qpres) &
                       + M8(3,4) * Q(i,j,k-1,qpres) &
                       - M8(6,4) * Q(i,j,k  ,qpres) &
                       - M8(6,3) * Q(i,j,k+1,qpres) &
                       - M8(6,2) * Q(i,j,k+2,qpres)
                M8p(4) = M8(4,1) * Q(i,j,k-4,qpres) &
                       + M8(4,2) * Q(i,j,k-3,qpres) &
                       + M8(4,3) * Q(i,j,k-2,qpres) &
                       + M8(4,4) * Q(i,j,k-1,qpres) &
                       - M8(5,4) * Q(i,j,k  ,qpres) &
                       - M8(5,3) * Q(i,j,k+1,qpres) &
                       - M8(5,2) * Q(i,j,k+2,qpres) &
                       - M8(5,1) * Q(i,j,k+3,qpres)
                M8p(5) = M8(5,1) * Q(i,j,k-4,qpres) &
                       + M8(5,2) * Q(i,j,k-3,qpres) &
                       + M8(5,3) * Q(i,j,k-2,qpres) &
                       + M8(5,4) * Q(i,j,k-1,qpres) &
                       - M8(4,4) * Q(i,j,k  ,qpres) &
                       - M8(4,3) * Q(i,j,k+1,qpres) &
                       - M8(4,2) * Q(i,j,k+2,qpres) &
                       - M8(4,1) * Q(i,j,k+3,qpres)
                M8p(6) = M8(6,2) * Q(i,j,k-3,qpres) &
                       + M8(6,3) * Q(i,j,k-2,qpres) &
                       + M8(6,4) * Q(i,j,k-1,qpres) &
                       - M8(3,4) * Q(i,j,k  ,qpres) &
                       - M8(3,3) * Q(i,j,k+1,qpres) &
                       - M8(3,2) * Q(i,j,k+2,qpres) &
                       - M8(3,1) * Q(i,j,k+3,qpres)
                M8p(7) = M8(7,3) * Q(i,j,k-2,qpres) &
                       + M8(7,4) * Q(i,j,k-1,qpres) &
                       - M8(2,4) * Q(i,j,k  ,qpres) &
                       - M8(2,3) * Q(i,j,k+1,qpres) &
                       - M8(2,2) * Q(i,j,k+2,qpres) &
                       - M8(2,1) * Q(i,j,k+3,qpres)
                M8p(8) = M8(8,4) * Q(i,j,k-1,qpres) &
                       - M8(1,4) * Q(i,j,k  ,qpres) &
                       - M8(1,3) * Q(i,j,k+1,qpres) &
                       - M8(1,2) * Q(i,j,k+2,qpres) &
                       - M8(1,1) * Q(i,j,k+3,qpres)
   !EXPAND             Hgz(i,j,k,iene) = dot_product(mmtmp, Q(i,j,k-4:k+3,qtemp)) &
   !EXPAND                  +           dot_product(     dpe(i,j,k-4:k+3), M8p)
                Hgz(i,j,k,iene) =  &
                   ( ( Q(i,j,k-4,qtemp)*mmtmp(1) + Q(i,j,k-3,qtemp)*mmtmp(2) &
                     + Q(i,j,k-2,qtemp)*mmtmp(3) + Q(i,j,k-1,qtemp)*mmtmp(4) &
                     + Q(i,j,k  ,qtemp)*mmtmp(5) + Q(i,j,k+1,qtemp)*mmtmp(6) &
                     + Q(i,j,k+2,qtemp)*mmtmp(7) + Q(i,j,k+3,qtemp)*mmtmp(8) ) &
                   + ( dpe(i,j,k-4)*M8p(1) + dpe(i,j,k-3)*M8p(2) &
                     + dpe(i,j,k-2)*M8p(3) + dpe(i,j,k-1)*M8p(4) &
                     + dpe(i,j,k  )*M8p(5) + dpe(i,j,k+1)*M8p(6) &
                     + dpe(i,j,k+2)*M8p(7) + dpe(i,j,k+3)*M8p(8) ) )

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

   !EXPAND                M8X = matmul(M8, Q(i,j,k-4:k+3,qxn))
                   M8X(1) = M8(1,1) * Q(i,j,k-4,qxn) &
                          + M8(1,2) * Q(i,j,k-3,qxn) &
                          + M8(1,3) * Q(i,j,k-2,qxn) &
                          + M8(1,4) * Q(i,j,k-1,qxn) &
                          - M8(8,4) * Q(i,j,k  ,qxn)
                   M8X(2) = M8(2,1) * Q(i,j,k-4,qxn) &
                          + M8(2,2) * Q(i,j,k-3,qxn) &
                          + M8(2,3) * Q(i,j,k-2,qxn) &
                          + M8(2,4) * Q(i,j,k-1,qxn) &
                          - M8(7,4) * Q(i,j,k  ,qxn) &
                          - M8(7,3) * Q(i,j,k+1,qxn)
                   M8X(3) = M8(3,1) * Q(i,j,k-4,qxn) &
                          + M8(3,2) * Q(i,j,k-3,qxn) &
                          + M8(3,3) * Q(i,j,k-2,qxn) &
                          + M8(3,4) * Q(i,j,k-1,qxn) &
                          - M8(6,4) * Q(i,j,k  ,qxn) &
                          - M8(6,3) * Q(i,j,k+1,qxn) &
                          - M8(6,2) * Q(i,j,k+2,qxn)
                   M8X(4) = M8(4,1) * Q(i,j,k-4,qxn) &
                          + M8(4,2) * Q(i,j,k-3,qxn) &
                          + M8(4,3) * Q(i,j,k-2,qxn) &
                          + M8(4,4) * Q(i,j,k-1,qxn) &
                          - M8(5,4) * Q(i,j,k  ,qxn) &
                          - M8(5,3) * Q(i,j,k+1,qxn) &
                          - M8(5,2) * Q(i,j,k+2,qxn) &
                          - M8(5,1) * Q(i,j,k+3,qxn)
                   M8X(5) = M8(5,1) * Q(i,j,k-4,qxn) &
                          + M8(5,2) * Q(i,j,k-3,qxn) &
                          + M8(5,3) * Q(i,j,k-2,qxn) &
                          + M8(5,4) * Q(i,j,k-1,qxn) &
                          - M8(4,4) * Q(i,j,k  ,qxn) &
                          - M8(4,3) * Q(i,j,k+1,qxn) &
                          - M8(4,2) * Q(i,j,k+2,qxn) &
                          - M8(4,1) * Q(i,j,k+3,qxn)
                   M8X(6) = M8(6,2) * Q(i,j,k-3,qxn) &
                          + M8(6,3) * Q(i,j,k-2,qxn) &
                          + M8(6,4) * Q(i,j,k-1,qxn) &
                          - M8(3,4) * Q(i,j,k  ,qxn) &
                          - M8(3,3) * Q(i,j,k+1,qxn) &
                          - M8(3,2) * Q(i,j,k+2,qxn) &
                          - M8(3,1) * Q(i,j,k+3,qxn)
                   M8X(7) = M8(7,3) * Q(i,j,k-2,qxn) &
                          + M8(7,4) * Q(i,j,k-1,qxn) &
                          - M8(2,4) * Q(i,j,k  ,qxn) &
                          - M8(2,3) * Q(i,j,k+1,qxn) &
                          - M8(2,2) * Q(i,j,k+2,qxn) &
                          - M8(2,1) * Q(i,j,k+3,qxn)
                   M8X(8) = M8(8,4) * Q(i,j,k-1,qxn) &
                          - M8(1,4) * Q(i,j,k  ,qxn) &
                          - M8(1,3) * Q(i,j,k+1,qxn) &
                          - M8(1,2) * Q(i,j,k+2,qxn) &
                          - M8(1,1) * Q(i,j,k+3,qxn)

   !EXPAND                Htmp(n) = dot_product(dpy(i,j,k-4:k+3,n), M8P) &
   !EXPAND                     +    dot_product(Ddiag(i,j,k-4:k+3,n), M8X)
                   Htmp(n) =  &
                      ( ( dpy(i,j,k-4,n)*M8P(1) + dpy(i,j,k-3,n)*M8P(2) &
                        + dpy(i,j,k-2,n)*M8P(3) + dpy(i,j,k-1,n)*M8P(4) &
                        + dpy(i,j,k  ,n)*M8P(5) + dpy(i,j,k+1,n)*M8P(6) &
                        + dpy(i,j,k+2,n)*M8P(7) + dpy(i,j,k+3,n)*M8P(8) ) &
                      + ( Ddiag(i,j,k-4,n)*M8X(1) + Ddiag(i,j,k-3,n)*M8X(2) &
                        + Ddiag(i,j,k-2,n)*M8X(3) + Ddiag(i,j,k-1,n)*M8X(4) &
                        + Ddiag(i,j,k  ,n)*M8X(5) + Ddiag(i,j,k+1,n)*M8X(6) &
                        + Ddiag(i,j,k+2,n)*M8X(7) + Ddiag(i,j,k+3,n)*M8X(8) ) )

   !EXPAND                Hgz(i,j,k,iene) = Hgz(i,j,k,iene) &
   !EXPAND                     +    dot_product(dxe(i,j,k-4:k+3,n), M8X)
                   Hgz(i,j,k,iene) = Hgz(i,j,k,iene)+ &
                      ( dxe(i,j,k-4,n)*M8X(1) + dxe(i,j,k-3,n)*M8X(2) &
                      + dxe(i,j,k-2,n)*M8X(3) + dxe(i,j,k-1,n)*M8X(4) &
                      + dxe(i,j,k  ,n)*M8X(5) + dxe(i,j,k+1,n)*M8X(6) &
                      + dxe(i,j,k+2,n)*M8X(7) + dxe(i,j,k+3,n)*M8X(8) )

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = (Q(i,j,k-1,qyn) + Q(i,j,k,qyn)) / 2.d0
                end do

                do n = 1, nspecies
                   Hgz(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = (Q(i,j,k-1,qhn) + Q(i,j,k,qhn)) / 2.d0
                   Hgz(i,j,k,iene) =  Hgz(i,j,k,iene) - Ytmp(n) * hhalf * Htot
                end do

                ! add z-direction Fdif
                if (k >= lo(3)+1) then
                   Uprime(i,j,k-1,imx) = Uprime(i,j,k-1,imx) + (Hgz(i,j,k,imx) - Hgz(i,j,k-1,imx)) * dx2inv(3)
                   Uprime(i,j,k-1,imy) = Uprime(i,j,k-1,imy) + (Hgz(i,j,k,imy) - Hgz(i,j,k-1,imy)) * dx2inv(3)
                   Uprime(i,j,k-1,imz) = Uprime(i,j,k-1,imz) + (Hgz(i,j,k,imz) - Hgz(i,j,k-1,imz)) * dx2inv(3)
                   Uprime(i,j,k-1,iene) = Uprime(i,j,k-1,iene) + (Hgz(i,j,k,iene) - Hgz(i,j,k-1,iene)) * dx2inv(3)
                   do n=1,nspecies
                      iryn = iry1+n-1
                      Uprime(i,j,k-1,iryn) = Uprime(i,j,k-1,iryn) + (Hgz(i,j,k,iryn) - Hgz(i,j,k-1,iryn)) * dx2inv(3)
                   end do
                end if
             end if
    ! ------- END z-direction -------
          end do
       end do
    end do
    !$omp end do nowait
    
    !$omp barrier

    ! add kinetic energy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Uprime(i,j,k,iene) = Uprime(i,j,k,iene) &
                  + Uprime(i,j,k,imx)*Q(i,j,k,qu) &
                  + Uprime(i,j,k,imy)*Q(i,j,k,qv) &
                  + Uprime(i,j,k,imz)*Q(i,j,k,qw)
          end do
       end do
    end do
    !$omp end do 
    
    !$omp end parallel

    deallocate(Hgx,Hgy,Hgz,dpy,dxe,dpe,vsp,vsm)

    !call destroy(bpt_diffterm)                !! ^^^^^^^^^^^^^^^^^^^^^^^ timer

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

    !$omp parallel private(i,j,k,n,un)
    !$omp do 
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

!             un = Q(i-4:i+4,j,k,qu)

!EXPAND             Fhyp(i,j,k,irho) = Fhyp(i,j,k,irho) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imx) ) 
             Uprime(i,j,k,irho) = Uprime(i,j,k,irho) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imx)-U(i-1,j,k,imx)) &
                + D8(2)*(U(i+2,j,k,imx)-U(i-2,j,k,imx)) &
                + D8(3)*(U(i+3,j,k,imx)-U(i-3,j,k,imx)) &
                + D8(4)*(U(i+4,j,k,imx)-U(i-4,j,k,imx)) )

!EXPAND             Fhyp(i,j,k,imx) = Fhyp(i,j,k,imx) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imx)*un+Q(i-4:i+4,j,k,qpres) )
             Uprime(i,j,k,imx) = Uprime(i,j,k,imx) - dxinv(1) * &
                ( D8(1)*((U(i+1,j,k,imx)*Q(i+1,j,k,qu)+Q(i+1,j,k,qpres))-(U(i-1,j,k,imx)*Q(i-1,j,k,qu)+Q(i-1,j,k,qpres))) &
                + D8(2)*((U(i+2,j,k,imx)*Q(i+2,j,k,qu)+Q(i+2,j,k,qpres))-(U(i-2,j,k,imx)*Q(i-2,j,k,qu)+Q(i-2,j,k,qpres))) &
                + D8(3)*((U(i+3,j,k,imx)*Q(i+3,j,k,qu)+Q(i+3,j,k,qpres))-(U(i-3,j,k,imx)*Q(i-3,j,k,qu)+Q(i-3,j,k,qpres))) &
                + D8(4)*((U(i+4,j,k,imx)*Q(i+4,j,k,qu)+Q(i+4,j,k,qpres))-(U(i-4,j,k,imx)*Q(i-4,j,k,qu)+Q(i-4,j,k,qpres))) )

!EXPAND             Fhyp(i,j,k,imy) = Fhyp(i,j,k,imy) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imy)*un ) 
             Uprime(i,j,k,imy) = Uprime(i,j,k,imy) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imy)*Q(i+1,j,k,qu)-U(i-1,j,k,imy)*Q(i-1,j,k,qu)) &
                + D8(2)*(U(i+2,j,k,imy)*Q(i+2,j,k,qu)-U(i-2,j,k,imy)*Q(i-2,j,k,qu)) &
                + D8(3)*(U(i+3,j,k,imy)*Q(i+3,j,k,qu)-U(i-3,j,k,imy)*Q(i-3,j,k,qu)) &
                + D8(4)*(U(i+4,j,k,imy)*Q(i+4,j,k,qu)-U(i-4,j,k,imy)*Q(i-4,j,k,qu)) )

!EXPAND             Fhyp(i,j,k,imz) = Fhyp(i,j,k,imz) - dxinv(1) * &
!EXPAND                  first_deriv_8( U(i-4:i+4,j,k,imz)*un ) 
             Uprime(i,j,k,imz) = Uprime(i,j,k,imz) - dxinv(1) * &
                ( D8(1)*(U(i+1,j,k,imz)*Q(i+1,j,k,qu)-U(i-1,j,k,imz)*Q(i-1,j,k,qu)) &
                + D8(2)*(U(i+2,j,k,imz)*Q(i+2,j,k,qu)-U(i-2,j,k,imz)*Q(i-2,j,k,qu)) &
                + D8(3)*(U(i+3,j,k,imz)*Q(i+3,j,k,qu)-U(i-3,j,k,imz)*Q(i-3,j,k,qu)) &
                + D8(4)*(U(i+4,j,k,imz)*Q(i+4,j,k,qu)-U(i-4,j,k,imz)*Q(i-4,j,k,qu)) )

!EXPAND             Fhyp(i,j,k,iene) = Fhyp(i,j,k,iene) - dxinv(1) * &
!EXPAND                  first_deriv_8( (U(i-4:i+4,j,k,iene)+Q(i-4:i+4,j,k,qpres))*un )
             Uprime(i,j,k,iene) = Uprime(i,j,k,iene) - dxinv(1) * &
                ( D8(1)*((U(i+1,j,k,iene)+Q(i+1,j,k,qpres))*Q(i+1,j,k,qu)-(U(i-1,j,k,iene)+Q(i-1,j,k,qpres))*Q(i-1,j,k,qu)) &
                + D8(2)*((U(i+2,j,k,iene)+Q(i+2,j,k,qpres))*Q(i+2,j,k,qu)-(U(i-2,j,k,iene)+Q(i-2,j,k,qpres))*Q(i-2,j,k,qu)) &
                + D8(3)*((U(i+3,j,k,iene)+Q(i+3,j,k,qpres))*Q(i+3,j,k,qu)-(U(i-3,j,k,iene)+Q(i-3,j,k,qpres))*Q(i-3,j,k,qu)) &
                + D8(4)*((U(i+4,j,k,iene)+Q(i+4,j,k,qpres))*Q(i+4,j,k,qu)-(U(i-4,j,k,iene)+Q(i-4,j,k,qpres))*Q(i-4,j,k,qu)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(1) * &
!EXPAND                     first_deriv_8( U(i-4:i+4,j,k,n)*un )
                Uprime(i,j,k,iryn) = Uprime(i,j,k,iryn) - dxinv(1) * &
                   ( D8(1)*(U(i+1,j,k,iryn)*Q(i+1,j,k,qu)-U(i-1,j,k,iryn)*Q(i-1,j,k,qu)) &
                   + D8(2)*(U(i+2,j,k,iryn)*Q(i+2,j,k,qu)-U(i-2,j,k,iryn)*Q(i-2,j,k,qu)) &
                   + D8(3)*(U(i+3,j,k,iryn)*Q(i+3,j,k,qu)-U(i-3,j,k,iryn)*Q(i-3,j,k,qu)) &
                   + D8(4)*(U(i+4,j,k,iryn)*Q(i+4,j,k,qu)-U(i-4,j,k,iryn)*Q(i-4,j,k,qu)) )
             end do

!EXPAND             Fhyp(i,j,k,irho)=Fhyp(i,j,k,irho) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imy) )
             Uprime(i,j,k,irho) = Uprime(i,j,k,irho) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imy)-U(i,j-1,k,imy)) &
                + D8(2)*(U(i,j+2,k,imy)-U(i,j-2,k,imy)) &
                + D8(3)*(U(i,j+3,k,imy)-U(i,j-3,k,imy)) &
                + D8(4)*(U(i,j+4,k,imy)-U(i,j-4,k,imy)) )

!EXPAND             Fhyp(i,j,k,imx)=Fhyp(i,j,k,imx) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imx)*un )
             Uprime(i,j,k,imx) = Uprime(i,j,k,imx) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imx)*Q(i,j+1,k,qv)-U(i,j-1,k,imx)*Q(i,j-1,k,qv)) &
                + D8(2)*(U(i,j+2,k,imx)*Q(i,j+2,k,qv)-U(i,j-2,k,imx)*Q(i,j-2,k,qv)) &
                + D8(3)*(U(i,j+3,k,imx)*Q(i,j+3,k,qv)-U(i,j-3,k,imx)*Q(i,j-3,k,qv)) &
                + D8(4)*(U(i,j+4,k,imx)*Q(i,j+4,k,qv)-U(i,j-4,k,imx)*Q(i,j-4,k,qv)) )

!EXPAND             Fhyp(i,j,k,imy)=Fhyp(i,j,k,imy) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imy)*un+Q(i,j-4:j+4,k,qpres) )
             Uprime(i,j,k,imy) = Uprime(i,j,k,imy) - dxinv(2) * &
                ( D8(1)*((U(i,j+1,k,imy)*Q(i,j+1,k,qv)+Q(i,j+1,k,qpres))-(U(i,j-1,k,imy)*Q(i,j-1,k,qv)+Q(i,j-1,k,qpres))) &
                + D8(2)*((U(i,j+2,k,imy)*Q(i,j+2,k,qv)+Q(i,j+2,k,qpres))-(U(i,j-2,k,imy)*Q(i,j-2,k,qv)+Q(i,j-2,k,qpres))) &
                + D8(3)*((U(i,j+3,k,imy)*Q(i,j+3,k,qv)+Q(i,j+3,k,qpres))-(U(i,j-3,k,imy)*Q(i,j-3,k,qv)+Q(i,j-3,k,qpres))) &
                + D8(4)*((U(i,j+4,k,imy)*Q(i,j+4,k,qv)+Q(i,j+4,k,qpres))-(U(i,j-4,k,imy)*Q(i,j-4,k,qv)+Q(i,j-4,k,qpres))) )

!EXPAND             Fhyp(i,j,k,imz)=Fhyp(i,j,k,imz) - dxinv(2) * &
!EXPAND                  first_deriv_8( U(i,j-4:j+4,k,imz)*un )
             Uprime(i,j,k,imz) = Uprime(i,j,k,imz) - dxinv(2) * &
                ( D8(1)*(U(i,j+1,k,imz)*Q(i,j+1,k,qv)-U(i,j-1,k,imz)*Q(i,j-1,k,qv)) &
                + D8(2)*(U(i,j+2,k,imz)*Q(i,j+2,k,qv)-U(i,j-2,k,imz)*Q(i,j-2,k,qv)) &
                + D8(3)*(U(i,j+3,k,imz)*Q(i,j+3,k,qv)-U(i,j-3,k,imz)*Q(i,j-3,k,qv)) &
                + D8(4)*(U(i,j+4,k,imz)*Q(i,j+4,k,qv)-U(i,j-4,k,imz)*Q(i,j-4,k,qv)) )

!EXPAND             Fhyp(i,j,k,iene)=Fhyp(i,j,k,iene) - dxinv(2) * &
!EXPAND                  first_deriv_8( (U(i,j-4:j+4,k,iene)+Q(i,j-4:j+4,k,qpres))*un )
             Uprime(i,j,k,iene) = Uprime(i,j,k,iene) - dxinv(2) * &
                ( D8(1)*((U(i,j+1,k,iene)+Q(i,j+1,k,qpres))*Q(i,j+1,k,qv)-(U(i,j-1,k,iene)+Q(i,j-1,k,qpres))*Q(i,j-1,k,qv)) &
                + D8(2)*((U(i,j+2,k,iene)+Q(i,j+2,k,qpres))*Q(i,j+2,k,qv)-(U(i,j-2,k,iene)+Q(i,j-2,k,qpres))*Q(i,j-2,k,qv)) &
                + D8(3)*((U(i,j+3,k,iene)+Q(i,j+3,k,qpres))*Q(i,j+3,k,qv)-(U(i,j-3,k,iene)+Q(i,j-3,k,qpres))*Q(i,j-3,k,qv)) &
                + D8(4)*((U(i,j+4,k,iene)+Q(i,j+4,k,qpres))*Q(i,j+4,k,qv)-(U(i,j-4,k,iene)+Q(i,j-4,k,qpres))*Q(i,j-4,k,qv)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(2) * &
!EXPAND                     first_deriv_8( U(i,j-4:j+4,k,n)*un )
                Uprime(i,j,k,iryn) = Uprime(i,j,k,iryn) - dxinv(2) * &
                   ( D8(1)*(U(i,j+1,k,iryn)*Q(i,j+1,k,qv)-U(i,j-1,k,iryn)*Q(i,j-1,k,qv)) &
                   + D8(2)*(U(i,j+2,k,iryn)*Q(i,j+2,k,qv)-U(i,j-2,k,iryn)*Q(i,j-2,k,qv)) &
                   + D8(3)*(U(i,j+3,k,iryn)*Q(i,j+3,k,qv)-U(i,j-3,k,iryn)*Q(i,j-3,k,qv)) &
                   + D8(4)*(U(i,j+4,k,iryn)*Q(i,j+4,k,qv)-U(i,j-4,k,iryn)*Q(i,j-4,k,qv)) )
             end do

!EXPAND             Fhyp(i,j,k,irho)=Fhyp(i,j,k,irho) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imz) )
             Uprime(i,j,k,irho) = Uprime(i,j,k,irho) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imz)-U(i,j,k-1,imz)) &
                + D8(2)*(U(i,j,k+2,imz)-U(i,j,k-2,imz)) &
                + D8(3)*(U(i,j,k+3,imz)-U(i,j,k-3,imz)) &
                + D8(4)*(U(i,j,k+4,imz)-U(i,j,k-4,imz)) )

!EXPAND             Fhyp(i,j,k,imx)=Fhyp(i,j,k,imx) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imx)*qwUN )
             Uprime(i,j,k,imx) = Uprime(i,j,k,imx) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imx)*Q(i,j,k+1,qw)-U(i,j,k-1,imx)*Q(i,j,k-1,qw)) &
                + D8(2)*(U(i,j,k+2,imx)*Q(i,j,k+2,qw)-U(i,j,k-2,imx)*Q(i,j,k-2,qw)) &
                + D8(3)*(U(i,j,k+3,imx)*Q(i,j,k+3,qw)-U(i,j,k-3,imx)*Q(i,j,k-3,qw)) &
                + D8(4)*(U(i,j,k+4,imx)*Q(i,j,k+4,qw)-U(i,j,k-4,imx)*Q(i,j,k-4,qw)) )

!EXPAND             Fhyp(i,j,k,imy)=Fhyp(i,j,k,imy) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imy)*qwUN )
             Uprime(i,j,k,imy) = Uprime(i,j,k,imy) - dxinv(3) * &
                ( D8(1)*(U(i,j,k+1,imy)*Q(i,j,k+1,qw)-U(i,j,k-1,imy)*Q(i,j,k-1,qw)) &
                + D8(2)*(U(i,j,k+2,imy)*Q(i,j,k+2,qw)-U(i,j,k-2,imy)*Q(i,j,k-2,qw)) &
                + D8(3)*(U(i,j,k+3,imy)*Q(i,j,k+3,qw)-U(i,j,k-3,imy)*Q(i,j,k-3,qw)) &
                + D8(4)*(U(i,j,k+4,imy)*Q(i,j,k+4,qw)-U(i,j,k-4,imy)*Q(i,j,k-4,qw)) )

!EXPAND             Fhyp(i,j,k,imz)=Fhyp(i,j,k,imz) - dxinv(3) * &
!EXPAND                  first_deriv_8( U(i,j,k-4:k+4,imz)*qwUN+Q(i,j,k-4:k+4,qpres) )
             Uprime(i,j,k,imz) = Uprime(i,j,k,imz) - dxinv(3) * &
                ( D8(1)*((U(i,j,k+1,imz)*Q(i,j,k+1,qw)+Q(i,j,k+1,qpres))-(U(i,j,k-1,imz)*Q(i,j,k-1,qw)+Q(i,j,k-1,qpres))) &
                + D8(2)*((U(i,j,k+2,imz)*Q(i,j,k+2,qw)+Q(i,j,k+2,qpres))-(U(i,j,k-2,imz)*Q(i,j,k-2,qw)+Q(i,j,k-2,qpres))) &
                + D8(3)*((U(i,j,k+3,imz)*Q(i,j,k+3,qw)+Q(i,j,k+3,qpres))-(U(i,j,k-3,imz)*Q(i,j,k-3,qw)+Q(i,j,k-3,qpres))) &
                + D8(4)*((U(i,j,k+4,imz)*Q(i,j,k+4,qw)+Q(i,j,k+4,qpres))-(U(i,j,k-4,imz)*Q(i,j,k-4,qw)+Q(i,j,k-4,qpres))) )

!EXPAND             Fhyp(i,j,k,iene)=Fhyp(i,j,k,iene) - dxinv(3) * &
!EXPAND                  first_deriv_8( (U(i,j,k-4:k+4,iene)+Q(i,j,k-4:k+4,qpres))*qwUN )
             Uprime(i,j,k,iene) = Uprime(i,j,k,iene) - dxinv(3) * &
                ( D8(1)*((U(i,j,k+1,iene)+Q(i,j,k+1,qpres))*Q(i,j,k+1,qw)-(U(i,j,k-1,iene)+Q(i,j,k-1,qpres))*Q(i,j,k-1,qw)) &
                + D8(2)*((U(i,j,k+2,iene)+Q(i,j,k+2,qpres))*Q(i,j,k+2,qw)-(U(i,j,k-2,iene)+Q(i,j,k-2,qpres))*Q(i,j,k-2,qw)) &
                + D8(3)*((U(i,j,k+3,iene)+Q(i,j,k+3,qpres))*Q(i,j,k+3,qw)-(U(i,j,k-3,iene)+Q(i,j,k-3,qpres))*Q(i,j,k-3,qw)) &
                + D8(4)*((U(i,j,k+4,iene)+Q(i,j,k+4,qpres))*Q(i,j,k+4,qw)-(U(i,j,k-4,iene)+Q(i,j,k-4,qpres))*Q(i,j,k-4,qw)) )

             do n = 1, nspecies
                iryn = iry1 + n - 1
!EXPAND                Fhyp(i,j,k,n) = Fhyp(i,j,k,n) - dxinv(3) * &
!EXPAND                     first_deriv_8(U(i,j,k-4:k+4,n)*qwUN)
                Uprime(i,j,k,iryn) = Uprime(i,j,k,iryn) - dxinv(3) * &
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
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,rho)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Unew(i,j,k,imx ) = Zero * Unew(i,j,k,imx ) + One * U(i,j,k,imx ) + dt * Uprime(i,j,k,imx )
             Unew(i,j,k,imy ) = Zero * Unew(i,j,k,imy ) + One * U(i,j,k,imy ) + dt * Uprime(i,j,k,imy )
             Unew(i,j,k,imz ) = Zero * Unew(i,j,k,imz ) + One * U(i,j,k,imz ) + dt * Uprime(i,j,k,imz )
             Unew(i,j,k,iene) = Zero * Unew(i,j,k,iene) + One * U(i,j,k,iene) + dt * Uprime(i,j,k,iene)
             rho = 0.d0
             do n = 1, nspecies
                iryn = iry1+n-1
                Unew(i,j,k,iryn) = Zero * Unew(i,j,k,iryn) + One * U(i,j,k,iryn) + dt * Uprime(i,j,k,iryn)
                rho = rho + Unew(i,j,k,iryn)
             end do
             Unew(i,j,k,irho) = rho
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    

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
