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

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable, dimension(:,:,:) :: vsp,vsm, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hg, dpy, dxe

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
    allocate(Fhyp(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons))
    allocate(mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(lam(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(Ddiag(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies))
    allocate( Fdif(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons))


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

    !$omp workshare
!EXPAND    Fdif = 0.d0
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fdif(i,j,k,irho) = 0.d0
             Fdif(i,j,k,imx ) = 0.d0
             Fdif(i,j,k,imy ) = 0.d0
             Fdif(i,j,k,imz ) = 0.d0
             Fdif(i,j,k,iene) = 0.d0
             do n=1,nspecies
                iryn = iry1+n-1
                Fdif(i,j,k,iryn) = 0.d0
             enddo
          enddo
       enddo
    enddo
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
          do i=lo(1),hi(1)
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
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)   
          do i=dlo(1),dhi(1)
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
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
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
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

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
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + &
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
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dmvywzdx
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dmuydx
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dmuzdx
          end do
       end do
    end do
    !$omp end do

    ! d()/dy
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
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
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dmvxdy
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dmuxwzdy
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dmvzdy
          end do
       end do
    end do
    !$omp end do 
    
    ! d()/dz
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
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
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + dmwxdz
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + dmwydz
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + dmuxvydz
          end do
       end do
    end do
    !$omp end do nowait

    !$omp end parallel

    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz)

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    !$omp parallel &
    !$omp private(i,j,k,n,qxn,qyn,qhn,Htot,Htmp,Ytmp,hhalf,M8p,M8X,mmtmp)

    !$omp workshare
    dpe = 0.d0
    !$omp end workshare

    !$OMP DO
    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
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
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1

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
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp, Q(i-4:i+3,j,k,qu))
             Hg(i,j,k,imx) =  &
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
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp, Q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imy) =  &
                ( Q(i-4,j,k,qv)*mmtmp(1) + Q(i-3,j,k,qv)*mmtmp(2) &
                + Q(i-2,j,k,qv)*mmtmp(3) + Q(i-1,j,k,qv)*mmtmp(4) &
                + Q(i  ,j,k,qv)*mmtmp(5) + Q(i+1,j,k,qv)*mmtmp(6) &
                + Q(i+2,j,k,qv)*mmtmp(7) + Q(i+3,j,k,qv)*mmtmp(8) )
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp, Q(i-4:i+3,j,k,qw))
             Hg(i,j,k,imz) =  &
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
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp, Q(i-4:i+3,j,k,qtemp)) &
!EXPAND                  &         + dot_product(     dpe(i-4:i+3,j,k), M8p)
             Hg(i,j,k,iene) =  &
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

!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
!EXPAND                     +    dot_product(dxe(i-4:i+3,j,k,n), M8X)
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i-4,j,k,n)*M8X(1) + dxe(i-3,j,k,n)*M8X(2) &
                   + dxe(i-2,j,k,n)*M8X(3) + dxe(i-1,j,k,n)*M8X(4) &
                   + dxe(i  ,j,k,n)*M8X(5) + dxe(i+1,j,k,n)*M8X(6) &
                   + dxe(i+2,j,k,n)*M8X(7) + dxe(i+3,j,k,n)*M8X(8) )

                Htot = Htot + Htmp(n)
                Ytmp(n) = (Q(i-1,j,k,qyn) + Q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (Q(i-1,j,k,qhn) + Q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add x-direction Fdif
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fdif(i,j,k,imx ) = Fdif(i,j,k,imx ) + (Hg(i+1,j,k,imx ) - Hg(i,j,k,imx )) * dx2inv(1)
             Fdif(i,j,k,imy ) = Fdif(i,j,k,imy ) + (Hg(i+1,j,k,imy ) - Hg(i,j,k,imy )) * dx2inv(1)
             Fdif(i,j,k,imz ) = Fdif(i,j,k,imz ) + (Hg(i+1,j,k,imz ) - Hg(i,j,k,imz )) * dx2inv(1)
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + (Hg(i+1,j,k,iene) - Hg(i,j,k,iene)) * dx2inv(1)
             do n=1,nspecies
                iryn = iry1+n-1
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + (Hg(i+1,j,k,iryn) - Hg(i,j,k,iryn)) * dx2inv(1)
             end do
          end do
       end do
    end do
    !$omp end do nowait
    ! ------- END x-direction -------

    !$omp barrier

    ! ------- BEGIN y-direction -------
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             
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
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp, Q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imx) =  &
                ( Q(i,j-4,k,qu)*mmtmp(1) + Q(i,j-3,k,qu)*mmtmp(2) &
                + Q(i,j-2,k,qu)*mmtmp(3) + Q(i,j-1,k,qu)*mmtmp(4) &
                + Q(i,j  ,k,qu)*mmtmp(5) + Q(i,j+1,k,qu)*mmtmp(6) &
                + Q(i,j+2,k,qu)*mmtmp(7) + Q(i,j+3,k,qu)*mmtmp(8) )
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp, Q(i,j-4:j+3,k,qw))
             Hg(i,j,k,imz) =  &
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
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp, Q(i,j-4:j+3,k,qv))
             Hg(i,j,k,imy) =  &
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
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp, Q(i,j-4:j+3,k,qtemp)) &
!EXPAND                  +           dot_product(     dpe(i,j-4:j+3,k), M8p)
             Hg(i,j,k,iene) =  &
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

!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
!EXPAND                     +    dot_product(dxe(i,j-4:j+3,k,n), M8X)
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i,j-4,k,n)*M8X(1) + dxe(i,j-3,k,n)*M8X(2) &
                   + dxe(i,j-2,k,n)*M8X(3) + dxe(i,j-1,k,n)*M8X(4) &
                   + dxe(i,j  ,k,n)*M8X(5) + dxe(i,j+1,k,n)*M8X(6) &
                   + dxe(i,j+2,k,n)*M8X(7) + dxe(i,j+3,k,n)*M8X(8) )

                Htot = Htot + Htmp(n)
                Ytmp(n) = (Q(i,j-1,k,qyn) + Q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (Q(i,j-1,k,qhn) + Q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add y-direction Fdif
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fdif(i,j,k,imx ) = Fdif(i,j,k,imx ) + (Hg(i,j+1,k,imx ) - Hg(i,j,k,imx )) * dx2inv(2)
             Fdif(i,j,k,imy ) = Fdif(i,j,k,imy ) + (Hg(i,j+1,k,imy ) - Hg(i,j,k,imy )) * dx2inv(2)
             Fdif(i,j,k,imz ) = Fdif(i,j,k,imz ) + (Hg(i,j+1,k,imz ) - Hg(i,j,k,imz )) * dx2inv(2)
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + (Hg(i,j+1,k,iene) - Hg(i,j,k,iene)) * dx2inv(2)
             do n=1,nspecies
                iryn = iry1+n-1
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + (Hg(i,j+1,k,iryn) - Hg(i,j,k,iryn)) * dx2inv(2)
             end do
          end do
       end do
    end do
    !$omp end do nowait
    ! ------- END y-direction -------

    !$omp barrier

    ! ------- BEGIN z-direction -------
    !$omp do
    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

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
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp, Q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imx) =  &
                ( Q(i,j,k-4,qu)*mmtmp(1) + Q(i,j,k-3,qu)*mmtmp(2) &
                + Q(i,j,k-2,qu)*mmtmp(3) + Q(i,j,k-1,qu)*mmtmp(4) &
                + Q(i,j,k  ,qu)*mmtmp(5) + Q(i,j,k+1,qu)*mmtmp(6) &
                + Q(i,j,k+2,qu)*mmtmp(7) + Q(i,j,k+3,qu)*mmtmp(8) )
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp, Q(i,j,k-4:k+3,qv))
             Hg(i,j,k,imy) =  &
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
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp, Q(i,j,k-4:k+3,qw))
             Hg(i,j,k,imz) =  &
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
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp, Q(i,j,k-4:k+3,qtemp)) &
!EXPAND                  +           dot_product(     dpe(i,j,k-4:k+3), M8p)
             Hg(i,j,k,iene) =  &
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

!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) &
!EXPAND                     +    dot_product(dxe(i,j,k-4:k+3,n), M8X)
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i,j,k-4,n)*M8X(1) + dxe(i,j,k-3,n)*M8X(2) &
                   + dxe(i,j,k-2,n)*M8X(3) + dxe(i,j,k-1,n)*M8X(4) &
                   + dxe(i,j,k  ,n)*M8X(5) + dxe(i,j,k+1,n)*M8X(6) &
                   + dxe(i,j,k+2,n)*M8X(7) + dxe(i,j,k+3,n)*M8X(8) )

                Htot = Htot + Htmp(n)
                Ytmp(n) = (Q(i,j,k-1,qyn) + Q(i,j,k,qyn)) / 2.d0
             end do

             do n = 1, nspecies
                Hg(i,j,k,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = (Q(i,j,k-1,qhn) + Q(i,j,k,qhn)) / 2.d0
                Hg(i,j,k,iene) =  Hg(i,j,k,iene) - Ytmp(n) * hhalf * Htot
             end do

          end do
       end do
    end do
    !$omp end do

    ! add z-direction Fdif
    !$omp do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             Fdif(i,j,k,imx) = Fdif(i,j,k,imx) + (Hg(i,j,k+1,imx) - Hg(i,j,k,imx)) * dx2inv(3)
             Fdif(i,j,k,imy) = Fdif(i,j,k,imy) + (Hg(i,j,k+1,imy) - Hg(i,j,k,imy)) * dx2inv(3)
             Fdif(i,j,k,imz) = Fdif(i,j,k,imz) + (Hg(i,j,k+1,imz) - Hg(i,j,k,imz)) * dx2inv(3)
             Fdif(i,j,k,iene) = Fdif(i,j,k,iene) + (Hg(i,j,k+1,iene) - Hg(i,j,k,iene)) * dx2inv(3)
             do n=1,nspecies
                iryn = iry1+n-1
                Fdif(i,j,k,iryn) = Fdif(i,j,k,iryn) + (Hg(i,j,k+1,iryn) - Hg(i,j,k,iryn)) * dx2inv(3)
             end do
          end do
       end do
    end do
    !$omp end do nowait
    ! ------- END z-direction -------
    
    !$omp barrier

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
    
    !$omp end parallel

    deallocate(Hg,dpy,dxe,dpe,vsp,vsm)

  end subroutine advance

end module advance_module
