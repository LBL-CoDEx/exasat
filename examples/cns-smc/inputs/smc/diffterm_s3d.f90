module advance_SMC

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

  public S3D_diffterm_1_2

contains

subroutine S3D_diffterm_1_2(lo,hi,ng,dx,Q,Fdif,mu,xi,qx,qy,qz, lam, Ddiag, nprim, ncons, nspecies)

    integer,          intent(in ) :: lo(3),hi(3),ng,ncons, nprim, nspecies
    double precision, intent(in ) :: dx(3)
    double precision, intent(in ) :: Q  (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(in ) :: mu (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: lam (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(in ) :: xi (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng)
    double precision, intent(out) :: qx (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(out) :: qy (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(out) :: qz (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    double precision, intent(out) :: Fdif(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    double precision, intent(in )  :: Ddiag(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nspecies)
 
    double precision, allocatable, dimension(:,:,:) :: vp, dpe, FE
    double precision, allocatable, dimension(:,:,:,:) :: dpy, FY
    ! Ddiag: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! NOT USING ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision ::rhoVc
    integer          :: i,j,k,n, qxn, qyn, qhn, idXn, iryn, qx1, iry1, qh1


    double precision, allocatable, dimension(:,:,:) :: vsm

    double precision :: dxinv(3), divu
    double precision :: dmvxdy,dmwxdz,dmvywzdx
    double precision :: dmuydx,dmwydz,dmuxwzdy
    double precision :: dmuzdx,dmvzdy,dmuxvydz
    double precision :: tauxx,tauyy,tauzz 


    double precision :: M8p(8), M8X(8), mmtmp(8)
    double precision, save, dimension(8,8) :: M8
    ! for 8th-order first derivatives                                                                                   
    double precision,dimension(4),parameter :: D8 = (/ 0.8d0, -0.2d0, 4.d0/105.d0, -1.d0/280.d0 /)

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
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qx(i,j,k,idXn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qx(i,j,k,idXn)*Q(i,j,k,qhn)
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

       !$omp do
    do n=1,nspecies    
       iryn = iry1+n-1

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
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qy(i,j,k,idXn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qy(i,j,k,idXn)*Q(i,j,k,qhn)
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

       !$omp do
    do n=1,nspecies    
       iryn = iry1+n-1

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
                idXn = idX1+n-1
                qhn = qh1+n-1
                FY(i,j,k,n) = Ddiag(i,j,k,n)*qz(i,j,k,idXn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                FE(i,j,k) = FE(i,j,k) + Ddiag(i,j,k,n)*qz(i,j,k,idXn)*Q(i,j,k,qhn)
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

  end subroutine S3D_diffterm_1_2

end module advance_SMC
