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

  subroutine advance (lo,hi,ng,dx,dt,a,b,Unew,U,Uold,eta,alam,courno)

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3),dt,a,b
    integer                       :: i,j,k

    double precision, intent(inout) :: Unew(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in   ) :: U(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in   ) :: Uold(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)

    double precision, intent(inout), optional :: courno

    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: dxinv(3)

    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:,:) :: Q,D,F
    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
    double precision :: divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz,txx,tyy,tzz
    double precision :: mechwork, uxy,uxz,vyz,wzx,wzy,vyx

    double precision, parameter :: OneThird   = 1.0d0/3.0d0
    double precision, parameter :: TwoThirds  = 2.0d0/3.0d0
    double precision, parameter :: FourThirds = 4.0d0/3.0d0

    double precision, parameter :: CENTER = -205.d0/72.d0
    double precision, parameter :: OFF1   =    8.d0/5.d0
    double precision, parameter :: OFF2   =   -0.2d0
    double precision, parameter :: OFF3   =    8.d0/315.d0
    double precision, parameter :: OFF4   =   -1.d0/560.d0

    allocate(Q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6))
    allocate(D(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5))
    allocate(F(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5))

    !$OMP PARALLEL DO PRIVATE(i,j,k,eint,rhoinv)
    do k = lo(3)-ng,hi(3)+ng
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng

             rhoinv     = 1.0d0/U(i,j,k,irho)
             Q(i,j,k,1) = U(i,j,k,irho)
             Q(i,j,k,qu) = U(i,j,k,imx)*rhoinv
             Q(i,j,k,qv) = U(i,j,k,imy)*rhoinv
             Q(i,j,k,qw) = U(i,j,k,imz)*rhoinv

             eint = U(i,j,k,iene)*rhoinv - 0.5d0*(Q(i,j,k,qu)**2 + Q(i,j,k,qv)**2 + Q(i,j,k,qw)**2)

             Q(i,j,k,qpres) = (GAMMA-1.d0)*eint*U(i,j,k,irho)
             Q(i,j,k,6) = eint/CV

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if ( present(courno) ) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,c,courx,coury,courz) REDUCTION(max:courmx,courmy,courmz)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)

                c     = sqrt(GAMMA*Q(i,j,k,qpres)/Q(i,j,k,1))
                courx = ( c+abs(Q(i,j,k,qu)) ) / dx(1)
                coury = ( c+abs(Q(i,j,k,qv)) ) / dx(2) 
                courz = ( c+abs(Q(i,j,k,qw)) ) / dx(3)

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

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = Q(i+1,j,k,qu)
             unp2 = Q(i+2,j,k,qu)
             unp3 = Q(i+3,j,k,qu)
             unp4 = Q(i+4,j,k,qu)

             unm1 = Q(i-1,j,k,qu)
             unm2 = Q(i-2,j,k,qu)
             unm3 = Q(i-3,j,k,qu)
             unm4 = Q(i-4,j,k,qu)

             F(i,j,k,irho)= - &
                   (ALP*(U(i+1,j,k,imx)-U(i-1,j,k,imx)) &
                  + BET*(U(i+2,j,k,imx)-U(i-2,j,k,imx)) &
                  + GAM*(U(i+3,j,k,imx)-U(i-3,j,k,imx)) &
                  + DEL*(U(i+4,j,k,imx)-U(i-4,j,k,imx)))*dxinv(1)

             F(i,j,k,imx)= - &
                   (ALP*(U(i+1,j,k,imx)*unp1-U(i-1,j,k,imx)*unm1 &
                  + (Q(i+1,j,k,qpres)-Q(i-1,j,k,qpres)))               &
                  + BET*(U(i+2,j,k,imx)*unp2-U(i-2,j,k,imx)*unm2 &
                  + (Q(i+2,j,k,qpres)-Q(i-2,j,k,qpres)))               &
                  + GAM*(U(i+3,j,k,imx)*unp3-U(i-3,j,k,imx)*unm3 &
                  + (Q(i+3,j,k,qpres)-Q(i-3,j,k,qpres)))               &
                  + DEL*(U(i+4,j,k,imx)*unp4-U(i-4,j,k,imx)*unm4 &
                  + (Q(i+4,j,k,qpres)-Q(i-4,j,k,qpres))))*dxinv(1)

             F(i,j,k,imy)= - &
                   (ALP*(U(i+1,j,k,imy)*unp1-U(i-1,j,k,imy)*unm1) &
                  + BET*(U(i+2,j,k,imy)*unp2-U(i-2,j,k,imy)*unm2) &
                  + GAM*(U(i+3,j,k,imy)*unp3-U(i-3,j,k,imy)*unm3) &
                  + DEL*(U(i+4,j,k,imy)*unp4-U(i-4,j,k,imy)*unm4))*dxinv(1)

             F(i,j,k,imz)= - &
                   (ALP*(U(i+1,j,k,imz)*unp1-U(i-1,j,k,imz)*unm1) &
                  + BET*(U(i+2,j,k,imz)*unp2-U(i-2,j,k,imz)*unm2) &
                  + GAM*(U(i+3,j,k,imz)*unp3-U(i-3,j,k,imz)*unm3) &
                  + DEL*(U(i+4,j,k,imz)*unp4-U(i-4,j,k,imz)*unm4))*dxinv(1)

             F(i,j,k,iene)= - &
                   (ALP*(U(i+1,j,k,iene)*unp1-U(i-1,j,k,iene)*unm1 &
                  + (Q(i+1,j,k,qpres)*unp1-Q(i-1,j,k,qpres)*unm1))       &
                  + BET*(U(i+2,j,k,iene)*unp2-U(i-2,j,k,iene)*unm2 &
                  + (Q(i+2,j,k,qpres)*unp2-Q(i-2,j,k,qpres)*unm2))       &
                  + GAM*(U(i+3,j,k,iene)*unp3-U(i-3,j,k,iene)*unm3 &
                  + (Q(i+3,j,k,qpres)*unp3-Q(i-3,j,k,qpres)*unm3))       &
                  + DEL*(U(i+4,j,k,iene)*unp4-U(i-4,j,k,iene)*unm4 &
                  + (Q(i+4,j,k,qpres)*unp4-Q(i-4,j,k,qpres)*unm4)))*dxinv(1) 
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = Q(i,j+1,k,qv)
             unp2 = Q(i,j+2,k,qv)
             unp3 = Q(i,j+3,k,qv)
             unp4 = Q(i,j+4,k,qv)

             unm1 = Q(i,j-1,k,qv)
             unm2 = Q(i,j-2,k,qv)
             unm3 = Q(i,j-3,k,qv)
             unm4 = Q(i,j-4,k,qv)

             F(i,j,k,irho)=F(i,j,k,irho) - &
                   (ALP*(U(i,j+1,k,imy)-U(i,j-1,k,imy)) &
                  + BET*(U(i,j+2,k,imy)-U(i,j-2,k,imy)) &
                  + GAM*(U(i,j+3,k,imy)-U(i,j-3,k,imy)) &
                  + DEL*(U(i,j+4,k,imy)-U(i,j-4,k,imy)))*dxinv(2)

             F(i,j,k,imx)=F(i,j,k,imx) - &
                   (ALP*(U(i,j+1,k,imx)*unp1-U(i,j-1,k,imx)*unm1) &
                  + BET*(U(i,j+2,k,imx)*unp2-U(i,j-2,k,imx)*unm2) &
                  + GAM*(U(i,j+3,k,imx)*unp3-U(i,j-3,k,imx)*unm3) &
                  + DEL*(U(i,j+4,k,imx)*unp4-U(i,j-4,k,imx)*unm4))*dxinv(2)

             F(i,j,k,imy)=F(i,j,k,imy) - &
                   (ALP*(U(i,j+1,k,imy)*unp1-U(i,j-1,k,imy)*unm1 &
                  + (Q(i,j+1,k,qpres)-Q(i,j-1,k,qpres)))               &
                  + BET*(U(i,j+2,k,imy)*unp2-U(i,j-2,k,imy)*unm2 &
                  + (Q(i,j+2,k,qpres)-Q(i,j-2,k,qpres)))               &
                  + GAM*(U(i,j+3,k,imy)*unp3-U(i,j-3,k,imy)*unm3 &
                  + (Q(i,j+3,k,qpres)-Q(i,j-3,k,qpres)))               &
                  + DEL*(U(i,j+4,k,imy)*unp4-U(i,j-4,k,imy)*unm4 &
                  + (Q(i,j+4,k,qpres)-Q(i,j-4,k,qpres))))*dxinv(2)

             F(i,j,k,imz)=F(i,j,k,imz) - &
                   (ALP*(U(i,j+1,k,imz)*unp1-U(i,j-1,k,imz)*unm1) &
                  + BET*(U(i,j+2,k,imz)*unp2-U(i,j-2,k,imz)*unm2) &
                  + GAM*(U(i,j+3,k,imz)*unp3-U(i,j-3,k,imz)*unm3) &
                  + DEL*(U(i,j+4,k,imz)*unp4-U(i,j-4,k,imz)*unm4))*dxinv(2)

             F(i,j,k,iene)=F(i,j,k,iene) - &
                   (ALP*(U(i,j+1,k,iene)*unp1-U(i,j-1,k,iene)*unm1 &
                  + (Q(i,j+1,k,qpres)*unp1-Q(i,j-1,k,qpres)*unm1))       &
                  + BET*(U(i,j+2,k,iene)*unp2-U(i,j-2,k,iene)*unm2 &
                  + (Q(i,j+2,k,qpres)*unp2-Q(i,j-2,k,qpres)*unm2))       &
                  + GAM*(U(i,j+3,k,iene)*unp3-U(i,j-3,k,iene)*unm3 &
                  + (Q(i,j+3,k,qpres)*unp3-Q(i,j-3,k,qpres)*unm3))       &
                  + DEL*(U(i,j+4,k,iene)*unp4-U(i,j-4,k,iene)*unm4 &
                  + (Q(i,j+4,k,qpres)*unp4-Q(i,j-4,k,qpres)*unm4)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             unp1 = Q(i,j,k+1,qw)
             unp2 = Q(i,j,k+2,qw)
             unp3 = Q(i,j,k+3,qw)
             unp4 = Q(i,j,k+4,qw)

             unm1 = Q(i,j,k-1,qw)
             unm2 = Q(i,j,k-2,qw)
             unm3 = Q(i,j,k-3,qw)
             unm4 = Q(i,j,k-4,qw)

             F(i,j,k,irho)=F(i,j,k,irho) - &
                   (ALP*(U(i,j,k+1,imz)-U(i,j,k-1,imz)) &
                  + BET*(U(i,j,k+2,imz)-U(i,j,k-2,imz)) &
                  + GAM*(U(i,j,k+3,imz)-U(i,j,k-3,imz)) &
                  + DEL*(U(i,j,k+4,imz)-U(i,j,k-4,imz)))*dxinv(3)

             F(i,j,k,imx)=F(i,j,k,imx) - &
                   (ALP*(U(i,j,k+1,imx)*unp1-U(i,j,k-1,imx)*unm1) &
                  + BET*(U(i,j,k+2,imx)*unp2-U(i,j,k-2,imx)*unm2) &
                  + GAM*(U(i,j,k+3,imx)*unp3-U(i,j,k-3,imx)*unm3) &
                  + DEL*(U(i,j,k+4,imx)*unp4-U(i,j,k-4,imx)*unm4))*dxinv(3)

             F(i,j,k,imy)=F(i,j,k,imy) - &
                   (ALP*(U(i,j,k+1,imy)*unp1-U(i,j,k-1,imy)*unm1) &
                  + BET*(U(i,j,k+2,imy)*unp2-U(i,j,k-2,imy)*unm2) &
                  + GAM*(U(i,j,k+3,imy)*unp3-U(i,j,k-3,imy)*unm3) &
                  + DEL*(U(i,j,k+4,imy)*unp4-U(i,j,k-4,imy)*unm4))*dxinv(3)

             F(i,j,k,imz)=F(i,j,k,imz) - &
                   (ALP*(U(i,j,k+1,imz)*unp1-U(i,j,k-1,imz)*unm1 &
                  + (Q(i,j,k+1,qpres)-Q(i,j,k-1,qpres)))               &
                  + BET*(U(i,j,k+2,imz)*unp2-U(i,j,k-2,imz)*unm2 &
                  + (Q(i,j,k+2,qpres)-Q(i,j,k-2,qpres)))               &
                  + GAM*(U(i,j,k+3,imz)*unp3-U(i,j,k-3,imz)*unm3 &
                  + (Q(i,j,k+3,qpres)-Q(i,j,k-3,qpres)))               &
                  + DEL*(U(i,j,k+4,imz)*unp4-U(i,j,k-4,imz)*unm4 &
                  + (Q(i,j,k+4,qpres)-Q(i,j,k-4,qpres))))*dxinv(3)

             F(i,j,k,iene)=F(i,j,k,iene) - &
                   (ALP*(U(i,j,k+1,iene)*unp1-U(i,j,k-1,iene)*unm1 &
                  + (Q(i,j,k+1,qpres)*unp1-Q(i,j,k-1,qpres)*unm1))       &
                  + BET*(U(i,j,k+2,iene)*unp2-U(i,j,k-2,iene)*unm2 &
                  + (Q(i,j,k+2,qpres)*unp2-Q(i,j,k-2,qpres)*unm2))       &
                  + GAM*(U(i,j,k+3,iene)*unp3-U(i,j,k-3,iene)*unm3 &
                  + (Q(i,j,k+3,qpres)*unp3-Q(i,j,k-3,qpres)*unm3))       &
                  + DEL*(U(i,j,k+4,iene)*unp4-U(i,j,k-4,iene)*unm4 &
                  + (Q(i,j,k+4,qpres)*unp4-Q(i,j,k-4,qpres)*unm4)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    D(:,:,:,irho) = 0.0d0

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1),hi(1)

             ux(i,j,k)= &
                   (ALP*(Q(i+1,j,k,qu)-Q(i-1,j,k,qu)) &
                  + BET*(Q(i+2,j,k,qu)-Q(i-2,j,k,qu)) &
                  + GAM*(Q(i+3,j,k,qu)-Q(i-3,j,k,qu)) &
                  + DEL*(Q(i+4,j,k,qu)-Q(i-4,j,k,qu)))*dxinv(1)

             vx(i,j,k)= &
                   (ALP*(Q(i+1,j,k,qv)-Q(i-1,j,k,qv)) &
                  + BET*(Q(i+2,j,k,qv)-Q(i-2,j,k,qv)) &
                  + GAM*(Q(i+3,j,k,qv)-Q(i-3,j,k,qv)) &
                  + DEL*(Q(i+4,j,k,qv)-Q(i-4,j,k,qv)))*dxinv(1)

             wx(i,j,k)= &
                   (ALP*(Q(i+1,j,k,qw)-Q(i-1,j,k,qw)) &
                  + BET*(Q(i+2,j,k,qw)-Q(i-2,j,k,qw)) &
                  + GAM*(Q(i+3,j,k,qw)-Q(i-3,j,k,qw)) &
                  + DEL*(Q(i+4,j,k,qw)-Q(i-4,j,k,qw)))*dxinv(1)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2),hi(2)   
          do i=lo(1)-ng,hi(1)+ng

             uy(i,j,k)= &
                   (ALP*(Q(i,j+1,k,qu)-Q(i,j-1,k,qu)) &
                  + BET*(Q(i,j+2,k,qu)-Q(i,j-2,k,qu)) &
                  + GAM*(Q(i,j+3,k,qu)-Q(i,j-3,k,qu)) &
                  + DEL*(Q(i,j+4,k,qu)-Q(i,j-4,k,qu)))*dxinv(2)

             vy(i,j,k)= &
                   (ALP*(Q(i,j+1,k,qv)-Q(i,j-1,k,qv)) &
                  + BET*(Q(i,j+2,k,qv)-Q(i,j-2,k,qv)) &
                  + GAM*(Q(i,j+3,k,qv)-Q(i,j-3,k,qv)) &
                  + DEL*(Q(i,j+4,k,qv)-Q(i,j-4,k,qv)))*dxinv(2)

             wy(i,j,k)= &
                   (ALP*(Q(i,j+1,k,qw)-Q(i,j-1,k,qw)) &
                  + BET*(Q(i,j+2,k,qw)-Q(i,j-2,k,qw)) &
                  + GAM*(Q(i,j+3,k,qw)-Q(i,j-3,k,qw)) &
                  + DEL*(Q(i,j+4,k,qw)-Q(i,j-4,k,qw)))*dxinv(2)
          enddo
       enddo
    enddo
    !$OMP END DO NOWAIT

    !$OMP DO
    do k=lo(3),hi(3)
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             uz(i,j,k)= &
                   (ALP*(Q(i,j,k+1,qu)-Q(i,j,k-1,qu)) &
                  + BET*(Q(i,j,k+2,qu)-Q(i,j,k-2,qu)) &
                  + GAM*(Q(i,j,k+3,qu)-Q(i,j,k-3,qu)) &
                  + DEL*(Q(i,j,k+4,qu)-Q(i,j,k-4,qu)))*dxinv(3)

             vz(i,j,k)= &
                   (ALP*(Q(i,j,k+1,qv)-Q(i,j,k-1,qv)) &
                  + BET*(Q(i,j,k+2,qv)-Q(i,j,k-2,qv)) &
                  + GAM*(Q(i,j,k+3,qv)-Q(i,j,k-3,qv)) &
                  + DEL*(Q(i,j,k+4,qv)-Q(i,j,k-4,qv)))*dxinv(3)

             wz(i,j,k)= &
                   (ALP*(Q(i,j,k+1,qw)-Q(i,j,k-1,qw)) &
                  + BET*(Q(i,j,k+2,qw)-Q(i,j,k-2,qw)) &
                  + GAM*(Q(i,j,k+3,qw)-Q(i,j,k-3,qw)) &
                  + DEL*(Q(i,j,k+4,qw)-Q(i,j,k-4,qw)))*dxinv(3)
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DO PRIVATE(i,j,k,uxx,uyy,uzz,vyx,wzx)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             uxx = (CENTER*Q(i,j,k,qu) &
                  + OFF1*(Q(i+1,j,k,qu)+Q(i-1,j,k,qu)) &
                  + OFF2*(Q(i+2,j,k,qu)+Q(i-2,j,k,qu)) &
                  + OFF3*(Q(i+3,j,k,qu)+Q(i-3,j,k,qu)) &
                  + OFF4*(Q(i+4,j,k,qu)+Q(i-4,j,k,qu)))*dxinv(1)**2

             uyy = (CENTER*Q(i,j,k,qu) &
                  + OFF1*(Q(i,j+1,k,qu)+Q(i,j-1,k,qu)) &
                  + OFF2*(Q(i,j+2,k,qu)+Q(i,j-2,k,qu)) &
                  + OFF3*(Q(i,j+3,k,qu)+Q(i,j-3,k,qu)) &
                  + OFF4*(Q(i,j+4,k,qu)+Q(i,j-4,k,qu)))*dxinv(2)**2

             uzz = (CENTER*Q(i,j,k,qu) &
                  + OFF1*(Q(i,j,k+1,qu)+Q(i,j,k-1,qu)) &
                  + OFF2*(Q(i,j,k+2,qu)+Q(i,j,k-2,qu)) &
                  + OFF3*(Q(i,j,k+3,qu)+Q(i,j,k-3,qu)) &
                  + OFF4*(Q(i,j,k+4,qu)+Q(i,j,k-4,qu)))*dxinv(3)**2

             vyx = (ALP*(vy(i+1,j,k)-vy(i-1,j,k)) &
                  + BET*(vy(i+2,j,k)-vy(i-2,j,k)) &
                  + GAM*(vy(i+3,j,k)-vy(i-3,j,k)) &
                  + DEL*(vy(i+4,j,k)-vy(i-4,j,k)))*dxinv(1)

             wzx = (ALP*(wz(i+1,j,k)-wz(i-1,j,k)) &
                  + BET*(wz(i+2,j,k)-wz(i-2,j,k)) &
                  + GAM*(wz(i+3,j,k)-wz(i-3,j,k)) &
                  + DEL*(wz(i+4,j,k)-wz(i-4,j,k)))*dxinv(1)

             D(i,j,k,imx) = eta*(FourThirds*uxx + uyy + uzz + OneThird*(vyx+wzx))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,vxx,vyy,vzz,uxy,wzy)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             vxx = (CENTER*Q(i,j,k,qv) &
                  + OFF1*(Q(i+1,j,k,qv)+Q(i-1,j,k,qv)) &
                  + OFF2*(Q(i+2,j,k,qv)+Q(i-2,j,k,qv)) &
                  + OFF3*(Q(i+3,j,k,qv)+Q(i-3,j,k,qv)) &
                  + OFF4*(Q(i+4,j,k,qv)+Q(i-4,j,k,qv)))*dxinv(1)**2

             vyy = (CENTER*Q(i,j,k,qv) &
                  + OFF1*(Q(i,j+1,k,qv)+Q(i,j-1,k,qv)) &
                  + OFF2*(Q(i,j+2,k,qv)+Q(i,j-2,k,qv)) &
                  + OFF3*(Q(i,j+3,k,qv)+Q(i,j-3,k,qv)) &
                  + OFF4*(Q(i,j+4,k,qv)+Q(i,j-4,k,qv)))*dxinv(2)**2

             vzz = (CENTER*Q(i,j,k,qv) &
                  + OFF1*(Q(i,j,k+1,qv)+Q(i,j,k-1,qv)) &
                  + OFF2*(Q(i,j,k+2,qv)+Q(i,j,k-2,qv)) &
                  + OFF3*(Q(i,j,k+3,qv)+Q(i,j,k-3,qv)) &
                  + OFF4*(Q(i,j,k+4,qv)+Q(i,j,k-4,qv)))*dxinv(3)**2

             uxy = (ALP*(ux(i,j+1,k)-ux(i,j-1,k)) &
                  + BET*(ux(i,j+2,k)-ux(i,j-2,k)) &
                  + GAM*(ux(i,j+3,k)-ux(i,j-3,k)) &
                  + DEL*(ux(i,j+4,k)-ux(i,j-4,k)))*dxinv(2)

             wzy = (ALP*(wz(i,j+1,k)-wz(i,j-1,k)) &
                  + BET*(wz(i,j+2,k)-wz(i,j-2,k)) &
                  + GAM*(wz(i,j+3,k)-wz(i,j-3,k)) &
                  + DEL*(wz(i,j+4,k)-wz(i,j-4,k)))*dxinv(2)

             D(i,j,k,imy) = eta*(vxx + FourThirds*vyy + vzz + OneThird*(uxy+wzy))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,wxx,wyy,wzz,uxz,vyz)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             wxx = (CENTER*Q(i,j,k,qw) &
                  + OFF1*(Q(i+1,j,k,qw)+Q(i-1,j,k,qw)) &
                  + OFF2*(Q(i+2,j,k,qw)+Q(i-2,j,k,qw)) &
                  + OFF3*(Q(i+3,j,k,qw)+Q(i-3,j,k,qw)) &
                  + OFF4*(Q(i+4,j,k,qw)+Q(i-4,j,k,qw)))*dxinv(1)**2

             wyy = (CENTER*Q(i,j,k,qw) &
                  + OFF1*(Q(i,j+1,k,qw)+Q(i,j-1,k,qw)) &
                  + OFF2*(Q(i,j+2,k,qw)+Q(i,j-2,k,qw)) &
                  + OFF3*(Q(i,j+3,k,qw)+Q(i,j-3,k,qw)) &
                  + OFF4*(Q(i,j+4,k,qw)+Q(i,j-4,k,qw)))*dxinv(2)**2

             wzz = (CENTER*Q(i,j,k,qw) &
                  + OFF1*(Q(i,j,k+1,qw)+Q(i,j,k-1,qw)) &
                  + OFF2*(Q(i,j,k+2,qw)+Q(i,j,k-2,qw)) &
                  + OFF3*(Q(i,j,k+3,qw)+Q(i,j,k-3,qw)) &
                  + OFF4*(Q(i,j,k+4,qw)+Q(i,j,k-4,qw)))*dxinv(3)**2

             uxz = (ALP*(ux(i,j,k+1)-ux(i,j,k-1)) &
                  + BET*(ux(i,j,k+2)-ux(i,j,k-2)) &
                  + GAM*(ux(i,j,k+3)-ux(i,j,k-3)) &
                  + DEL*(ux(i,j,k+4)-ux(i,j,k-4)))*dxinv(3)

             vyz = (ALP*(vy(i,j,k+1)-vy(i,j,k-1)) &
                  + BET*(vy(i,j,k+2)-vy(i,j,k-2)) &
                  + GAM*(vy(i,j,k+3)-vy(i,j,k-3)) &
                  + DEL*(vy(i,j,k+4)-vy(i,j,k-4)))*dxinv(3)

             D(i,j,k,imz) = eta*(wxx + wyy + FourThirds*wzz + OneThird*(uxz+vyz))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,txx,tyy,tzz) &
    !$OMP PRIVATE(divu,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,mechwork)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             txx = (CENTER*Q(i,j,k,6) &
                  + OFF1*(Q(i+1,j,k,6)+Q(i-1,j,k,6)) &
                  + OFF2*(Q(i+2,j,k,6)+Q(i-2,j,k,6)) &
                  + OFF3*(Q(i+3,j,k,6)+Q(i-3,j,k,6)) &
                  + OFF4*(Q(i+4,j,k,6)+Q(i-4,j,k,6)))*dxinv(1)**2

             tyy = (CENTER*Q(i,j,k,6) &
                  + OFF1*(Q(i,j+1,k,6)+Q(i,j-1,k,6)) &
                  + OFF2*(Q(i,j+2,k,6)+Q(i,j-2,k,6)) &
                  + OFF3*(Q(i,j+3,k,6)+Q(i,j-3,k,6)) &
                  + OFF4*(Q(i,j+4,k,6)+Q(i,j-4,k,6)))*dxinv(2)**2

             tzz = (CENTER*Q(i,j,k,6) &
                  + OFF1*(Q(i,j,k+1,6)+Q(i,j,k-1,6)) &
                  + OFF2*(Q(i,j,k+2,6)+Q(i,j,k-2,6)) &
                  + OFF3*(Q(i,j,k+3,6)+Q(i,j,k-3,6)) &
                  + OFF4*(Q(i,j,k+4,6)+Q(i,j,k-4,6)))*dxinv(3)**2

             divu  = TwoThirds*(ux(i,j,k)+vy(i,j,k)+wz(i,j,k))
             tauxx = 2.d0*ux(i,j,k) - divu
             tauyy = 2.d0*vy(i,j,k) - divu
             tauzz = 2.d0*wz(i,j,k) - divu
             tauxy = uy(i,j,k)+vx(i,j,k)
             tauxz = uz(i,j,k)+wx(i,j,k)
             tauyz = vz(i,j,k)+wy(i,j,k)

             mechwork = tauxx*ux(i,j,k) + &
                        tauyy*vy(i,j,k) + &
                        tauzz*wz(i,j,k) + tauxy**2+tauxz**2+tauyz**2

             mechwork = eta*mechwork &
                  + D(i,j,k,imx)*Q(i,j,k,qu) &
                  + D(i,j,k,imy)*Q(i,j,k,qv) &
                  + D(i,j,k,imz)*Q(i,j,k,qw)

             D(i,j,k,iene) = alam*(txx+tyy+tzz) + mechwork
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(ux)
    deallocate(uy)
    deallocate(uz)
    deallocate(vx)
    deallocate(vy)
    deallocate(vz)
    deallocate(wx)
    deallocate(wy)
    deallocate(wz)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             Unew(i,j,k,irho) = a * U(i,j,k,irho) + b * (Uold(i,j,k,irho) + &
                                dt * (D(i,j,k,irho) + F(i,j,k,irho)))
             Unew(i,j,k,imx)  = a * U(i,j,k,imx) + b * (Uold(i,j,k,imx) + &
                                dt * (D(i,j,k,imx) + F(i,j,k,imx)))
             Unew(i,j,k,imy)  = a * U(i,j,k,imy) + b * (Uold(i,j,k,imy) + &
                                dt * (D(i,j,k,imy) + F(i,j,k,imy)))
             Unew(i,j,k,imz)  = a * U(i,j,k,imz) + b * (Uold(i,j,k,imz) + &
                                dt * (D(i,j,k,imz) + F(i,j,k,imz)))
             Unew(i,j,k,iene) = a * U(i,j,k,iene) + b * (Uold(i,j,k,iene) + &
                                dt * (D(i,j,k,iene) + F(i,j,k,iene)))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(Q)
    deallocate(D)
    deallocate(F)

  end subroutine advance 

end module advance_module

