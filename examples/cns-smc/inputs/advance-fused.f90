 module advance_module

  use bl_error_module
  use multifab_module
  use bl_prof_module

  implicit none

  private
  !
  ! These index constants are shared with the initial data routine.
  !
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


  subroutine advance (U,dt,dx,cfl,eta,alam)

    use bl_prof_module

    type(multifab),   intent(inout) :: U
    double precision, intent(out  ) :: dt
    double precision, intent(in   ) :: dx(U%dim), cfl, eta, alam

    integer          :: lo(U%dim), hi(U%dim), i, j, k, m, n, nc, ng
    double precision :: courno, courno_proc, a, b
    type(layout)     :: la
    type(multifab)   :: Unew, Q

    double precision, pointer, dimension(:,:,:,:) :: up, unp, un2p, qp
    !
    ! Some arithmetic constants.
    !
    double precision, parameter :: Zero          = 0.d0
    double precision, parameter :: One           = 1.d0
    double precision, parameter :: OneThird      = 1.d0/3.d0
    double precision, parameter :: TwoThirds     = 2.d0/3.d0
    double precision, parameter :: OneQuarter    = 1.d0/4.d0
    double precision, parameter :: ThreeQuarters = 3.d0/4.d0

    type(bl_prof_timer), save :: bpt_advance

    call build(bpt_advance, "bpt_advance")

    nc = ncomp(U)
    ng = nghost(U)
    la = get_layout(U)
    !
    ! Sync U prior to calculating D & F.
    !

    call multifab_fill_boundary(U)

    call multifab_build(Q,     la, nc+1, ng)
    call multifab_build(Unew,  la, nc,   ng)
    call multifab_build(Unew2, la, nc,   ng)
    !
    ! Calculate primitive variables based on U.
    !
    ! Also calculate courno so we can set "dt".
    !
    courno_proc = 1.0d-50


    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng,courno=courno_proc)
    end do

    call parallel_reduce(courno, courno_proc, MPI_MAX)

    dt = cfl / courno

    if ( parallel_IOProcessor() ) then
       print*, "dt,courno", dt, courno
    end if

    !
    ! Calculate D and F at time N.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       unp => dataptr(Unew,n)
       up  => dataptr(U,   n)
       qp  => dataptr(Q,   n)
       fp  => dataptr(F,   n)
       dp  => dataptr(D,   n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       a = Zero
       b = One

       call updateU(lo,hi,ng,dx,dt,a,b,unp,up,unp,qp,ETA,ALAM)
    end do

    !
    ! Sync U^1/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate primitive variables based on U^1/3.
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng)
    end do
    !
    ! Calculate D and F at time N+1/3.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       un2p => dataptr(Unew2,n)
       unp  => dataptr(Unew, n)
       up   => dataptr(U,    n)
       qp   => dataptr(Q,    n)
       fp   => dataptr(F,    n)
       dp   => dataptr(D,    n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       a = ThreeQuarters
       b = OneQuarter

       call updateU(lo,hi,ng,dx,dt,a,b,un2p,unp,up,qp,ETA,ALAM)
    end do
    !
    ! Sync U^2/3 prior to calculating D & F.
    !
    call multifab_fill_boundary(Unew)
    !
    ! Calculate primitive variables based on U^2/3.
    !
    do n=1,nboxes(Q)
       if ( remote(Q,n) ) cycle

       up => dataptr(Unew,n)
       qp => dataptr(Q,   n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim(lo,hi,up,qp,dx,ng)
    end do
    !
    ! Calculate D and F at time N+2/3.
    !
    do n=1,nboxes(U)
       if ( remote(U,n) ) cycle

       un2p => dataptr(Unew2,n)
       unp  => dataptr(Unew, n)
       up   => dataptr(U,    n)
       qp   => dataptr(Q,    n)
       fp   => dataptr(F,    n)
       dp   => dataptr(D,    n)

       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       a = OneThird
       b = TwoThirds

       call updateU(lo,hi,ng,dx,dt,a,b,up,un2p,up,qp,ETA,ALAM)
    end do

    call destroy(Unew)
    call destroy(Q)
    call destroy(F)
    call destroy(D)

    call destroy(bpt_advance)

  end subroutine advance



  subroutine ctoprim (lo,hi,u,q,dx,ng,courno)

    use bl_prof_module

    integer,          intent(in ) :: lo(3), hi(3), ng
    double precision, intent(in ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,5)
    double precision, intent(in ) :: dx(3)
    double precision, intent(out) :: q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,6)

    double precision, intent(inout), optional :: courno

    integer          :: i, j, k
    double precision :: c, eint, courx, coury, courz, courmx, courmy, courmz, rhoinv
    double precision :: dx1inv, dx2inv, dx3inv, CVinv

    double precision, parameter :: GAMMA = 1.4d0
    double precision, parameter :: CV    = 8.3333333333d6

    type(bl_prof_timer), save :: bpt_ctoprim, bpt_ctoprim_loop1, bpt_ctoprim_loop2

    call build(bpt_ctoprim, "bpt_ctoprim")

    CVinv = 1.0d0 / CV

    courmx = -Huge(courmx)
    courmy = -Huge(courmy)
    courmz = -Huge(courmz)

    dx1inv = 1.0d0 / dx(1)
    dx2inv = 1.0d0 / dx(2)
    dx3inv = 1.0d0 / dx(3)

    call build(bpt_ctoprim_loop1, "bpt_ctoprim_loop1")
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
             q(i,j,k,6) = eint * CVinv

             if ( present(courno) .and. &
                  lo(1) <= i .and. i <= hi(1) .and. &
                  lo(2) <= j .and. j <= hi(2) .and. &
                  lo(3) <= k .and. k <= hi(3) ) then

                c     = sqrt(GAMMA*q(i,j,k,5)/q(i,j,k,1))
                courx = ( c+abs(q(i,j,k,2)) ) * dx1inv
                coury = ( c+abs(q(i,j,k,3)) ) * dx2inv
                courz = ( c+abs(q(i,j,k,4)) ) * dx3inv

                courmx = max( courmx, courx )
                courmy = max( courmy, coury )
                courmz = max( courmz, courz )

             end if

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_ctoprim_loop1)

    !
    ! Compute running max of Courant number over grids.
    !
    if ( present(courno) ) then
       courno = max( courmx, courmy, courmz , courno )
    end if

    call destroy(bpt_ctoprim)

  end subroutine ctoprim


  subroutine updateU (lo,hi,ng,dx,dt,a,b,unew,u,uold,q,eta,alam)

    use bl_prof_module

    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: dx(3),dt,a,b
    double precision, intent(out) :: unew(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) :: u(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) :: uold(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
    double precision, intent(in ) :: q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
    double precision, intent(in ) :: eta, alam

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz

    double precision :: unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4
    double precision :: flux_irho,flux_imx,flux_imy,flux_imz,flux_iene
    double precision :: difflux_imx,difflux_imy,difflux_imz,difflux_iene
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

    type(bl_prof_timer), save :: bpt_updateU, bpt_updateU_loop1, bpt_updateU_loop2


    call build(bpt_updateU, "bpt_updateU")

    allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))
    allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng))

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    !
    ! Compute first derivative terms over entire region including ghost zones
    !

    call build(bpt_updateU_loop1, "bpt_updateU_loop1")
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             if (i >= lo(1) .and. i <= hi(1)) then
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
             end if

             if (j >= lo(2) .and. j <= hi(2)) then
                uy(i,j,k)= &
                      (ALP*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                     + BET*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                     + GAM*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                     + DEL*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))*dxinv(2)

                vy(i,j,k)= &
                      (ALP*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                     + BET*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                     + GAM*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                     + DEL*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))*dxinv(2)

                wy(i,j,k)= &
                      (ALP*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                     + BET*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                     + GAM*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                     + DEL*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))*dxinv(2)
             end if

             if (k >= lo(3) .and. k <= hi(3)) then
                uz(i,j,k)= &
                      (ALP*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                     + BET*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                     + GAM*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                     + DEL*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))*dxinv(3)

                vz(i,j,k)= &
                      (ALP*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                     + BET*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                     + GAM*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                     + DEL*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))*dxinv(3)

                wz(i,j,k)= &
                      (ALP*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                     + BET*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                     + GAM*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                     + DEL*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))*dxinv(3)
             end if
          enddo
       enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call destroy(bpt_updateU_loop1)

    call build(bpt_updateU_loop2, "bpt_updateU_loop2")
    !$OMP PARALLEL DO PRIVATE(i,j,k,unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4) &
    !$OMP PRIVATE(uxx,uyy,uzz,vyx,wzx,vxx,vyy,vzz,uxy,wzy,wxx,wyy,wzz,uxz,vyz) &
    !$OMP PRIVATE(txx,tyy,tzz,divu,tauxx,tauyy,tauzz,tauxy,tauxz,tauyz,mechwork)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             !
             ! Compute diffterm first
             !

             uxx = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i+1,j,k,qu)+q(i-1,j,k,qu)) &
                  + OFF2*(q(i+2,j,k,qu)+q(i-2,j,k,qu)) &
                  + OFF3*(q(i+3,j,k,qu)+q(i-3,j,k,qu)) &
                  + OFF4*(q(i+4,j,k,qu)+q(i-4,j,k,qu)))*dxinv(1)**2

             uyy = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j+1,k,qu)+q(i,j-1,k,qu)) &
                  + OFF2*(q(i,j+2,k,qu)+q(i,j-2,k,qu)) &
                  + OFF3*(q(i,j+3,k,qu)+q(i,j-3,k,qu)) &
                  + OFF4*(q(i,j+4,k,qu)+q(i,j-4,k,qu)))*dxinv(2)**2

             uzz = (CENTER*q(i,j,k,qu) &
                  + OFF1*(q(i,j,k+1,qu)+q(i,j,k-1,qu)) &
                  + OFF2*(q(i,j,k+2,qu)+q(i,j,k-2,qu)) &
                  + OFF3*(q(i,j,k+3,qu)+q(i,j,k-3,qu)) &
                  + OFF4*(q(i,j,k+4,qu)+q(i,j,k-4,qu)))*dxinv(3)**2

             vyx = (ALP*(vy(i+1,j,k)-vy(i-1,j,k)) &
                  + BET*(vy(i+2,j,k)-vy(i-2,j,k)) &
                  + GAM*(vy(i+3,j,k)-vy(i-3,j,k)) &
                  + DEL*(vy(i+4,j,k)-vy(i-4,j,k)))*dxinv(1)

             wzx = (ALP*(wz(i+1,j,k)-wz(i-1,j,k)) &
                  + BET*(wz(i+2,j,k)-wz(i-2,j,k)) &
                  + GAM*(wz(i+3,j,k)-wz(i-3,j,k)) &
                  + DEL*(wz(i+4,j,k)-wz(i-4,j,k)))*dxinv(1)

             difflux_imx = eta*(FourThirds*uxx + uyy + uzz + OneThird*(vyx+wzx))

             vxx = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i+1,j,k,qv)+q(i-1,j,k,qv)) &
                  + OFF2*(q(i+2,j,k,qv)+q(i-2,j,k,qv)) &
                  + OFF3*(q(i+3,j,k,qv)+q(i-3,j,k,qv)) &
                  + OFF4*(q(i+4,j,k,qv)+q(i-4,j,k,qv)))*dxinv(1)**2

             vyy = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j+1,k,qv)+q(i,j-1,k,qv)) &
                  + OFF2*(q(i,j+2,k,qv)+q(i,j-2,k,qv)) &
                  + OFF3*(q(i,j+3,k,qv)+q(i,j-3,k,qv)) &
                  + OFF4*(q(i,j+4,k,qv)+q(i,j-4,k,qv)))*dxinv(2)**2

             vzz = (CENTER*q(i,j,k,qv) &
                  + OFF1*(q(i,j,k+1,qv)+q(i,j,k-1,qv)) &
                  + OFF2*(q(i,j,k+2,qv)+q(i,j,k-2,qv)) &
                  + OFF3*(q(i,j,k+3,qv)+q(i,j,k-3,qv)) &
                  + OFF4*(q(i,j,k+4,qv)+q(i,j,k-4,qv)))*dxinv(3)**2

             uxy = (ALP*(ux(i,j+1,k)-ux(i,j-1,k)) &
                  + BET*(ux(i,j+2,k)-ux(i,j-2,k)) &
                  + GAM*(ux(i,j+3,k)-ux(i,j-3,k)) &
                  + DEL*(ux(i,j+4,k)-ux(i,j-4,k)))*dxinv(2)

             wzy = (ALP*(wz(i,j+1,k)-wz(i,j-1,k)) &
                  + BET*(wz(i,j+2,k)-wz(i,j-2,k)) &
                  + GAM*(wz(i,j+3,k)-wz(i,j-3,k)) &
                  + DEL*(wz(i,j+4,k)-wz(i,j-4,k)))*dxinv(2)

             difflux_imy = eta*(vxx + FourThirds*vyy + vzz + OneThird*(uxy+wzy))

             wxx = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i+1,j,k,qw)+q(i-1,j,k,qw)) &
                  + OFF2*(q(i+2,j,k,qw)+q(i-2,j,k,qw)) &
                  + OFF3*(q(i+3,j,k,qw)+q(i-3,j,k,qw)) &
                  + OFF4*(q(i+4,j,k,qw)+q(i-4,j,k,qw)))*dxinv(1)**2

             wyy = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j+1,k,qw)+q(i,j-1,k,qw)) &
                  + OFF2*(q(i,j+2,k,qw)+q(i,j-2,k,qw)) &
                  + OFF3*(q(i,j+3,k,qw)+q(i,j-3,k,qw)) &
                  + OFF4*(q(i,j+4,k,qw)+q(i,j-4,k,qw)))*dxinv(2)**2

             wzz = (CENTER*q(i,j,k,qw) &
                  + OFF1*(q(i,j,k+1,qw)+q(i,j,k-1,qw)) &
                  + OFF2*(q(i,j,k+2,qw)+q(i,j,k-2,qw)) &
                  + OFF3*(q(i,j,k+3,qw)+q(i,j,k-3,qw)) &
                  + OFF4*(q(i,j,k+4,qw)+q(i,j,k-4,qw)))*dxinv(3)**2

             uxz = (ALP*(ux(i,j,k+1)-ux(i,j,k-1)) &
                  + BET*(ux(i,j,k+2)-ux(i,j,k-2)) &
                  + GAM*(ux(i,j,k+3)-ux(i,j,k-3)) &
                  + DEL*(ux(i,j,k+4)-ux(i,j,k-4)))*dxinv(3)

             vyz = (ALP*(vy(i,j,k+1)-vy(i,j,k-1)) &
                  + BET*(vy(i,j,k+2)-vy(i,j,k-2)) &
                  + GAM*(vy(i,j,k+3)-vy(i,j,k-3)) &
                  + DEL*(vy(i,j,k+4)-vy(i,j,k-4)))*dxinv(3)

             difflux_imz = eta*(wxx + wyy + FourThirds*wzz + OneThird*(uxz+vyz))

             txx = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i+1,j,k,6)+q(i-1,j,k,6)) &
                  + OFF2*(q(i+2,j,k,6)+q(i-2,j,k,6)) &
                  + OFF3*(q(i+3,j,k,6)+q(i-3,j,k,6)) &
                  + OFF4*(q(i+4,j,k,6)+q(i-4,j,k,6)))*dxinv(1)**2

             tyy = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j+1,k,6)+q(i,j-1,k,6)) &
                  + OFF2*(q(i,j+2,k,6)+q(i,j-2,k,6)) &
                  + OFF3*(q(i,j+3,k,6)+q(i,j-3,k,6)) &
                  + OFF4*(q(i,j+4,k,6)+q(i,j-4,k,6)))*dxinv(2)**2

             tzz = (CENTER*q(i,j,k,6) &
                  + OFF1*(q(i,j,k+1,6)+q(i,j,k-1,6)) &
                  + OFF2*(q(i,j,k+2,6)+q(i,j,k-2,6)) &
                  + OFF3*(q(i,j,k+3,6)+q(i,j,k-3,6)) &
                  + OFF4*(q(i,j,k+4,6)+q(i,j,k-4,6)))*dxinv(3)**2

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
                  + difflux_imx*q(i,j,k,qu) &
                  + difflux_imy*q(i,j,k,qv) &
                  + difflux_imz*q(i,j,k,qw)

             difflux_iene = alam*(txx+tyy+tzz) + mechwork

             !
             ! Compute hypterm second
             !

             unp1 = q(i+1,j,k,qu)
             unp2 = q(i+2,j,k,qu)
             unp3 = q(i+3,j,k,qu)
             unp4 = q(i+4,j,k,qu)

             unm1 = q(i-1,j,k,qu)
             unm2 = q(i-2,j,k,qu)
             unm3 = q(i-3,j,k,qu)
             unm4 = q(i-4,j,k,qu)

             flux_irho = - &
                   (ALP*(u(i+1,j,k,imx)-u(i-1,j,k,imx)) &
                  + BET*(u(i+2,j,k,imx)-u(i-2,j,k,imx)) &
                  + GAM*(u(i+3,j,k,imx)-u(i-3,j,k,imx)) &
                  + DEL*(u(i+4,j,k,imx)-u(i-4,j,k,imx)))*dxinv(1)

             flux_imx= - &
                   (ALP*(u(i+1,j,k,imx)*unp1-u(i-1,j,k,imx)*unm1 &
                  + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))               &
                  + BET*(u(i+2,j,k,imx)*unp2-u(i-2,j,k,imx)*unm2 &
                  + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))               &
                  + GAM*(u(i+3,j,k,imx)*unp3-u(i-3,j,k,imx)*unm3 &
                  + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))               &
                  + DEL*(u(i+4,j,k,imx)*unp4-u(i-4,j,k,imx)*unm4 &
                  + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres))))*dxinv(1)

             flux_imy= - &
                   (ALP*(u(i+1,j,k,imy)*unp1-u(i-1,j,k,imy)*unm1) &
                  + BET*(u(i+2,j,k,imy)*unp2-u(i-2,j,k,imy)*unm2) &
                  + GAM*(u(i+3,j,k,imy)*unp3-u(i-3,j,k,imy)*unm3) &
                  + DEL*(u(i+4,j,k,imy)*unp4-u(i-4,j,k,imy)*unm4))*dxinv(1)

             flux_imz= - &
                   (ALP*(u(i+1,j,k,imz)*unp1-u(i-1,j,k,imz)*unm1) &
                  + BET*(u(i+2,j,k,imz)*unp2-u(i-2,j,k,imz)*unm2) &
                  + GAM*(u(i+3,j,k,imz)*unp3-u(i-3,j,k,imz)*unm3) &
                  + DEL*(u(i+4,j,k,imz)*unp4-u(i-4,j,k,imz)*unm4))*dxinv(1)

             flux_iene= - &
                   (ALP*(u(i+1,j,k,iene)*unp1-u(i-1,j,k,iene)*unm1 &
                  + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))       &
                  + BET*(u(i+2,j,k,iene)*unp2-u(i-2,j,k,iene)*unm2 &
                  + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))       &
                  + GAM*(u(i+3,j,k,iene)*unp3-u(i-3,j,k,iene)*unm3 &
                  + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))       &
                  + DEL*(u(i+4,j,k,iene)*unp4-u(i-4,j,k,iene)*unm4 &
                  + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4)))*dxinv(1) 

             unp1 = q(i,j+1,k,qv)
             unp2 = q(i,j+2,k,qv)
             unp3 = q(i,j+3,k,qv)
             unp4 = q(i,j+4,k,qv)

             unm1 = q(i,j-1,k,qv)
             unm2 = q(i,j-2,k,qv)
             unm3 = q(i,j-3,k,qv)
             unm4 = q(i,j-4,k,qv)

             flux_irho=flux_irho - &
                   (ALP*(u(i,j+1,k,imy)-u(i,j-1,k,imy)) &
                  + BET*(u(i,j+2,k,imy)-u(i,j-2,k,imy)) &
                  + GAM*(u(i,j+3,k,imy)-u(i,j-3,k,imy)) &
                  + DEL*(u(i,j+4,k,imy)-u(i,j-4,k,imy)))*dxinv(2)

             flux_imx=flux_imx - &
                   (ALP*(u(i,j+1,k,imx)*unp1-u(i,j-1,k,imx)*unm1) &
                  + BET*(u(i,j+2,k,imx)*unp2-u(i,j-2,k,imx)*unm2) &
                  + GAM*(u(i,j+3,k,imx)*unp3-u(i,j-3,k,imx)*unm3) &
                  + DEL*(u(i,j+4,k,imx)*unp4-u(i,j-4,k,imx)*unm4))*dxinv(2)

             flux_imy=flux_imy - &
                   (ALP*(u(i,j+1,k,imy)*unp1-u(i,j-1,k,imy)*unm1 &
                  + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))               &
                  + BET*(u(i,j+2,k,imy)*unp2-u(i,j-2,k,imy)*unm2 &
                  + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))               &
                  + GAM*(u(i,j+3,k,imy)*unp3-u(i,j-3,k,imy)*unm3 &
                  + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))               &
                  + DEL*(u(i,j+4,k,imy)*unp4-u(i,j-4,k,imy)*unm4 &
                  + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres))))*dxinv(2)

             flux_imz=flux_imz - &
                   (ALP*(u(i,j+1,k,imz)*unp1-u(i,j-1,k,imz)*unm1) &
                  + BET*(u(i,j+2,k,imz)*unp2-u(i,j-2,k,imz)*unm2) &
                  + GAM*(u(i,j+3,k,imz)*unp3-u(i,j-3,k,imz)*unm3) &
                  + DEL*(u(i,j+4,k,imz)*unp4-u(i,j-4,k,imz)*unm4))*dxinv(2)

             flux_iene=flux_iene - &
                   (ALP*(u(i,j+1,k,iene)*unp1-u(i,j-1,k,iene)*unm1 &
                  + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))       &
                  + BET*(u(i,j+2,k,iene)*unp2-u(i,j-2,k,iene)*unm2 &
                  + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))       &
                  + GAM*(u(i,j+3,k,iene)*unp3-u(i,j-3,k,iene)*unm3 &
                  + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))       &
                  + DEL*(u(i,j+4,k,iene)*unp4-u(i,j-4,k,iene)*unm4 &
                  + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4)))*dxinv(2)

             unp1 = q(i,j,k+1,qw)
             unp2 = q(i,j,k+2,qw)
             unp3 = q(i,j,k+3,qw)
             unp4 = q(i,j,k+4,qw)

             unm1 = q(i,j,k-1,qw)
             unm2 = q(i,j,k-2,qw)
             unm3 = q(i,j,k-3,qw)
             unm4 = q(i,j,k-4,qw)

             flux_irho=flux_irho - &
                   (ALP*(u(i,j,k+1,imz)-u(i,j,k-1,imz)) &
                  + BET*(u(i,j,k+2,imz)-u(i,j,k-2,imz)) &
                  + GAM*(u(i,j,k+3,imz)-u(i,j,k-3,imz)) &
                  + DEL*(u(i,j,k+4,imz)-u(i,j,k-4,imz)))*dxinv(3)

             flux_imx=flux_imx - &
                   (ALP*(u(i,j,k+1,imx)*unp1-u(i,j,k-1,imx)*unm1) &
                  + BET*(u(i,j,k+2,imx)*unp2-u(i,j,k-2,imx)*unm2) &
                  + GAM*(u(i,j,k+3,imx)*unp3-u(i,j,k-3,imx)*unm3) &
                  + DEL*(u(i,j,k+4,imx)*unp4-u(i,j,k-4,imx)*unm4))*dxinv(3)

             flux_imy=flux_imy - &
                   (ALP*(u(i,j,k+1,imy)*unp1-u(i,j,k-1,imy)*unm1) &
                  + BET*(u(i,j,k+2,imy)*unp2-u(i,j,k-2,imy)*unm2) &
                  + GAM*(u(i,j,k+3,imy)*unp3-u(i,j,k-3,imy)*unm3) &
                  + DEL*(u(i,j,k+4,imy)*unp4-u(i,j,k-4,imy)*unm4))*dxinv(3)

             flux_imz=flux_imz - &
                   (ALP*(u(i,j,k+1,imz)*unp1-u(i,j,k-1,imz)*unm1 &
                  + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))               &
                  + BET*(u(i,j,k+2,imz)*unp2-u(i,j,k-2,imz)*unm2 &
                  + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))               &
                  + GAM*(u(i,j,k+3,imz)*unp3-u(i,j,k-3,imz)*unm3 &
                  + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))               &
                  + DEL*(u(i,j,k+4,imz)*unp4-u(i,j,k-4,imz)*unm4 &
                  + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres))))*dxinv(3)

             flux_iene=flux_iene - &
                   (ALP*(u(i,j,k+1,iene)*unp1-u(i,j,k-1,iene)*unm1 &
                  + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))       &
                  + BET*(u(i,j,k+2,iene)*unp2-u(i,j,k-2,iene)*unm2 &
                  + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))       &
                  + GAM*(u(i,j,k+3,iene)*unp3-u(i,j,k-3,iene)*unm3 &
                  + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))       &
                  + DEL*(u(i,j,k+4,iene)*unp4-u(i,j,k-4,iene)*unm4 &
                  + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4)))*dxinv(3)

             !
             ! Update U: If Unew is the same as U, this turns the update
             ! into a Gauss-Seidel type update where some newer values of U
             ! are used with older values of U when computing flux terms.
             !
             unew(i,j,k,irho) = a * uold(i,j,k,irho) + &
                  b * (u(i,j,k,irho) + dt * flux_irho)
             unew(i,j,k,imx) = a * uold(i,j,k,imx) + &
                  b * (u(i,j,k,imx) + dt * (difflux_imx + flux_imx))
             unew(i,j,k,imy) = a * uold(i,j,k,imy) + &
                  b * (u(i,j,k,imy) + dt * (difflux_imy + flux_imy))
             unew(i,j,k,imz) = a * uold(i,j,k,imz) + &
                  b * (u(i,j,k,imz) + dt * (difflux_imz + flux_imz))
             unew(i,j,k,iene) = a * uold(i,j,k,iene) + &
                  b * (u(i,j,k,iene) + dt * (difflux_iene + flux_iene))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call destroy(bpt_updateU_loop2)

    call destroy(bpt_updateU)

  end subroutine updateU

end module advance_module

