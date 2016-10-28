module mg_kernels_mod

    implicit none

contains

    subroutine gsrb_3d(color,u,f,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3), color
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bx(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: by(lo(1):, lo(2):, lo(3):)
        double precision, intent(in)  :: bz(lo(1):, lo(2):, lo(3):)
        double precision, intent(out) :: u(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)

        integer :: i, j, k, offset

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                offset = modulo(j+k+color, 2)
                do i = lo(1)+offset,hi(1),2
                    u(i,j,k) = ( bx(i+1,j,  k  )*u(i+1,j,  k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                               + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                               + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                               - dx**2*f(i,j,k)) / &
                               ( bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                               + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ) )
                end do
            end do
        end do

    end subroutine gsrb_3d

    subroutine restriction_3d(crse,fine,crse_lo,crse_hi,fine_lo,fine_hi)

        implicit none

        integer, intent(in) :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: fine(fine_lo(1):,fine_lo(2):,fine_lo(3):)
        double precision, intent(out) :: crse(crse_lo(1):,crse_lo(2):,crse_lo(3):)

        integer :: i,j,k

        do k = crse_lo(3), crse_hi(3)
            do j = crse_lo(2), crse_hi(2)
                do i = crse_lo(1), crse_hi(1)
                    crse(i,j,k) = (fine(2*i  ,2*j  ,2*k) + fine(2*i  ,2*j  ,2*k+1) &
                                 + fine(2*i+1,2*j  ,2*k) + fine(2*i+1,2*j  ,2*k+1) &
                                 + fine(2*i  ,2*j+1,2*k) + fine(2*i  ,2*j+1,2*k+1) &
                                 + fine(2*i+1,2*j+1,2*k) + fine(2*i+1,2*j+1,2*k+1)) / 8.d0
                end do
            end do
        end do

    end subroutine restriction_3d

    subroutine pc_prolongation_3d(fine,crse,ng,fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)

        implicit none

        integer, intent(in) :: ng, prolongation_type
        integer, intent(in) :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:,crse_lo(3)-ng:)
        double precision, intent(out) :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:,fine_lo(3)-ng:)

        integer :: i, j, k, ii, jj, kk

        do k = fine_lo(3),fine_hi(3)
            do j = fine_lo(2),fine_hi(2)
                do i = fine_lo(1),fine_hi(1)
                    ii = i/2
                    jj = j/2
                    kk = k/2
                    fine(i,j,k) = crse(ii,jj,kk) + fine(i,j,k)
                end do
            end do
        end do

    end subroutine pc_prolongation_3d

    subroutine trilinear_prolongation_3d(fine,crse,ng,fine_lo,fine_hi,crse_lo,crse_hi,prolongation_type)

        implicit none

        integer, intent(in) :: ng, prolongation_type
        integer, intent(in) :: crse_lo(3), crse_hi(3), fine_lo(3), fine_hi(3)
        double precision, intent(in)  :: crse(crse_lo(1)-ng:,crse_lo(2)-ng:,crse_lo(3)-ng:)
        double precision, intent(out) :: fine(fine_lo(1)-ng:,fine_lo(2)-ng:,fine_lo(3)-ng:)

        integer :: i, j, k, ii, jj, kk, si, sj, sk
        double precision, parameter :: c1 = 27.d0/64.d0
        double precision, parameter :: c2 = 9.d0/64.d0
        double precision, parameter :: c3 = 3.d0/64.d0
        double precision, parameter :: c4 = 1.d0/64.d0

        do k = fine_lo(3),fine_hi(3)
            do j = fine_lo(2),fine_hi(2)
                do i = fine_lo(1),fine_hi(1)
                    ii = i/2
                    jj = j/2
                    kk = k/2
                    si = 2*mod(i,2)-1 ! -1 if even, +1 if odd
                    sj = 2*mod(j,2)-1 ! -1 if even, +1 if odd
                    sk = 2*mod(k,2)-1 ! -1 if even, +1 if odd
                    fine(i,j,k) = c1*crse(ii,jj,kk) + c4*crse(ii+si,jj+sj,kk+sk) + fine(i,j,k) &
                                + c2*(crse(ii+si,jj   ,kk   ) + crse(ii   ,jj+sj,kk   ) + crse(ii   ,jj   ,kk+sk)) &
                                + c3*(crse(ii   ,jj+sj,kk+sk) + crse(ii+si,jj   ,kk+sk) + crse(ii+si,jj+sj,kk   ))
                end do
            end do
        end do

    end subroutine trilinear_prolongation_3d

    subroutine applyop_3d(g,u,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: bx(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: by(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: bz(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(out) :: g(lo(1):,lo(2):,lo(3):)

        integer :: i, j, k
        double precision :: lapl

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    g(i,j,k) = ( bx(i+1,j  ,k  )*u(i+1,j  ,k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                               + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                               + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                               -(bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                               + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ))*u(i,j,k) ) / dx**2
                end do
            end do
        end do

    end subroutine applyop_3d

    subroutine residual_3d(res,u,f,bx,by,bz,ng,lo,hi,dx)

        implicit none

        integer         , intent(in)  :: ng ! Will always be 1 for us
        integer         , intent(in)  :: lo(3), hi(3)
        double precision, intent(in)  :: dx
        double precision, intent(in)  :: f(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: bx(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: by(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: bz(lo(1):,lo(2):,lo(3):)
        double precision, intent(in)  :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
        double precision, intent(out) :: res(lo(1):,lo(2):,lo(3):)

        integer :: i, j, k
        double precision :: lapl

        do k = lo(3),hi(3)
            do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                    lapl = ( bx(i+1,j  ,k  )*u(i+1,j  ,k  ) + bx(i,j,k)*u(i-1,j  ,k  ) &
                           + by(i  ,j+1,k  )*u(i  ,j+1,k  ) + by(i,j,k)*u(i  ,j-1,k  ) &
                           + bz(i  ,j  ,k+1)*u(i  ,j  ,k+1) + bz(i,j,k)*u(i  ,j  ,k-1) &
                           -(bx(i+1,j,k) + by(i,j+1,k) + bz(i,j,k+1) &
                           + bx(i  ,j,k) + by(i,j  ,k) + bz(i,j,k  ))*u(i,j,k) ) / dx**2
                    res(i,j,k) = f(i,j,k) - lapl
                end do
            end do
        end do

    end subroutine residual_3d


end module mg_kernels_mod

