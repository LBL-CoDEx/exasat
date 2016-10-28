subroutine div(lo,hi,ng,alpha,beta,gamma,u,x,y,z)

  integer,          intent(in ) :: lo(3),hi(3),ng
  double precision, intent(in ) :: alpha,beta,gamma
  integer                       :: i,j,k

  double precision, intent(out  ) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision, intent(in   ) :: x(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision, intent(in   ) :: y(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
  double precision, intent(in   ) :: z(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)

  !$OMP PARALLEL DO PRIVATE(i,j,k,eint,rhoinv)
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           u(i,j,k) = alpha * (x(i+1,j,k) - x(i-1,j,k)) + &
                      beta  * (y(i,j+1,k) - y(i,j-1,k)) + &
                      gamma * (z(i,j,k+1) - z(i,j,k-1))
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine div
