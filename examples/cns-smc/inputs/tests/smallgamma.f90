subroutine calc_gamma( k, nx, ny, nz, avmolwt, n_spec, yspecies, molwt_c)
    !---------------------------------------------                                                       
    ! compute gamma = cp/cv = cp/(cp-R)                                                                  
    !---------------------------------------------                                                         
    implicit none
    integer nx, ny, nz, n_spec, m,k
    
    real, intent(inout), dimension(k:nx,k:ny,k:nz):: avmolwt 
    real, intent(in) , dimension(k:nx,k:ny,k:nz, n_spec) :: yspecies
    real, intent(in), dimension(n_spec) :: molwt_c

    m = n_spec
    avmolwt(:,:,:) = yspecies(:,:,:,m) * molwt_c(m)

    do m = 1, n_spec-1
       avmolwt(:,:,:) = yspecies(:,:,:,m) * molwt_c(m)
    enddo

    return
  end subroutine calc_gamma
