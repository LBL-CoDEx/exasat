subroutine calc_gamma( gamma, cpmix,k, mixMW, nx, ny, nz, pressure, rho, Ru, temp, avmolwt, n_spec, yspecies, molwt_c)
    !---------------------------------------------                                                       
    ! compute gamma = cp/cv = cp/(cp-R)                                                                  
    !---------------------------------------------                                                         
    implicit none
    integer nx, ny, nz, n_spec, m,k
    
    real, intent(in), dimension(2:nx,2:ny,2:nz) :: cpmix, mixMW, rho, Ru, temp
    real, intent(out), dimension(2:nx,2:ny,2:nz) :: gamma, pressure
    real, intent(inout), dimension(k:nx,k:ny,k:nz):: avmolwt 
    real, intent(in) , dimension(k:nx,k:ny,k:nz, n_spec) :: yspecies
    real, intent(in), dimension(n_spec) :: molwt_c

!!$    gamma = cpmix/(cpmix-Ru/mixMW)                                                                      
    gamma = cpmix/(cpmix-gamma)

!!$    pressure = rho*Ru*temp/mixMW                                                                                                                                                           
    pressure = rho*Ru*temp*avmolwt

    m = n_spec
    avmolwt(:,:,:) = yspecies(:,:,:,m) * molwt_c(m)

    do m = 1, n_spec-1
       avmolwt(:,:,:) = avmolwt(:,:,:) + yspecies(:,:,:,m) * molwt_c(m)
    enddo


    return
  end subroutine calc_gamma
