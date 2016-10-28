!21 species
subroutine chemterm_3d(lo,hi,ng,Q,Uprime, ncons, nprim, nspecies)

  integer, intent(in)::nspecies, ncons, nprim
! up is UPrime that has no ghost cells
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: Q (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: Uprime(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer, parameter :: qy1 = 8
    integer :: iwrk, i,j,k, n, iry1
    double precision :: Yt(nspecies), wdot(nspecies), molecular_weight(nspecies), rwrk

    double precision, save ::k_f_save(84)
    double precision, save ::Kc_save(84)

    double precision, parameter:: imw(21) = (/ 1.0 / 2.015940, &    !/*H2 !
                                1.0 / 1.007970, &    !/*H !
                                1.0 / 15.999400, &   !/*O !
                                1.0 / 31.998800, &   !/*O2 !
                                1.0 / 17.007370, &   !/*OH !
                                1.0 / 18.015340, &   !/*H2O !
                                1.0 / 33.006770, &   !/*HO2 !
                                1.0 / 14.027090, &   !/*CH2 !
                                1.0 / 14.027090, &   !/*CH2(S) !
                                1.0 / 15.035060, &   !/*CH3 !
                                1.0 / 16.043030, &   !/*CH4 !
                                1.0 / 28.010550, &   !/*CO !
                                1.0 / 44.009950, &   !/*CO2 !
                                1.0 / 29.018520, &   !/*HCO !
                                1.0 / 30.026490, &   !/*CH2O !
                                1.0 / 31.034460, &   !/*CH3O !
                                1.0 / 28.054180, &   !/*C2H4 !
                                1.0 / 29.062150, &   !/*C2H5 !
                                1.0 / 30.070120, &   !/*C2H6 !
                                1.0 / 28.013400, &   !/*N2 !
                                1.0 / 39.948000 /)  !/*AR !

    double precision:: sc(21)  !/*temporary storage !
    double precision:: qdot 

    double precision:: mixture                  !/*mixture concentration !
    double precision:: g_RT(21)                 !/*Gibbs free energy !
    double precision:: Kc                       !/*equilibrium constant !
    double precision:: k_f                      !/*forward reaction rate !
    double precision:: k_r                      !/*reverse reaction rate !
    double precision:: q_f                      !/*forward progress rate !
    double precision:: q_r                      !/*reverse progress rate !
    double precision:: phi_f                    !/*forward phase space factor !
    double precision:: phi_r                    !/*reverse phase space factor !
    double precision:: alpha                    !/*enhancement !
    double precision:: redP                     !/*reduced pressure !
    double precision:: logPred                  !/*log of above !
    double precision:: F                        !/*fallof rate enhancement !

    double precision:: F_troe                   !/*TROE integer::ermediate !
    double precision:: logFcent                 !/*TROE integer::ermediate !
    double precision:: troe                     !/*TROE integer::ermediate !
    double precision:: troe_c                   !/*TROE integer::ermediate !
    double precision:: troe_n                   !/*TROE integer::ermediate !

    double precision:: tsc(5) != ( log(T), T, T*T, T*T*T, T*T*T*T )  !/*temperature cache !

    double precision:: invT != 1.0 / tsc(2) 

    !/*reference concentration: P_atm / (RT) in inverse mol/m^3 !
    double precision:: refC != 101325 / 8.31451 / T 
   !/*temperature !
    double precision:: T != tsc(2), invT = 1.0 / T 
    integer:: qyn, qrho, qtemp

    double precision:: rho, T_save

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             do n = 1, nspecies
                qyn = n + qy1 - 1
                Yt(n) = Q(i,j,k,qyn)
             end do 

             rho = Q(i,j,k,qrho)
             T =   Q(i,j,k,qtemp)

             !/*See Eq 8 with an extra 1e6 so c goes to SI !
             do n =1, nspecies
                sc(n) = 1e6 * rho * Yt(n)* imw(n)
             end do

             tsc(1) = log(T)
             tsc(2) = T
             tsc(3) = T*T
             tsc(4) = T*T*T
             tsc(5) = T*T*T*T   !/*temperature cache !
             
             invT = 1.0 / T
             refC = 101325 / 8.31451 / T

    
             !/*compute the mixture concentration !
             mixture = 0.0 
             do n = 1,nspecies
                mixture = mixture +  sc(n) 
             end do

             !/*species 0: H2 !
             do n= 1, nspecies
                g_RT(n) = &
                     -9.179351730000000e+02 * invT &
                     +1.661320882000000e+00 &
                     -2.344331120000000e+00 * tsc(1) &
                     -3.990260375000000e-03 * tsc(2) &
                     +3.246358500000000e-06 * tsc(3) &
                     -1.679767450000000e-09 * tsc(4) &
                     +3.688058805000000e-13 * tsc(5) 
             end do
             
             
             !/*zero out wdot !
             do n = 1, nspecies
                wdot(n) = 0.0 
             end do
             
    if (T .ne. T_save) then
        T_save = T 

        k_f_save(1) = 1e-12 * 5e+17*exp(-1*tsc(1)) 
        k_f_save(2) = 1e-06 * 50000*exp(2.67*tsc(1)-3165.2328279116868543*invT) 
        k_f_save(3) = 1e-06 * 2e+13 
        k_f_save(4) = 1e-06 * 8e+13 
        k_f_save(5) = 1e-06 * 1.5e+13 
        k_f_save(6) = 1e-06 * 8.43e+13 
        k_f_save(7) = 1e-06 * 1.02e+09*exp(1.5*tsc(1)-4327.6633259205891591*invT) 
        k_f_save(8) = 1e-12 * 6.02e+14*exp(-1509.6499974141590883*invT) 
        k_f_save(9) = 1e-06 * 3e+13 
        k_f_save(10) = 1e-06 * 3e+13 
        k_f_save(11) = 1e-06 * 3.9e+13*exp(-1781.3869969487077469*invT) 
        k_f_save(12) = 1e-06 * 1.92e+07*exp(1.83*tsc(1)-110.7076664770383303*invT) 
        k_f_save(13) = 1e-06 * 1.32e+14 
        k_f_save(14) = 1e-06 * 8.98e+07*exp(1.92*tsc(1)-2863.3028284288548093*invT) 
        k_f_save(15) = 1e-06 * 2.5e+12*exp(-24053.756625465601246*invT) 
        k_f_save(16) = 1e-06 * 1e+14*exp(-20128.666632188789663*invT) 
        k_f_save(17) = 1e-12 * 2.8e+18*exp(-0.86*tsc(1)) 
        k_f_save(18) = 1e-12 * 3e+20*exp(-1.72*tsc(1)) 
        k_f_save(19) = 1e-12 * 9.38e+18*exp(-0.76*tsc(1)) 
        k_f_save(20) = 1e-12 * 3.75e+20*exp(-1.72*tsc(1)) 
        k_f_save(21) = 1e-12 * 7e+17*exp(-0.8*tsc(1)) 
        k_f_save(22) = 1e-06 * 8.3e+13*exp(-7252.8618042434254676*invT) 
        k_f_save(23) = 1e-12 * 1e+18*exp(-1*tsc(1)) 
        k_f_save(24) = 1e-12 * 9e+16*exp(-0.6*tsc(1)) 
        k_f_save(25) = 1e-12 * 6e+19*exp(-1.25*tsc(1)) 
        k_f_save(26) = 1e-12 * 5.5e+20*exp(-2*tsc(1)) 
        k_f_save(27) = 1e-12 * 2.2e+22*exp(-2*tsc(1)) 
        k_f_save(28) = 1e-06 * 2.8e+13*exp(-537.43539907944057177*invT) 
        k_f_save(29) = 1e-06 * 1.34e+14*exp(-319.54258278599701271*invT) 
        k_f_save(30) = 1e-06 * 2.5e+16*exp(-0.8*tsc(1)) 
        k_f_save(31) = 1e-06 * 1.27e+16*exp(-0.63*tsc(1)-192.73198300320765952*invT) 
        k_f_save(32) = 1e-06 * 6.6e+08*exp(1.62*tsc(1)-5454.8686573231616421*invT) 
        k_f_save(33) = 1e-06 * 1.09e+12*exp(0.48*tsc(1)+130.83633310922709825*invT) 
        k_f_save(34) = 1e-06 * 7.34e+13 
        k_f_save(35) = 1e-06 * 5.4e+11*exp(0.454*tsc(1)-1308.3633310922712099*invT) 
        k_f_save(36) = 1e-06 * 2.3e+10*exp(1.05*tsc(1)-1648.0345805104568626*invT) 
        k_f_save(37) = 1e-06 * 3.2e+13 
        k_f_save(38) = 1e-06 * 1.08e+12*exp(0.454*tsc(1)-915.85433176458980142*invT) 
        k_f_save(39) = 1e-06 * 5.21e+17*exp(-0.99*tsc(1)-795.08233197145705162*invT) 
        k_f_save(40) = 1e-06 * 1.15e+08*exp(1.9*tsc(1)-3789.2214935095394139*invT) 
        k_f_save(41) = 1e-06 * 4.3e+07*exp(1.5*tsc(1)-40056.046598055690993*invT) 
        k_f_save(42) = 1e-06 * 2.16e+08*exp(1.51*tsc(1)-1726.0331637101885462*invT) 
        k_f_save(43) = 1e-06 * 35700*exp(2.4*tsc(1)+1061.7871648479585929*invT) 
        k_f_save(44) = 1e-06 * 2.9e+13*exp(+251.60833290235981963*invT) 
        k_f_save(45) = 1e-06 * 2e+13 
        k_f_save(46) = 1e-06 * 3e+13 
        k_f_save(47) = 1e-06 * 5.6e+07*exp(1.6*tsc(1)-2727.4343286615808211*invT) 
        k_f_save(48) = 1e-06 * 2.501e+13 
        k_f_save(49) = 1e-06 * 1e+08*exp(1.6*tsc(1)-1570.0359973107254064*invT) 
        k_f_save(50) = 1e-06 * 4.76e+07*exp(1.228*tsc(1)-35.225166606330375885*invT) 
        k_f_save(51) = 1e-06 * 5e+13 
        k_f_save(52) = 1e-06 * 3.43e+09*exp(1.18*tsc(1)+224.93784961470970529*invT) 
        k_f_save(53) = 1e-06 * 3.54e+06*exp(2.12*tsc(1)-437.79849925010609013*invT) 
        k_f_save(54) = 1e-06 * 2e+13 
        k_f_save(55) = 1e-06 * 1e+12 
        k_f_save(56) = 1e-06 * 2e+13 
        k_f_save(57) = 1e-06 * 1.5e+14*exp(-11875.913312991384373*invT) 
        k_f_save(58) = 1e-06 * 1.32e+13*exp(-754.82499870707954415*invT) 
        k_f_save(59) = 1e-06 * 500000*exp(2*tsc(1)-3638.256493768123164*invT) 
        k_f_save(60) = 1e-06 * 4e+13 
        k_f_save(61) = 1e-06 * 2.46e+06*exp(2*tsc(1)-4161.6018262050320118*invT) 
        k_f_save(62) = 1e-06 * 1.5e+13*exp(-301.92999948283181766*invT) 
        k_f_save(63) = 1e-06 * 9e+12*exp(-301.92999948283181766*invT) 
        k_f_save(64) = 1e-06 * 2.8e+13 
        k_f_save(65) = 1e-06 * 1.2e+13 
        k_f_save(66) = 1e-06 * 7e+13 
        k_f_save(67) = 1e-06 * 3e+13 
        k_f_save(68) = 1e-06 * 1.2e+13*exp(+286.83349950869023814*invT) 
        k_f_save(69) = 1e-06 * 1.6e+13*exp(+286.83349950869023814*invT) 
        k_f_save(70) = 1e-06 * 9e+12 
        k_f_save(71) = 1e-06 * 7e+12 
        k_f_save(72) = 1e-06 * 1.4e+13 
        k_f_save(73) = 1e-06 * 2.675e+13*exp(-14492.639975175927248*invT) 
        k_f_save(74) = 1e-06 * 3.6e+10*exp(-4498.7569922941938785*invT) 
        k_f_save(75) = 1e-06 * 2.12e+16*exp(-0.97*tsc(1)-311.99433279892622295*invT) 
        k_f_save(76) = 1e-06 * 4.99e+12*exp(0.1*tsc(1)-5334.096657530029006*invT) 
        k_f_save(77) = 1e-06 * 2.648e+13 
        k_f_save(78) = 1e-06 * 3320*exp(2.81*tsc(1)-2948.8496616156576238*invT) 
        k_f_save(79) = 1e-06 * 6.14e+06*exp(1.74*tsc(1)-5258.6141576593208811*invT) 
        k_f_save(80) = 1e-06 * 2.244e+18*exp(-1*tsc(1)-8554.6833186802341515*invT) 
        k_f_save(81) = 1e-06 * 1.87e+17*exp(-1*tsc(1)-8554.6833186802341515*invT) 
        k_f_save(82) = 1e-06 * 7.6e+12*exp(-201.28666632188787844*invT) 
        k_f_save(83) = 1e-06 * 4.28e-13*exp(7.6*tsc(1)+1776.3548302906606295*invT) 
        k_f_save(84) = 1e-06 * 8.4e+11*exp(-1949.9645799932889076*invT) 

        Kc_save(1) = 1.0 / (refC) * exp((g_RT(3) + g_RT(2)) - (g_RT(5))) 
        Kc_save(2) = exp((g_RT(3) + g_RT(1)) - (g_RT(2) + g_RT(5))) 
        Kc_save(3) = exp((g_RT(3) + g_RT(7)) - (g_RT(5) + g_RT(4))) 
        Kc_save(4) = exp((g_RT(3) + g_RT(8)) - (g_RT(2) + g_RT(14))) 
        Kc_save(5) = exp((g_RT(3) + g_RT(9)) - (g_RT(2) + g_RT(14))) 
        Kc_save(6) = exp((g_RT(3) + g_RT(10)) - (g_RT(2) + g_RT(15))) 
        Kc_save(7) = exp((g_RT(3) + g_RT(11)) - (g_RT(5) + g_RT(10))) 
        Kc_save(8) = 1.0 / (refC) * exp((g_RT(3) + g_RT(12)) - (g_RT(13))) 
        Kc_save(9) = exp((g_RT(3) + g_RT(14)) - (g_RT(5) + g_RT(12))) 
        Kc_save(10) = exp((g_RT(3) + g_RT(14)) - (g_RT(2) + g_RT(13))) 
        Kc_save(11) = exp((g_RT(3) + g_RT(15)) - (g_RT(5) + g_RT(14))) 
        Kc_save(12) = exp((g_RT(3) + g_RT(17)) - (g_RT(10) + g_RT(14))) 
        Kc_save(13) = exp((g_RT(3) + g_RT(18)) - (g_RT(10) + g_RT(15))) 
        Kc_save(14) = exp((g_RT(3) + g_RT(19)) - (g_RT(5) + g_RT(18))) 
        Kc_save(15) = exp((g_RT(4) + g_RT(12)) - (g_RT(3) + g_RT(13))) 
        Kc_save(16) = exp((g_RT(4) + g_RT(15)) - (g_RT(7) + g_RT(14))) 
        Kc_save(17) = 1.0 / (refC) * exp((g_RT(2) + g_RT(4)) - (g_RT(7))) 
        Kc_save(18) = 1.0 / (refC) * exp((g_RT(2) + 2 * g_RT(4)) - (g_RT(7) + g_RT(4))) 
        Kc_save(19) = 1.0 / (refC) * exp((g_RT(2) + g_RT(4) + g_RT(6)) - (g_RT(7) + g_RT(6))) 
        Kc_save(20) = 1.0 / (refC) * exp((g_RT(2) + g_RT(4) + g_RT(20)) - (g_RT(7) + g_RT(20))) 
        Kc_save(21) = 1.0 / (refC) * exp((g_RT(2) + g_RT(4) + g_RT(21)) - (g_RT(7) + g_RT(21))) 
        Kc_save(22) = exp((g_RT(2) + g_RT(4)) - (g_RT(3) + g_RT(5))) 
        Kc_save(23) = 1.0 / (refC) * exp((2 * g_RT(2)) - (g_RT(1))) 
        Kc_save(24) = 1.0 / (refC) * exp((2 * g_RT(2) + g_RT(1)) - (2 * g_RT(1))) 
        Kc_save(25) = 1.0 / (refC) * exp((2 * g_RT(2) + g_RT(6)) - (g_RT(1) + g_RT(6))) 
        Kc_save(26) = 1.0 / (refC) * exp((2 * g_RT(2) + g_RT(13)) - (g_RT(1) + g_RT(13))) 
        Kc_save(27) = 1.0 / (refC) * exp((g_RT(2) + g_RT(5)) - (g_RT(6))) 
        Kc_save(28) = exp((g_RT(2) + g_RT(7)) - (g_RT(4) + g_RT(1))) 
        Kc_save(29) = exp((g_RT(2) + g_RT(7)) - (2 * g_RT(5))) 
        Kc_save(30) = 1.0 / (refC) * exp((g_RT(2) + g_RT(8)) - (g_RT(10))) 
        Kc_save(31) = 1.0 / (refC) * exp((g_RT(2) + g_RT(10)) - (g_RT(11))) 
        Kc_save(32) = exp((g_RT(2) + g_RT(11)) - (g_RT(10) + g_RT(1))) 
        Kc_save(33) = 1.0 / (refC) * exp((g_RT(2) + g_RT(14)) - (g_RT(15))) 
        Kc_save(34) = exp((g_RT(2) + g_RT(14)) - (g_RT(1) + g_RT(12))) 
        Kc_save(35) = 1.0 / (refC) * exp((g_RT(2) + g_RT(15)) - (g_RT(16))) 
        Kc_save(36) = exp((g_RT(2) + g_RT(15)) - (g_RT(14) + g_RT(1))) 
        Kc_save(37) = exp((g_RT(2) + g_RT(16)) - (g_RT(5) + g_RT(10))) 
        Kc_save(38) = 1.0 / (refC) * exp((g_RT(2) + g_RT(17)) - (g_RT(18))) 
        Kc_save(39) = 1.0 / (refC) * exp((g_RT(2) + g_RT(18)) - (g_RT(19))) 
        Kc_save(40) = exp((g_RT(2) + g_RT(19)) - (g_RT(18) + g_RT(1))) 
        Kc_save(41) = 1.0 / (refC) * exp((g_RT(1) + g_RT(12)) - (g_RT(15))) 
        Kc_save(42) = exp((g_RT(5) + g_RT(1)) - (g_RT(2) + g_RT(6))) 
        Kc_save(43) = exp((2 * g_RT(5)) - (g_RT(3) + g_RT(6))) 
        Kc_save(44) = exp((g_RT(5) + g_RT(7)) - (g_RT(4) + g_RT(6))) 
        Kc_save(45) = exp((g_RT(5) + g_RT(8)) - (g_RT(2) + g_RT(15))) 
        Kc_save(46) = exp((g_RT(5) + g_RT(9)) - (g_RT(2) + g_RT(15))) 
        Kc_save(47) = exp((g_RT(5) + g_RT(10)) - (g_RT(8) + g_RT(6))) 
        Kc_save(48) = exp((g_RT(5) + g_RT(10)) - (g_RT(9) + g_RT(6))) 
        Kc_save(49) = exp((g_RT(5) + g_RT(11)) - (g_RT(10) + g_RT(6))) 
        Kc_save(50) = exp((g_RT(5) + g_RT(12)) - (g_RT(2) + g_RT(13))) 
        Kc_save(51) = exp((g_RT(5) + g_RT(14)) - (g_RT(6) + g_RT(12))) 
        Kc_save(52) = exp((g_RT(5) + g_RT(15)) - (g_RT(14) + g_RT(6))) 
        Kc_save(53) = exp((g_RT(5) + g_RT(19)) - (g_RT(18) + g_RT(6))) 
        Kc_save(54) = exp((g_RT(7) + g_RT(8)) - (g_RT(5) + g_RT(15))) 
        Kc_save(55) = exp((g_RT(7) + g_RT(10)) - (g_RT(4) + g_RT(11))) 
        Kc_save(56) = exp((g_RT(7) + g_RT(10)) - (g_RT(5) + g_RT(16))) 
        Kc_save(57) = exp((g_RT(7) + g_RT(12)) - (g_RT(5) + g_RT(13))) 
        Kc_save(58) = exp((g_RT(8) + g_RT(4)) - (g_RT(5) + g_RT(14))) 
        Kc_save(59) = exp((g_RT(8) + g_RT(1)) - (g_RT(2) + g_RT(10))) 
        Kc_save(60) = exp((g_RT(8) + g_RT(10)) - (g_RT(2) + g_RT(17))) 
        Kc_save(61) = exp((g_RT(8) + g_RT(11)) - (2 * g_RT(10))) 
        Kc_save(62) = exp((g_RT(9) + g_RT(20)) - (g_RT(8) + g_RT(20))) 
        Kc_save(63) = exp((g_RT(9) + g_RT(21)) - (g_RT(8) + g_RT(21))) 
        Kc_save(64) = refC * exp((g_RT(9) + g_RT(4)) - (g_RT(2) + g_RT(5) + g_RT(12))) 
        Kc_save(65) = exp((g_RT(9) + g_RT(4)) - (g_RT(12) + g_RT(6))) 
        Kc_save(66) = exp((g_RT(9) + g_RT(1)) - (g_RT(10) + g_RT(2))) 
        Kc_save(67) = exp((g_RT(9) + g_RT(6)) - (g_RT(8) + g_RT(6))) 
        Kc_save(68) = exp((g_RT(9) + g_RT(10)) - (g_RT(2) + g_RT(17))) 
        Kc_save(69) = exp((g_RT(9) + g_RT(11)) - (2 * g_RT(10))) 
        Kc_save(70) = exp((g_RT(9) + g_RT(12)) - (g_RT(8) + g_RT(12))) 
        Kc_save(71) = exp((g_RT(9) + g_RT(13)) - (g_RT(8) + g_RT(13))) 
        Kc_save(72) = exp((g_RT(9) + g_RT(13)) - (g_RT(12) + g_RT(15))) 
        Kc_save(73) = exp((g_RT(10) + g_RT(4)) - (g_RT(3) + g_RT(16))) 
        Kc_save(74) = exp((g_RT(10) + g_RT(4)) - (g_RT(5) + g_RT(15))) 
        Kc_save(75) = 1.0 / (refC) * exp((2 * g_RT(10)) - (g_RT(19))) 
        Kc_save(76) = exp((2 * g_RT(10)) - (g_RT(2) + g_RT(18))) 
        Kc_save(77) = exp((g_RT(10) + g_RT(14)) - (g_RT(11) + g_RT(12))) 
        Kc_save(78) = exp((g_RT(10) + g_RT(15)) - (g_RT(14) + g_RT(11))) 
        Kc_save(79) = exp((g_RT(10) + g_RT(19)) - (g_RT(18) + g_RT(11))) 
        Kc_save(80) = refC * exp((g_RT(14) + g_RT(6)) - (g_RT(2) + g_RT(12) + g_RT(6))) 
        Kc_save(81) = refC * exp((g_RT(14)) - (g_RT(2) + g_RT(12))) 
        Kc_save(82) = exp((g_RT(14) + g_RT(4)) - (g_RT(7) + g_RT(12))) 
        Kc_save(83) = exp((g_RT(16) + g_RT(4)) - (g_RT(7) + g_RT(15))) 
        Kc_save(84) = exp((g_RT(18) + g_RT(4)) - (g_RT(7) + g_RT(17))) 
    end if


    !/*reaction 1: O + H + M <=> OH + M !
    phi_f = sc(3)*sc(2) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = alpha * k_f_save(1) 
    q_f = phi_f * k_f 
    phi_r = sc(5) 
    Kc = Kc_save(1) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 

    !/*reaction 2: O + H2 <=> H + OH !
    phi_f = sc(3)*sc(1) 
    k_f = k_f_save(2) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(5) 
    Kc = Kc_save(2) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 

    !/*reaction 3: O + HO2 <=> OH + O2 !
    phi_f = sc(3)*sc(7) 
    k_f = k_f_save(3) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(4) 
    Kc = Kc_save(3) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(4) =    wdot(4) + 1 * qdot 

    !/*reaction 4: O + CH2 <=> H + HCO !
    phi_f = sc(3)*sc(8) 
    k_f = k_f_save(4) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(14) 
    Kc = Kc_save(4) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 5: O + CH2(S) <=> H + HCO !
    phi_f = sc(3)*sc(9) 
    k_f = k_f_save(5) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(14) 
    Kc = Kc_save(5) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 6: O + CH3 <=> H + CH2O !
    phi_f = sc(3)*sc(10) 
    k_f = k_f_save(6) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(15) 
    Kc = Kc_save(6) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 7: O + CH4 <=> OH + CH3 !
    phi_f = sc(3)*sc(11) 
    k_f = k_f_save(7) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(10) 
    Kc = Kc_save(7) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(11) =    wdot(11) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 

    !/*reaction 8: O + CO + M <=> CO2 + M !
    phi_f = sc(3)*sc(12) 
    alpha = mixture + sc(1) + 5*sc(4) + 5*sc(6) + sc(11) + 0.5*sc(12) + 2.5*sc(13) + 2*sc(19) -0.5*sc(21) 
    k_f = alpha * k_f_save(8) 
    q_f = phi_f * k_f 
    phi_r = sc(13) 
    Kc = Kc_save(8) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 9: O + HCO <=> OH + CO !
    phi_f = sc(3)*sc(14) 
    k_f = k_f_save(9) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(12) 
    Kc = Kc_save(9) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 10: O + HCO <=> H + CO2 !
    phi_f = sc(3)*sc(14) 
    k_f = k_f_save(10) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(13) 
    Kc = Kc_save(10) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 11: O + CH2O <=> OH + HCO !
    phi_f = sc(3)*sc(15) 
    k_f = k_f_save(11) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(14) 
    Kc = Kc_save(11) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 12: O + C2H4 <=> CH3 + HCO !
    phi_f = sc(3)*sc(17) 
    k_f = k_f_save(12) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(14) 
    Kc = Kc_save(12) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(17) =    wdot(17) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 13: O + C2H5 <=> CH3 + CH2O !
    phi_f = sc(3)*sc(18) 
    k_f = k_f_save(13) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(15) 
    Kc = Kc_save(13) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(18) =    wdot(18) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 14: O + C2H6 <=> OH + C2H5 !
    phi_f = sc(3)*sc(19) 
    k_f = k_f_save(14) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(18) 
    Kc = Kc_save(14) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(3) =    wdot(3) - 1 * qdot 
    wdot(19) =    wdot(19) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 

    !/*reaction 15: O2 + CO <=> O + CO2 !
    phi_f = sc(4)*sc(12) 
    k_f = k_f_save(15) 
    q_f = phi_f * k_f 
    phi_r = sc(3)*sc(13) 
    Kc = Kc_save(15) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 16: O2 + CH2O <=> HO2 + HCO !
    phi_f = sc(4)*sc(15) 
    k_f = k_f_save(16) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(14) 
    Kc = Kc_save(16) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 17: H + O2 + M <=> HO2 + M !
    phi_f = sc(2)*sc(4) 
    alpha = mixture -1*sc(4) -1*sc(6) -0.25*sc(12) + 0.5*sc(13) + 0.5*sc(19) -1*sc(20) -1*sc(21) 
    k_f = alpha * k_f_save(17) 
    q_f = phi_f * k_f 
    phi_r = sc(7) 
    Kc = Kc_save(17) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 

    !/*reaction 18: H + 2 O2 <=> HO2 + O2 !
    phi_f = sc(2)*sc(4)*sc(4) 
    k_f = k_f_save(18) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(4) 
    Kc = Kc_save(18) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 2 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(4) =    wdot(4) + 1 * qdot 

    !/*reaction 19: H + O2 + H2O <=> HO2 + H2O !
    phi_f = sc(2)*sc(4)*sc(6) 
    k_f = k_f_save(19) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(6) 
    Kc = Kc_save(19) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 20: H + O2 + N2 <=> HO2 + N2 !
    phi_f = sc(2)*sc(4)*sc(20) 
    k_f = k_f_save(20) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(20) 
    Kc = Kc_save(20) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(20) =    wdot(20) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(20) =    wdot(20) + 1 * qdot 

    !/*reaction 21: H + O2 + AR <=> HO2 + AR !
    phi_f = sc(2)*sc(4)*sc(21) 
    k_f = k_f_save(21) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(21) 
    Kc = Kc_save(21) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(21) =    wdot(21) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(21) =    wdot(21) + 1 * qdot 

    !/*reaction 22: H + O2 <=> O + OH !
    phi_f = sc(2)*sc(4) 
    k_f = k_f_save(22) 
    q_f = phi_f * k_f 
    phi_r = sc(3)*sc(5) 
    Kc = Kc_save(22) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 

    !/*reaction 23: 2 H + M <=> H2 + M !
    phi_f = sc(2)*sc(2) 
    alpha = mixture -1*sc(1) -1*sc(6) + sc(11) -1*sc(13) + 2*sc(19) -0.37*sc(21) 
    k_f = alpha * k_f_save(23) 
    q_f = phi_f * k_f 
    phi_r = sc(1) 
    Kc = Kc_save(23) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 2 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 24: 2 H + H2 <=> 2 H2 !
    phi_f = sc(2)*sc(2)*sc(1) 
    k_f = k_f_save(24) 
    q_f = phi_f * k_f 
    phi_r = sc(1)*sc(1) 
    Kc = Kc_save(24) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 2 * qdot 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(1) =    wdot(1) + 2 * qdot 

    !/*reaction 25: 2 H + H2O <=> H2 + H2O !
    phi_f = sc(2)*sc(2)*sc(6) 
    k_f = k_f_save(25) 
    q_f = phi_f * k_f 
    phi_r = sc(1)*sc(6) 
    Kc = Kc_save(25) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 2 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 26: 2 H + CO2 <=> H2 + CO2 !
    phi_f = sc(2)*sc(2)*sc(13) 
    k_f = k_f_save(26) 
    q_f = phi_f * k_f 
    phi_r = sc(1)*sc(13) 
    Kc = Kc_save(26) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 2 * qdot 
    wdot(13) =    wdot(13) - 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 27: H + OH + M <=> H2O + M !
    phi_f = sc(2)*sc(5) 
    alpha = mixture -0.27*sc(1) + 2.65*sc(6) + sc(11) + 2*sc(19) -0.62*sc(21) 
    k_f = alpha * k_f_save(27) 
    q_f = phi_f * k_f 
    phi_r = sc(6) 
    Kc = Kc_save(27) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 28: H + HO2 <=> O2 + H2 !
    phi_f = sc(2)*sc(7) 
    k_f = k_f_save(28) 
    q_f = phi_f * k_f 
    phi_r = sc(4)*sc(1) 
    Kc = Kc_save(28) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(4) =    wdot(4) + 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 29: H + HO2 <=> 2 OH !
    phi_f = sc(2)*sc(7) 
    k_f = k_f_save(29) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(5) 
    Kc = Kc_save(29) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(5) =    wdot(5) + 2 * qdot 

    !/*reaction 30: H + CH2 (+M) <=> CH3 (+M) !
    phi_f = sc(2)*sc(8) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(30) 
    redP = 1e-12 * alpha / k_f * 3.2e+27*exp(-3.14*tsc(1)-618.95649893980521483*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.32*exp(T/(-78)))+ (0.68*exp(T/(-1995)))+ (exp(-5590/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(10) 
    Kc = Kc_save(30) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 

    !/*reaction 31: H + CH3 (+M) <=> CH4 (+M) !
    phi_f = sc(2)*sc(10) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(31) 
    redP = 1e-12 * alpha / k_f * 2.477e+32*exp(-4.76*tsc(1)-1227.8486645635159675*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.217*exp(T/(-74)))+ (0.783*exp(T/(-2941)))+ (exp(-6964/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(11) 
    Kc = Kc_save(31) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(11) =    wdot(11) + 1 * qdot 

    !/*reaction 32: H + CH4 <=> CH3 + H2 !
    phi_f = sc(2)*sc(11) 
    k_f = k_f_save(32) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(1) 
    Kc = Kc_save(32) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(11) =    wdot(11) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 33: H + HCO (+M) <=> CH2O (+M) !
    phi_f = sc(2)*sc(14) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(33) 
    redP = 1e-12 * alpha / k_f * 1.35e+24*exp(-2.57*tsc(1)-717.08374877172548167*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.2176*exp(T/(-271)))+ (0.7824*exp(T/(-2755)))+ (exp(-6570/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(15) 
    Kc = Kc_save(33) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 34: H + HCO <=> H2 + CO !
    phi_f = sc(2)*sc(14) 
    k_f = k_f_save(34) 
    q_f = phi_f * k_f 
    phi_r = sc(1)*sc(12) 
    Kc = Kc_save(34) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 35: H + CH2O (+M) <=> CH3O (+M) !
    phi_f = sc(2)*sc(15) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) 
    k_f = k_f_save(35) 
    redP = 1e-12 * alpha / k_f * 2.2e+30*exp(-4.8*tsc(1)-2797.8846618742413739*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.242*exp(T/(-94)))+ (0.758*exp(T/(-1555)))+ (exp(-4200/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(16) 
    Kc = Kc_save(35) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(16) =    wdot(16) + 1 * qdot 

    !/*reaction 36: H + CH2O <=> HCO + H2 !
    phi_f = sc(2)*sc(15) 
    k_f = k_f_save(36) 
    q_f = phi_f * k_f 
    phi_r = sc(14)*sc(1) 
    Kc = Kc_save(36) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 37: H + CH3O <=> OH + CH3 !
    phi_f = sc(2)*sc(16) 
    k_f = k_f_save(37) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(10) 
    Kc = Kc_save(37) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(16) =    wdot(16) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 

    !/*reaction 38: H + C2H4 (+M) <=> C2H5 (+M) !
    phi_f = sc(2)*sc(17) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(38) 
    redP = 1e-12 * alpha / k_f * 1.2e+32*exp(-7.62*tsc(1)-3507.4201606588962932*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.0247*exp(T/(-210)))+ (0.9753*exp(T/(-984)))+ (exp(-4374/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(18) 
    Kc = Kc_save(38) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(17) =    wdot(17) - 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 

    !/*reaction 39: H + C2H5 (+M) <=> C2H6 (+M) !
    phi_f = sc(2)*sc(18) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(39) 
    redP = 1e-12 * alpha / k_f * 1.99e+32*exp(-7.08*tsc(1)-3364.0034109045514015*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.1578*exp(T/(-125)))+ (0.8422*exp(T/(-2219)))+ (exp(-6882/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(19) 
    Kc = Kc_save(39) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(18) =    wdot(18) - 1 * qdot 
    wdot(19) =    wdot(19) + 1 * qdot 

    !/*reaction 40: H + C2H6 <=> C2H5 + H2 !
    phi_f = sc(2)*sc(19) 
    k_f = k_f_save(40) 
    q_f = phi_f * k_f 
    phi_r = sc(18)*sc(1) 
    Kc = Kc_save(40) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(2) =    wdot(2) - 1 * qdot 
    wdot(19) =    wdot(19) - 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 41: H2 + CO (+M) <=> CH2O (+M) !
    phi_f = sc(1)*sc(12) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(41) 
    redP = 1e-12 * alpha / k_f * 5.07e+27*exp(-3.42*tsc(1)-42446.325760628104035*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.068*exp(T/(-197)))+ (0.932*exp(T/(-1540)))+ (exp(-10300/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(15) 
    Kc = Kc_save(41) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 42: OH + H2 <=> H + H2O !
    phi_f = sc(5)*sc(1) 
    k_f = k_f_save(42) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(6) 
    Kc = Kc_save(42) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 43: 2 OH <=> O + H2O !
    phi_f = sc(5)*sc(5) 
    k_f = k_f_save(43) 
    q_f = phi_f * k_f 
    phi_r = sc(3)*sc(6) 
    Kc = Kc_save(43) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 2 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 44: OH + HO2 <=> O2 + H2O !
    phi_f = sc(5)*sc(7) 
    k_f = k_f_save(44) 
    q_f = phi_f * k_f 
    phi_r = sc(4)*sc(6) 
    Kc = Kc_save(44) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(4) =    wdot(4) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 45: OH + CH2 <=> H + CH2O !
    phi_f = sc(5)*sc(8) 
    k_f = k_f_save(45) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(15) 
    Kc = Kc_save(45) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 46: OH + CH2(S) <=> H + CH2O !
    phi_f = sc(5)*sc(9) 
    k_f = k_f_save(46) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(15) 
    Kc = Kc_save(46) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 47: OH + CH3 <=> CH2 + H2O !
    phi_f = sc(5)*sc(10) 
    k_f = k_f_save(47) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(6) 
    Kc = Kc_save(47) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 48: OH + CH3 <=> CH2(S) + H2O !
    phi_f = sc(5)*sc(10) 
    k_f = k_f_save(48) 
    q_f = phi_f * k_f 
    phi_r = sc(9)*sc(6) 
    Kc = Kc_save(48) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(9) =    wdot(9) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 49: OH + CH4 <=> CH3 + H2O !
    phi_f = sc(5)*sc(11) 
    k_f = k_f_save(49) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(6) 
    Kc = Kc_save(49) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(11) =    wdot(11) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 50: OH + CO <=> H + CO2 !
    phi_f = sc(5)*sc(12) 
    k_f = k_f_save(50) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(13) 
    Kc = Kc_save(50) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 51: OH + HCO <=> H2O + CO !
    phi_f = sc(5)*sc(14) 
    k_f = k_f_save(51) 
    q_f = phi_f * k_f 
    phi_r = sc(6)*sc(12) 
    Kc = Kc_save(51) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 52: OH + CH2O <=> HCO + H2O !
    phi_f = sc(5)*sc(15) 
    k_f = k_f_save(52) 
    q_f = phi_f * k_f 
    phi_r = sc(14)*sc(6) 
    Kc = Kc_save(52) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 53: OH + C2H6 <=> C2H5 + H2O !
    phi_f = sc(5)*sc(19) 
    k_f = k_f_save(53) 
    q_f = phi_f * k_f 
    phi_r = sc(18)*sc(6) 
    Kc = Kc_save(53) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(19) =    wdot(19) - 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 54: HO2 + CH2 <=> OH + CH2O !
    phi_f = sc(7)*sc(8) 
    k_f = k_f_save(54) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(15) 
    Kc = Kc_save(54) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 55: HO2 + CH3 <=> O2 + CH4 !
    phi_f = sc(7)*sc(10) 
    k_f = k_f_save(55) 
    q_f = phi_f * k_f 
    phi_r = sc(4)*sc(11) 
    Kc = Kc_save(55) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(4) =    wdot(4) + 1 * qdot 
    wdot(11) =    wdot(11) + 1 * qdot 

    !/*reaction 56: HO2 + CH3 <=> OH + CH3O !
    phi_f = sc(7)*sc(10) 
    k_f = k_f_save(56) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(16) 
    Kc = Kc_save(56) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(16) =    wdot(16) + 1 * qdot 

    !/*reaction 57: HO2 + CO <=> OH + CO2 !
    phi_f = sc(7)*sc(12) 
    k_f = k_f_save(57) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(13) 
    Kc = Kc_save(57) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(7) =    wdot(7) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 58: CH2 + O2 <=> OH + HCO !
    phi_f = sc(8)*sc(4) 
    k_f = k_f_save(58) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(14) 
    Kc = Kc_save(58) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 

    !/*reaction 59: CH2 + H2 <=> H + CH3 !
    phi_f = sc(8)*sc(1) 
    k_f = k_f_save(59) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(10) 
    Kc = Kc_save(59) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 

    !/*reaction 60: CH2 + CH3 <=> H + C2H4 !
    phi_f = sc(8)*sc(10) 
    k_f = k_f_save(60) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(17) 
    Kc = Kc_save(60) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(17) =    wdot(17) + 1 * qdot 

    !/*reaction 61: CH2 + CH4 <=> 2 CH3 !
    phi_f = sc(8)*sc(11) 
    k_f = k_f_save(61) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(10) 
    Kc = Kc_save(61) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(11) =    wdot(11) - 1 * qdot 
    wdot(10) =    wdot(10) + 2 * qdot 

    !/*reaction 62: CH2(S) + N2 <=> CH2 + N2 !
    phi_f = sc(9)*sc(20) 
    k_f = k_f_save(62) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(20) 
    Kc = Kc_save(62) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(20) =    wdot(20) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(20) =    wdot(20) + 1 * qdot 

    !/*reaction 63: CH2(S) + AR <=> CH2 + AR !
    phi_f = sc(9)*sc(21) 
    k_f = k_f_save(63) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(21) 
    Kc = Kc_save(63) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(21) =    wdot(21) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(21) =    wdot(21) + 1 * qdot 

    !/*reaction 64: CH2(S) + O2 <=> H + OH + CO !
    phi_f = sc(9)*sc(4) 
    k_f = k_f_save(64) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(5)*sc(12) 
    Kc = Kc_save(64) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 65: CH2(S) + O2 <=> CO + H2O !
    phi_f = sc(9)*sc(4) 
    k_f = k_f_save(65) 
    q_f = phi_f * k_f 
    phi_r = sc(12)*sc(6) 
    Kc = Kc_save(65) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 66: CH2(S) + H2 <=> CH3 + H !
    phi_f = sc(9)*sc(1) 
    k_f = k_f_save(66) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(2) 
    Kc = Kc_save(66) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(1) =    wdot(1) - 1 * qdot 
    wdot(10) =    wdot(10) + 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 

    !/*reaction 67: CH2(S) + H2O <=> CH2 + H2O !
    phi_f = sc(9)*sc(6) 
    k_f = k_f_save(67) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(6) 
    Kc = Kc_save(67) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 68: CH2(S) + CH3 <=> H + C2H4 !
    phi_f = sc(9)*sc(10) 
    k_f = k_f_save(68) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(17) 
    Kc = Kc_save(68) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(17) =    wdot(17) + 1 * qdot 

    !/*reaction 69: CH2(S) + CH4 <=> 2 CH3 !
    phi_f = sc(9)*sc(11) 
    k_f = k_f_save(69) 
    q_f = phi_f * k_f 
    phi_r = sc(10)*sc(10) 
    Kc = Kc_save(69) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(11) =    wdot(11) - 1 * qdot 
    wdot(10) =    wdot(10) + 2 * qdot 

    !/*reaction 70: CH2(S) + CO <=> CH2 + CO !
    phi_f = sc(9)*sc(12) 
    k_f = k_f_save(70) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(12) 
    Kc = Kc_save(70) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(12) =    wdot(12) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 71: CH2(S) + CO2 <=> CH2 + CO2 !
    phi_f = sc(9)*sc(13) 
    k_f = k_f_save(71) 
    q_f = phi_f * k_f 
    phi_r = sc(8)*sc(13) 
    Kc = Kc_save(71) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(13) =    wdot(13) - 1 * qdot 
    wdot(8) =    wdot(8) + 1 * qdot 
    wdot(13) =    wdot(13) + 1 * qdot 

    !/*reaction 72: CH2(S) + CO2 <=> CO + CH2O !
    phi_f = sc(9)*sc(13) 
    k_f = k_f_save(72) 
    q_f = phi_f * k_f 
    phi_r = sc(12)*sc(15) 
    Kc = Kc_save(72) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(9) =    wdot(9) - 1 * qdot 
    wdot(13) =    wdot(13) - 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 73: CH3 + O2 <=> O + CH3O !
    phi_f = sc(10)*sc(4) 
    k_f = k_f_save(73) 
    q_f = phi_f * k_f 
    phi_r = sc(3)*sc(16) 
    Kc = Kc_save(73) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 
    wdot(16) =    wdot(16) + 1 * qdot 

    !/*reaction 74: CH3 + O2 <=> OH + CH2O !
    phi_f = sc(10)*sc(4) 
    k_f = k_f_save(74) 
    q_f = phi_f * k_f 
    phi_r = sc(5)*sc(15) 
    Kc = Kc_save(74) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(5) =    wdot(5) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 75: 2 CH3 (+M) <=> C2H6 (+M) !
    phi_f = sc(10)*sc(10) 
    alpha = mixture + sc(1) + 5*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) -0.3*sc(21) 
    k_f = k_f_save(75) 
    redP = 1e-12 * alpha / k_f * 1.77e+32*exp(-9.67*tsc(1)-3130.0076613053565779*invT) 
    F = redP / (1 + redP) 
    logPred = log10(redP) 
    logFcent = log10((0.4675*exp(T/(-151)))+ (0.5325*exp(T/(-1038)))+ (exp(-4970/T))) 
    troe_c = -.4 - .67 * logFcent 
    troe_n = .75 - 1.27 * logFcent 
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
    F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
    F =    F * F_troe 
    k_f =    k_f * F 
    q_f = phi_f * k_f 
    phi_r = sc(19) 
    Kc = Kc_save(75) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 2 * qdot 
    wdot(19) =    wdot(19) + 1 * qdot 

    !/*reaction 76: 2 CH3 <=> H + C2H5 !
    phi_f = sc(10)*sc(10) 
    k_f = k_f_save(76) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(18) 
    Kc = Kc_save(76) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 2 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 

    !/*reaction 77: CH3 + HCO <=> CH4 + CO !
    phi_f = sc(10)*sc(14) 
    k_f = k_f_save(77) 
    q_f = phi_f * k_f 
    phi_r = sc(11)*sc(12) 
    Kc = Kc_save(77) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(11) =    wdot(11) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 78: CH3 + CH2O <=> HCO + CH4 !
    phi_f = sc(10)*sc(15) 
    k_f = k_f_save(78) 
    q_f = phi_f * k_f 
    phi_r = sc(14)*sc(11) 
    Kc = Kc_save(78) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(15) =    wdot(15) - 1 * qdot 
    wdot(14) =    wdot(14) + 1 * qdot 
    wdot(11) =    wdot(11) + 1 * qdot 

    !/*reaction 79: CH3 + C2H6 <=> C2H5 + CH4 !
    phi_f = sc(10)*sc(19) 
    k_f = k_f_save(79) 
    q_f = phi_f * k_f 
    phi_r = sc(18)*sc(11) 
    Kc = Kc_save(79) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(10) =    wdot(10) - 1 * qdot 
    wdot(19) =    wdot(19) - 1 * qdot 
    wdot(18) =    wdot(18) + 1 * qdot 
    wdot(11) =    wdot(11) + 1 * qdot 

    !/*reaction 80: HCO + H2O <=> H + CO + H2O !
    phi_f = sc(14)*sc(6) 
    k_f = k_f_save(80) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(12)*sc(6) 
    Kc = Kc_save(80) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 81: HCO + M <=> H + CO + M !
    phi_f = sc(14) 
    alpha = mixture + sc(1) -1*sc(6) + sc(11) + 0.5*sc(12) + sc(13) + 2*sc(19) 
    k_f = alpha * k_f_save(81) 
    q_f = phi_f * k_f 
    phi_r = sc(2)*sc(12) 
    Kc = Kc_save(81) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(2) =    wdot(2) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 82: HCO + O2 <=> HO2 + CO !
    phi_f = sc(14)*sc(4) 
    k_f = k_f_save(82) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(12) 
    Kc = Kc_save(82) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(14) =    wdot(14) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(12) =    wdot(12) + 1 * qdot 

    !/*reaction 83: CH3O + O2 <=> HO2 + CH2O !
    phi_f = sc(16)*sc(4) 
    k_f = k_f_save(83) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(15) 
    Kc = Kc_save(83) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(16) =    wdot(16) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(15) =    wdot(15) + 1 * qdot 

    !/*reaction 84: C2H5 + O2 <=> HO2 + C2H4 !
    phi_f = sc(18)*sc(4) 
    k_f = k_f_save(84) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(17) 
    Kc = Kc_save(84) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(18) =    wdot(18) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(17) =    wdot(17) + 1 * qdot 

    !/*zero out wdot !
    do n = 1, nspecies
       wdot(n) = wdot(n) * 1.0e-6
    end do
    
    do n = 1, nspecies
       qyn = n + iry1 - 1
       Uprime(i,j,k, qyn) = Uprime(i,j,k,qyn) + wdot(n) * molecular_weight(n)
    end do


!    call ckwyr(q(i,j,k,qrho), q(i,j,k,qtemp), Yt, iwrk, rwrk, wdot)
!    Uprime(i,j,k,iry1:) = Uprime(i,j,k,iry1:) + wdot * molecular_weight
             
 end do
end do
end do
!$omp end parallel do 

end subroutine chemterm_3d


