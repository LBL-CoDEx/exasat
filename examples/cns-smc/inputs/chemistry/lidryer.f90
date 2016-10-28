!9 species
subroutine chemterm_3d(lo,hi,ng,Q,Uprime, ncons, nprim, nspecies)
  integer, intent(in)::ncons, nprim, nspecies

! up is UPrime that has no ghost cells
    integer,          intent(in ) :: lo(3),hi(3),ng
    double precision, intent(in ) :: Q (-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,nprim)
    double precision, intent(inout) :: Uprime(    lo(1):hi(1)   ,    lo(2):hi(2)   ,    lo(3):hi(3)   ,ncons)

    integer, parameter :: qy1 = 8
    integer :: iwrk, i,j,k, n
    double precision :: Yt(nspecies), wdot(nspecies), molecular_weight(nspecies), rwrk

    double precision, save ::k_f_save(21) !nreactions)
    double precision, save ::Kc_save(21) !nreactions)

    double precision, parameter:: imw(9) = (/ 1.0 / 2.01594,&   !/*H2 !
                               1.0 / 31.9988,&   !/*O2 !
                               1.0 / 18.01534,&  !/*H2O !
                               1.0 / 1.00797,&   !/*H !
                               1.0 / 15.9994,&   !/*O !
                               1.0 / 17.00737,&  !/*OH !
                               1.0 / 33.00677,&  !/*HO2 !
                               1.0 / 34.01474,&  !/*H2O2 !
                               1.0 / 28.0134 /)  !/*N2 !

    double precision:: sc(9)  !/*temporary storage !
    double precision:: qdot 

    double precision:: mixture                  !/*mixture concentration !
    double precision:: g_RT(9)                 !/*Gibbs free energy !
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

!             Yt = q(i,j,k,qy1:qy1+nspecies-1)
             rho = Q(i,j,k,qrho)
             T =   Q(i,j,k,qtemp)

             !sc code goes here
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
                mixture = mixture + sc(n) 
             end do

             !gibbs code goes here
             !!!!!!!!!!!!!!!!!!!!!
             !gibbs code goes here
             do n=1, nspecies
                g_RT(n) = &
                     -9.179351730000000e+02 * invT &
                     +1.661320882000000e+00 &
                     -2.344331120000000e+00 * tsc(1) &
                     -3.990260375000000e-03 * tsc(2) &
                     +3.246358500000000e-06 * tsc(3) &
                     -1.679767450000000e-09 * tsc(4) &
                     +3.688058805000000e-13 * tsc(5) 
             end do

             !!!!!!!!!!!!!!!!!!!!!

             !/*zero out wdot !
             do n = 1,nspecies
                wdot(n) = 0.0 
             end do
             
             if (T .ne. T_save) then
                T_save = T 

                !!!!!!!!!!!!!!!!!!!!!
                !production rate goes here
                !!!!!!!!!!!!!!!!!!!!!

                k_f_save(1) = 1e-06 * 3.547e+15*exp(-0.406*tsc(1)-8352.8934356925419706*invT) 
                k_f_save(2) = 1e-06 * 50800*exp(2.67*tsc(1)-3165.2328279116868543*invT) 
                k_f_save(3) = 1e-06 * 2.16e+08*exp(1.51*tsc(1)-1726.0331637101885462*invT) 
                k_f_save(4) = 1e-06 * 2.97e+06*exp(2.02*tsc(1)-6743.1033217832446098*invT) 
                k_f_save(5) = 1e-06 * 4.577e+19*exp(-1.4*tsc(1)-52525.75557669664704*invT) 
                k_f_save(6) = 1e-12 * 6.165e+15*exp(-0.5*tsc(1)) 
                k_f_save(7) = 1e-12 * 4.714e+18*exp(-1*tsc(1)) 
                k_f_save(8) = 1e-12 * 3.8e+22*exp(-2*tsc(1)) 
                k_f_save(9) = 1e-06 * 1.475e+12*exp(0.6*tsc(1)) 
                k_f_save(10) = 1e-06 * 1.66e+13*exp(-414.14731595728432012*invT) 
                k_f_save(11) = 1e-06 * 7.079e+13*exp(-148.44891641239232172*invT) 
                k_f_save(12) = 1e-06 * 3.25e+13 
                k_f_save(13) = 1e-06 * 2.89e+13*exp(+250.09868290494571852*invT) 
                k_f_save(14) = 1e-06 * 4.2e+14*exp(-6029.5420896721516328*invT) 
                k_f_save(15) = 1e-06 * 1.3e+11*exp(+819.89091359562974048*invT) 
                k_f_save(16) = 1 * 2.951e+14*exp(-24370.783124922574643*invT) 
                k_f_save(17) = 1e-06 * 2.41e+13*exp(-1997.7701632447369775*invT) 
                k_f_save(18) = 1e-06 * 4.82e+13*exp(-4000.5724931475210724*invT) 
                k_f_save(19) = 1e-06 * 9.55e+06*exp(2*tsc(1)-1997.7701632447369775*invT) 
                k_f_save(20) = 1e-06 * 1e+12 
                k_f_save(21) = 1e-06 * 5.8e+14*exp(-4809.2416750957054319*invT) 

                Kc_save(1) = exp((g_RT(4) + g_RT(2)) - (g_RT(5) + g_RT(6))) 
                Kc_save(2) = exp((g_RT(5) + g_RT(1)) - (g_RT(4) + g_RT(6))) 
                Kc_save(3) = exp((g_RT(1) + g_RT(6)) - (g_RT(3) + g_RT(4))) 
                Kc_save(4) = exp((g_RT(5) + g_RT(3)) - (g_RT(6) + g_RT(6))) 
                Kc_save(5) = refC * exp((g_RT(1)) - (g_RT(4) + g_RT(4))) 
                Kc_save(6) = 1.0 / (refC) * exp((g_RT(5) + g_RT(5)) - (g_RT(2))) 
                Kc_save(7) = 1.0 / (refC) * exp((g_RT(5) + g_RT(4)) - (g_RT(6))) 
                Kc_save(8) = 1.0 / (refC) * exp((g_RT(4) + g_RT(6)) - (g_RT(3))) 
                Kc_save(9) = 1.0 / (refC) * exp((g_RT(4) + g_RT(2)) - (g_RT(7))) 
                Kc_save(10) = exp((g_RT(7) + g_RT(4)) - (g_RT(1) + g_RT(2))) 
                Kc_save(11) = exp((g_RT(7) + g_RT(4)) - (g_RT(6) + g_RT(6))) 
                Kc_save(12) = exp((g_RT(7) + g_RT(5)) - (g_RT(2) + g_RT(6))) 
                Kc_save(13) = exp((g_RT(7) + g_RT(6)) - (g_RT(3) + g_RT(2))) 
                Kc_save(14) = exp((g_RT(7) + g_RT(7)) - (g_RT(8) + g_RT(2))) 
                Kc_save(15) = exp((g_RT(7) + g_RT(7)) - (g_RT(8) + g_RT(2))) 
                Kc_save(16) = refC * exp((g_RT(8)) - (g_RT(6) + g_RT(6))) 
                Kc_save(17) = exp((g_RT(8) + g_RT(4)) - (g_RT(3) + g_RT(6))) 
                Kc_save(18) = exp((g_RT(8) + g_RT(4)) - (g_RT(7) + g_RT(1))) 
                Kc_save(19) = exp((g_RT(8) + g_RT(5)) - (g_RT(6) + g_RT(7))) 
                Kc_save(20) = exp((g_RT(8) + g_RT(6)) - (g_RT(7) + g_RT(3))) 
                Kc_save(21) = exp((g_RT(8) + g_RT(6)) - (g_RT(7) + g_RT(3))) 
             end if

             !/*reaction 1: H + O2 <=> O + OH !
             phi_f = sc(4)*sc(2) 
             k_f = k_f_save(1) 
             q_f = phi_f * k_f 
             phi_r = sc(5)*sc(6) 
             Kc = Kc_save(1) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(2) =    wdot(2) - 1 * qdot 
             wdot(5) =    wdot(5) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             
             !/*reaction 2: O + H2 <=> H + OH !
             phi_f = sc(5)*sc(1) 
             k_f = k_f_save(2) 
             q_f = phi_f * k_f 
             phi_r = sc(4)*sc(6) 
             Kc = Kc_save(2) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(1) =    wdot(1) - 1 * qdot 
             wdot(4) =    wdot(4) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             
             !/*reaction 3: H2 + OH <=> H2O + H !
             phi_f = sc(1)*sc(6) 
             k_f = k_f_save(3) 
             q_f = phi_f * k_f 
             phi_r = sc(3)*sc(4) 
             Kc = Kc_save(3) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(1) =    wdot(1) - 1 * qdot 
             wdot(6) =    wdot(6) - 1 * qdot 
             wdot(3) =    wdot(3) + 1 * qdot 
             wdot(4) =    wdot(4) + 1 * qdot 
             
             !/*reaction 4: O + H2O <=> OH + OH !
             phi_f = sc(5)*sc(3) 
             k_f = k_f_save(4) 
             q_f = phi_f * k_f 
             phi_r = sc(6)*sc(6) 
             Kc = Kc_save(4) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(3) =    wdot(3) - 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             
             !/*reaction 5: H2 + M <=> H + H + M !
             phi_f = sc(1) 
             alpha = mixture + 1.5*sc(1) + 11*sc(3) 
             k_f = alpha * k_f_save(5) 
             q_f = phi_f * k_f 
             phi_r = sc(4)*sc(4) 
             Kc = Kc_save(5) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(1) =    wdot(1) - 1 * qdot 
             wdot(4) =    wdot(4) + 1 * qdot 
             wdot(4) =    wdot(4) + 1 * qdot 
             
             !/*reaction 6: O + O + M <=> O2 + M !
             phi_f = sc(5)*sc(5) 
             alpha = mixture + 1.5*sc(1) + 11*sc(3) 
             k_f = alpha * k_f_save(6) 
             q_f = phi_f * k_f 
             phi_r = sc(2) 
             Kc = Kc_save(6) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 
             
             !/*reaction 7: O + H + M <=> OH + M !
             phi_f = sc(5)*sc(4) 
             alpha = mixture + 1.5*sc(1) + 11*sc(3) 
             k_f = alpha * k_f_save(7) 
             q_f = phi_f * k_f 
             phi_r = sc(6) 
             Kc = Kc_save(7) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             
             !/*reaction 8: H + OH + M <=> H2O + M !
             phi_f = sc(4)*sc(6) 
             alpha = mixture + 1.5*sc(1) + 11*sc(3) 
             k_f = alpha * k_f_save(8) 
             q_f = phi_f * k_f 
             phi_r = sc(3) 
             Kc = Kc_save(8) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(6) =    wdot(6) - 1 * qdot 
             wdot(3) =    wdot(3) + 1 * qdot 
             
             !/*reaction 9: H + O2 (+M) <=> HO2 (+M) !
             phi_f = sc(4)*sc(2) 
             alpha = mixture + sc(1) + 10*sc(3) -0.22*sc(2) 
             k_f = k_f_save(9) 
             redP = 1e-12 * alpha / k_f * 6.366e+20*exp(-1.72*tsc(1)-264.08810621431689469*invT) 
             F = redP / (1 + redP) 
             logPred = log10(redP) 
             logFcent = log10((0.2*exp(T/(-e-30)))+ (0.8*exp(T/(-e+30)))) 
             troe_c = -.4 - .67 * logFcent 
             troe_n = .75 - 1.27 * logFcent 
             troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
             F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
             F =    F * F_troe 
             k_f =    k_f * F 
             q_f = phi_f * k_f 
             phi_r = sc(7) 
             Kc = Kc_save(9) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(2) =    wdot(2) - 1 * qdot 
             wdot(7) =    wdot(7) + 1 * qdot 
             
             !/*reaction 10: HO2 + H <=> H2 + O2 !
             phi_f = sc(7)*sc(4) 
             k_f = k_f_save(10) 
             q_f = phi_f * k_f 
             phi_r = sc(1)*sc(2) 
             Kc = Kc_save(10) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(1) =    wdot(1) + 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 

             !/*reaction 11: HO2 + H <=> OH + OH !
             phi_f = sc(7)*sc(4) 
             k_f = k_f_save(11) 
             q_f = phi_f * k_f 
             phi_r = sc(6)*sc(6) 
             Kc = Kc_save(11) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(4) =    wdot(4) - 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 

             !/*reaction 12: HO2 + O <=> O2 + OH !
             phi_f = sc(7)*sc(5) 
             k_f = k_f_save(12) 
             q_f = phi_f * k_f 
             phi_r = sc(2)*sc(6) 
             Kc = Kc_save(12) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(5) =    wdot(5) - 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             
             !/*reaction 13: HO2 + OH <=> H2O + O2 !
             phi_f = sc(7)*sc(6) 
             k_f = k_f_save(13) 
             q_f = phi_f * k_f 
             phi_r = sc(3)*sc(2) 
             Kc = Kc_save(13) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(6) =    wdot(6) - 1 * qdot 
             wdot(3) =    wdot(3) + 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 

             !/*reaction 14: HO2 + HO2 <=> H2O2 + O2 !
             phi_f = sc(7)*sc(7) 
             k_f = k_f_save(14) 
             q_f = phi_f * k_f 
             phi_r = sc(8)*sc(2) 
             Kc = Kc_save(14) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(8) =    wdot(8) + 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 
             
             !/*reaction 15: HO2 + HO2 <=> H2O2 + O2 !
             phi_f = sc(7)*sc(7) 
             k_f = k_f_save(15) 
             q_f = phi_f * k_f 
             phi_r = sc(8)*sc(2) 
             Kc = Kc_save(15) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(7) =    wdot(7) - 1 * qdot 
             wdot(8) =    wdot(8) + 1 * qdot 
             wdot(2) =    wdot(2) + 1 * qdot 
             
             !/*reaction 16: H2O2 (+M) <=> OH + OH (+M) !
             phi_f = sc(8) 
             alpha = mixture + 1.5*sc(1) + 11*sc(3) 
             k_f = k_f_save(16) 
             redP = 1e-6 * alpha / k_f * 1.202e+17*exp(-22896.358294114746968*invT) 
             F = redP / (1 + redP) 
             logPred = log10(redP) 
             logFcent = log10((0.5*exp(T/(-e-30)))+ (0.5*exp(T/(-1e+30)))) 
             troe_c = -.4 - .67 * logFcent 
             troe_n = .75 - 1.27 * logFcent 
             troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred)) 
             F_troe = pow(10, logFcent / (1.0 + troe*troe)) 
             F =    F * F_troe 
             k_f =    k_f * F 
             q_f = phi_f * k_f 
             phi_r = sc(6)*sc(6) 
             Kc = Kc_save(16) 
             k_r = k_f / Kc 
             q_r = phi_r * k_r 
             qdot = q_f - q_r 
             wdot(8) =    wdot(8) - 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 
             wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 17: H2O2 + H <=> H2O + OH !
    phi_f = sc(8)*sc(4) 
    k_f = k_f_save(17) 
    q_f = phi_f * k_f 
    phi_r = sc(3)*sc(6) 
    Kc = Kc_save(17) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 

    !/*reaction 18: H2O2 + H <=> HO2 + H2 !
    phi_f = sc(8)*sc(4) 
    k_f = k_f_save(18) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(1) 
    Kc = Kc_save(18) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(4) =    wdot(4) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(1) =    wdot(1) + 1 * qdot 

    !/*reaction 19: H2O2 + O <=> OH + HO2 !
    phi_f = sc(8)*sc(5) 
    k_f = k_f_save(19) 
    q_f = phi_f * k_f 
    phi_r = sc(6)*sc(7) 
    Kc = Kc_save(19) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(5) =    wdot(5) - 1 * qdot 
    wdot(6) =    wdot(6) + 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 

    !/*reaction 20: H2O2 + OH <=> HO2 + H2O !
    phi_f = sc(8)*sc(6) 
    k_f = k_f_save(20) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(3) 
    Kc = Kc_save(20) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 

    !/*reaction 21: H2O2 + OH <=> HO2 + H2O !
    phi_f = sc(8)*sc(6) 
    k_f = k_f_save(21) 
    q_f = phi_f * k_f 
    phi_r = sc(7)*sc(3) 
    Kc = Kc_save(21) 
    k_r = k_f / Kc 
    q_r = phi_r * k_r 
    qdot = q_f - q_r 
    wdot(8) =    wdot(8) - 1 * qdot 
    wdot(6) =    wdot(6) - 1 * qdot 
    wdot(7) =    wdot(7) + 1 * qdot 
    wdot(3) =    wdot(3) + 1 * qdot 


    !/*convert to chemkin units */
    do n=1, nspecies 
        wdot(n) = wdot(n) *1.0e-6
    end do

             !             call ckwyr(q(i,j,k,qrho), q(i,j,k,qtemp), Yt, iwrk, rwrk, wdot)
             !up(i,j,k,iry1:) = up(i,j,k,iry1:) + wdot * molecular_weight
    do n = 1, nspecies
       qyn = n + iry1 - 1
       Uprime(i,j,k, qyn) = Uprime(i,j,k,qyn) + wdot(n) * molecular_weight(n)
    end do
    
         end do
       end do
    end do
!$omp end parallel do 

end subroutine chemterm_3d


