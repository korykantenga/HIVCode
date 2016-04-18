!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Algorithm to Solve for Value Functions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE solvevfun

    USE mod_globalvar
    IMPLICIT NONE

    !Solve for Ex-Post Long-Term Value Functions given W_s^bbeta(pphi,x), W_s^iiota(pphi,x)

    !Ex-Post LT
    !call vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap,bprob)

    !Ex-Ante LT
    !call vf_exante_l(wfun_bl,wfun_il,vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !Ex-Post ST

    !Ex-Ante ST + Policy Function

END SUBROUTINE solvevfun

!------------------------------------------------------------------------------
! Solve Ex-Post Long-Term Value Functions
!------------------------------------------------------------------------------
! d = bbeta or iiota, z = l
SUBROUTINE vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap,bprob)

    !Input  : W_s^bbeta(pphi,x),W_s^iiota(pphi,x),bprob (probability mass func for discretised b)
    !Output : V_l^bbeta(pphi,varpphi,x)
    !dim(pphi) = pmap, dim(x) = bmap

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER bmap,pmap
    DOUBLE PRECISION, INTENT(in) :: wfun_b(pmap,bmap),wfun_i(pmap,bmap),bprob(bmap)
    DOUBLE PRECISION, INTENT(out) :: vfun_bl(2,2,bmap),vfun_il(2,2,bmap)
    DOUBLE PRECISION uut
    DOUBLE PRECISION vhelp1(bmap),chelp1,vhelp2(bmap),chelp2
    DOUBLE PRECISION vhelp3(bmap),chelp3,chelp4

    !DOT_PRODUCT(bprob,)
    !V_l^bbeta(0,0,x)
    CALL util(uut,qincome-mcost)
    vhelp1 = aalphapr*AAIDS + (1.0-xxi)*(aalphapr*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:))) + &
        xxi*((1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)))
    chelp1 = 1.0-bbeta*(1.0-xxi)*((1.0-aalphapr)**2)

    vfun_bl(1,1,:) = (uut + ppref_o + bbeta*vhelp1)/(chelp1)

    !V_l^bbeta(1,0,x)
    CALL util(uut,qincome)
    vhelp2 = (1.0-ggamma_p)*aalpha*AAIDS + &
        (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        aalphapr*ggamma_p*DOT_PRODUCT(bprob,wfun_b(pmap,:)) + &
        aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*DOT_PRODUCT(bprob,wfun_b(1,:))) + &
        xxi*((1.0-aalpha)*(1.0-ggamma_p)*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        ggamma_p*DOT_PRODUCT(bprob,wfun_b(pmap,:)))
    chelp2 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(2,1,:) = (uut + ppref_o + bbeta*vhelp2)/(chelp2)

    !V_l^bbeta(0,1,x)
    CALL util(uut,qincome-mcost)
    vhelp3 = aalphapr*AAIDS + &
        (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        aalpha*(1.0-aalphapr)*(1.0-ggamma_p)*DOT_PRODUCT(bprob,wfun_b(1,:))) + &
        xxi*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:))
    chelp3 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(1,2,:) = (uut + ppref_o + bbeta*vhelp3)/(chelp3)

    !V_l^bbeta(1,1,x)
    CALL util(uut,qincome)
    chelp4 = 1.0 - bbeta*(1.0-xxi)

    vfun_bl(2,2,:) = (uut + lpref_o + xxi*bbeta*DOT_PRODUCT(bprob,wfun_b(pmap,:)))/chelp4

    !---------------------------------------------------------------------

    !V_l^iiota(0,0,x)
    CALL util(uut,qincome-mcost)
    vhelp1 = aalphapr*AAIDS + &
        (1.0-xxi)*(aalphapr*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:))+ &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
        ((1.0-aalphapr)**2)*vfun_bl(1,1,:)) + &
        xxi*((1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))))
    chelp1 = 1.0-iiota*(1.0-xxi)*(1.0-eeta)*((1.0-aalphapr)**2)

    vfun_il(1,1,:) = (uut + ppref_y + iiota*vhelp1)/(chelp1)

    !V_l^iiota(1,0,x)
    CALL util(uut,qincome)
    vhelp2 = (1.0-ggamma_p)*aalpha*AAIDS + &
        (1.0-xxi)*eeta*(1.0-aalphapr)*ggamma_p*vfun_bl(2,1,:) + &
        (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*(eeta*vfun_bl(1,1,:) + &
        (1.0-eeta)*vfun_il(1,1,:)) + &
        aalphapr*ggamma_p*(eeta*DOT_PRODUCT(bprob,wfun_b(pmap,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(pmap,:))) + &
        aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*(DOT_PRODUCT(bprob,wfun_b(1,:))*eeta + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:)))) + &
        xxi*((1.0-aalpha)*(1.0-ggamma_p)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
        ggamma_p*(eeta*DOT_PRODUCT(bprob,wfun_b(pmap,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(pmap,:))))
    chelp2 = 1.0 - iiota*(1.0-xxi)*(1.0-eeta)*((1.0-aalphapr)*ggamma_p)

    vfun_il(2,1,:) = (uut + ppref_y + iiota*vhelp2)/(chelp2)

    !V_l^iiota(0,1,x)
    CALL util(uut,qincome-mcost)
    vhelp3 = aalphapr*AAIDS + (1.0-xxi)*((1.0-aalphapr)*ggamma_p*eeta*vfun_bl(1,2,:) + &
        (1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*(eeta*vfun_bl(1,1,:)+ &
        (1.0-eeta)*vfun_il(1,1,:)) + &
        aalpha*(1.0-aalphapr)*(1.0-ggamma_p)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:)))) + &
        xxi*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:)))
    chelp3 = 1.0 - iiota*(1.0-xxi)*(1.0-eeta)*(1.0-aalphapr)*ggamma_p

    vfun_il(1,2,:) = (uut + ppref_y + iiota*vhelp3)/(chelp3)

    !V_l^iiota(1,1,x)
    CALL util(uut,qincome)
    chelp4 = 1.0 - iiota*(1.0-eeta)*(1.0-xxi)

    vfun_il(2,2,:) = (uut + lpref_y + iiota*eeta*(1.0-xxi)*vfun_bl(2,2,:) + &
        xxi*iiota*(eeta*DOT_PRODUCT(bprob,wfun_b(pmap,:))+(1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(pmap,:))))/chelp4

END SUBROUTINE vf_expost_l

!------------------------------------------------------------------------------
! Solve Ex-Ante Long-Term Value Functions
!------------------------------------------------------------------------------
SUBROUTINE vf_exante_l(wfun_bl,wfun_il,vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !Input  : V(pphi,varpphi,x), W_s^d(pphi,x)
    !Output : W_l^d(pphi,x)

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER bmap, pmap
    DOUBLE PRECISION, INTENT(in) :: wfun_b(pmap,bmap),wfun_i(pmap,bmap)
    DOUBLE PRECISION, INTENT(in) :: vfun_bl(2,2,bmap),vfun_il(2,2,bmap)
    DOUBLE PRECISION, INTENT(out) :: wfun_bl(2,bmap),wfun_il(2,bmap)
    DOUBLE PRECISION vhelp1(bmap,2),vhelp2(bmap,2)

    !W_l^bbeta(pphi,x)

    vhelp1(:,1) = vfun_bl(1,2,:)-wfun_b(1,:) !HIV+,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_bl(1,1,:)-wfun_b(1,:) !HIV+,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_bl(1,:) = nnu_l*MAXVAL(vhelp1,2) + (1.0-nnu_l)*MAXVAL(vhelp2,2) + wfun_b(1,:)

    vhelp1(:,1) = vfun_bl(2,2,:)-wfun_b(pmap,:) !HIV-,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_bl(2,1,:)-wfun_b(pmap,:) !HIV-,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_bl(2,:) = nnu_l*MAXVAL(vhelp1,2) + (1.0-nnu_l)*MAXVAL(vhelp2,2) + wfun_b(pmap,:)

    !-----------------------------------------------------------------------------------
    !W_l^iiota(pphi,x)

    vhelp1(:,1) = vfun_il(1,2,:)-wfun_i(1,:) !HIV+,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_il(1,1,:)-wfun_i(1,:) !HIV+,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_il(1,:) = nnu_l*MAXVAL(vhelp1,2) + (1.0-nnu_l)*MAXVAL(vhelp2,2) + wfun_i(1,:)

    vhelp1(:,1) = vfun_il(2,2,:)-wfun_i(pmap,:) !HIV-,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_il(2,1,:)-wfun_i(pmap,:) !HIV-,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_il(2,:) = nnu_l*MAXVAL(vhelp1,2) + (1.0-nnu_l)*MAXVAL(vhelp2,2) + wfun_i(pmap,:)


END SUBROUTINE vf_exante_l

!------------------------------------------------------------------------------
! Solve Ex-Post Short-Term Value Functions
!------------------------------------------------------------------------------
SUBROUTINE vf_expost_s(vfa_b,vfb_b,vfp_b,vfa_i,vfb_i,vfp_i, &
    wfun_bl,wfun_il,wfun_b,wfun_i,pmap,bmap,bprob,pphigrid,bgrid)

    ! Inputs : Ex-Ante Value Functions
    ! Outputs: Ex-Post Short Term Value Functions

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER bmap,pmap, iCount, bCount
    DOUBLE PRECISION, INTENT(in) :: bprob(bmap),pphigrid(pmap),bgrid(bmap)
    DOUBLE PRECISION, INTENT(in) :: wfun_b(pmap,bmap),wfun_i(pmap,bmap)
    DOUBLE PRECISION, INTENT(in) :: wfun_bl(2,bmap),wfun_il(2,bmap)

    DOUBLE PRECISION, INTENT(out) :: vfa_b(pmap),vfa_i(pmap)
    DOUBLE PRECISION, INTENT(out) :: vfb_b(pmap,bmap),vfb_i(pmap,bmap)
    DOUBLE PRECISION, INTENT(out) :: vfp_b(pmap),vfp_i(pmap)

    DOUBLE PRECISION uut,probhelp,pphitd(pmap)
    DOUBLE PRECISION pphia(pmap),pphib(pmap),pphip(pmap)        !Updated phi'
    DOUBLE PRECISION wfun_b_apr(pmap,bmap),wfun_i_apr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_A(phi)
    DOUBLE PRECISION wfun_b_bpr(pmap,bmap),wfun_i_bpr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_B(phi)
    DOUBLE PRECISION wfun_b_ppr(pmap,bmap),wfun_i_ppr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_P(phi)

    call PPhiPRIME(pphia,pphib,pphip,pphigrid,pmap)

    !Interpolate Continuation Values
    !$OMP PARALLEL DO
    do iCount = 1,bmap
        !W_s^bbeta(phi',x)
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphia,wfun_b_apr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphib,wfun_b_bpr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphip,wfun_b_ppr(:,iCount))
        !W_s^iiota(phi',x)
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphia,wfun_i_apr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphib,wfun_i_bpr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphip,wfun_i_ppr(:,iCount))
    end do
    !$OMP END PARALELL DO

    !Abstinence, BBETA, PPhi_A(pphipr,pphi)
    CALL util(uut,qincome)
    vfa_b = uut + (1.0-pphigrid)*aalpha*bbeta*AAIDS + bbeta*(1.0-((1.0-pphigrid)*aalpha))*( &
        (1.0-eepsilon)*MATMUL(wfun_b_apr,bprob) + &
        eepsilon*(pphigrid*MATMUL(wfun_b_apr,bprob) + &
        (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:))))
    CALL util(uut,qincome-mcost)
    vfa_b(1) = uut + bbeta*aalphapr*AAIDS + bbeta*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) !Replace for HIV+ aware

    !Abstinence, IIOTA, PPhi_A(pphipr,pphi)
    CALL util(uut,qincome)
    vfa_i = uut + (1.0-pphigrid)*aalpha*iiota*AAIDS + iiota*(1.0-((1.0-pphigrid)*aalpha))*( &
        (1.0-eepsilon)*(eeta*MATMUL(wfun_b_apr,bprob) + (1.0-eeta)*MATMUL(wfun_i_apr,bprob)) + &
        eepsilon*(pphigrid*(eeta*MATMUL(wfun_b_apr,bprob) + (1.0-eeta)*MATMUL(wfun_i_apr,bprob)) + &
        (1.0-pphigrid)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:)))))
    CALL util(uut,qincome-mcost)
    vfa_i(1) = uut + iiota*aalphapr*AAIDS + iiota*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) !Replace for HIV+ aware

    !Protected Sex, BBETA, PPhi_P(pphipr,pphi)
    CALL util(uut,qincome)
    CALL PR_p(probhelp,pphigrid)
    pphitd = pphigrid - pphigrid*(1.0-nnu_p)*(1.0-ggamma_p)
    vfp_b  = uut + ppref_o + bbeta*(1.0-probhelp)*AAIDS + bbeta*(1.0-mmu)*probhelp*( &
        (1.0-eepsilon)*MATMUL(wfun_b_ppr,bprob) + eepsilon*( &
        pphigrid*MATMUL(wfun_b_ppr,bprob) + &
        (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:)))) + &
        bbeta*mmu*probhelp*(pphitd*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
        (1.0-pphitd)*DOT_PRODUCT(bprob,wfun_bl(1,:)))
    CALL util(uut,qincome-mcost)
    vfp_b(1) = uut + ppref_o + bbeta*aalphapr*AAIDS + &
        bbeta*(1.0-mmu)*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        bbeta*mmu*DOT_PRODUCT(bprob,wfun_bl(1,:))

    !Protected Sex, IIOTA, PPhi_P(pphipr,pphi)
    CALL util(uut,qincome)
    !CALL PR_p(probhelp,pphigrid)
    !pphitd = pphigrid - pphigrid*(1.0-nnu_p)*(1.0-ggamma_p)
    vfp_i  = uut + ppref_y + bbeta*(1.0-probhelp)*AAIDS + bbeta*(1.0-mmu)*probhelp*( &
        (1.0-eepsilon)*(eeta*MATMUL(wfun_b_ppr,bprob) + (1.0-eeta)*MATMUL(wfun_i_ppr,bprob)) + &
        eepsilon*( &
        pphigrid*(eeta*MATMUL(wfun_b_ppr,bprob) + (1.0-eeta)*MATMUL(wfun_i_ppr,bprob)) + &
        (1.0-pphigrid)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))))) + &
        bbeta*mmu*probhelp*(pphitd*(eeta*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(2,:))) + &
        (1.0-pphitd)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:))))
    CALL util(uut,qincome-mcost)
    vfp_i(1) = uut + ppref_o + bbeta*aalphapr*AAIDS + &
        bbeta*(1.0-mmu)*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
        bbeta*mmu*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:)))

    !Unprotected Sex
    CALL PR_b(probhelp,pphigrid)
    pphitd = pphigrid - pphigrid*(1.0-nnu_b)*(1.0-ggamma_b)

    do bCount = 1,bmap

        !Unprotected Sex, BBETA, PPhi_B(pphipr,pphi)
        CALL util(uut,qincome)
        vfb_b(:,bCount) = uut + bgrid(bCount) + bbeta*(1.0-probhelp)*AAIDS + &
            bbeta*(1.0-mmu)*probhelp*( &
            (1.0-eepsilon)*MATMUL(wfun_b_ppr,bprob) + eepsilon*( &
            pphigrid*MATMUL(wfun_b_ppr,bprob) + &
            (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:)))) + &
            bbeta*mmu*probhelp*(pphitd*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
            (1.0-pphitd)*DOT_PRODUCT(bprob,wfun_bl(1,:)))
        CALL util(uut,qincome-mcost)
        vfb_b(1,bCount) = uut + bgrid(bCount) + bbeta*aalphapr*AAIDS + &
            bbeta*(1.0-mmu)*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
            bbeta*mmu*DOT_PRODUCT(bprob,wfun_bl(1,:))

        !Unprotected Sex, IIOTA, PPhi_B(pphipr,pphi)
        CALL util(uut,qincome)
        vfb_i(:,bCount)  = uut + bgrid(bCount) + bbeta*(1.0-probhelp)*AAIDS + &
            bbeta*(1.0-mmu)*probhelp*( &
            (1.0-eepsilon)*(eeta*MATMUL(wfun_b_ppr,bprob) + (1.0-eeta)*MATMUL(wfun_i_ppr,bprob)) + &
            eepsilon*( &
            pphigrid*(eeta*MATMUL(wfun_b_ppr,bprob) + (1.0-eeta)*MATMUL(wfun_i_ppr,bprob)) + &
            (1.0-pphigrid)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))))) + &
            bbeta*mmu*probhelp*(pphitd*(eeta*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(2,:))) + &
            (1.0-pphitd)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:))))
        CALL util(uut,qincome-mcost)
        vfb_i(1,bCount) = uut + bgrid(bCount) + bbeta*aalphapr*AAIDS + &
            bbeta*(1.0-mmu)*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
            bbeta*mmu*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:)))

    end do

END SUBROUTINE vf_expost_s

!------------------------------------------------------------------------------
! Solve Ex-Ante Short-Term Value Functions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Solve Policy Functions using FOC and Corner Solution
!------------------------------------------------------------------------------
SUBROUTINE solvepolicyfun_p(policy_p,vfun_a,vfun_p,smap)
    !Input: Vectorized Value Functions

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER smap
    DOUBLE PRECISION, INTENT(in) :: vfun_a(smap),vfun_p(smap)
    DOUBLE PRECISION, INTENT(out) :: policy_p(smap)
    DOUBLE PRECISION vhelp(smap,2),phelp1(smap),phelp2(smap)
    DOUBLE PRECISION ppower

    ppower = 1.0/kkappa_p

    vhelp(:,1) = vfun_p - vfun_a
    CALL zerov(vhelp(:,2),smap)

    phelp1 = MAXVAL(vhelp,2)

    phelp2 = (1.0/(oomega_p*(kkappa_p+1.0)))*(phelp1**ppower)

    policy_p = phelp2/(1.0+phelp2)


END SUBROUTINE solvepolicyfun_p

SUBROUTINE solvepolicyfun_b(policy_b,vfun_a,vfun_b,smap)
    !Input: Vectorized Value Functions

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER smap
    DOUBLE PRECISION, INTENT(in) :: vfun_a(smap),vfun_b(smap)
    DOUBLE PRECISION, INTENT(out) :: policy_b(smap)
    DOUBLE PRECISION vhelp(smap,2),phelp1(smap),phelp2(smap)
    DOUBLE PRECISION ppower

    ppower = 1.0/kkappa_b

    vhelp(:,1) = vfun_b - vfun_a
    CALL zerov(vhelp(:,2),smap)

    phelp1 = MAXVAL(vhelp,2)

    phelp2 = (1.0/(oomega_b*(kkappa_b+1.0)))*(phelp1**ppower)

    policy_b = phelp2/(1.0+phelp2)


END SUBROUTINE solvepolicyfun_b
