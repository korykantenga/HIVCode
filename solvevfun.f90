!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Algorithm to Solve for Value Functions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE solvevfun

    USE mod_globalvar
    IMPLICIT NONE

    !Solve for Ex-Post Long-Term Value Functions given W_s^bbeta(pphi,x), W_s^iiota(pphi,x)

    !-------------------------
    !Value Functions
    !-------------------------

    !-----------------
    !Ex-Post LT
    !-----------------
    !call vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap,bprob)

    !-----------------
    !Ex-Ante LT
    !-----------------
    !call vf_exante_l(wfun_bl,wfun_il,vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !-----------------
    !Ex-Post ST
    !-----------------
    !call vf_expost_s(vfa_b,vfb_b,vfp_b,vfa_i,vfb_i,vfp_i, &
    !                 wfun_bl,wfun_il,wfun_b,wfun_i,pmap,bmap,bprob,pphigrid,bgrid)

    !-------------------------
    !Policy Functions
    !-------------------------
    !-----------------
    !Protected, Old
    !-----------------
    !call solvepolicyfun_p(policyp_bbeta,vfa_b,vfp_b,pmap)
    !-----------------
    !Protected, Young
    !-----------------
    !call solvepolicyfun_p(policyp_iiota,vfa_i,vfp_i,pmap)
    !-----------------
    !Unprotected, Old
    !-----------------
    !call vec2VF(vecfun_b,vfb_b,pmap,bmap)
    !call solvepolicyfun_b(policy_b,vfa_b,vecfun_b,pmap*bmap)
    !-----------------
    !Unprotected, Young
    !-----------------

    !Ex-Ante ST


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

    DOUBLE PRECISION uut,uutpos,probhelp,pphitd(pmap)
    DOUBLE PRECISION pphia(pmap),pphib(pmap),pphip(pmap)        !Updated phi'
    DOUBLE PRECISION wfun_b_apr(pmap,bmap),wfun_i_apr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_A(phi)
    DOUBLE PRECISION wfun_b_bpr(pmap,bmap),wfun_i_bpr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_B(phi)
    DOUBLE PRECISION wfun_b_ppr(pmap,bmap),wfun_i_ppr(pmap,bmap)!Interpolation W(phi',x),phi'=PPhi_P(phi)

    CALL PPhiPRIME(pphia,pphib,pphip,pphigrid,pmap)
    CALL util(uut,qincome)
    CALL util(uutpos,qincome-mcost)

    !Interpolate Continuation Values
    !$OMP PARALLEL DO
    do iCount = 1,bmap
        !W_s^bbeta(phi',x)
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphia,wfun_b_apr(:,iCount))
        !print *,wfun_b_apr(:,iCount)
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphib,wfun_b_bpr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_b(:,iCount),pmap,pphip,wfun_b_ppr(:,iCount))
        !W_s^iiota(phi',x)
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphia,wfun_i_apr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphib,wfun_i_bpr(:,iCount))
        call LinInterp(pmap,pphigrid,wfun_i(:,iCount),pmap,pphip,wfun_i_ppr(:,iCount))
    end do
    !$OMP END PARALELL DO

    !Abstinence, BBETA, PPhi_A(pphipr,pphi)
    vfa_b = uut + (1.0-pphigrid)*aalpha*bbeta*AAIDS + bbeta*(1.0-((1.0-pphigrid)*aalpha))*( &
        (1.0-eepsilon)*MATMUL(wfun_b_apr,bprob) + &
        eepsilon*(pphigrid*MATMUL(wfun_b_apr,bprob) + &
        (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:))))
    vfa_b(1) = uutpos + bbeta*aalphapr*AAIDS + bbeta*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) !Replace for HIV+ aware

    !Abstinence, IIOTA, PPhi_A(pphipr,pphi)
    vfa_i = uut + (1.0-pphigrid)*aalpha*iiota*AAIDS + iiota*(1.0-((1.0-pphigrid)*aalpha))*( &
        (1.0-eepsilon)*(eeta*MATMUL(wfun_b_apr,bprob) + (1.0-eeta)*MATMUL(wfun_i_apr,bprob)) + &
        eepsilon*(pphigrid*(eeta*MATMUL(wfun_b_apr,bprob) + (1.0-eeta)*MATMUL(wfun_i_apr,bprob)) + &
        (1.0-pphigrid)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:)))))
    vfa_i(1) = uutpos + iiota*aalphapr*AAIDS + iiota*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) !Replace for HIV+ aware

    !Protected Sex, BBETA, PPhi_P(pphipr,pphi)
    CALL PR_p(probhelp,pphigrid)
    pphitd = pphigrid - pphigrid*(1.0-nnu_p)*(1.0-ggamma_p)
    vfp_b  = uut + ppref_o + bbeta*(1.0-probhelp)*AAIDS + bbeta*(1.0-mmu)*probhelp*( &
        (1.0-eepsilon)*MATMUL(wfun_b_ppr,bprob) + eepsilon*( &
        pphigrid*MATMUL(wfun_b_ppr,bprob) + &
        (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:)))) + &
        bbeta*mmu*probhelp*(pphitd*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
        (1.0-pphitd)*DOT_PRODUCT(bprob,wfun_bl(1,:)))
    vfp_b(1) = uutpos + ppref_o + bbeta*aalphapr*AAIDS + &
        bbeta*(1.0-mmu)*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        bbeta*mmu*DOT_PRODUCT(bprob,wfun_bl(1,:))

    !Protected Sex, IIOTA, PPhi_P(pphipr,pphi)
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
    vfp_i(1) = uutpos + ppref_o + bbeta*aalphapr*AAIDS + &
        bbeta*(1.0-mmu)*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
        bbeta*mmu*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
        (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:)))

    !Unprotected Sex
    CALL PR_b(probhelp,pphigrid)
    pphitd = pphigrid - pphigrid*(1.0-nnu_b)*(1.0-ggamma_b)

    !$OMP PARALLEL DO
    do bCount = 1,bmap

        !Unprotected Sex, BBETA, PPhi_B(pphipr,pphi)
        vfb_b(:,bCount) = uut + bgrid(bCount) + bbeta*(1.0-probhelp)*AAIDS + &
            bbeta*(1.0-mmu)*probhelp*( &
            (1.0-eepsilon)*MATMUL(wfun_b_ppr,bprob) + eepsilon*( &
            pphigrid*MATMUL(wfun_b_ppr,bprob) + &
            (1.0-pphigrid)*DOT_PRODUCT(bprob,wfun_b(1,:)))) + &
            bbeta*mmu*probhelp*(pphitd*DOT_PRODUCT(bprob,wfun_bl(2,:)) + &
            (1.0-pphitd)*DOT_PRODUCT(bprob,wfun_bl(1,:)))
        vfb_b(1,bCount) = uutpos + bgrid(bCount) + bbeta*aalphapr*AAIDS + &
            bbeta*(1.0-mmu)*(1.0-aalphapr)*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
            bbeta*mmu*DOT_PRODUCT(bprob,wfun_bl(1,:))

        !Unprotected Sex, IIOTA, PPhi_B(pphipr,pphi)
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
        vfb_i(1,bCount) = uutpos + bgrid(bCount) + bbeta*aalphapr*AAIDS + &
            bbeta*(1.0-mmu)*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_b(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_i(1,:))) + &
            bbeta*mmu*(1.0-aalphapr)*(eeta*DOT_PRODUCT(bprob,wfun_bl(1,:)) + &
            (1.0-eeta)*DOT_PRODUCT(bprob,wfun_il(1,:)))

    end do
    !$OMP END PARALLEL DO

END SUBROUTINE vf_expost_s
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Solve Policy Functions using FOC and Corner Solution
!------------------------------------------------------------------------------
!--------------------------------
! Vectorize Value Function
!--------------------------------
SUBROUTINE vec1VF(vecVF,vfun,dim1,dim2)

    INTEGER(4) dim1,dim2,ic1,xc
    DOUBLE PRECISION vecVF(dim1*dim2),vfun(dim1)

    xc = 1
    do ic1 = 1,dim1
        do ic2 = 1,dim2
            vecVF(xc) = vfun(ic1)
            xc = xc + 1
        end do
    end do

END SUBROUTINE vec1VF
!--------------------------------
SUBROUTINE vec2VF(vecVF,vfun,dim1,dim2)

    INTEGER(4) dim1,dim2, ic1,ic2,xc
    DOUBLE PRECISION vecVF(dim1*dim2),vfun(dim1,dim2)

    xc = 1
    do ic1 = 1,dim1
        do ic2 = 1,dim2
            vecVF(xc) = vfun(ic1,ic2)
            xc = xc + 1
        end do
    end do

END SUBROUTINE vec2VF

!--------------------------------


SUBROUTINE solvepolicyfun(policy_p,policy_b,vfun_a,vfun_b,vfun_p,smap)


    USE mod_globalvar
    IMPLICIT NONE

    INTEGER smap, ix
    DOUBLE PRECISION, INTENT(in) :: vfun_a(smap),vfun_p(smap),vfun_b(smap)
    DOUBLE PRECISION, INTENT(out) :: policy_p(smap),policy_b(smap)
    DOUBLE PRECISION phelp1,phelp2,help1(smap),help2(smap)
    DOUBLE PRECISION ppower,llambda(smap),constant,aa,bb,cc,discrm

    ppower   = 1.0/kkappa
    constant = 1.0/(oomega_p*oomega_b*(kkappa+1.0))

    !Binding Constraint

    !Interior Solution and Binding Constraint
    do ix = 1,smap

        aa = constant
        bb = constant*(vfun_p(ix) - vfun_a(ix) + vfun_b(ix) - vfun_a(ix))
        cc = (vfun_p(ix)-vfun_a(ix))*(vfun_b(ix)-vfun_a(ix)) - 1.0
        discrm = (bb**2) - 4*aa*cc
        llambda(ix) = (bb - sqrt(discrm))/(2.0*constant)
        print *, llambda(ix)

        if (llambda(ix).GT.0.0) then

            help1(ix) = (vfun_p(ix) - vfun_a(ix) - llambda(ix))
            help2(ix) = (vfun_b(ix) - vfun_a(ix) - llambda(ix))
            help1(ix) = (help1(ix)/(oomega_p*(kkappa+1.0)))**ppower
            help2(ix) = (help2(ix)/(oomega_b*(kkappa+1.0)))**ppower
            help1(ix) = help1(ix)/(1.0+help1(ix))
            help2(ix) = help2(ix)/(1.0+help2(ix))


        else

            phelp1 = MAX(vfun_p(ix) - vfun_a(ix),0.0)
            phelp2 = (1.0/(oomega_p*(kkappa + 1.0)))*(phelp1**ppower)

            policy_p(ix) = phelp2/(1.0+phelp2)

            phelp1 = MAX(vfun_b(ix) - vfun_a(ix),0.0)
            phelp2 = (1.0/(oomega_b*(kkappa + 1.0)))*(phelp1**ppower)

            policy_b(ix) = phelp2/(1.0+phelp2)

        end if
    end do

            print *, help1
            print *, help2



END SUBROUTINE solvepolicyfun

!------------------------------------------------------------------------------
! Solve Ex-Ante Short-Term Value Functions
!------------------------------------------------------------------------------


