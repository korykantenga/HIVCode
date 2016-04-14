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
    !vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !Ex-Ante LT
    !vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !Ex-Post ST

    !Ex-Ante ST + Policy Function

END SUBROUTINE solvevfun

!------------------------------------------------------------------------------
! Solve Ex-Post Long-Term Value Functions
!------------------------------------------------------------------------------
! d = bbeta or iiota, z = l
SUBROUTINE vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

    !Input  : W_s^bbeta(pphi,x), W_s^iiota(pphi,x)
    !Output : V_l^bbeta(pphi,varpphi,x)
    !dim(pphi) = pmap, dim(x) = bmap

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER bmap,pmap
    DOUBLE PRECISION, INTENT(in) :: wfun_b(pmap,bmap),wfun_i(pmap,bmap)
    DOUBLE PRECISION, INTENT(out) :: vfun_bl(2,2,bmap),vfun_il(2,2,bmap)
    DOUBLE PRECISION uut
    DOUBLE PRECISION vhelp1(bmap),chelp1,vhelp2(bmap),chelp2
    DOUBLE PRECISION vhelp3(bmap),chelp3,chelp4


    !V_l^bbeta(0,0,x)
    CALL util(uut,qincome-mcost)
    vhelp1 = aalphapr*AAIDS + (1.0-xxi)*(aalphapr*(1.0-aalphapr)*wfun_b(1,:)) + &
        xxi*((1.0-aalphapr)*wfun_b(1,:))
    chelp1 = 1.0-bbeta*(1.0-xxi)*((1.0-aalphapr)**2)

    vfun_bl(1,1,:) = (uut + ppref_o + bbeta*vhelp1)/(chelp1)

    !V_l^bbeta(1,0,x)
    CALL util(uut,qincome)
    vhelp2 = (1.0-ggamma_p)*aalpha*AAIDS + (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        aalphapr*ggamma_p*wfun_b(pmap,:) + aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*wfun_b(1,:)) + &
        xxi*((1.0-aalpha)*(1.0-ggamma_p)*wfun_b(1,:) + ggamma_p*wfun_b(pmap,:))
    chelp2 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(2,1,:) = (uut + ppref_o + bbeta*vhelp2)/(chelp2)

    !V_l^bbeta(0,1,x)
    CALL util(uut,qincome-mcost)
    vhelp3 = aalphapr*AAIDS + (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        aalpha*(1.0-aalphapr)*(1.0-ggamma_p)*wfun_b(1,:)) + &
        xxi*(1.0-aalphapr)*wfun_b(1,:)
    chelp3 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(1,2,:) = (uut + ppref_o + bbeta*vhelp3)/(chelp3)

    !V_l^bbeta(1,1,x)
    CALL util(uut,qincome)
    chelp4 = 1.0 - bbeta*(1.0-xxi)

    vfun_bl(2,2,:) = (uut + lpref_o + xxi*bbeta*wfun_b(pmap,:))/chelp4

    !---------------------------------------------------------------------

    !V_l^iiota(0,0,x)
    CALL util(uut,qincome-mcost)
    vhelp1 = aalphapr*AAIDS + &
        (1.0-xxi)*(aalphapr*(1.0-aalphapr)*(eeta*wfun_b(1,:)+(1.0-eeta)*wfun_i(1,:)) + &
        ((1.0-aalphapr)**2)*vfun_bl(1,1,:)) + &
        xxi*((1.0-aalphapr)*(eeta*wfun_b(1,:)+(1.0-eeta)*wfun_i(1,:)))
    chelp1 = 1.0-iiota*(1.0-xxi)*(1.0-eeta)*((1.0-aalphapr)**2)

    vfun_il(1,1,:) = (uut + ppref_y + iiota*vhelp1)/(chelp1)

    !V_l^iiota(1,0,x)
    CALL util(uut,qincome)
    vhelp2 = (1.0-ggamma_p)*aalpha*AAIDS + (1.0-xxi)*eeta*(1.0-aalphapr)*ggamma_p*vfun_bl(2,1,:) + &
        (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*(eeta*vfun_bl(1,1,:) + &
        (1.0-eeta)*vfun_il(1,1,:)) + &
        aalphapr*ggamma_p*(eeta*wfun_b(pmap,:)+(1.0-eeta)*wfun_i(pmap,:)) + &
        aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*(wfun_b(1,:)*eeta+(1.0-eeta)*wfun_i(1,:))) + &
        xxi*((1.0-aalpha)*(1.0-ggamma_p)*(eeta*wfun_b(1,:)+(1.0-eeta)*wfun_i(1,:)) + &
        ggamma_p*(eeta*wfun_b(pmap,:)+(1.0-eeta)*wfun_i(pmap,:)))
    chelp2 = 1.0 - iiota*(1.0-xxi)*(1.0-eeta)*((1.0-aalphapr)*ggamma_p)

    vfun_il(2,1,:) = (uut + ppref_y + iiota*vhelp2)/(chelp2)

    !V_l^iiota(0,1,x)
    CALL util(uut,qincome-mcost)
    vhelp3 = aalphapr*AAIDS + (1.0-xxi)*((1.0-aalphapr)*ggamma_p*eeta*vfun_bl(1,2,:) + &
        (1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*(eeta*vfun_bl(1,1,:)+ &
        (1.0-eeta)*vfun_il(1,1,:)) + &
        aalpha*(1.0-aalphapr)*(1.0-ggamma_p)*(eeta*wfun_b(1,:)+(1.0-eeta)*wfun_i(1,:))) + &
        xxi*(1.0-aalphapr)*(eeta*wfun_b(1,:)+(1.0-eeta)*wfun_i(1,:))
    chelp3 = 1.0 - iiota*(1.0-xxi)*(1.0-eeta)*(1.0-aalphapr)*ggamma_p

    vfun_il(1,2,:) = (uut + ppref_y + iiota*vhelp3)/(chelp3)

    !V_l^iiota(1,1,x)
    CALL util(uut,qincome)
    chelp4 = 1.0 - iiota*(1.0-eeta)*(1.0-xxi)

    vfun_il(2,2,:) = (uut + lpref_y + iiota*eeta*(1.0-xxi)*vfun_bl(2,2,:) + &
        xxi*iiota*(eeta*wfun_b(pmap,:)+(1.0-eeta)*wfun_i(pmap,:)))/chelp4

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

    wfun_bl(1,:) = nnu_l*MAXVAL(vhelp1,2) + (1-nnu_l)*MAXVAL(vhelp2,2) + wfun_b(1,:)

    vhelp1(:,1) = vfun_bl(2,2,:)-wfun_b(pmap,:) !HIV-,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_bl(2,1,:)-wfun_b(pmap,:) !HIV-,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_bl(2,:) = nnu_l*MAXVAL(vhelp1,2) + (1-nnu_l)*MAXVAL(vhelp2,2) + wfun_b(pmap,:)

    !-----------------------------------------------------------------------------------
    !W_l^iiota(pphi,x)

    vhelp1(:,1) = vfun_il(1,2,:)-wfun_i(1,:) !HIV+,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_il(1,1,:)-wfun_i(1,:) !HIV+,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_il(1,:) = nnu_l*MAXVAL(vhelp1,2) + (1-nnu_l)*MAXVAL(vhelp2,2) + wfun_i(1,:)

    vhelp1(:,1) = vfun_il(2,2,:)-wfun_i(pmap,:) !HIV-,HIV- match
    CALL zerov(vhelp1(:,2),bmap)

    vhelp2(:,1) = vfun_il(2,1,:)-wfun_i(pmap,:) !HIV-,HIV+ match
    CALL zerov(vhelp2(:,2),bmap)

    wfun_il(2,:) = nnu_l*MAXVAL(vhelp1,2) + (1-nnu_l)*MAXVAL(vhelp2,2) + wfun_i(pmap,:)


END SUBROUTINE vf_exante_l

!------------------------------------------------------------------------------
! Solve Ex-Post Short-Term Value Functions
!------------------------------------------------------------------------------

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
