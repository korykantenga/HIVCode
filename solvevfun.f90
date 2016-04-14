!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Algorithm to Solve for Value Functions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE solvevfun

    USE mod_globalvar
    IMPLICIT NONE

    !Solve for Ex-Post Long-Term Value Functions given W_s^bbeta(pphi,x), W_s^iiota(pphi,x)
    !vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)

END SUBROUTINE solvevfun

!------------------------------------------------------------------------------
! Solve Ex-Post Long-Term Value Functions
!------------------------------------------------------------------------------
! d = bbeta or iiota, z = l
SUBROUTINE vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)
    !Input: W_s^bbeta(pphi,x), W_s^iiota(pphi,x)
    !Output: V_l^bbeta(pphi,varpphi,x)
    !dim(pphi)=pmap,dim(x)=bmap

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

    vfun_bl(1,1,:) = (uut + ppref + bbeta*vhelp1)/(chelp1)

    !V_l^bbeta(1,0,x)
    CALL util(uut,qincome)
    vhelp2 = (1.0-ggamma_p)*aalpha*AAIDS + (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        aalphapr*ggamma_p*wfun_b(pmap,:) + aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*wfun_b(1,:)) + &
        xxi*((1.0-aalpha)*(1.0-ggamma_p)*wfun_b(1,:) + ggamma_p*wfun_b(pmap,:))
    chelp2 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(2,1,:) = (uut + ppref + bbeta*vhelp2)/(chelp2)

    !V_l^bbeta(0,1,x)
    CALL util(uut,qincome-mcost)
    vhelp3 = aalphapr*AAIDS + (1.0-xxi)*((1.0-ggamma_p)*(1.0-aalpha)*(1.0-aalphapr)*vfun_bl(1,1,:) + &
        (1.0-aalphapr)*(aalpha*(1.0-ggamma_p))*wfun_b(pmap,:) + aalphapr*(1.0-aalpha)*(1.0-ggamma_p)*wfun_b(1,:)) + &
        xxi*(1.0-aalphapr)*wfun_b(pmap,:)
    chelp3 = 1.0 - bbeta*(1.0-xxi)*((1.0-aalphapr)*ggamma_p)

    vfun_bl(1,2,:) = (uut + ppref + bbeta*vhelp3)/(chelp3)

    !V_l^bbeta(1,1,x)
    CALL util(uut,qincome)
    chelp4 = 1.0 - bbeta*(1.0-xxi)

    vfun_bl(2,2,:) = (uut + lpref + xxi*bbeta*wfun_b(pmap,:))/chelp4

END SUBROUTINE vf_expost_l

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

!------------------------------------------------------------------------------
! Abstinence Updating Rule
!------------------------------------------------------------------------------
SUBROUTINE PPhi_A(pphipr,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,pphipr

    pphipr = pphi*(1.0/(pphi+(1.0-pphi)*(1.0-aalpha)))

END SUBROUTINE PPhi_A
!------------------------------------------------------------------------------
! Unprotected Sex Updating Rule
!------------------------------------------------------------------------------
SUBROUTINE PPhi_B(pphipr,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,pphipr,h1,h2,h3

    h1 = nnu_b + ggamma_b*(1.0-nnu_b)
    h2 = pphi*(nnu_b + ggamma_b*(1.0-nnu_b) + (1.0-aalpha)*(1.0-ggamma_b)*(1.0-nnu_b))
    h3 = (1.0-aalpha)*(1.0-pphi)

    pphipr = pphi*(h1/(h2+h3))

END SUBROUTINE PPhi_B
!------------------------------------------------------------------------------
! Protected Sex Updating Rule
!------------------------------------------------------------------------------
SUBROUTINE PPhi_P(pphipr,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,pphipr,h1,h2,h3

    h1 = nnu_p + ggamma_p*(1.0-nnu_p)
    h2 = pphi*(nnu_p + ggamma_p*(1.0-nnu_p) + (1.0-aalpha)*(1.0-ggamma_p)*(1.0-nnu_p))
    h3 = (1.0-aalpha)*(1.0-pphi)

    pphipr = pphi*(h1/(h2+h3))

END SUBROUTINE PPhi_P
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Transition Related Probabilities
!------------------------------------------------------------------------------
SUBROUTINE PR_b(prob,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*((1.0-pphi)+(1.0-nnu_b)*(1.0-ggamma_b))

END SUBROUTINE PR_b
!------------------------------------------------------------------------------
SUBROUTINE PR_p(prob,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*((1.0-pphi)+(1.0-nnu_p)*(1.0-ggamma_p))

END SUBROUTINE PR_p
!------------------------------------------------------------------------------
SUBROUTINE PR_a(prob,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*(1.0-pphi)

END SUBROUTINE PR_a
!------------------------------------------------------------------------------
