!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Algorithm to Solve for Value Functions
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

SUBROUTINE solvevfun

    USE mod_globalvar
    IMPLICIT NONE



END SUBROUTINE solvevfun

!------------------------------------------------------------------------------
! Solve Ex-Post Long-Term Value Functions
!------------------------------------------------------------------------------
! d = bbeta, z = l
SUBROUTINE vf_expost_bl(vfun,wfun,smap)
!Input W_s^bbeta(0,x)

    USE mod_globalvar
    IMPLICIT NONE

    INTEGER smap
    DOUBLE PRECISION, INTENT(in) :: wfun(smap)
    DOUBLE PRECISION, INTENT(out) :: vfun(2,2,smap)
    DOUBLE PRECISION

    vfun(1,1,:) =


END SUBROUTINE vf_expost_bl
! d = iiota, z = l
SUBROUTINE vf_expost_il(vfun,wfunb,wfuni,smap)
!Input W_s^bbeta(0,x), W_s^iiota(0,x)

END SUBROUTINE vf_expost_il

!------------------------------------------------------------------------------
! Solve Policy Functions using FOC and Corner Solution
!------------------------------------------------------------------------------
SUBROUTINE solvepolicyfun_p(policy_p,vfun_a,vfun_p,smap)

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
