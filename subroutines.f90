!------------------------------------------------------------------------------
! Fill Vector with Zeros
!------------------------------------------------------------------------------
SUBROUTINE zerov(x,n)

    implicit none

    INTEGER(4) n,zcounter
    REAL(8) x(n)

    zeroloop: do zcounter = 1,n
        x(zcounter) = 0.0
    end do zeroloop

END SUBROUTINE zerov

!------------------------------------------------------------------------------
! Utility Function
!------------------------------------------------------------------------------
SUBROUTINE util(u,x)

    implicit none

    REAL(8) x,u

    u = log(x)

END SUBROUTINE util


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

