MODULE mod_globalvar

    IMPLICIT NONE

    ! All parameters are quarterly

    DOUBLE PRECISION aalpha   ! Pr(AIDS|HIV+,status unaware)
    DOUBLE PRECISION aalphapr ! Pr(AIDS|HIV+,status aware)

    DOUBLE PRECISION bbeta    ! Old discount factor net of mortality
    DOUBLE PRECISION iiota    ! Young discount factor net of mortarlity

    DOUBLE PRECISION ggamma_p ! Protected sex transmission rate
    DOUBLE PRECISION ggamma_b ! Unprotected sex transmission rate

    DOUBLE PRECISION ppref    ! Preference for protected sex
    DOUBLE PRECISION lpref    ! Preference for unprotected sex in LT relationship

    DOUBLE PRECISION nnu_b    ! Prevalence rate in unprotected sex market
    DOUBLE PRECISION nnu_p    ! Prevalence rate in protected sex market
    DOUBLE PRECISION nnu_l    ! Prevalence rate in long-term relationships

    DOUBLE PRECISION mmu      ! Long-term match rate
    DOUBLE PRECISION eeta     ! Ageing rate
    DOUBLE PRECISION xxi      ! Exogenous breakup rate
    DOUBLE PRECISION eepsilon ! Testing rate

    DOUBLE PRECISION rrho     ! Exogenous Birth Rate
    DOUBLE PRECISION ddelta   ! Exogenous death rate
    DOUBLE PRECISION ddelta1  ! Exogenous death rate for HIV+, status aware
    DOUBLE PRECISION ddelta2  ! Exogenous death rate for AIDS

    DOUBLE PRECISION oomega_p ! Cost scale parameter
    DOUBLE PRECISION kkappa_p ! Cost exponential parameter
    DOUBLE PRECISION oomega_b ! Cost scale parameter
    DOUBLE PRECISION kkappa_b ! Cost exponential parameter

    DOUBLE PRECISION AAIDS    ! Continuation value in AIDS stages

    DOUBLE PRECISION qincome  ! Quarterly income
    DOUBLE PRECISION mcost    ! Cost to HIV+ status-aware agents

END MODULE mod_globalvar
