program main

    USE mod_globalvar
    implicit none

    aalpha   = 0.025
    aalphapr = 0.01

    rrho    = 0.156
    ddelta  = 0.006
    ddelta1 = 0.025
    ddelta2 = 0.125

    bbeta   = (1-ddelta)*0.995
    iiota   = (1-ddelta)*0.8955
    eeta    = 0.06

    ggamma_b = 0.94
    ggamma_p = 0.98

    xxi      = 0.03
    qincome  = 500.0
    mcost    = 100.0
    ppref    = 7.5
    lpref    = 7.0
    AAIDS    = 5.0

    oomega_p = 0.5
    kkappa_p  = 0.2

    oomega_b = 0.5
    kkappa_b = 0.2

    mmu       = 0.3
    eepsilon  = 0.5

    call solvevfun

end program main
