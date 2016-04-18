program main

    USE mod_globalvar
    implicit none

    INTEGER bmap,pmap
    DOUBLE PRECISION bprob(50)
    DOUBLE PRECISION vfun_bl(2,2,50),vfun_il(2,2,50),wfun_b(50,50),wfun_i(50,50)
    DOUBLE PRECISION wfun_bl(2,50),wfun_il(2,50)
    DOUBLE PRECISION vfa_b(50),vfa_i(50)
    DOUBLE PRECISION vfb_b(50,50),vfb_i(50,50)
    DOUBLE PRECISION vfp_b(50),vfp_i(50)
    DOUBLE PRECISION pphigrid(50),bgrid(50)
    INTEGER iCount1,iCount2

    bmap = 50
    pmap = 50

    CALL phiGrid(pphigrid,0.0,50)
    CALL bprefGrid(bgrid,bprob,50)

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
    ppref_o  = 9.0
    ppref_y  = 6.5
    lpref_o  = 8.0
    lpref_y  = 7.0

    AAIDS    = 5.0

    oomega_p = 0.5
    kkappa_p  = 0.2

    oomega_b = 0.5
    kkappa_b = 0.2

    mmu       = 0.3
    eepsilon  = 0.5

    do iCount1 = 1,50
        do iCount2 = 1,50
            wfun_b(iCount1,iCount2) = 0.0
            wfun_i(iCount1,iCount2) = 0.0
        end do
    end do

    call vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,50,50,bprob)

    call vf_exante_l(wfun_bl,wfun_il,vfun_bl,vfun_il,wfun_b,wfun_i,50,50)

    call vf_expost_s(vfa_b,vfb_b,vfp_b,vfa_i,vfb_i,vfp_i, &
        wfun_bl,wfun_il,wfun_b,wfun_i,pmap,bmap,bprob,pphigrid,bgrid)




end program main
