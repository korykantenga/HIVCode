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

    DOUBLE PRECISION vecfun_b(50*50),vecfun_p(50*50),vecfun_a(50*50)
    DOUBLE PRECISION policyb_bbeta(50*50),policyp_bbeta(50*50)


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
    ppref_o  = 8.0
    ppref_y  = 6.5
    lpref_o  = 8.0
    lpref_y  = 7.0

    AAIDS    = 5.0

    oomega_p = 0.49
    oomega_b = 0.5
    kkappa   = 0.2

    mmu       = 0.3
    eepsilon  = 0.5

    !Initialise W_s^bbeta(pphi,x),W_s^iiota(pphi,x)
    do iCount1 = 1,50
        do iCount2 = 1,50
            wfun_b(iCount1,iCount2) = 0.0
            wfun_i(iCount1,iCount2) = 0.0
        end do
    end do

    !Updating Value Functions
    call vf_expost_l(vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap,bprob)
    call vf_exante_l(wfun_bl,wfun_il,vfun_bl,vfun_il,wfun_b,wfun_i,pmap,bmap)
    call vf_expost_s(vfa_b,vfb_b,vfp_b,vfa_i,vfb_i,vfp_i, &
        wfun_bl,wfun_il,wfun_b,wfun_i,pmap,bmap,bprob,pphigrid,bgrid)

    call vec2VF(vecfun_b,vfb_b,pmap,bmap)
    call vec1VF(vecfun_a,vfa_b,pmap,bmap)
    call vec1VF(vecfun_p,vfp_b,pmap,bmap)

    call solvepolicyfun(policyp_bbeta,policyb_bbeta,vecfun_a,vecfun_b,vecfun_p,pmap*bmap)

end program main
