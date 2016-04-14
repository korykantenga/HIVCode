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
