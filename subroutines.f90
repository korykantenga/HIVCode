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
! Construct PPHI' GRID
!------------------------------------------------------------------------------
SUBROUTINE PPhiPRIME(pphia,pphib,pphip,pphigrid,pmap)

    USE mod_globalvar

    IMPLICIT NONE


    INTEGER(4) pmap, iCount
    DOUBLE PRECISION, INTENT(in) :: pphigrid(1:pmap)
    DOUBLE PRECISION, INTENT(out):: pphia(1:pmap),pphib(1:pmap),pphip(1:pmap)

    !$OMP PARALLEL DO
    do iCount = 1,pmap
        call PPhi_A(pphia(iCount),pphigrid(iCount))
        call PPhi_B(pphib(iCount),pphigrid(iCount))
        call PPhi_P(pphip(iCount),pphigrid(iCount))
    end do
    !$OMP END PARALLEL DO

END SUBROUTINE PPhiPRIME
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
! Transition Related Probabilities
!------------------------------------------------------------------------------
SUBROUTINE PR_b(prob,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*((1.0-pphi)+pphi*(1.0-nnu_b)*(1.0-ggamma_b))

END SUBROUTINE PR_b
!------------------------------------------------------------------------------
SUBROUTINE PR_p(prob,pphi)
!Prob(No Aids Symptoms|Protected Sex)
    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*((1.0-pphi) + pphi*(1.0-nnu_p)*(1.0-ggamma_p))

END SUBROUTINE PR_p
!------------------------------------------------------------------------------
SUBROUTINE PR_a(prob,pphi)

    USE mod_globalvar
    IMPLICIT NONE

    DOUBLE PRECISION pphi,prob

    prob = 1.0 - aalpha*(1.0-pphi)

END SUBROUTINE PR_a
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Linear Interpolation of Value Function
!------------------------------------------------------------------------------
SUBROUTINE LinInterp (n,x,y,ni,xi,yi)
    !this does linear interpolation of (x,y) at points xi
    !requires x to be sorted in ascending order
    !extrapolates out of range
    IMPLICIT NONE
    INTEGER, INTENT(in)        :: n,ni
    real(8), INTENT(in)    :: x(:),y(:),xi(:)
    real(8), INTENT(out)    :: yi(:)
    real(8), DIMENSION(ni)    :: xL,xH,yL,yH
    INTEGER                    :: i,locL(ni)


    DO i = 1,ni

        LocL(i) = MAXLOC(x,1,MASK=xi(i)>x)

        IF (xi(i)<=x(1)) THEN
            LocL(i) = 1
        END IF

        IF (LocL(i)>=n) THEN
            LocL(i) = n-1
        END IF


        xL(i) = x(locL(i))
        xH(i) = x(locL(i)+1)
        yL(i) = y(locL(i))
        yH(i) = y(locL(i)+1)

        yi(i) = yL(i) + (xi(i)-xL(i))*((yH(i)-yL(i))/(xH(i)-xL(i)))

    END DO

END SUBROUTINE LinInterp


    SUBROUTINE BiLinInterp1(nx,x,ny,y,f,xi,yi,fi)
        !this does linear interpolation of f(x,y) at points (xi,yi)
        !requires x to be sorted in ascending order
        !
        !extrapolates out of range
        IMPLICIT NONE
        INTEGER, INTENT(in)        :: nx,ny
        real(8), INTENT(in)    :: x(:),y(:),f(:,:),xi,yi
        real(8), INTENT(out)    :: fi
        real(8)                :: xL,xH,yL,yH,fLL,fHH,fLH,fHL,dxdy
        INTEGER                    :: xlocL,ylocL

        xlocL = MAXLOC(x,1,MASK=xi>x)
        ylocL = MAXLOC(y,1,MASK=yi>y)

        IF (xi<=x(1)) THEN
            xlocL = 1
        END IF

        IF (xLocL>=nx) THEN
            xLocL = nx-1
        END IF

        IF (yi<=y(1)) THEN
            ylocL = 1
        END IF

        IF (yLocL>=ny) THEN
            yLocL = ny-1
        END IF

        xL  = x(xlocL)
        xH  = x(xlocL +1)
        yL  = y(ylocL)
        yH  = y(ylocL +1)
        fLL = f(xlocL,ylocL)
        fLH = f(xlocL,ylocL+1)
        fHL = f(xlocL+1,ylocL)
        fHH = f(xlocL+1,ylocL+1)

        dxdy = (xH-xL)*(yH-yL)
        fi = fLL*(xH-xi)*(yH-yi)/(dxdy) + fHL*(xi-xL)*(yH-yi)/(dxdy) + fLH*(xH-xi)*(yi-yL)/(dxdy) + fHH*(xi-xL)*(yi-yL)/(dxdy)


    END SUBROUTINE BiLinInterp1

