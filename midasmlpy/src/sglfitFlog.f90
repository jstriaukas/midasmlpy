! ------------------------------------------------------------------
! sglfitF.f90: block coordinate descent for sg-LASSO logistic regression.
! ------------------------------------------------------------------
SUBROUTINE sglfit(gamma,ngroups,gindex,nobs,nvars,&
     x, y, pf,dfmax, pmax, nlam,flmin, ulam, eps, peps,&
     isd, intr, maxit, nalam, b0, beta, ibeta, nbeta, alam,&
     npass,jerr)

  
    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    REAL*8, INTENT(IN) :: gamma
    INTEGER, INTENT(IN) :: ngroups
    INTEGER*4, INTENT(IN) :: gindex(:)
    INTEGER, INTENT(IN) :: nobs
    INTEGER, INTENT(IN) :: nvars, dfmax, pmax, nlam
    INTEGER, INTENT(IN) :: isd, intr, maxit
    INTEGER, INTENT(OUT) :: nalam, npass, jerr
    INTEGER*4, INTENT(OUT) :: nbeta(nlam), ibeta(pmax)
    REAL*8, INTENT(IN) :: flmin, eps, peps
    REAL*8, INTENT(IN) :: x(:, :), y(:)
    REAL*8, INTENT(IN) :: pf(:)
    REAL*8, INTENT(OUT) :: b0(nlam), beta(pmax, nlam)
    REAL*8, INTENT(OUT) :: alam(nlam)
    REAL*8, INTENT(IN) :: ulam(:)
    REAL*8  ulam_(nlam)

    INTEGER j, l, nk, ierr
    INTEGER, DIMENSION(:), ALLOCATABLE :: ju
    REAL*8, DIMENSION(:), ALLOCATABLE :: xmean, xnorm, maj
    REAL*8, DIMENSION(:), ALLOCATABLE :: pf_
    REAL*8 maxlam, tmp

    nalam = 0
    b0 = 0.D0
    beta = 0.D0
    ibeta = 0
    nbeta = 0
    alam = 0.D0
    npass = 0
    jerr = 0
    ulam_ = ulam

   
    ALLOCATE(ju(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(xmean(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(xnorm(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    ALLOCATE(maj(1:nvars), STAT=ierr)
    jerr = jerr + ierr
    IF (jerr /= 0) RETURN
    CALL chkvars(nobs, nvars, x, ju)

    ALLOCATE(pf_(SIZE(pf)))

    IF (MAXVAL(pf) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF
    pf_ = MAX(0.0D0, pf)

    CALL standard(nobs, nvars, x, ju, isd, intr, xmean, xnorm, maj)
    ! -------------------- COMPUTE LAMBDA --------------------- !
    IF (ulam(1) .EQ. -1.0D0) THEN
        CALL maxlambda(nvars, nobs, x, y, gamma, gindex, ngroups, pf, maxlam)
        ulam_(1) = maxlam
        DO j = 2, nlam
            tmp = LOG(maxlam) + (LOG(maxlam*flmin) - LOG(maxlam)) * (j - 1) / (nlam - 1)
            ulam_(j) = EXP(tmp)
        END DO
    END IF

    ! ------------- CHOOSE MODEL BASED ON GAMMA ------------- !
    IF (gamma == 1.0D0) THEN
        CALL lassologisticfitpathF(nobs, nvars, x, y, ju, pf,&
               pmax, nlam, flmin, ulam_, eps, maxit, nalam, b0, beta,&
               nbeta, alam)
    ELSE
        CALL sgllogisticfitpathF(ngroups, gindex, nobs, nvars,&
               x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam_, eps, peps,& 
               nalam, b0, beta, ibeta, nbeta, alam, npass, intr)
    END IF
    IF (jerr > 0) RETURN
      
    DO l = 1, nalam
        nk = nbeta(l)
        IF (isd == 1) THEN
            DO j = 1, nk
                beta(j,l) = beta(j,l)/xnorm(ibeta(j))
            END DO
        END IF
        IF (intr == 1) THEN
            b0(l)=b0(l)-DOT_PRODUCT(beta(1:nk,l),xmean(ibeta(1:nk)))
        END IF
    END DO
      
    DEALLOCATE(ju,xmean,xnorm,maj)
    RETURN
END SUBROUTINE sglfit
      
      
SUBROUTINE sgllogisticfitpathF( ngroups, gindex, nobs, nvars, &
       x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, &
       nalam, b0, beta, ibeta, nbeta, alam, npass, intr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nobs, nvars, dfmax, pmax, nlam, maxit, ngroups
    INTEGER, INTENT(IN) :: gindex(ngroups)
    REAL*8, INTENT(IN) :: eps, peps, ulam(nlam), pf(nvars), x(nobs,nvars), y(nobs)
    INTEGER, INTENT(INOUT) :: ju(nvars), npass, nalam
    INTEGER, INTENT(INOUT) :: ibeta(pmax), nbeta(nlam)
    REAL*8, INTENT(INOUT) :: b0(nlam), beta(pmax, nlam), alam(nlam)
    INTEGER, INTENT(IN) :: intr
    REAL*8 :: flmin

    ! Local variables
    INTEGER l, m, k, ni, nj, gstart, gend
    REAL*8 d, newb, grad, hess, al, sumw
    REAL*8, ALLOCATABLE :: r(:), pvec(:), xb(:)

    ALLOCATE(r(nobs), pvec(nobs), xb(nobs))
    r = 0.0D0
    pvec = 0.0D0
    xb = 0.0D0

    ! Initialization
    nalam = 0
    b0 = 0.0D0
    beta = 0.0D0
    ibeta = 0
    nbeta = 0
    alam = 0.0D0
    npass = 0

    ! Begin path over lambda values
    DO l = 1, nlam
        al = ulam(l)
        ni = 0
        ! Initialize predictions and residuals
        xb = MATMUL(x, beta(:, l)) + b0(l)
        pvec = 1.0D0 / (1.0D0 + EXP(-xb))
        r = y - pvec  ! residuals for logistic regression

        ! Iterate over variables for coordinate descent
        DO k = 1, ngroups
            gstart = gindex(k) + 1
            IF (k > 1) gstart = gindex(k-1) + 1
            gend = gindex(k)

            ! Calculate the gradient and hessian for the group
            grad = 0.0D0
            hess = 0.0D0
            DO m = gstart, gend
                IF (ju(m) == 1) THEN
                    sumw = SUM(x(:, m) * (pvec * (1.0D0 - pvec)))
                    grad = SUM(x(:, m) * r) - pf(m) * beta(m, l)
                    hess = sumw + pf(m)
                    newb = beta(m, l) + grad / hess
                    d = newb - beta(m, l)
                    IF (ABS(d) > eps) THEN
                        beta(m, l) = newb
                        xb = xb + d * x(:, m)  ! Update predictions
                        pvec = 1.0D0 / (1.0D0 + EXP(-xb))  ! Update probabilities
                        r = y - pvec  ! Update residuals
                    END IF
                END IF
            END DO
        END DO
        ! Update intercept
        IF (intr == 1) THEN
            sumw = SUM(pvec * (1.0D0 - pvec))
            grad = SUM(r)
            hess = sumw
            newb = b0(l) + grad / hess
            d = newb - b0(l)
            IF (ABS(d) > eps) THEN
                b0(l) = newb
            END IF
        END IF

        alam(l) = al
        nbeta(l) = ni
        nalam = l
        npass = npass + 1
        IF (MAXVAL(ABS(r)) < peps) EXIT  ! Convergence check based on residuals
    END DO

    DEALLOCATE(r, pvec, xb)

END SUBROUTINE sgllogisticfitpathF


SUBROUTINE lassologisticfitpathF( nobs, nvars, x, y, ju, pf,&
      pmax, nlam, flmin, ulam, eps, maxit, nalam, b0,&
       beta, nbeta, alam)

    IMPLICIT NONE
    ! -------- INPUT VARIABLES -------- !
    INTEGER, INTENT(IN) :: nobs, nvars, pmax, nlam, maxit
    INTEGER, INTENT(INOUT) :: nalam
    INTEGER, INTENT(INOUT) :: ju(nvars), nbeta(nlam)
    REAL*8, INTENT(IN) :: eps, ulam(nlam), pf(nvars)
    REAL*8, INTENT(IN) :: x(nobs, nvars), y(nobs), maj(nvars)
    REAL*8, INTENT(INOUT) :: b0(nlam), beta(pmax, nlam), alam(nlam)
    REAL*8, INTENT(IN) :: flmin

    ! -------- LOCAL DECLARATIONS -------- !
    INTEGER l, k, ni, nj, mnl
    REAL*8 d, dif, u, v, al, grad, hess, bnorm
    REAL*8, DIMENSION(:), ALLOCATABLE :: b, oldbeta, r, pvec, xb
    REAL*8 sumw

    ! -------- ALLOCATE VARIABLES -------- !
    ALLOCATE(b(0:nvars), oldbeta(0:nvars), r(nobs), pvec(nobs), xb(nobs))

    b = 0.0D0
    oldbeta = 0.0D0
    r = y
    pvec = 0.0D0
    xb = 0.0D0

    mnl = MIN(nlam, 6)  ! Control loop variable for maximum iterations

    ! ---------------- INITIALIZATION ---------------- !
    DO l = 1, nlam
        al = ulam(l)
        b = 0.0D0
        oldbeta = 0.0D0
        xb = MATMUL(x, beta(:, l)) + b0(l)
        pvec = 1.0D0 / (1.0D0 + EXP(-xb))
        r = y - pvec  ! residuals for logistic regression
        nj = 0
        ni = 0

        ! ------------------ OUTER LOOP -------------------- !
        DO WHILE (ni < maxit AND nj < pmax)
            nj = nj + 1
            ! ---------------- GRADIENT AND HESSIAN CALCULATION ---------------- !
            FORALL (k = 1:nvars)
                sumw = SUM(pvec * (1.0 - pvec) * x(:, k)**2)
                grad = DOT_PRODUCT(x(:, k), r) - pf(k) * b(k)
                hess = sumw + pf(k)
                ! ------------------ COORDINATE DESCENT ---------------- !
                u = b(k) + grad / hess
                v = ABS(u) - al / hess
                IF (v > 0.0D0) THEN
                    tmp = SIGN(v, u)
                ELSE
                    tmp = 0.0D0
                END IF
                d = tmp - b(k)
                IF (ABS(d) > eps) THEN
                    b(k) = tmp
                    xb = xb + d * x(:, k)  ! Update predictions
                    pvec = 1.0D0 / (1.0D0 + EXP(-xb))  ! Update probabilities
                    r = y - pvec  ! Update residuals
                    ni = ni + 1  ! Update iterations count
                    IF (ni > pmax) EXIT
                END IF
            END FORALL
        END DO
        ! ----------- FINAL UPDATE & SAVE RESULTS ------------ !
        beta(1:nvars,l) = b
        nbeta(l) = nvars
        b0(l) = b(0)
        alam(l) = al
        nalam = l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
        IF (ni >= maxit) THEN
            jerr = -1000 - l
            EXIT
        END IF
    END DO

    DEALLOCATE(b, oldbeta, r, pvec, xb)

END SUBROUTINE lassologisticfitpathF

      SUBROUTINE prox_sgl(gstart, gend, nvars, nobs, x, r,&
           b, al, gamma, pf, peps, gw)
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER gstart, gend, nvars, nobs
      Real*8 x(nobs,nvars), r(nobs), b(nvars), al
      Real*8 gamma, pf(nvars), peps, gw
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER :: g
      Real*8 u, scl, tmp, maxg, normg, d 
      Real*8 bold(nvars), vg
      Real*8, Parameter :: big = 9.9D30
      ! -------- BEGIN PROGRAM -------- !
      Do
       bold(gstart:gend) = b(gstart:gend)
        !--------- LASSO PART ----------!
        Do g = gstart, gend
          u = b(g) + DOT_PRODUCT(x(:,g),r)/nobs
          !S(.) map
          tmp=max(u-al*gamma*pf(g),0.0D0)-max(-u-al*gamma*pf(g),0.0D0)
          b(g) = tmp
        End Do
         !--------- g-LASSO PART ----------!
         ! L2 norm of b_g
         normg = NORM2(b(gstart:gend))
         ! Initialize storage vars
         maxg = 0.0D0
         vg = gw * al * (1.0D0-gamma)/normg
         If (normg .EQ. 0.0D0) Then
           vg = big
         End If
         Do g = gstart, gend
           scl = 1.0D0 - pf(g) * vg
           scl = MAX(scl, 0.0D0)
           ! l_2,1 norm map
           tmp = scl*b(g)
           d = tmp - bold(g)
           r = r - x(:,g)*d
           maxg = MAX(maxg, ABS(d))
           b(g) = tmp
         End Do
      !--------- CHECK CONVERGENCE ----------!
        If (maxg < peps) Exit
      !Exit
      End Do
      End SUBROUTINE prox_sgl

      SUBROUTINE standard(nobs, nvars, x, ju, isd, intr,&
           xmean, xnorm, maj)     
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER nobs, nvars, isd, intr, ju(nvars)
      Real*8 x(nobs, nvars), xmean(nvars)
      Real*8 xnorm(nvars), maj(nvars)
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER j
      Real*8 xmsq, xvar
      ! -------- BEGIN PROGRAM -------- !
      IF (intr == 0) THEN
      DO j = 1, nvars
      IF (ju(j) == 1) THEN
      xmean(j) = 0.0D0
      maj(j) = DOT_PRODUCT(x(:,j),x(:,j))/nobs
      IF (isd == 1) THEN
      xmsq = (SUM(x(:,j))/nobs)**2
      xvar = maj(j) - xmsq
      xnorm(j) = SQRT(xvar)
      x(:,j) = x(:,j)/xnorm(j)
      maj(j) = 1.0D0 + xmsq/xvar
      END IF
      END IF
      END DO
      ELSE                   
      DO j = 1, nvars                                  
      IF (ju(j) == 1) THEN                         
      xmean(j) = SUM(x(:,j))/nobs  ! MEAN                        
      x(:,j) = x(:,j) - xmean(j)    
      maj(j) = DOT_PRODUCT(x(:,j),x(:,j))/nobs                                             
      IF (isd == 1) THEN
      xnorm(j) = SQRT(maj(j))  ! STANDARD DEVIATION              
      x(:,j) = x(:,j)/xnorm(j)
      maj(j) = 1.0D0
      END IF                                                        
      END IF                                     
      END DO  
      END IF                           
      END SUBROUTINE standard
      
      ! ---------------------------------------------------------- !
      SUBROUTINE chkvars(nobs, nvars, x, ju)
      IMPLICIT NONE
      ! -------- INPUT VARIABLES -------- !
      INTEGER :: nobs,nvars,ju(nvars)
      Real*8 :: x(nobs,nvars)
      ! -------- LOCAL DECLARATIONS -------- !
      INTEGER :: i,j
      Real*8 :: t
      ! -------- BEGIN PROGRAM -------- ! 
      DO j = 1, nvars
      ju(j) = 0
      t = x(1,j)
      DO i = 2, nobs
      IF (x(i,j) /= t) THEN
      ju(j) = 1
      EXIT
      END IF
      END DO
      END DO
      END SUBROUTINE chkvars





    SUBROUTINE maxlambda(nvars, nobs, x, y, gamma, gindex, ngroups, pf, maxlam)

        IMPLICIT NONE
        ! -------- INPUT VARIABLES -------- !
        INTEGER :: nvars, nobs, ngroups, gindex(ngroups)
        DOUBLE PRECISION :: x(nobs,nvars), y(nobs), pf(nvars), gamma, maxlam
        ! -------- LOCAL DECLARATIONS -------- !
        INTEGER :: k, c, nzvars
        INTEGER :: gstart, gend, gs, gj
        DOUBLE PRECISION :: gw, xy(nvars), r(nobs)
        DOUBLE PRECISION :: wmaxg(ngroups), lb, rb

        ! -------- BEGIN PROGRAM -------- !
        c = 0
        r = y
        nzvars = 0
        DO k = 1, nvars
            IF (pf(k) .EQ. 0.0D0) THEN
            nzvars = nzvars + 1
            END IF
        END DO
        IF (nzvars .NE. 0) THEN
            !CALL rnz(nvars, nobs, nzvars, y, x, r, pf)
        END IF
        xy = MATMUL(TRANSPOSE(x),r)/nobs


        IF (gamma .EQ. 1.0D0) THEN
            maxlam = MAXVAL(ABS(xy))
        ELSE
            DO k = 1, ngroups
                gend = gindex(k)
                IF (k == 1) THEN
                    gstart = 1
                ELSE
                    gstart = gindex(k-1) + 1
                END IF
                gs = gend - gstart + 1
                gw = 0.0D0
                DO gj = gstart, gend
                    gw = gw + pf(gj)
                END DO
                gw = SQRT(gw)
                IF (gw == 0.0D0) THEN
                    wmaxg(k) = 0.0D0
                ELSE
                    IF (gamma .EQ. 0.0D0) THEN
                        !rb = NORM2(xy(gstart:gend))
                        rb = SQRT(DOT_PRODUCT(xy(gstart:gend), xy(gstart:gend)))
                        wmaxg(k) = rb/gw
                    ELSE
                        lb = 0.0D0
                        rb = MAXVAL(ABS(xy(gstart:gend)))/gamma
                        CALL solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy, nvars)
                        wmaxg(k) = rb
                    END IF
                END IF
            END DO
            maxlam = MAXVAL(wmaxg)
        END IF
        !--- ADD SMALL NUMBER TO ENSURE b = 0 @ maxlam (DUE TO NUMERICAL IMPRESSION)
        maxlam = maxlam + 1E-5

    END SUBROUTINE maxlambda


    SUBROUTINE  solvewmaxg(gstart, gend, gamma, lb, rb, gw, pf, xy, nvars)
        IMPLICIT NONE
        ! -------- INPUT VARIABLES -------- !
        INTEGER :: gstart, gend, nvars
        DOUBLE PRECISION :: gamma, lb, rb, gw, pf(nvars), xy(nvars)
        ! -------- LOCAL DECLARATIONS -------- !
        INTEGER :: stopflag, indexi
        DOUBLE PRECISION ::  tol = 1E-13, mp, fl, fm, fr, tmpl, tmpm, tmpr

        stopflag = 0
        DO
            mp = 0.5 * (lb + rb)
            fl = 0.0D0
            fm = 0.0D0
            fr = 0.0D0
            tmpl = 0.0D0
            tmpm = 0.0D0
            tmpr = 0.0D0
            DO indexi =  gstart, gend
                tmpl = ABS(xy(indexi)) - gamma * lb * pf(indexi)
                tmpm = ABS(xy(indexi)) - gamma * mp * pf(indexi)
                tmpr = ABS(xy(indexi)) - gamma * rb * pf(indexi)
                IF (tmpl > 0.0D0) THEN
                    fl = fl + tmpl * tmpl
                END IF
                IF (tmpm > 0.0D0) THEN
                    fm = fm + tmpm * tmpm
                END IF
                IF (tmpr > 0.0D0) THEN
                    fr = fr + tmpr * tmpr
                END IF
            END DO
            fl = fl - (1.0D0 - gamma) * (1.0D0 - gamma) * lb * lb * gw * gw
            fm = fm - (1.0D0 - gamma) * (1.0D0 - gamma) * mp * mp * gw * gw
            fr = fr - (1.0D0 - gamma) * (1.0D0 - gamma) * rb * rb * gw * gw
            IF (fl * fm < 0.0D0) THEN
                IF (ABS(lb - mp) > tol) THEN
                    rb = mp
                ELSE
                    stopflag = 1
                END IF
            ELSE
                IF (fm * fr < 0.0D0) THEN
                    IF (ABS(mp - rb) > tol) THEN
                        lb = mp
                    ELSE
                        stopflag = 1
                    END IF
                ELSE
                    stopflag = 1
                END IF
            END IF
           IF (stopflag .EQ. 1) EXIT
        END DO
        rb = mp

    END SUBROUTINE solvewmaxg
