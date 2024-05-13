MODULE sgl_subfuns

   USE spmatmul
   IMPLICIT NONE

   CONTAINS

   SUBROUTINE strong_rule (is_in_E_set, ga, pf, tlam, alsparse)
      IMPLICIT NONE
      INTEGER :: g, k
      INTEGER, DIMENSION (:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: pf
      DOUBLE PRECISION, INTENT(in) :: tlam, alsparse
      DOUBLE PRECISION :: z
      k = SIZE(is_in_E_set)
      z = tlam * (1 - alsparse)

      DO g = 1, k
         IF (is_in_E_set(g) == 1) CYCLE
         IF (ga(g) > pf(g) * z) is_in_E_set(g) = 1
      ENDDO
      RETURN
   END SUBROUTINE strong_rule

   SUBROUTINE kkt_check(is_in_E_set, violation, bn, ix, iy, vl, pf,pfl1,&
        lam1ma, bs, lama, ga, nvars)
      IMPLICIT NONE
      INTEGER :: g, startix, endix, nvars
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) :: bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: vl
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

      DO g = 1, bn
         IF (is_in_E_set(g) == 1) CYCLE
         startix = ix(g)
         endix = iy(g)
         ALLOCATE(s(bs(g)))
         s = vl(startix:endix)
         CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
         snorm = SQRT(DOT_PRODUCT(s,s))
         ga(g) = snorm
         IF(ga(g) > pf(g) * lam1ma) THEN
            is_in_E_set(g) = 1
            violation = 1
         ENDIF
         DEALLOCATE(s)
      ENDDO
      RETURN
   END SUBROUTINE kkt_check


   SUBROUTINE update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, pfl1, lam1ma, x,&
         isDifZero, nobs, r, gamg, maxDif,nvars, lb, ub)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, t_for_sg, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k

      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      oldb = b(startix:endix)
      s = MATMUL(r, x(:, startix:endix))/nobs
      s = s*t_for_sg + b(startix:endix)
      CALL softthresh(s, lama*t_for_sg*pfl1, bsg)
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lam1ma * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(lb, b(k)), ub)
         ENDDO
      ELSE
         b(startix:endix) = 0.0D0
      ENDIF
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF (ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         r = r - MATMUL(x(:,startix:endix), dd)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE update_step


   SUBROUTINE strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,pfl1,lam1ma,bs,&
         lama,ga,is_in_S_set,x,r,nobs,nvars,vl)
      IMPLICIT NONE
      INTEGER, INTENT(in)::nobs
      INTEGER, INTENT(in)::nvars
      DOUBLE PRECISION,INTENT(in):: x(nobs, nvars)
      DOUBLE PRECISION, INTENT(in):: r(nobs)
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) ::bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

      violation = 0
      DO g = 1, bn
         IF (is_in_S_set(g) == 1) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            s = MATMUL(r, x(:,startix:endix)) / nobs
            vl(startix:endix) = s
            CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ga(g) = snorm
            DEALLOCATE(s)
            IF (is_in_E_set(g) == 1) CYCLE
            IF (ga(g) > pf(g) * lam1ma) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE strong_kkt_check

   SUBROUTINE sp_update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, pfl1, lam1ma, x,&
         xidx, xcptr, nnz, isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: bsg, nobs, nvars, nnz
      INTEGER, INTENT(in) :: startix, endix
      DOUBLE PRECISION :: gamg
      DOUBLE PRECISION, INTENT(inout) :: maxDif
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb, s, dd
      DOUBLE PRECISION, DIMENSION (0:nvars), INTENT(inout) :: b
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: r
      DOUBLE PRECISION :: snorm, tea
      DOUBLE PRECISION, INTENT(in) :: lama, t_for_sg, pfg, lam1ma, lb, ub
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k

      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      s = 0.0D0
      oldb = b(startix:endix)
      ! print *, oldb

      CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
      s = s * t_for_sg / nobs + b(startix:endix)
      CALL softthresh(s, lama * t_for_sg * pfl1, bsg)
      snorm = SQRT(DOT_PRODUCT(s,s))
      tea = snorm - t_for_sg * lam1ma * pfg
      IF (tea > 0.0D0) THEN
         b(startix:endix) = s * tea / snorm
         DO k = startix, endix
            b(k) = MIN(MAX(b(k), lb), ub)
         ENDDO
      ELSE
         b(startix:endix) = 0.0D0
      ENDIF
      ALLOCATE(dd(bsg))
      dd = b(startix:endix) - oldb
      IF(ANY(dd .ne. 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         CALL ymspax(x, xidx, xcptr, nobs, nvars, nnz, dd, r, startix, endix, bsg)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE sp_update_step


   SUBROUTINE sp_strong_kkt_check(is_in_E_set,violation,bn,ix,iy,pf,pfl1,lam1ma,bs,&
         lama,ga,is_in_S_set,x,xidx,xcptr,nnz,r,nobs,nvars,vl)

      IMPLICIT NONE
      INTEGER, INTENT(in) :: nobs, nvars, nnz
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      INTEGER, INTENT(in) :: xidx(nnz)
      INTEGER, INTENT(in) :: xcptr(nvars + 1)
      DOUBLE PRECISION, INTENT(in):: r(nobs)
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: vl
      INTEGER :: g, startix, endix
      INTEGER, INTENT(in) :: bn
      INTEGER, INTENT(in) ::bs(bn)
      INTEGER, INTENT(in) :: ix(bn), iy(bn)
      INTEGER, DIMENSION(:), INTENT(inout) :: is_in_E_set
      INTEGER, DIMENSION(:), INTENT(in) :: is_in_S_set
      DOUBLE PRECISION, DIMENSION (:), INTENT(inout) :: ga
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s
      DOUBLE PRECISION :: snorm
      DOUBLE PRECISION, INTENT(in) :: pf(bn)
      DOUBLE PRECISION, INTENT(in) :: pfl1(nvars)
      INTEGER, INTENT(inout) :: violation
      DOUBLE PRECISION, INTENT(in) :: lam1ma, lama

      violation = 0
      DO g = 1, bn
         IF(is_in_S_set(g) == 1) THEN
            startix = ix(g)
            endix = iy(g)
            ALLOCATE(s(bs(g)))
            s = 0.0D0
            CALL spatx(x, xidx, xcptr, nobs, nvars, nnz, r, s, startix, endix)
            vl(startix:endix) = s / nobs
            CALL softthresh(s, lama * pfl1(startix:endix), bs(g))
            snorm = SQRT(dot_PRODUCT(s,s))
            ! print *, "kkt snorm = ", snorm
            ga(g) = snorm
            DEALLOCATE(s)
            IF(is_in_E_set(g) == 1) CYCLE
            IF(ga(g) > pf(g) * lam1ma) THEN
               is_in_E_set(g) = 1
               violation = 1
            ENDIF
         ENDIF
      ENDDO
      RETURN
   END SUBROUTINE sp_strong_kkt_check

END MODULE sgl_subfuns