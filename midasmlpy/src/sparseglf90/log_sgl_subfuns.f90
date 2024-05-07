MODULE log_sgl_subfuns

  USE spmatmul
  IMPLICIT NONE

CONTAINS


  SUBROUTINE log_update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, pfl1, lam1ma,&
       x, y, isDifZero, nobs, r, gamg, maxDif, nvars, lb, ub)
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
      DOUBLE PRECISION, DIMENSION (:), INTENT(in) :: x(nobs,nvars)
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, INTENT(in) :: y(nobs)
      INTEGER, INTENT(inout) :: isDifZero
      INTEGER :: k

      ALLOCATE(s(bsg))
      ALLOCATE(oldb(bsg))
      isDifZero = 0
      oldb = b(startix:endix)
      s = MATMUL(y/(1.0D0+exp(r)), x(:, startix:endix))/nobs
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
      IF (ANY(ABS(dd) > 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         r = r + y * MATMUL(x(:,startix:endix), dd)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE log_update_step



   SUBROUTINE log_sp_update_step(bsg, startix, endix, b, lama, t_for_sg, pfg, pfl1,&
        lam1ma, x, y, xidx, xcptr, nnz, isDifZero, nobs, r, gamg, maxDif,&
        nvars, lb, ub)

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
      DOUBLE PRECISION, INTENT(in) :: pfl1(bsg)
      DOUBLE PRECISION, INTENT(in) :: x(nnz)
      DOUBLE PRECISION, INTENT(in) :: y(nobs)
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
      CALL spatx(x,xidx, xcptr, nobs, nvars, nnz, y/(1.0D0+exp(r)), s, startix, endix)
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
      IF(ANY(ABS(dd) > 0.0D0)) THEN
         maxDif = MAX(maxDif, gamg**2 * DOT_PRODUCT(dd,dd))
         CALL ypbspax(x, y, xidx, xcptr, nobs, nvars, nnz, dd, r, startix, endix, bsg)
         isDifZero = 1
      ENDIF
      DEALLOCATE(s, oldb, dd)
      RETURN
   END SUBROUTINE log_sp_update_step



END MODULE log_sgl_subfuns
