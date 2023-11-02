MODULE spmatmul

   IMPLICIT NONE

   CONTAINS
   ! ----
   ! These functions perform Ax=y and Atx = y for A in CSC form
   ! A is given by (a, ridx, cptr)
   ! a is the value, size nnz
   ! ridx is the row index of each nz
   ! cptr is of length ncol(A)+1 where cptr[j] is the first entry of a in column j, last entry == nnz
   ! to slice into columns cj:ck, you need a[cptr[cj]:(cptrr[ck+1]-1)] and
   ! ridx[cptr[cj]:(cptr[ck+1]-1)] and cptr[cj:ck])
     ! We don't actually need Ax->y ever, we need y <- y-A[,cj:ck] x (ymspax)
     ! and we need y <- y + b*A[,cj:ck]x (ypbspax)
   SUBROUTINE ymspax (a, ridx, cptr, n, p, nnz, x, y, cj, ck, lx)
      IMPLICIT NONE
      INTEGER :: n, p, nnz, cj, ck, lx
      DOUBLE PRECISION, INTENT(in) :: a(nnz)
      INTEGER, INTENT(in) :: ridx(nnz)
      INTEGER, INTENT(in) :: cptr(p+1)
      DOUBLE PRECISION, INTENT(in) :: x(lx)
      DOUBLE PRECISION, INTENT(inout) :: y(n)

      INTEGER :: i, j, k

      DO i = cj, ck
         k = cptr(i + 1) - 1
         DO j = cptr(i), k
            y(ridx(j)) = y(ridx(j)) - x(i - cj + 1) * a(j)
            ! r = r - MATMUL(x(:,startix:endix), dd) (a->x, x->dd, y->r)
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE ymspax

    SUBROUTINE ypbspax (a, b, ridx, cptr, n, p, nnz, x, y, cj, ck, lx)
      IMPLICIT NONE
      INTEGER :: n, p, nnz, cj, ck, lx
      DOUBLE PRECISION, INTENT(in) :: a(nnz), b(n)
      INTEGER, INTENT(in) :: ridx(nnz)
      INTEGER, INTENT(in) :: cptr(p+1)
      DOUBLE PRECISION, INTENT(in) :: x(lx)
      DOUBLE PRECISION, INTENT(inout) :: y(n)

      INTEGER :: i, j, k

      DO i = cj, ck
         k = cptr(i + 1) - 1
         DO j = cptr(i), k
            y(ridx(j)) = y(ridx(j)) + x(i - cj + 1) * a(j) * b(ridx(j))
            ! r = r + y * MATMUL(x(:,startix:endix), dd)  (a->x, b->y, x->dd, y->r)
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE ypbspax

   SUBROUTINE spatx (a, ridx, cptr, n, p, nnz, x, y, cj, ck)
      IMPLICIT NONE
      INTEGER :: n, p, nnz, cj, ck
      DOUBLE PRECISION, INTENT(in) :: a(nnz)
      INTEGER, INTENT(in) :: ridx(nnz)
      INTEGER, INTENT(in) :: cptr(p+1)
      DOUBLE PRECISION, INTENT(in) :: x(n)
      DOUBLE PRECISION, INTENT(inout) :: y(ck - cj + 1)

      INTEGER :: i, j, k
      y = 0.0D0

      DO i = cj, ck
         k = i - cj + 1
         DO j = cptr(i), (cptr(i + 1) - 1)
            y(k) = y(k) + x(ridx(j)) * a(j)
            ! s = MATMUL(r, x(:,startix:endix)) (a->x, x->r, y->s)
         ENDDO
      ENDDO
      RETURN
   END SUBROUTINE spatx

   SUBROUTINE softthresh(vec, thresh, n)
      IMPLICIT NONE
      INTEGER :: n, it
      DOUBLE PRECISION :: sg
      DOUBLE PRECISION :: vec(n)
      DOUBLE PRECISION :: thresh(n)
      DO it=1,n
         sg = vec(it)
         vec(it) = SIGN(MAX(ABS(sg) - thresh(it), 0.0D0), sg)
      ENDDO
      RETURN
   END SUBROUTINE softthresh

END MODULE spmatmul
