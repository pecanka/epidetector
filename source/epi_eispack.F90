MODULE EISPACK

  IMPLICIT NONE
  
 CONTAINS

SUBROUTINE bakvec ( n, t, e, m, z, ierr )

!*****************************************************************************80
!
!! BAKVEC determines eigenvectors by reversing the FIGI transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a nonsymmetric tridiagonal
!    matrix by back transforming those of the corresponding symmetric
!    matrix determined by FIGI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) T(N,3), contains the nonsymmetric matrix.  Its
!    subdiagonal is stored in the positions 2:N of the first column,
!    its diagonal in positions 1:N of the second column,
!    and its superdiagonal in positions 1:N-1 of the third column.
!    T(1,1) and T(N,3) are arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the symmetric matrix.  E(1) is arbitrary.
!    On output, the contents of E have been destroyed.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M), contains the eigenvectors.
!    On output, they have been transformed as requested.
!
!    Output, INTEGER ( kind = 4 ) IERR, an error flag.
!    0, for normal RETURN,
!    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
!    In this case, the symmetric matrix is not similar
!    to the original matrix, and the eigenvectors
!    cannot be found by this program.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) e(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) j
  REAL ( kind = 8 ) t(n,3)
  REAL ( kind = 8 ) z(n,m)

  ierr = 0

  IF( m == 0 ) then
    RETURN
  ENDIF

  e(1) = 1.0D+00
  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    IF( e(i) == 0.0D+00 ) then
      IF( t(i,1) /= 0.0D+00 .or. t(i-1,3) /= 0.0D+00 ) then
        ierr = 2 * n + i
        RETURN
      ENDIF
      e(i) = 1.0D+00
    else
      e(i) = e(i-1) * e(i) / t(i-1,3)
    ENDIF
  ENDDO

  DO j = 1, m
    z(2:n,j) = z(2:n,j) * e(2:n)
  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE balanc ( n, a, low, igh, scale )

!*****************************************************************************80
!
!! BALANC balances a real matrix before eigenvalue calculations.
!
!  Discussion:
!
!    This SUBROUTINE balances a real matrix and isolates eigenvalues
!    whenever possible.
!
!    Suppose that the principal submatrix in rows LOW through IGH
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!
!      SCALE(J) = P(J),    J = 1,...,LOW-1,
!               = D(J,J),  J = LOW,...,IGH,
!               = P(J)     J = IGH+1,...,N.
!
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is RETURNed for LOW if IGH is zero formally.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) A(N,N), the N by N matrix.  On output,
!    the matrix has been balanced.
!
!    Output, INTEGER ( kind = 4 ) LOW, IGH, indicate that A(I,J) is equal to
!    zero if
!    (1) I is greater than J and
!    (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!    Output, REAL ( kind = 8 ) SCALE(N), contains information determining the
!    permutations and scaling factors used.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b2
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) iexc
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  logical noconv
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) radix
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale(n)

  radix = 16.0D+00

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 CONTINUE

  scale(m) = j

  IF( j /= m ) then

    DO i = 1, l
      call r8_swap ( a(i,j), a(i,m) )
    ENDDO

    DO i = k, n
      call r8_swap ( a(j,i), a(m,i) )
    ENDDO

  ENDIF

!50 CONTINUE

  IF( iexc == 2 ) then
    go to 130
  ENDIF
!
!  Search for rows isolating an eigenvalue and push them down.
!
!80 CONTINUE

  IF( l == 1 ) then
    low = k
    igh = l
    RETURN
  ENDIF

  l = l - 1

100 CONTINUE

  DO j = l, 1, -1

     DO i = 1, l
       IF( i /= j ) then
         IF( a(j,i) /= 0.0D+00 ) then
           go to 120
         ENDIF
       ENDIF
     ENDDO

     m = l
     iexc = 1
     go to 20

120  CONTINUE

  ENDDO

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 CONTINUE

  k = k + 1

140 CONTINUE

  DO j = k, l

    DO i = k, l
      IF( i /= j ) then
        IF( a(i,j) /= 0.0D+00 ) then
          go to 170
        ENDIF
      ENDIF
    ENDDO

    m = k
    iexc = 2
    go to 20

170 CONTINUE

  ENDDO
!
!  Balance the submatrix in rows K to L.
!
  scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
  noconv = .true.

  DO WHILE ( noconv )

    noconv = .false.

    DO i = k, l

      c = 0.0D+00
      r = 0.0D+00

      DO j = k, l
        IF( j /= i ) then
          c = c + abs ( a(j,i) )
          r = r + abs ( a(i,j) )
        ENDIF
      ENDDO
!
!  Guard against zero C or R due to underflow.
!
      IF( c /= 0.0D+00 .and. r /= 0.0D+00 ) then

        g = r / radix
        f = 1.0D+00
        s = c + r

        DO WHILE ( c < g )
          f = f * radix
          c = c * b2
        ENDDO

        g = r * radix

        DO WHILE ( g <= c )
          f = f / radix
          c = c / b2
        ENDDO
!
!  Balance.
!
        IF( ( c + r ) / f < 0.95D+00 * s ) then

          g = 1.0D+00 / f
          scale(i) = scale(i) * f
          noconv = .true.

          a(i,k:n) = a(i,k:n) * g
          a(1:l,i) = a(1:l,i) * f

        ENDIF

      ENDIF

    ENDDO

  ENDDO

  low = k
  igh = l

  RETURN

END SUBROUTINE

SUBROUTINE balbak ( n, low, igh, scale, m, z )

!*****************************************************************************80
!
!! BALBAK determines eigenvectors by undoing the BALANC transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a real general matrix by
!    back transforming those of the corresponding balanced matrix
!    determined by BALANC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Parlett and Reinsch,
!    Numerische Mathematik,
!    Volume 13, pages 293-304, 1969.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, column indices determined by BALANC.
!
!    Input, REAL ( kind = 8 ) SCALE(N), contains information determining
!    the permutations and scaling factors used by BALANC.
!
!    Input, INTEGER ( kind = 4 ) M, the number of columns of Z to be
!    back-transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M), contains the real and imaginary 
!    parts of the eigenvectors, which, on RETURN, have been back-transformed.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) low
  !REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale(n)
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) z(n,m)

  IF( m <= 0 ) then
    RETURN
  ENDIF

  IF( igh /= low ) then
    DO i = low, igh
      z(i,1:m) = z(i,1:m) * scale(i)
    ENDDO
  ENDIF

  DO ii = 1, n

    i = ii

    IF( i < low .or. igh < i ) then

      IF( i < low ) then
        i = low - ii
      ENDIF

      k = int ( scale(i) )

      IF( k /= i ) then

        DO j = 1, m
          t      = z(i,j)
          z(i,j) = z(k,j)
          z(k,j) = t
        ENDDO

      ENDIF

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE bandr ( n, mb, a, d, e, e2, matz, z )

!*****************************************************************************80
!
!! BANDR reduces a symmetric band matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    This SUBROUTINE reduces a real symmetric band matrix
!    to a symmetric tridiagonal matrix using and optionally
!    accumulating orthogonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) MB, is the (half) band width of the matrix,
!    defined as the number of adjacent diagonals, including the principal
!    diagonal, required to specify the non-zero portion of the
!    lower triangle of the matrix.
!
!    Input/output, REAL ( kind = 8 ) A(N,MB).  On input, contains the lower
!    triangle of the symmetric band input matrix stored as an N by MB array.
!    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
!    column, its next subdiagonal in the last N+2-MB positions of the second
!    column, further subdiagonals similarly, and finally its principal diagonal
!    in the N positions of the last column.  Contents of storages not part of
!    the matrix are arbitrary.  On output, A has been destroyed, except for
!    its last two columns which contain a copy of the tridiagonal matrix.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, REAL ( kind = 8 ) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
!    Input, logical MATZ, should be set to TRUE if the transformation matrix is
!    to be accumulated, and to FALSE otherwise.
!
!    Output, REAL ( kind = 8 ) Z(N,N), the orthogonal transformation matrix
!    produced in the reduction if MATZ has been set to TRUE.  Otherwise, Z is
!    not referenced.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mb
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,mb)
  REAL ( kind = 8 ) b1
  REAL ( kind = 8 ) b2
  REAL ( kind = 8 ) c2
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) dmin
  REAL ( kind = 8 ) dminrt
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f1
  REAL ( kind = 8 ) f2
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) i1
  INTEGER ( kind = 4 ) i2
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) j1
  INTEGER ( kind = 4 ) j2
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) kr
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) m1
  logical matz
  INTEGER ( kind = 4 ) maxl
  INTEGER ( kind = 4 ) maxr
  INTEGER ( kind = 4 ) mr
  INTEGER ( kind = 4 ) r
  INTEGER ( kind = 4 ) r1
  REAL ( kind = 8 ) s2
  REAL ( kind = 8 ) u
  INTEGER ( kind = 4 ) ugl
  REAL ( kind = 8 ) z(n,n)

  dmin = epsilon ( dmin )
  dminrt = sqrt ( dmin )
!
!  Initialize the diagonal scaling matrix.
!
  d(1:n) = 1.0D+00

  IF( matz ) then

    z(1:n,1:n) = 0.0D+00

    DO i = 1, n
      z(i,i) = 1.0D+00
    ENDDO

  ENDIF

  m1 = mb - 1

  IF( m1 < 1 ) then
    d(1:n) = a(1:n,mb)
    e(1:n) = 0.0D+00
    e2(1:n) = 0.0D+00
    RETURN
  ENDIF

  IF( m1 /= 1 ) then

    DO k = 1, n - 2

      maxr = min ( m1, n - k )

      DO r1 = 2, maxr

        r = maxr + 2 - r1
        kr = k + r
        mr = mb - r
        g = a(kr,mr)
        a(kr-1,1) = a(kr-1,mr+1)
        ugl = k

        DO j = kr, n, m1

          j1 = j - 1
          j2 = j1 - 1

          IF( g == 0.0D+00 ) then
            exit
          ENDIF

          b1 = a(j1,1) / g
          b2 = b1 * d(j1) / d(j)
          s2 = 1.0D+00 / ( 1.0D+00 + b1 * b2 )

          IF( s2 < 0.5D+00 ) then

            b1 = g / a(j1,1)
            b2 = b1 * d(j) / d(j1)
            c2 = 1.0D+00 - s2
            d(j1) = c2 * d(j1)
            d(j) = c2 * d(j)
            f1 = 2.0D+00 * a(j,m1)
            f2 = b1 * a(j1,mb)
            a(j,m1) = -b2 * ( b1 * a(j,m1) - a(j,mb) ) - f2 + a(j,m1)
            a(j1,mb) = b2 * ( b2 * a(j,mb) + f1 ) + a(j1,mb)
            a(j,mb) = b1 * ( f2 - f1 ) + a(j,mb)

            DO l = ugl, j2
              i2 = mb - j + l
              u = a(j1,i2+1) + b2 * a(j,i2)
              a(j,i2) = -b1 * a(j1,i2+1) + a(j,i2)
              a(j1,i2+1) = u
            ENDDO

            ugl = j
            a(j1,1) = a(j1,1) + b2 * g

            IF( j /= n ) then

              maxl = min ( m1, n - j1 )

              DO l = 2, maxl
                i1 = j1 + l
                i2 = mb - l
                u = a(i1,i2) + b2 * a(i1,i2+1)
                a(i1,i2+1) = -b1 * a(i1,i2) + a(i1,i2+1)
                a(i1,i2) = u
              ENDDO

              i1 = j + m1

              IF( i1 <= n ) then
                g = b2 * a(i1,1)
              ENDIF

            ENDIF

            IF( matz ) then

              DO l = 1, n
                u = z(l,j1) + b2 * z(l,j)
                z(l,j) = -b1 * z(l,j1) + z(l,j)
                z(l,j1) = u
              ENDDO

            ENDIF

          else

            u = d(j1)
            d(j1) = s2 * d(j)
            d(j) = s2 * u
            f1 = 2.0D+00 * a(j,m1)
            f2 = b1 * a(j,mb)
            u = b1 * ( f2 - f1 ) + a(j1,mb)
            a(j,m1) = b2 * ( b1 * a(j,m1) - a(j1,mb) ) + f2 - a(j,m1)
            a(j1,mb) = b2 * ( b2 * a(j1,mb) + f1 ) + a(j,mb)
            a(j,mb) = u

            DO l = ugl, j2
              i2 = mb - j + l
              u = b2 * a(j1,i2+1) + a(j,i2)
              a(j,i2) = -a(j1,i2+1) + b1 * a(j,i2)
              a(j1,i2+1) = u
            ENDDO

            ugl = j
            a(j1,1) = b2 * a(j1,1) + g
 
            IF( j /= n ) then

              maxl = min ( m1, n - j1 )

              DO l = 2, maxl
                i1 = j1 + l
                i2 = mb - l
                u = b2 * a(i1,i2) + a(i1,i2+1)
                a(i1,i2+1) = -a(i1,i2) + b1 * a(i1,i2+1)
                a(i1,i2) = u
              ENDDO

              i1 = j + m1

              IF( i1 <= n ) then
                g = a(i1,1)
                a(i1,1) = b1 * a(i1,1)
              ENDIF

            ENDIF

            IF( matz ) then

              DO l = 1, n
                u = b2 * z(l,j1) + z(l,j)
                z(l,j) = -z(l,j1) + b1 * z(l,j)
                z(l,j1) = u
              ENDDO

            ENDIF

          ENDIF

        ENDDO

      ENDDO
!
!  Rescale to avoid underflow or overflow.
!
      IF( mod ( k, 64 ) == 0 ) then

        DO j = k, n

          IF( d(j) < dmin ) then

            maxl = max ( 1, mb + 1 - j )

            a(j,maxl:m1) = dminrt * a(j,maxl:m1)

            IF( j /= n ) then

              maxl = min ( m1, n - j )

              DO l = 1, maxl
                i1 = j + l
                i2 = mb - l
                a(i1,i2) = dminrt * a(i1,i2)
              ENDDO

            ENDIF

            IF( matz ) then
              z(1:n,j) = dminrt * z(1:n,j)
            ENDIF

            a(j,mb) = dmin * a(j,mb)
            d(j) = d(j) / dmin

          ENDIF

        ENDDO

      ENDIF

    ENDDO

  ENDIF
!
!  Form square root of scaling matrix.
!
  e(2:n) = sqrt ( d(2:n) )

  IF( matz ) then

    DO k = 2, n
      z(1:n,k) = z(1:n,k) * e(k)
    ENDDO

  ENDIF

  u = 1.0D+00

  DO j = 2, n
    a(j,m1) = u * e(j) * a(j,m1)
    u = e(j)
    e2(j) = a(j,m1)**2
    a(j,mb) = d(j) * a(j,mb)
    d(j) = a(j,mb)
    e(j) = a(j,m1)
  ENDDO

  d(1) = a(1,mb)
  e(1) = 0.0D+00
  e2(1) = 0.0D+00

  RETURN

END SUBROUTINE

SUBROUTINE bandv ( n, mbw, a, e21, m, w, z, ierr )

!*****************************************************************************80
!
!! BANDV finds eigenvectors from eigenvalues, for a real symmetric band matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvectors of a real symmetric
!    band matrix corresponding to specified eigenvalues, using inverse
!    iteration.  The SUBROUTINE may also be used to solve systems
!    of linear equations with a symmetric or non-symmetric band
!    coefficient matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) MBW, the number of columns of the array A used
!    to store the band matrix.  If the matrix is symmetric, MBW is its (half)
!    band width, denoted MB and defined as the number of adjacent
!    diagonals, including the principal diagonal, required to
!    specify the non-zero portion of the lower triangle of the
!    matrix.  If the SUBROUTINE is being used to solve systems
!    of linear equations and the coefficient matrix is not
!    symmetric, it must however have the same number of adjacent
!    diagonals above the main diagonal as below, and in this
!    case, MBW=2*MB-1.
!
!    Input, REAL ( kind = 8 ) A(N,MBW), the lower triangle of the symmetric
!    band input matrix stored as an N by MB array.  Its lowest subdiagonal is
!    stored in the last N+1-MB positions of the first column, its next
!    subdiagonal in the last N+2-MB positions of the second column, further
!    subdiagonals similarly, and finally its principal diagonal in the N
!    positions of column MB.  If the SUBROUTINE is being used to solve systems
!    of linear equations, and the coefficient matrix is not symmetric, A is
!    N by 2*MB-1 instead, with lower triangle as above and with its first
!    superdiagonal stored in the first N-1 positions of column MB+1, its
!    second superdiagonal in the first N-2 positions of column MB+2, further
!    superdiagonals similarly, and finally its highest superdiagonal in
!    the first N+1-MB positions of the last column.  Contents of storages
!    not part of the matrix are arbitrary.
!
!    Input, REAL ( kind = 8 ) E21, specifies the ordering of the eigenvalues
!    and contains 0.0 if the eigenvalues are in ascending order, or 2.0 if
!    the eigenvalues are in descending order.  If the SUBROUTINE is being used
!    to solve systems of linear equations, E21 should be set to 1.0
!    if the coefficient matrix is symmetric and to -1.0 if not.
!
!    Input, INTEGER ( kind = 4 ) M, the number of specified eigenvalues or the
!    number of systems of linear equations.
!
!    Input, REAL ( kind = 8 ) W(M), contains the M eigenvalues in ascending or
!    descending order.  If the SUBROUTINE is being used to solve systems of
!    linear equations (A-W(1:M)*I) * X(1:M) = B(1:M), where I is the identity
!    matrix, W should be set accordingly.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the constant matrix
!    columns B(1:M), if the SUBROUTINE is used to solve systems of linear
!    equations.  On output, the associated set of orthogonal eigenvectors.
!    Any vector which fails to converge is set to zero.  If the
!    routine is used to solve systems of linear equations,
!    Z contains the solution matrix columns X(1:M).
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
!    converge, or if the R-th system of linear equations is nearly singular.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mbw
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,mbw)
  REAL ( kind = 8 ) e21
  REAL ( kind = 8 ) eps2
  REAL ( kind = 8 ) eps3
  REAL ( kind = 8 ) eps4
  INTEGER ( kind = 4 ) group
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ij
  INTEGER ( kind = 4 ) ij1
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) kj
  INTEGER ( kind = 4 ) kj1
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) m1
  INTEGER ( kind = 4 ) m21
  INTEGER ( kind = 4 ) maxj
  INTEGER ( kind = 4 ) maxk
  INTEGER ( kind = 4 ) mb
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) order
  !REAL ( kind = 8 ) pythag
  INTEGER ( kind = 4 ) r
  REAL ( kind = 8 ) rv(n*(2*mbw-1))
  REAL ( kind = 8 ) rv6(n)
  REAL ( kind = 8 ) u
  REAL ( kind = 8 ) uk
  REAL ( kind = 8 ) v
  REAL ( kind = 8 ) w(m)
  REAL ( kind = 8 ) x0
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) xu
  REAL ( kind = 8 ) z(n,m)

  ierr = 0

  IF( m == 0 ) then
    RETURN
  ENDIF

  x0 = 0.0D+00

  IF( e21 < 0.0D+00 ) then
    mb = ( mbw + 1 ) / 2
  else
    mb = mbw
  ENDIF

  m1 = mb - 1
  m21 = m1 + mb
  order = 1.0D+00 - abs ( e21 )
!
!  Find vectors by inverse iteration.
!
  DO r = 1, m

    its = 1
    x1 = w(r)

    IF( r /= 1 ) then
      go to 100
    ENDIF
!
!  Compute norm of matrix.
!
    norm = 0.0D+00

    DO j = 1, mb

      jj = mb + 1 - j
      kj = jj + m1
      ij = 1
      v = 0.0D+00

      DO i = jj, n

        v = v + abs ( a(i,j) )

        IF( e21 < 0.0D+00 ) then
          v = v + abs ( a(ij,kj) )
          ij = ij + 1
        ENDIF

      ENDDO

      norm = max ( norm, v )

    ENDDO

    IF( e21 < 0.0D+00 ) then
      norm = 0.5D+00 * norm
    ENDIF
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by eps3,
!  EPS4 is taken very small to avoid overflow.
!
    IF( norm == 0.0D+00 ) then
      norm = 1.0D+00
    ENDIF

    eps2 = 0.001D+00 * norm * abs ( order)
    eps3 = abs ( norm ) * epsilon ( norm )
    uk = n
    uk = sqrt ( uk )
    eps4 = uk * eps3

80  CONTINUE

    group = 0
    go to 120
!
!  Look for close or coincident roots.
!
100 CONTINUE

    IF( eps2 <= abs ( x1 - x0 ) ) then
      go to 80
    ENDIF

    group = group + 1

    IF( order * ( x1 - x0 ) <= 0.0D+00 ) then
      x1 = x0 + order * eps3
    ENDIF
!
!  Expand matrix, subtract eigenvalue, and initialize vector.
!
120 CONTINUE

    DO i = 1, n

      ij = i + min ( 0, i - m1 ) * n
      kj = ij + mb * n
      ij1 = kj + m1 * n

      IF( m1 == 0 ) then
        go to 180
      ENDIF

      DO j = 1, m1

        IF( ij <= m1 ) then
          IF( ij <= 0 ) then
            rv(ij1) = 0.0D+00
            ij1 = ij1 + n
          ENDIF
        else
          rv(ij) = a(i,j)
        ENDIF

        ij = ij + n
        ii = i + j

        IF( ii <= n ) then

          jj = mb - j

          IF( e21 < 0.0D+00 ) then
            ii = i
            jj = mb + j
          ENDIF

          rv(kj) = a(ii,jj)
          kj = kj + n

        ENDIF

      ENDDO

180   CONTINUE

      rv(ij) = a(i,mb) - x1
      rv6(i) = eps4
      IF( order == 0.0D+00 ) then
        rv6(i) = z(i,r)
      ENDIF

    ENDDO

    IF( m1 /= 0 ) then
!
!  Elimination with interchanges.
!
      DO i = 1, n

        ii = i + 1
        maxk = min ( i + m1 - 1, n )
        maxj = min ( n - i, m21 - 2 ) * n

        DO k = i, maxk

          kj1 = k
          j = kj1 + n
          jj = j + maxj

          DO kj = j, jj, n
            rv(kj1) = rv(kj)
            kj1 = kj
          ENDDO

          rv(kj1) = 0.0D+00

        ENDDO

        IF( i /= n ) then

          u = 0.0D+00
          maxk = min ( i + m1, n )
          maxj = min ( n - ii, m21 - 2 ) * n

          DO j = i, maxk
            IF( abs ( u ) <= abs ( rv(j) ) ) then
              u = rv(j)
              k = j
            ENDIF
          ENDDO

          j = i + n
          jj = j + maxj

          IF( k /= i ) then

            kj = k

            DO ij = i, jj, n
              call r8_swap ( rv(ij), rv(kj) )
              kj = kj + n
            ENDDO

            IF( order == 0.0D+00 ) then
              call r8_swap ( rv6(i), rv6(k) )
            ENDIF

          ENDIF

          IF( u /= 0.0D+00 ) then

            DO k = ii, maxk

              v = rv(k) / u
              kj = k

              DO ij = j, jj, n
                kj = kj + n
                rv(kj) = rv(kj) - v * rv(ij)
              ENDDO

              IF( order == 0.0D+00 ) then
                rv6(k) = rv6(k) - v * rv6(i)
              ENDIF

            ENDDO

          ENDIF

        ENDIF

      ENDDO

    ENDIF
!
!  Back substitution.
!
600 CONTINUE

    DO ii = 1, n

      i = n + 1 - ii
      maxj = min ( ii, m21 )

      IF( maxj /= 1 ) then

        ij1 = i
        j = ij1 + n
        jj = j + ( maxj - 2 ) * n

        DO ij = j, jj, n
          ij1 = ij1 + 1
          rv6(i) = rv6(i) - rv(ij) * rv6(ij1)
        ENDDO

      ENDIF

      v = rv(i)
!
!  Error: nearly singular linear system.
!
      IF( abs ( v ) < eps3 ) then
        IF( order == 0.0D+00 ) then
          ierr = -r
        ENDIF
        v = sign ( eps3, v )
      ENDIF

      rv6(i) = rv6(i) / v

    ENDDO

    xu = 1.0D+00

    IF( order == 0.0D+00 ) then
      go to 870
    ENDIF
!
!  Orthogonalize with respect to previous members of group.
!
    DO jj = 1, group

      j = r - group - 1 + jj

      xu = dot_product ( rv6(1:n), z(1:n,j) )

      rv6(1:n) = rv6(1:n) - xu * z(1:n,j)

    ENDDO

    norm = sum ( abs ( rv6(1:n) ) )
!
!  Choose a new starting vector.
!
    IF( norm < 0.1D+00 ) then

      IF( its < n ) then
        its = its + 1
        xu = eps4 / ( uk + 1.0D+00 )
        rv6(1) = eps4
        rv6(2:n) = xu
        rv6(its) = rv6(its) - eps4 * uk
        go to 600
      else
        ierr = -r
        xu = 0.0D+00
        go to 870
      ENDIF

    ENDIF
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
    u = 0.0D+00
    DO i = 1, n
      u = pythag ( u, rv6(i) )
    ENDDO

    xu = 1.0D+00 / u

870 CONTINUE

    z(1:n,r) = rv6(1:n) * xu

    x0 = x1

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE bisect ( n, eps1, d, e, e2, lb, ub, mm, m, w, ind, ierr )

!*****************************************************************************80
!
!! BISECT computes some eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvalues of a real symmetric
!    tridiagonal matrix which lie in a specified interval, using bisection.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) EPS1, is an absolute error tolerance for
!    the computed eigenvalues.  If the input EPS1 is non-positive, it is reset
!    for each submatrix to a default value, namely, minus the product of the
!    relative machine precision and the 1-norm of the submatrix.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the input matrix.
!
!    Input, REAL ( kind = 8 ) E(N), contains in E(2:N) the subdiagonal elements
!    of the matrix.  E(1) is arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E2(N).  On input, the squares of the
!    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of 
!    E2, corresponding to elements of E regarded as negligible, have been
!    replaced by zero, causing the matrix to split into a direct sum of
!    submatrices.  E2(1) is also set to zero.
!
!    Input, REAL ( kind = 8 ) LB, UB, define the interval to be searched for
!    eigenvalues.  If LB is not less than UB, no eigenvalues will be found.
!
!    Input, INTEGER ( kind = 4 ) MM, an upper bound for the number of
!    eigenvalues in the interval.  Warning: if more than MM eigenvalues are
!    determined to lie in the interval, an error RETURN is made with no
!    eigenvalues found.
!
!    Output, INTEGER ( kind = 4 ) M, the number of eigenvalues determined to lie
!    in (LB,UB).
!
!    Output, REAL ( kind = 8 ) W(M), the eigenvalues in ascending order.
!
!    Output, INTEGER ( kind = 4 ) IND(MM), contains in its first M positions
!    the submatrix indices associated with the corresponding eigenvalues in W:
!    1 for eigenvalues belonging to the first submatrix from the top, 2 for
!    those belonging to the second submatrix, and so on.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    3*N+1, if M exceeds MM.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) eps1
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ind(mm)
  INTEGER ( kind = 4 ) isturm
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) lb
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) m1
  INTEGER ( kind = 4 ) m2
  INTEGER ( kind = 4 ) p
  INTEGER ( kind = 4 ) q
  INTEGER ( kind = 4 ) r
  REAL ( kind = 8 ) rv4(n)
  REAL ( kind = 8 ) rv5(n)
  INTEGER ( kind = 4 ) s
  REAL ( kind = 8 ) t1
  REAL ( kind = 8 ) t2
  INTEGER ( kind = 4 ) tag
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) u
  REAL ( kind = 8 ) ub
  REAL ( kind = 8 ) v
  REAL ( kind = 8 ) w(mm)
  REAL ( kind = 8 ) x0
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) xu

  ierr = 0
  s = 0
  tag = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries.
!
  e2(1) = 0.0D+00

  DO i = 2, n

    tst1 = abs ( d(i) ) + abs ( d(i-1) )
    tst2 = tst1 + abs ( e(i) )

    IF( tst2 <= tst1 ) then
      e2(i) = 0.0D+00
    ENDIF

  ENDDO
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 CONTINUE

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 CONTINUE

  m = m - s

  IF( mm < m ) then
    go to 980
  ENDIF

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining
!  interval by the Gerschgorin bounds.
!
100 CONTINUE

  IF( r == m ) then
    go to 1001
  ENDIF

  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0D+00

  DO q = p, n

    x1 = u
    u = 0.0D+00
    v = 0.0D+00

    IF( q /= n ) then
      u = abs ( e(q+1) )
      v = e2(q+1)
    ENDIF

    xu = min ( d(q) - ( x1 + u ), xu )
    x0 = max ( d(q) + ( x1 + u ), x0 )

    IF( v == 0.0D+00 ) then
      exit
    ENDIF

  ENDDO

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
  IF( eps1 <= 0.0D+00 ) then
    eps1 = -x1
  ENDIF

  IF( p /= q ) then
    go to 180
  ENDIF
!
!  Check for an isolated root within interval.
!
  IF( d(p) < t1 .or. t2 <= d(p) ) then
    go to 940
  ENDIF

  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900

  180 CONTINUE

  x1 = x1 * ( q - p + 1 )
  lb = max ( t1, xu - x1 )
  ub = min ( t2, x0 + x1 )
  x1 = lb
  isturm = 3
  go to 320

  200 CONTINUE

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

  220 CONTINUE

  m2 = s
  IF( m2 < m1 ) then
    go to 940
  ENDIF
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5
  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for the K-th eigenvalue.
!
  k = m2

250 CONTINUE

     xu = lb

     DO ii = m1, k
       i = m1 + k - ii
       IF( xu < rv4(i) ) then
         xu = rv4(i)
         go to 280
       ENDIF
     ENDDO

  280 CONTINUE

   x0 = min ( x0, rv5(k) )
!
!  Next bisection step.
!
  300    CONTINUE

     x1 = ( xu + x0 ) * 0.5D+00

     IF( ( x0 - xu ) <= abs ( eps1 ) ) then
       go to 420
     ENDIF

     tst1 = 2.0D+00 * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + ( x0 - xu )
     IF( tst2 == tst1 ) then
       go to 420
     ENDIF
!
!  Sturm sequence.
!
320  CONTINUE

     s = p - 1
     u = 1.0D+00

     DO i = p, q

        IF( u == 0.0D+00 ) then
          v = abs ( e(i) ) / epsilon ( v )
          IF( e2(i) == 0.0D+00 ) then
            v = 0.0D+00
          ENDIF
        else
          v = e2(i) / u
        ENDIF

        u = d(i) - x1 - v
        IF( u < 0.0D+00 ) then
          s = s + 1
        ENDIF

     ENDDO

     go to (60,80,200,220,360), isturm
!
!  Refine intervals.
!
  360 CONTINUE

     IF( k <= s ) then
       go to 400
     ENDIF

     xu = x1

     IF( s < m1 ) then
       rv4(m1) = x1
       go to 300
     ENDIF

!  380 CONTINUE

     rv4(s+1) = x1

     IF( x1 < rv5(s) ) then
       rv5(s) = x1
     ENDIF

     go to 300
400  CONTINUE
     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420 CONTINUE

  rv5(k) = x1
  k = k - 1

  IF( m1 <= k ) then
    go to 250
  ENDIF
!
!  Order eigenvalues tagged with their submatrix associations.
!
900 CONTINUE

  s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  DO l = 1, r

    IF( j <= s ) then

      IF( m2 < k ) then
        exit
      ENDIF

      IF( w(l) <= rv5(k) ) then
        j = j + 1
        cycle
      ENDIF

      DO ii = j, s
        i = l + s - ii
        w(i+1) = w(i)
        ind(i+1) = ind(i)
      ENDDO

    ENDIF

    w(l) = rv5(k)
    ind(l) = tag
    k = k + 1

  ENDDO

940 CONTINUE

  IF( q < n ) then
    go to 100
  ENDIF

  go to 1001
!
!  Set error: underestimate of number of eigenvalues in interval.
!
980 CONTINUE

  ierr = 3 * n + 1

 1001 CONTINUE

  lb = t1
  ub = t2

  RETURN

END SUBROUTINE

SUBROUTINE bqr ( n, mb, a, t, r, ierr )

!*****************************************************************************80
!
!! BQR finds the smallest eigenvalue of a real symmetric band matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalue of smallest magnitude of a real
!    symmetric band matrix using the QR algorithm with shifts of origin.
!    Consecutive calls can be made to find further eigenvalues.
!
!    Note that for a subsequent call, N should be replaced by N-1, but
!    MB should not be altered even when it exceeds the current N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) MB, the (half) band width of the matrix,
!    defined as the number of adjacent diagonals, including the principal
!    diagonal, required to specify the non-zero portion of the
!    lower triangle of the matrix.
!
!    Input/output, REAL ( kind = 8 ) A(N,MB).  On input, A contains the lower
!    triangle of the symmetric band input matrix stored as an N by MB array.
!    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
!    column, its next subdiagonal in the last N+2-MB positions of the
!    second column, further subdiagonals similarly, and finally its principal
!    diagonal in the N positions of the last column.  Contents of storages
!    not part of the matrix are arbitrary.  On a subsequent call, its output
!    contents from the previous call should be passed.  On output, A contains
!    the transformed band matrix.  The matrix A+T*I derived from the output
!    parameters is similar to the input A+T*I to within rounding errors.
!    Its last row and column are null as long as IERR is zero.
!
!    Input/output, REAL ( kind = 8 ) T.  On input, T specifies the shift (of
!    eigenvalues) applied to the diagonal of A in forming the input matrix.
!    What is actually determined is the eigenvalue nearest to T of A+T*I, where
!    I is the identity matrix.  On a subsequent call, the output value of T
!    from the previous call should be passed if the next nearest eigenvalue
!    is sought.  On output, T contains the computed eigenvalue of A+T*I,
!    as long as IERR is zero.
!
!    Input/output, REAL ( kind = 8 ) R.  On input for the first call, R should 
!    be specified as zero, and as its output value from the previous call
!    on a subsequent call.  It is used to determine when the last row and
!    column of the transformed band matrix can be regarded as negligible.
!    On output, R contains the maximum of its input value and the norm of the
!    last column of the input matrix A.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, normal RETURN.
!    N, if the eigenvalue has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mb
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,mb)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ik
  INTEGER ( kind = 4 ) imult
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jk
  INTEGER ( kind = 4 ) jm
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) kj
  INTEGER ( kind = 4 ) kj1
  INTEGER ( kind = 4 ) kk
  INTEGER ( kind = 4 ) km
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) m1
  INTEGER ( kind = 4 ) m2
  INTEGER ( kind = 4 ) m21
  INTEGER ( kind = 4 ) m3
  INTEGER ( kind = 4 ) m31
  INTEGER ( kind = 4 ) m4
  INTEGER ( kind = 4 ) mk
  INTEGER ( kind = 4 ) mn
  INTEGER ( kind = 4 ) mz
  INTEGER ( kind = 4 ) ni
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) q
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) rv(2*mb*mb+4*mb-3)
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2

  ierr = 0
  m1 = min ( mb, n )
  m = m1 - 1
  m2 = m + m
  m21 = m2 + 1
  m3 = m21 + m
  m31 = m3 + 1
  m4 = m31 + m2
  mn = m + n
  mz = mb - m1
  its = 0
!
!  Test for convergence.
!
40 CONTINUE

  g = a(n,mb)

  IF( m == 0 ) then
    go to 360
  ENDIF

  f = 0.0D+00
  DO k = 1, m
    mk = k + mz
    f = f + abs ( a(n,mk) )
  ENDDO

  IF( its == 0 .and. r < f ) then
    r = f
  ENDIF

  tst1 = r
  tst2 = tst1 + f

  IF( tst2 <= tst1 ) then
    go to 360
  ENDIF

  IF( 30 <= its ) then
    ierr = n
    RETURN
  ENDIF

  its = its + 1
!
!  Form shift from bottom 2 by 2 minor.
!
  IF( f <= 0.25D+00 * r .or. 5 <= its ) then

    f = a(n,mb-1)

    IF( f /= 0.0D+00 ) then
      q = ( a(n-1,mb) - g ) / ( 2.0D+00 * f )
      s = pythag ( q, 1.0D+00 )
      g = g - f / ( q + sign ( s, q ) )
    ENDIF

    t = t + g

    a(1:n,mb) = a(1:n,mb) - g

  ENDIF

  rv(m31:m4) = 0.0D+00

  DO ii = 1, mn

     i = ii - m
     ni = n - ii

     IF( ni < 0 ) then
       go to 230
     ENDIF
!
!  Form column of shifted matrix A-G*I.
!
     l = max ( 1, 2-i )

     rv(1:m3) = 0.0D+00

     DO k = l, m1
       km = k + m
       mk = k + mz
       rv(km) = a(ii,mk)
     ENDDO

     ll = min ( m, ni )

     DO k = 1, ll
       km = k + m21
       ik = ii + k
       mk = mb - k
       rv(km) = a(ik,mk)
     ENDDO
!
!  Pre-multiply with Householder reflections.
!
     ll = m2
     imult = 0
!
!  Multiplication procedure.
!
140  CONTINUE

     kj = m4 - m1

     DO j = 1, ll

        kj = kj + m1
        jm = j + m3

        IF( rv(jm) /= 0.0D+00 ) then

          f = 0.0D+00

          DO k = 1, m1
            kj = kj + 1
            jk = j + k - 1
            f = f + rv(kj) * rv(jk)
          ENDDO

          f = f / rv(jm)
          kj = kj - m1

          DO k = 1, m1
            kj = kj + 1
            jk = j + k - 1
            rv(jk) = rv(jk) - rv(kj) * f
          ENDDO

          kj = kj - m1

        ENDIF

     ENDDO

     IF( imult /= 0 ) then
       go to 280
     ENDIF
!
!  Householder reflection.
!
     f = rv(m21)
     s = 0.0D+00
     rv(m4) = 0.0D+00
     scale = sum ( abs ( rv(m21:m3) ) )

     IF( scale == 0.0D+00 ) then
       go to 210
     ENDIF

     DO k = m21, m3
       s = s + ( rv(k) / scale )**2
     ENDDO

     s = scale * scale * s
     g = - sign ( sqrt ( s ), f )
     rv(m21) = g
     rv(m4) = s - f * g
     kj = m4 + m2 * m1 + 1
     rv(kj) = f - g

     DO k = 2, m1
       kj = kj + 1
       km = k + m2
       rv(kj) = rv(km)
     ENDDO
!
!  Save column of triangular factor R.
!
210  CONTINUE

     DO k = l, m1
       km = k + m
       mk = k + mz
       a(ii,mk) = rv(km)
     ENDDO

230  CONTINUE

     l = max ( 1, m1+1-i )
     IF( i <= 0 ) then
       go to 300
     ENDIF
!
!  Perform additional steps.
!
     rv(1:m21) = 0.0D+00
     ll = min ( m1, ni + m1 )
!
!  Get row of triangular factor R.
!
     DO kk = 1, ll
       k = kk - 1
       km = k + m1
       ik = i + k
       mk = mb - k
       rv(km) = a(ik,mk)
     ENDDO
!
!  Post-multiply with Householder reflections.
!
     ll = m1
     imult = 1
     go to 140
!
!  Store column of new a matrix.
!
280  CONTINUE

     DO k = l, m1
       mk = k + mz
       a(i,mk) = rv(k)
     ENDDO
!
!  Update Householder reflections.
!
300  CONTINUE

     IF( 1 < l ) then
       l = l - 1
     ENDIF

     kj1 = m4 + l * m1

     DO j = l, m2

       jm = j + m3
       rv(jm) = rv(jm+1)

       DO k = 1, m1
         kj1 = kj1 + 1
         kj = kj1 - m1
         rv(kj) = rv(kj1)
       ENDDO

     ENDDO

  ENDDO

  go to 40
!
!  Convergence.
!
360 CONTINUE

  t = t + g
  a(1:n,mb) = a(1:n,mb) - g

  DO k = 1, m1
    mk = k + mz
    a(n,mk) = 0.0D+00
  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE cbabk2 ( n, low, igh, scale, m, zr, zi )

!*****************************************************************************80
!
!! CBABK2 finds eigenvectors by undoing the CBAL transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    balanced matrix determined by CBAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, values determined by CBAL.
!
!    Input, REAL ( kind = 8 ) SCALE(N), information determining the permutations
!    and scaling factors used by CBAL.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) ZR(N,M), ZI(N,M).  On input, the real
!    and imaginary parts, respectively, of the eigenvectors to be back
!    transformed in their first M columns.  On output, the transformed
!    eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) low
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale(n)
  REAL ( kind = 8 ) zi(n,m)
  REAL ( kind = 8 ) zr(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  IF( igh /= low ) then

    DO i = low, igh

      s = scale(i)

      zr(i,1:m) = zr(i,1:m) * s
      zi(i,1:m) = zi(i,1:m) * s

    ENDDO

  ENDIF

  DO ii = 1, n

    i = ii

    IF( i < low .or. igh < i ) then

      IF( i < low ) then
        i = low - ii
      ENDIF

      k = INT(scale(i))

      IF( k /= i ) then

        DO j = 1, m
          call r8_swap ( zr(i,j), zr(k,j) )
          call r8_swap ( zi(i,j), zi(k,j) )
        ENDDO

      ENDIF

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE cbal ( n, ar, ai, low, igh, scale )

!*****************************************************************************80
!
!! CBAL balances a complex matrix before eigenvalue calculations.
!
!  Discussion:
!
!    This SUBROUTINE balances a complex matrix and isolates
!    eigenvalues whenever possible.
!
!    Suppose that the principal submatrix in rows low through igh
!    has been balanced, that P(J) denotes the index interchanged
!    with J during the permutation step, and that the elements
!    of the diagonal matrix used are denoted by D(I,J).  Then
!      SCALE(J) = P(J),    for J = 1,...,LOW-1
!               = D(J,J)       J = LOW,...,IGH
!               = P(J)         J = IGH+1,...,N.
!    The order in which the interchanges are made is N to IGH+1,
!    then 1 to LOW-1.
!
!    Note that 1 is RETURNed for IGH if IGH is zero formally.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real and
!    imaginary parts of the complex matrix to be balanced.  On output,
!    the real and imaginary parts of the balanced matrix.
!
!    Output, INTEGER ( kind = 4 ) LOW, IGH, are values such that AR(I,J)
!    and AI(I,J) are zero if I is greater than J and either J=1,...,LOW-1 or
!    I=IGH+1,...,N.
!
!    Output, REAL ( kind = 8 ) SCALE(N), information determining the
!    permutations and scaling factors used.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) b2
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) iexc
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  logical noconv
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) radix
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale(n)

  radix = 16.0D+00

  iexc = 0
  j = 0
  m = 0

  b2 = radix * radix
  k = 1
  l = n
  go to 100

20 CONTINUE

  scale(m) = j

  IF( j /= m ) then

    DO i = 1, l
      call r8_swap ( ar(i,j), ar(i,m) )
      call r8_swap ( ai(i,j), ai(i,m) )
    ENDDO

    DO i = k, n
      call r8_swap ( ar(j,i), ar(m,i) )
      call r8_swap ( ai(j,i), ai(m,i) )
    ENDDO

  ENDIF

  IF( iexc == 2 ) then
    go to 130
  ENDIF
!
!  Search for rows isolating an eigenvalue and push them down.
!
!80 CONTINUE

  IF( l == 1 ) then
    go to 280
  ENDIF

  l = l - 1

100 CONTINUE

  DO jj = 1, l

     j = l + 1 - jj

     DO i = 1, l
       IF( i /= j ) then
         IF( ar(j,i) /= 0.0D+00 .or. ai(j,i) /= 0.0D+00 ) then
           go to 120
         ENDIF
       ENDIF
     ENDDO

     m = l
     iexc = 1
     go to 20

120  CONTINUE

  ENDDO

  go to 140
!
!  Search for columns isolating an eigenvalue and push them left.
!
130 CONTINUE

  k = k + 1

140 CONTINUE

   DO j = k, l

     DO i = k, l
       IF( i /= j ) then
         IF( ar(i,j) /= 0.0D+00 .or. ai(i,j) /= 0.0D+00 ) then
           go to 170
         ENDIF
       ENDIF
     ENDDO

     m = k
     iexc = 2
     go to 20

170  CONTINUE

  ENDDO
!
!  Now balance the submatrix in rows k to l.
!
  scale(k:l) = 1.0D+00
!
!  Iterative loop for norm reduction.
!
190 CONTINUE

  noconv = .false.

  DO i = k, l

    c = 0.0D+00
    r = 0.0D+00

    DO j = k, l
      IF( j /= i ) then
        c = c + abs ( ar(j,i) ) + abs ( ai(j,i) )
        r = r + abs ( ar(i,j) ) + abs ( ai(i,j) )
      ENDIF
    ENDDO
!
!  Guard against zero C or R due to underflow.
!
     IF( c == 0.0D+00 .or. r == 0.0D+00 ) then
       go to 270
     ENDIF

     g = r / radix
     f = 1.0D+00
     s = c + r

     DO WHILE ( c < g )
       f = f * radix
       c = c * b2
     ENDDO

     g = r * radix

     DO WHILE  ( g <= c )
       f = f / radix
       c = c / b2
     ENDDO
!
!  Now balance.
!
     IF( ( c + r ) / f < 0.95D+00 * s ) then

       g = 1.0D+00 / f
       scale(i) = scale(i) * f
       noconv = .true.

       ar(i,k:n) = ar(i,k:n) * g
       ai(i,k:n) = ai(i,k:n) * g

       ar(1:l,i) = ar(1:l,i) * f
       ai(1:l,i) = ai(1:l,i) * f

     ENDIF

270  CONTINUE

  ENDDO

  IF( noconv ) then
    go to 190
  ENDIF

  280 CONTINUE

  low = k
  igh = l

  RETURN

END SUBROUTINE

SUBROUTINE cdiv ( ar, ai, br, bi, cr, ci )

!*****************************************************************************80
!
!! CDIV emulates complex division, using real arithmetic.
!
!  Discussion:
!
!    This routine performs complex division:
!
!      (CR,CI) = (AR,AI) / (BR,BI)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, REAL ( kind = 8 ) AR, AI, the real and imaginary parts of
!    the numerator.
!
!    Input, REAL ( kind = 8 ) BR, BI, the real and imaginary parts of
!    the denominator.
!
!    Output, REAL ( kind = 8 ) CR, CI, the real and imaginary parts of
!    the result.
!
  IMPLICIT NONE

  REAL ( kind = 8 ) ai
  REAL ( kind = 8 ) ais
  REAL ( kind = 8 ) ar
  REAL ( kind = 8 ) ars
  REAL ( kind = 8 ) bi
  REAL ( kind = 8 ) bis
  REAL ( kind = 8 ) br
  REAL ( kind = 8 ) brs
  REAL ( kind = 8 ) ci
  REAL ( kind = 8 ) cr
  REAL ( kind = 8 ) s

  s = abs ( br ) + abs ( bi )

  ars = ar / s
  ais = ai / s
  brs = br / s
  bis = bi / s

  s = brs**2 + bis**2
  cr = ( ars * brs + ais * bis ) / s
  ci = ( ais * brs - ars * bis ) / s

  RETURN

END SUBROUTINE

SUBROUTINE cg ( n, ar, ai, wr, wi, matz, zr, zi, ierr )

!*****************************************************************************80
!
!! CG gets eigenvalues and eigenvectors of a complex general matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of EISPACK SUBROUTINEs
!    to find the eigenvalues and eigenvectors (if desired)
!    of a complex general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real and
!    imaginary parts of the complex matrix.  On output, AR and AI
!    have been overwritten by other information.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts
!    of the eigenvalues.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is 0 if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are to be computed.
!
!    Output, REAL ( kind = 8 ) ZR(N,N), ZI(N,N), the real and imaginary parts,
!    respectively, of the eigenvectors, if MATZ is not zero.
!
!    Output, INTEGER ( kind = 4 ) IERR, an error completion code described in 
!    the documentation for COMQR and COMQR2.  The normal completion code 
!    is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  REAL ( kind = 8 ) fv3(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) is1
  INTEGER ( kind = 4 ) is2
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) zi(n,n)
  REAL ( kind = 8 ) zr(n,n)

  call cbal ( n, ar, ai, is1, is2, fv1 )

  call corth ( n, is1, is2, ar, ai, fv2, fv3 )

  IF( matz == 0 ) then

    call comqr ( n, is1, is2, ar, ai, wr, wi, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

  else

    call comqr2 ( n, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CG - Fatal error!'
      write ( *, '(a)' ) '  Nonzero error RETURN from COMQR2.'
      RETURN
    ENDIF

    call cbabk2 ( n, is1, is2, fv1, n, zr, zi )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE ch ( n, ar, ai, w, matz, zr, zi, ierr )

!*****************************************************************************80
!
!! CH gets eigenvalues and eigenvectors of a complex Hermitian matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of SUBROUTINEs from the
!    EISPACK eigensystem package to find the eigenvalues and eigenvectors
!    of a complex hermitian matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real and
!    imaginary parts of the complex matrix.  On output, AR and AI
!    have been overwritten by other information.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is 0 if only eigenvalues are desired, and
!    nonzero if both eigenvalues and eigenvectors are to be computed.
!
!    Output, REAL ( kind = 8 ) ZR(N,N), ZI(N,N), the real and imaginary parts,
!    respectively, of the eigenvectors, if MATZ is not zero.
!
!    Output, INTEGER ( kind = 4 ) IERR, an error completion code described in 
!    the documentation for TQLRAT and TQL2.  The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) fm1(2,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) zi(n,n)
  REAL ( kind = 8 ) zr(n,n)

  call htridi ( n, ar, ai, w, fv1, fv2, fm1 )

  IF( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

  else

    zr(1:n,1:n) = 0.0D+00

    DO i = 1, n
      zr(i,i) = 1.0D+00
    ENDDO

    call tql2 ( n, w, fv1, zr, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

    call htribk ( n, ar, ai, fm1, n, zr, zi )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE cinvit ( n, ar, ai, wr, wi, select, mm, m, zr, zi, ierr )

!*****************************************************************************80
!
!! CINVIT gets eigenvectors from eigenvalues, for a complex Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvectors of a complex upper
!    Hessenberg matrix corresponding to specified eigenvalues,
!    using inverse iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) AR(N,N), AI(N,N), the real and imaginary parts of
!    the complex Hessenberg matrix.
!
!    Input/output, REAL ( kind = 8 ) WR(N), WI(N).  On input, the real and
!    imaginary parts of the eigenvalues of the matrix.  The eigenvalues must
!    be stored in a manner identical to that of SUBROUTINE COMLR, which
!    recognizes possible splitting of the matrix.  On output, WR may have been
!    altered since close eigenvalues are perturbed slightly in searching for
!    independent eigenvectors.
!
!    Input, logical SELECT(N), specifies the eigenvectors to be found.  The
!    eigenvector corresponding to the J-th eigenvalue is specified by
!    setting SELECT(J) to TRUE.
!
!    Input, INTEGER ( kind = 4 ) MM, an upper bound for the number of 
!    eigenvectors to be found.
!
!    Output, INTEGER ( kind = 4 ) M, the number of eigenvectors actually found.
!
!    Output, REAL ( kind = 8 ) ZR(N,MM), ZI(N,MM), the real and imaginary parts
!    of the eigenvectors.  The eigenvectors are normalized so that the
!    component of largest magnitude is 1.
!    Any vector which fails the acceptance test is set to zero.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    -(2*N+1), if more than MM eigenvectors have been specified,
!    -K, if the iteration corresponding to the K-th value fails,
!    -(N+K), if both error situations occur.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) eps3
  REAL ( kind = 8 ) growto
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  REAL ( kind = 8 ) ilambd
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) km1
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) normv
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) rlambd
  REAL ( kind = 8 ) rm1(n,n)
  REAL ( kind = 8 ) rm2(n,n)
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) rv2(n)
  INTEGER ( kind = 4 ) s
  logical select(n)
  INTEGER ( kind = 4 ) uk
  REAL ( kind = 8 ) ukroot
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) zi(n,mm)
  REAL ( kind = 8 ) zr(n,mm)

  ierr = 0
  uk = 0
  s = 1

  DO k = 1, n

    IF( .not. select(k) ) then
      cycle
    ENDIF

    IF( mm < s ) then
      go to 1000
    ENDIF

    IF( k <= uk ) then
      go to 200
    ENDIF
!
!  Check for possible splitting.
!
     DO uk = k, n - 1

       IF( ar(uk+1,uk) == 0.0D+00 .and. ai(uk+1,uk) == 0.0D+00 ) then
         exit
       ENDIF

     ENDDO
!
!  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
!
     norm = 0.0D+00
     mp = 1

     DO i = 1, uk

       x = 0.0D+00
       DO j = mp, uk
         x = x + pythag ( ar(i,j), ai(i,j) )
       ENDDO

       norm = max ( norm, x )
       mp = i

     ENDDO
!
!  EPS3 replaces zero pivot in decomposition
!  and close roots are modified by EPS3.
!
     IF( norm == 0.0D+00 ) then
       norm = 1.0D+00
     ENDIF

     eps3 = abs ( norm ) * epsilon ( eps3 )
!
!  GROWTO is the criterion for growth.
!
     ukroot = uk
     ukroot = sqrt ( ukroot )
     growto = 0.1D+00 / ukroot

200  CONTINUE

     rlambd = wr(k)
     ilambd = wi(k)

     IF( k == 1 ) then
       go to 280
     ENDIF

     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue.
!
220  CONTINUE

     rlambd = rlambd + eps3

240  CONTINUE

     DO ii = 1, km1
        i = k - ii
        IF( select(i) .and. abs ( wr(i)-rlambd) < eps3 .and. &
            abs ( wi(i)-ilambd) < eps3 ) then
          go to 220
        ENDIF
     ENDDO

     wr(k) = rlambd
!
!  Form upper Hessenberg (ar,ai)-(rlambd,ilambd) * I
!  and initial complex vector.
!
280  CONTINUE

     mp = 1

     DO i = 1, uk

        DO j = mp, uk
          rm1(i,j) = ar(i,j)
          rm2(i,j) = ai(i,j)
        ENDDO

        rm1(i,i) = rm1(i,i) - rlambd
        rm2(i,i) = rm2(i,i) - ilambd
        mp = i
        rv1(i) = eps3

     ENDDO
!
!  Triangular decomposition with interchanges, replacing zero pivots by eps3.
!
     DO i = 2, uk

        mp = i - 1

        IF( pythag ( rm1(i,mp), rm2(i,mp) ) > &
             pythag ( rm1(mp,mp),rm2(mp,mp) ) ) then

          DO j = mp, uk
            call r8_swap ( rm1(i,j), rm1(mp,j) )
            call r8_swap ( rm2(i,j), rm2(mp,j) )
          ENDDO

        ENDIF

        IF( rm1(mp,mp) == 0.0D+00 .and. rm2(mp,mp) == 0.0D+00 ) then
          rm1(mp,mp) = eps3
        ENDIF

        call cdiv ( rm1(i,mp), rm2(i,mp), rm1(mp,mp), rm2(mp,mp), x, y )

        IF( x /= 0.0D+00 .or. y /= 0.0D+00 ) then

          DO j = i, uk
            rm1(i,j) = rm1(i,j) - x * rm1(mp,j) + y * rm2(mp,j)
            rm2(i,j) = rm2(i,j) - x * rm2(mp,j) - y * rm1(mp,j)
          ENDDO

        ENDIF

     ENDDO

     IF( rm1(uk,uk) == 0.0D+00 .and. rm2(uk,uk) == 0.0D+00 ) then
       rm1(uk,uk) = eps3
     ENDIF

     its = 0
!
!  Back substitution.
!
  660   CONTINUE

    DO ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = 0.0D+00

        DO j = i + 1, uk
          x = x - rm1(i,j) * rv1(j) + rm2(i,j) * rv2(j)
          y = y - rm1(i,j) * rv2(j) - rm2(i,j) * rv1(j)
        ENDDO

        call cdiv ( x, y, rm1(i,i), rm2(i,i), rv1(i), rv2(i) )

     ENDDO
!
!  Acceptance test for eigenvector and normalization.
!
     its = its + 1
     norm = 0.0D+00
     normv = 0.0D+00

     DO i = 1, uk
        x = pythag ( rv1(i), rv2(i) )
        IF( normv < x ) then
          normv = x
          j = i
        ENDIF
        norm = norm + x
     ENDDO

     IF( norm < growto ) then
       go to 840
     ENDIF
!
!  Accept vector.
!
     x = rv1(j)
     y = rv2(j)

     DO i = 1, uk
       call cdiv ( rv1(i), rv2(i), x, y, zr(i,s), zi(i,s) )
     ENDDO

     IF( uk == n ) then
       go to 940
     ENDIF

     j = uk + 1
     go to 900
!
!  Choose a new starting vector.
!
  840    CONTINUE

     IF( its < uk ) then

       x = ukroot
       y = eps3 / ( x + 1.0D+00 )

       rv1(1) = eps3
       rv1(2:uk) = y

       j = uk - its + 1
       rv1(j) = rv1(j) - eps3 * x
       go to 660

     ENDIF
!
!  Error: unaccepted eigenvector.
!
!  880    CONTINUE

     j = 1
     ierr = -k
!
!  Set remaining vector components to zero.
!
900    CONTINUE

       zr(j:n,s) = 0.0D+00
       zi(j:n,s) = 0.0D+00

940    CONTINUE

       s = s + 1

  ENDDO

  go to 1001
!
!  Set error: underestimate of eigenvector space required.
!
 1000 CONTINUE

  IF( ierr /= 0 ) then
    ierr = ierr - n
  ENDIF

  IF( ierr == 0 ) then
    ierr = - ( 2 * n + 1 )
  ENDIF

 1001 CONTINUE
  m = s - 1
  RETURN

END SUBROUTINE

SUBROUTINE combak ( n, low, igh, ar, ai, int, m, zr, zi )

!*****************************************************************************80
!
!! COMBAK determines eigenvectors by undoing the COMHES transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by COMHES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = to the order
!    of the matrix.
!
!    Input, REAL ( kind = 8 ) AR(N,IGH), AI(N,IGH), the multipliers which
!    were used in the reduction by COMHES in their lower triangles below
!    the subdiagonal.
!
!    Input, INTEGER ( kind = 4 ) INT(IGH), information on the rows and
!    columns interchanged in the reduction by COMHES.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) ZR(N,M), ZI(N,M).  On input, the real
!    and imaginary parts of the eigenvectors to be back transformed.  On
!    output, the real and imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,igh)
  REAL ( kind = 8 ) ar(n,igh)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) int(igh)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) zi(n,m)
  REAL ( kind = 8 ) zr(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  la = igh - 1

  IF( igh - 1 < low + 1 ) then
    RETURN
  ENDIF

  DO mm = low + 1, la

     mp = low + igh - mm

     DO i = mp + 1, igh

        xr = ar(i,mp-1)
        xi = ai(i,mp-1)

        IF( xr /= 0.0D+00 .or. xi /= 0.0D+00 ) then
          zr(i,1:m) = zr(i,1:m) + xr * zr(mp,1:m) - xi * zi(mp,1:m)
          zi(i,1:m) = zi(i,1:m) + xr * zi(mp,1:m) + xi * zr(mp,1:m)
       ENDIF

     ENDDO

     i = int(mp)

     IF( i /= mp ) then

       DO j = 1, m
         call r8_swap ( zr(i,j), zr(mp,j) )
         call r8_swap ( zi(i,j), zi(mp,j) )
       ENDDO

     ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE comhes ( n, low, igh, ar, ai, int )

!*****************************************************************************80
!
!! COMHES transforms a complex general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a complex general matrix, this SUBROUTINE
!    reduces a submatrix situated in rows and columns
!    LOW through IGH to upper Hessenberg form by
!    stabilized elementary similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real and
!    imaginary parts of the complex input matrix.  On output, the real and
!    imaginary parts of the Hessenberg matrix.  The multipliers which were
!    used in the reduction are stored in the remaining triangles under the
!    Hessenberg matrix.
!
!    Output, INTEGER ( kind = 4 ) INT(IGH), information on the rows and columns
!    interchanged in the reduction.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) int(igh)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr

  la = igh - 1

  DO m = low + 1, la

     xr = 0.0D+00
     xi = 0.0D+00
     i = m

     DO j = m, igh

       IF( abs ( ar(j,m-1) ) + abs ( ai(j,m-1) ) > &
         abs ( xr ) + abs ( xi ) ) then
         xr = ar(j,m-1)
         xi = ai(j,m-1)
         i = j
       ENDIF

     ENDDO

     int(m) = i
!
!  Interchange rows and columns of AR and AI.
!
     IF( i /= m ) then

       DO j = m - 1, n
         call r8_swap ( ar(i,j), ar(m,j) )
         call r8_swap ( ai(i,j), ai(m,j) )
       ENDDO

       DO j = 1, igh
         call r8_swap ( ar(j,i), ar(j,m) )
         call r8_swap ( ai(j,i), ai(j,m) )
       ENDDO

     ENDIF

    IF( xr /= 0.0D+00 .or. xi /= 0.0D+00 ) then

      DO i = m + 1, igh

        yr = ar(i,m-1)
        yi = ai(i,m-1)

        IF( yr /= 0.0D+00 .or. yi /= 0.0D+00 ) then

          call cdiv ( yr, yi, xr, xi, yr, yi )
          ar(i,m-1) = yr
          ai(i,m-1) = yi

          DO j = m, n
            ar(i,j) = ar(i,j) - yr * ar(m,j) + yi * ai(m,j)
            ai(i,j) = ai(i,j) - yr * ai(m,j) - yi * ar(m,j)
          ENDDO

          ar(1:igh,m) = ar(1:igh,m) + yr * ar(1:igh,i) - yi * ai(1:igh,i)
          ai(1:igh,m) = ai(1:igh,m) + yr * ai(1:igh,i) + yi * ar(1:igh,i)

        ENDIF

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE comlr ( n, low, igh, hr, hi, wr, wi, ierr )

!*****************************************************************************80
!
!! COMLR gets all eigenvalues of a complex upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a complex upper Hessenberg
!    matrix by the modified LR method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) HR(N,N), HI(N,N).  On input, the real and
!    imaginary parts of the complex upper Hessenberg matrix.  Their lower
!    triangles below the subdiagonal contain the multipliers which were used
!    in the reduction by COMHES if performed.  On output, the upper Hessenberg
!    portions of HR and HI have been destroyed.  Therefore, they must be
!    saved before calling COMLR if subsequent calculation of eigenvectors
!    is to be performed.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  If an error exit is made, the eigenvalues should be correct
!    for indices IERR+1,...,N.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!      eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm1
  REAL ( kind = 8 ) hi(n,n)
  REAL ( kind = 8 ) hr(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mm
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) sr
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr
  REAL ( kind = 8 ) zzi
  REAL ( kind = 8 ) zzr

  ierr = 0
!
!  Store roots isolated by CBAL.
!
  DO i = 1, n
    IF( i < low .or. igh < i ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    ENDIF
  ENDDO

  en = igh
  tr = 0.0D+00
  ti = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalue.
!
  220 CONTINUE

  IF( en < low ) then
    RETURN
  ENDIF

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 CONTINUE

  DO ll = low, en

     l = en + low - ll

     IF( l == low ) then
       exit
     ENDIF

     tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
       + abs ( hi(l,l) )
     tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )

     IF( tst2 == tst1 ) then
       exit
     ENDIF

  ENDDO
!
!  Form shift.
!
!300 CONTINUE

  IF( l == en ) then
    go to 660
  ENDIF

  IF( itn == 0 ) then
    ierr = en
    RETURN
  ENDIF

  IF( its == 10 .or. its == 20 ) then
    go to 320
  ENDIF

  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)

  IF( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
    go to 340
  ENDIF

  yr = ( hr(enm1,enm1) - sr) / 2.0D+00
  yi = ( hi(enm1,enm1) - si) / 2.0D+00
  call csroot ( yr**2-yi**2+xr, 2.0D+00*yr*yi+xi, zzr, zzi )

  IF( yr * zzr + yi * zzi < 0.0D+00 ) then
    zzr = -zzr
    zzi = -zzi
  ENDIF

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
  320 CONTINUE

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = abs ( hi(en,enm1) ) + abs ( hi(enm1,en-2) )

  340 CONTINUE

  DO i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  ENDDO

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  xr = abs ( hr(enm1,enm1) ) + abs ( hi(enm1,enm1) )
  yr = abs ( hr(en,enm1) ) + abs ( hi(en,enm1) )
  zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

  DO mm = l, enm1
    m = enm1 + l - mm
    IF( m == l ) then
      exit
    ENDIF
    yi = yr
    yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
    xi = zzr
    zzr = xr
    xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
    tst1 = zzr / yi * (zzr + xr + xi)
    tst2 = tst1 + yr
    IF( tst2 == tst1 ) then
      exit
    ENDIF
  ENDDO
!
!  Triangular decomposition H=L*R.
!
  DO i = m + 1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)

     IF( abs ( xr ) + abs ( xi ) >= abs ( yr ) + abs ( yi ) ) then
       go to 460
     ENDIF
!
!  Interchange rows of HR and HI.
!
     DO j = i - 1, en
       call r8_swap ( hr(i-1,j), hr(i,j) )
       call r8_swap ( hi(i-1,j), hi(i,j) )
     ENDDO

     call cdiv ( xr, xi, yr, yi, zzr, zzi )
     wr(i) = 1.0D+00
     go to 480

460 CONTINUE

     call cdiv ( yr, yi, xr, xi, zzr, zzi )
     wr(i) = -1.0D+00

480  CONTINUE

     hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     DO j = i, en
        hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
        hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     ENDDO

  ENDDO
!
!  Composition R*L=H.
!
  DO j = m + 1, en

    xr = hr(j,j-1)
    xi = hi(j,j-1)
    hr(j,j-1) = 0.0D+00
    hi(j,j-1) = 0.0D+00
!
!  Interchange columns of HR and HI, if necessary.
!
    IF( 0.0D+00 < wr(j) ) then

      DO i = l, j
        call r8_swap ( hr(i,j-1), hr(i,j) )
        call r8_swap ( hi(i,j-1), hi(i,j) )
      ENDDO

    ENDIF

    DO i = l, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    ENDDO

  ENDDO

  go to 240
!
!  A root found.
!
  660 CONTINUE

  wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220

END SUBROUTINE

SUBROUTINE comlr2 ( n, low, igh, int, hr, hi, wr, wi, zr, zi, ierr )

!*****************************************************************************80
!
!! COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues and eigenvectors of a complex
!    upper Hessenberg matrix by the modified LR method.  The eigenvectors
!    of a complex general matrix can also be found if COMHES has been used
!    to reduce this general matrix to Hessenberg form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input, INTEGER ( kind = 4 ) INT(IGH), information on the rows and columns
!    interchanged in the reduction by COMHES, if performed.  If the
!    eigenvectors of the Hessenberg matrix are desired, set INT(J)=J for these
!    elements.
!
!    Input/output, REAL ( kind = 8 ) HR(N,N), HI(N,N).  On input, the real
!    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
!    triangles below the subdiagonal contain the multipliers which were used in
!    the reduction by COMHES, if performed.  If the eigenvectors of the
!    Hessenberg matrix are desired, these elements must be set to zero.  On
!    output, the upper Hessenberg portions of HR and HI have been destroyed,
!    but the location HR(1,1) contains the norm of the triangularized matrix.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  If an error exit is made, the eigenvalues should be
!    correct for indices IERR+1,...,N.
!
!    Output, REAL ( kind = 8 ) ZR(N,N), ZI(N,N), the real and imaginary parts
!    of the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
!    is made, none of the eigenvectors has been found.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!      eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm1
  REAL ( kind = 8 ) hi(n,n)
  REAL ( kind = 8 ) hr(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) iend
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) int(igh)
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) sr
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr
  REAL ( kind = 8 ) zi(n,n)
  REAL ( kind = 8 ) zr(n,n)
  REAL ( kind = 8 ) zzi
  REAL ( kind = 8 ) zzr

  ierr = 0
!
!  Initialize the eigenvector matrix.
!
  zr(1:n,1:n) = 0.0D+00

  DO i = 1, n
    zr(i,i) = 1.0D+00
  ENDDO

  zi(1:n,1:n) = 0.0D+00
!
!  Form the matrix of accumulated transformations from the information left
!  by COMHES.
!
  iend = igh - low - 1

  DO ii = 1, iend

    i = igh - ii

    DO k = i + 1, igh
      zr(k,i) = hr(k,i-1)
      zi(k,i) = hi(k,i-1)
    ENDDO

    j = int(i)

    IF( i /= j ) then

      DO k = i, igh
        zr(i,k) = zr(j,k)
        zi(i,k) = zi(j,k)
        zr(j,k) = 0.0D+00
        zi(j,k) = 0.0D+00
      ENDDO

      zr(j,i) = 1.0D+00

    ENDIF

  ENDDO
!
!  Store roots isolated by CBAL.
!
  DO i = 1, n
    IF( i < low .or. igh < i ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    ENDIF
  ENDDO

  en = igh
  tr = 0.0D+00
  ti = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalue.
!
220 CONTINUE

  IF( en < low ) then
    go to 680
  ENDIF

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 CONTINUE

  DO ll = low, en

     l = en + low - ll

     IF( l == low ) then
       exit
     ENDIF

     tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
       + abs ( hi(l,l) )
     tst2 = tst1 + abs ( hr(l,l-1) ) + abs ( hi(l,l-1) )

     IF( tst2 == tst1 ) then
       exit
     ENDIF

  ENDDO
!
!  Form shift.
!
  IF( l == en ) then
    go to 660
  ENDIF

  IF( itn == 0 ) then
    ierr = en
    RETURN
  ENDIF

  IF( its == 10 .or. its == 20 ) then
    go to 320
  ENDIF

  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1) - hi(enm1,en) * hi(en,enm1)
  xi = hr(enm1,en) * hi(en,enm1) + hi(enm1,en) * hr(en,enm1)

  IF( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
    go to 340
  ENDIF

  yr = (hr(enm1,enm1) - sr) / 2.0D+00
  yi = (hi(enm1,enm1) - si) / 2.0D+00
  call csroot ( yr**2-yi**2+xr, 2.0D+00*yr*yi+xi, zzr, zzi )

  IF( yr * zzr + yi * zzi < 0.0D+00 ) then
    zzr = -zzr
    zzi = -zzi
  ENDIF

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
  320 CONTINUE

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = abs ( hi(en,enm1) ) + abs ( hi(enm1,en-2) )

  340 CONTINUE

  DO i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  ENDDO

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  xr = abs ( hr(enm1,enm1) ) + abs ( hi(enm1,enm1) )
  yr = abs ( hr(en,enm1) ) + abs ( hi(en,enm1) )
  zzr = abs ( hr(en,en) ) + abs ( hi(en,en) )

  DO mm = l, enm1
     m = enm1 + l - mm
     IF( m == l ) then
       exit
     ENDIF
     yi = yr
     yr = abs ( hr(m,m-1) ) + abs ( hi(m,m-1) )
     xi = zzr
     zzr = xr
     xr = abs ( hr(m-1,m-1) ) + abs ( hi(m-1,m-1) )
     tst1 = zzr / yi * (zzr + xr + xi)
     tst2 = tst1 + yr
     IF( tst2 == tst1 ) then
       exit
     ENDIF
  ENDDO
!
!  Triangular decomposition H=L*R.
!
  DO i = m + 1, en

     xr = hr(i-1,i-1)
     xi = hi(i-1,i-1)
     yr = hr(i,i-1)
     yi = hi(i,i-1)
!
!  Interchange rows of HR and HI.
!
     IF( abs ( xr ) + abs ( xi) < abs ( yr ) + abs ( yi ) ) then

       DO j = i - 1, n
         call r8_swap ( hr(i-1,j), hr(i,j) )
         call r8_swap ( hi(i-1,j), hi(i,j) )
       ENDDO

       call cdiv ( xr, xi, yr, yi, zzr, zzi )
       wr(i) = 1.0D+00

     else

       call cdiv ( yr, yi, xr, xi, zzr, zzi )
       wr(i) = -1.0D+00

     ENDIF

     hr(i,i-1) = zzr
     hi(i,i-1) = zzi

     DO j = i, n
       hr(i,j) = hr(i,j) - zzr * hr(i-1,j) + zzi * hi(i-1,j)
       hi(i,j) = hi(i,j) - zzr * hi(i-1,j) - zzi * hr(i-1,j)
     ENDDO

  ENDDO
!
!  Composition R*L=H.
!
  DO j = m + 1, en

     xr = hr(j,j-1)
     xi = hi(j,j-1)
     hr(j,j-1) = 0.0D+00
     hi(j,j-1) = 0.0D+00
!
!  Interchange columns of HR, HI, ZR, and ZI.
!
     IF( 0.0D+00 < wr(j) ) then

       DO i = 1, j
         call r8_swap ( hr(i,j-1), hr(i,j) )
         call r8_swap ( hi(i,j-1), hi(i,j) )
       ENDDO

       DO i = low, igh
         call r8_swap ( zr(i,j-1), zr(i,j) )
         call r8_swap ( zi(i,j-1), zi(i,j) )
       ENDDO

    ENDIF

    DO i = 1, j
      hr(i,j-1) = hr(i,j-1) + xr * hr(i,j) - xi * hi(i,j)
      hi(i,j-1) = hi(i,j-1) + xr * hi(i,j) + xi * hr(i,j)
    ENDDO
!
!  Accumulate transformations.
!
    DO i = low, igh
      zr(i,j-1) = zr(i,j-1) + xr * zr(i,j) - xi * zi(i,j)
      zi(i,j-1) = zi(i,j-1) + xr * zi(i,j) + xi * zr(i,j)
    ENDDO

  ENDDO

  go to 240
!
!  A root found.
!
  660 CONTINUE

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.
!  Backsubstitute to find vectors of upper triangular form.
!
  680 CONTINUE

  norm = 0.0D+00

  DO i = 1, n
    DO j = i, n
      tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
      norm = max ( norm, tr )
    ENDDO
  ENDDO

  hr(1,1) = norm
  IF( n == 1 ) then
    RETURN
  ENDIF

  IF( norm == 0.0D+00 ) then
    RETURN
  ENDIF

  DO nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = 1.0D+00
     hi(en,en) = 0.0D+00
     enm1 = en - 1

     DO ii = 1, enm1

        i = en - ii
        zzr = 0.0D+00
        zzi = 0.0D+00

        DO j = i + 1, en
          zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
          zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        ENDDO

        yr = xr - wr(i)
        yi = xi - wi(i)

        IF( yr == 0.0D+00 .and. yi == 0.0D+00 ) then

          tst1 = norm
          yr = tst1

          do
            yr = 0.01D+00 * yr
            tst2 = norm + yr
            IF( tst2 <=  tst1 ) then
              exit
            ENDIF
          ENDDO

        ENDIF

        call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
!
!  Overflow control.
!
        tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

        IF( tr /= 0.0D+00 ) then

          tst1 = tr
          tst2 = tst1 + 1.0D+00 / tst1

          IF( tst2 <= tst1 ) then

            hr(i:en,en) = hr(i:en,en) / tr
            hi(i:en,en) = hi(i:en,en) / tr

          ENDIF

        ENDIF

      ENDDO

  ENDDO
!
!  End backsubstitution.
!
  enm1 = n - 1
!
!  Vectors of isolated roots.
!
  DO i = 1, n - 1

    IF( i < low .or. igh < i ) then

      zr(i,i+1:n) = hr(i,i+1:n)
      zi(i,i+1:n) = hi(i,i+1:n)

    ENDIF

  ENDDO
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  DO jj = low, n - 1

    j = n + low - jj
    m = min ( j, igh )

    DO i = low, igh
      zzr = 0.0D+00
      zzi = 0.0D+00
      DO k = low, m
        zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
        zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
      ENDDO
      zr(i,j) = zzr
      zi(i,j) = zzi
    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE comqr ( n, low, igh, hr, hi, wr, wi, ierr )

!*****************************************************************************80
!
!! COMQR gets eigenvalues of a complex upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a complex
!    upper Hessenberg matrix by the QR method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) HR(N,N), HI(N,N).  On input, the real
!    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
!    triangles below the subdiagonal contain information about the unitary
!    transformations used in the reduction by CORTH, if performed.  On output,
!    the upper Hessenberg portions of HR and HI have been destroyed.
!    Therefore, they must be saved before calling COMQR if subsequent
!    calculation of eigenvectors is to be performed.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  If an error exit is made, the eigenvalues should be
!    correct for indices IERR+1,...,N.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!       eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm1
  REAL ( kind = 8 ) hi(n,n)
  REAL ( kind = 8 ) hr(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  REAL ( kind = 8 ) norm
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) sr
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr
  REAL ( kind = 8 ) zzi
  REAL ( kind = 8 ) zzr

  ierr = 0
!
!  Create real subdiagonal elements.
!
  l = low + 1

  DO i = l, igh

     ll = min ( i + 1, igh )

     IF( hi(i,i-1) /= 0.0D+00 ) then

     norm = pythag ( hr(i,i-1), hi(i,i-1) )
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = 0.0D+00

     DO j = i, igh
       si = yr * hi(i,j) - yi * hr(i,j)
       hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
       hi(i,j) = si
     ENDDO

     DO j = low, ll
       si = yr * hi(j,i) + yi * hr(j,i)
       hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
       hi(j,i) = si
     ENDDO

    ENDIF

  ENDDO
!
!  Store roots isolated by CBAL.
!
  DO i = 1, n
    IF( i < low .or. igh < i ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    ENDIF
  ENDDO

  en = igh
  tr = 0.0D+00
  ti = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalue.
!
  220 CONTINUE

  IF( en < low ) then
    RETURN
  ENDIF

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
  240 CONTINUE

  DO ll = low, en
    l = en + low - ll
    IF( l == low ) then
      exit
    ENDIF
    tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
      + abs ( hi(l,l) )
    tst2 = tst1 + abs ( hr(l,l-1) )
    IF( tst2 == tst1 ) then
      exit
    ENDIF
  ENDDO
!
!  Form shift.
!
  IF( l == en ) then
    go to 660
  ENDIF

  IF( itn == 0 ) then
    go to 1000
  ENDIF

  IF( its == 10 .or. its == 20 ) then
    go to 320
  ENDIF

  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  IF( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
    go to 340
  ENDIF

  yr = (hr(enm1,enm1) - sr) / 2.0D+00
  yi = (hi(enm1,enm1) - si) / 2.0D+00

  call csroot ( yr**2-yi**2+xr, 2.0D+00*yr*yi+xi, zzr, zzi )

  IF( yr * zzr + yi * zzi < 0.0D+00 ) then
    zzr = -zzr
    zzi = -zzi
  ENDIF

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
320 CONTINUE

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = 0.0D+00

340 CONTINUE

  DO i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  ENDDO

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows).
!
  DO i = l + 1, en

     sr = hr(i,i-1)
     hr(i,i-1) = 0.0D+00
     norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = 0.0D+00
     hi(i,i-1) = sr / norm

     DO j = i, en
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
    ENDDO

  ENDDO

  si = hi(en,en)

  IF( si /= 0.0D+00 ) then
    norm = pythag ( hr(en,en), si )
    sr = hr(en,en) / norm
    si = si / norm
    hr(en,en) = norm
    hi(en,en) = 0.0D+00
  ENDIF
!
!  Inverse operation (columns).
!
  DO j = l + 1, en

     xr = wr(j-1)
     xi = wi(j-1)

     DO i = l, j

        yr = hr(i,j-1)
        yi = 0.0D+00
        zzr = hr(i,j)
        zzi = hi(i,j)
        IF( i /= j ) then
          yi = hi(i,j-1)
          hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
        ENDIF
        hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
        hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
        hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     ENDDO

  ENDDO

  IF( si /= 0.0D+00 ) then

    DO i = l, en
      yr = hr(i,en)
      yi = hi(i,en)
      hr(i,en) = sr * yr - si * yi
      hi(i,en) = sr * yi + si * yr
    ENDDO

  ENDIF

  go to 240
!
!  A root found.
!
660 CONTINUE

  wr(en) = hr(en,en) + tr
  wi(en) = hi(en,en) + ti
  en = enm1
  go to 220
!
!  Set error: all eigenvalues have not converged after 30*n iterations.
!
1000 CONTINUE

  ierr = en
  RETURN

END SUBROUTINE

SUBROUTINE comqr2 ( n, low, igh, ortr, orti, hr, hi, wr, wi, zr, zi, ierr )

!*****************************************************************************80
!
!! COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues and eigenvectors
!    of a complex upper Hessenberg matrix by the QR
!    method.  The eigenvectors of a complex general matrix
!    can also be found if CORTH has been used to reduce
!    this general matrix to Hessenberg form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) ORTR(N), ORTI(N).  On input, information 
!    about the unitary transformations used in the reduction by CORTH, if 
!    performed.  If the eigenvectors of the Hessenberg matrix are desired, set 
!    ORTR(J) and ORTI(J) to 0.0D+00 for these elements.  On output, these arrays
!    have been overwritten.
!
!    Input/output, REAL ( kind = 8 ) HR(N,N), HI(N,N).  On input, the real and 
!    imaginary parts of the complex upper Hessenberg matrix.  Their lower 
!    triangles below the subdiagonal contain further information about the
!    transformations which were used in the reduction by CORTH, if performed.
!    If the eigenvectors of the Hessenberg matrix are desired, these elements
!    may be arbitrary.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  If an error exit is made, the eigenvalues should be
!    correct for indices IERR+1,...,N.
!
!    Output, REAL ( kind = 8 ) ZR(N,N), ZI(N,N), the real and imaginary parts of
!    the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
!    is made, none of the eigenvectors has been found.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!      eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm1
  REAL ( kind = 8 ) hi(n,n)
  REAL ( kind = 8 ) hr(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) iend
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) orti(igh)
  REAL ( kind = 8 ) ortr(igh)
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) sr
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr
  REAL ( kind = 8 ) zi(n,n)
  REAL ( kind = 8 ) zr(n,n)
  REAL ( kind = 8 ) zzi
  REAL ( kind = 8 ) zzr

  ierr = 0
!
!  Initialize eigenvector matrix.
!
  zr(1:n,1:n) = 0.0D+00

  DO i = 1, n
    zr(i,i) = 1.0D+00
  ENDDO

  zi(1:n,1:n) = 0.0D+00
!
!  Form the matrix of accumulated transformations from the information
!  left by CORTH.
!
  iend = igh - low - 1
!  IF( iend ) 180, 150, 105 
  IF( iend < 0 ) THEN 
    GO TO 180 
  ELSE IF( iend == 0 ) THEN
    GO TO 150
  ELSE 
    GO TO 105
  ENDIF

105 CONTINUE

  DO ii = 1, iend

     i = igh - ii

     IF( ortr(i) == 0.0D+00 .and. orti(i) == 0.0D+00 ) then
       go to 140
     ENDIF

     IF( hr(i,i-1) == 0.0D+00 .and. hi(i,i-1) == 0.0D+00 ) then
       go to 140
     ENDIF
!
!  Norm below is negative of H formed in CORTH.
!
     norm = hr(i,i-1) * ortr(i) + hi(i,i-1) * orti(i)

     DO k = i + 1, igh
       ortr(k) = hr(k,i-1)
       orti(k) = hi(k,i-1)
     ENDDO

     DO j = i, igh

        sr = 0.0D+00
        si = 0.0D+00

        DO k = i, igh
          sr = sr + ortr(k) * zr(k,j) + orti(k) * zi(k,j)
          si = si + ortr(k) * zi(k,j) - orti(k) * zr(k,j)
        ENDDO

        sr = sr / norm
        si = si / norm

        DO k = i, igh
          zr(k,j) = zr(k,j) + sr * ortr(k) - si * orti(k)
          zi(k,j) = zi(k,j) + sr * orti(k) + si * ortr(k)
        ENDDO

      ENDDO

140 CONTINUE

  ENDDO
!
!  Create real subdiagonal elements.
!
150 CONTINUE

  l = low + 1

  DO i = l, igh

     ll = min ( i + 1, igh )

     IF( hi(i,i-1) == 0.0D+00 ) then
       cycle
     ENDIF

     norm = pythag ( hr(i,i-1), hi(i,i-1) )
     yr = hr(i,i-1) / norm
     yi = hi(i,i-1) / norm
     hr(i,i-1) = norm
     hi(i,i-1) = 0.0D+00

     DO j = i, n
       si = yr * hi(i,j) - yi * hr(i,j)
       hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
       hi(i,j) = si
     ENDDO

     DO j = 1, ll
       si = yr * hi(j,i) + yi * hr(j,i)
       hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
       hi(j,i) = si
     ENDDO

     DO j = low, igh
       si = yr * zi(j,i) + yi * zr(j,i)
       zr(j,i) = yr * zr(j,i) - yi * zi(j,i)
       zi(j,i) = si
     ENDDO

  ENDDO
!
!  Store roots isolated by CBAL.
!
180 CONTINUE

  DO i = 1, n
    IF( i < low .or. igh < i ) then
      wr(i) = hr(i,i)
      wi(i) = hi(i,i)
    ENDIF
  ENDDO

  en = igh
  tr = 0.0D+00
  ti = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalue.
!
220 CONTINUE

  IF( en < low ) then
    go to 680
  ENDIF

  its = 0
  enm1 = en - 1
!
!  Look for single small sub-diagonal element.
!
240 CONTINUE

  DO ll = low, en
    l = en + low - ll
    IF( l == low ) then
      exit
    ENDIF
    tst1 = abs ( hr(l-1,l-1) ) + abs ( hi(l-1,l-1) ) + abs ( hr(l,l) ) &
      + abs ( hi(l,l) )
    tst2 = tst1 + abs ( hr(l,l-1) )
    IF( tst2 == tst1 ) then
      exit
    ENDIF
  ENDDO
!
!  Form shift.
!
  IF( l == en ) then
    go to 660
  ENDIF

  IF( itn == 0 ) then
    go to 1000
  ENDIF

  IF( its == 10 .or. its == 20 ) then
    go to 320
  ENDIF

  sr = hr(en,en)
  si = hi(en,en)
  xr = hr(enm1,en) * hr(en,enm1)
  xi = hi(enm1,en) * hr(en,enm1)
  IF( xr == 0.0D+00 .and. xi == 0.0D+00 ) then
    go to 340
  ENDIF

  yr = ( hr(enm1,enm1) - sr ) / 2.0D+00
  yi = ( hi(enm1,enm1) - si ) / 2.0D+00

  call csroot ( yr**2-yi**2+xr, 2.0D+00*yr*yi+xi, zzr, zzi )

  IF( yr * zzr + yi * zzi < 0.0D+00 ) then
    zzr = -zzr
    zzi = -zzi
  ENDIF

  call cdiv ( xr, xi, yr+zzr, yi+zzi, xr, xi )
  sr = sr - xr
  si = si - xi
  go to 340
!
!  Form exceptional shift.
!
320 CONTINUE

  sr = abs ( hr(en,enm1) ) + abs ( hr(enm1,en-2) )
  si = 0.0D+00

340 CONTINUE

  DO i = low, en
    hr(i,i) = hr(i,i) - sr
    hi(i,i) = hi(i,i) - si
  ENDDO

  tr = tr + sr
  ti = ti + si
  its = its + 1
  itn = itn - 1
!
!  Reduce to triangle (rows).
!
  DO i = l + 1, en

     sr = hr(i,i-1)
     hr(i,i-1) = 0.0D+00
     norm = pythag ( pythag ( hr(i-1,i-1), hi(i-1,i-1) ), sr )
     xr = hr(i-1,i-1) / norm
     wr(i-1) = xr
     xi = hi(i-1,i-1) / norm
     wi(i-1) = xi
     hr(i-1,i-1) = norm
     hi(i-1,i-1) = 0.0D+00
     hi(i,i-1) = sr / norm

     DO j = i, n
        yr = hr(i-1,j)
        yi = hi(i-1,j)
        zzr = hr(i,j)
        zzi = hi(i,j)
        hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
        hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
        hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
        hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
     ENDDO

  ENDDO

  si = hi(en,en)

  IF( si /= 0.0D+00 ) then

    norm = pythag ( hr(en,en), si )
    sr = hr(en,en) / norm
    si = si / norm
    hr(en,en) = norm
    hi(en,en) = 0.0D+00

    DO j = en + 1, n
      yr = hr(en,j)
      yi = hi(en,j)
      hr(en,j) = sr * yr + si * yi
      hi(en,j) = sr * yi - si * yr
    ENDDO

  ENDIF
!
!  Inverse operation (columns).
!
  DO j = l + 1, en

     xr = wr(j-1)
     xi = wi(j-1)

     DO i = 1, j

       yr = hr(i,j-1)
       yi = 0.0D+00
       zzr = hr(i,j)
       zzi = hi(i,j)

       IF( i /= j ) then
         yi = hi(i,j-1)
         hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
       ENDIF

       hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
       hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
       hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi

     ENDDO

     DO i = low, igh
       yr = zr(i,j-1)
       yi = zi(i,j-1)
       zzr = zr(i,j)
       zzi = zi(i,j)
       zr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
       zi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
       zr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
       zi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
     ENDDO

  ENDDO

  IF( si /= 0.0D+00 ) then

    DO i = 1, en
      yr = hr(i,en)
      yi = hi(i,en)
      hr(i,en) = sr * yr - si * yi
      hi(i,en) = sr * yi + si * yr
    ENDDO

    DO i = low, igh
      yr = zr(i,en)
      yi = zi(i,en)
      zr(i,en) = sr * yr - si * yi
      zi(i,en) = sr * yi + si * yr
    ENDDO

  ENDIF

  go to 240
!
!  A root found.
!
660 CONTINUE

  hr(en,en) = hr(en,en) + tr
  wr(en) = hr(en,en)
  hi(en,en) = hi(en,en) + ti
  wi(en) = hi(en,en)
  en = enm1
  go to 220
!
!  All roots found.
!  Backsubstitute to find vectors of upper triangular form.
!
680 CONTINUE

  norm = 0.0D+00

  DO i = 1, n
    DO j = i, n
      tr = abs ( hr(i,j) ) + abs ( hi(i,j) )
      norm = max ( norm, tr )
    ENDDO
  ENDDO

  IF( n == 1 ) then
    RETURN
  ENDIF

  IF( norm == 0.0D+00 ) then
    RETURN
  ENDIF

  DO nn = 2, n

     en = n + 2 - nn
     xr = wr(en)
     xi = wi(en)
     hr(en,en) = 1.0D+00
     hi(en,en) = 0.0D+00
     enm1 = en - 1

     DO ii = 1, enm1

        i = en - ii
        zzr = 0.0D+00
        zzi = 0.0D+00

        DO j = i + 1, en
          zzr = zzr + hr(i,j) * hr(j,en) - hi(i,j) * hi(j,en)
          zzi = zzi + hr(i,j) * hi(j,en) + hi(i,j) * hr(j,en)
        ENDDO

        yr = xr - wr(i)
        yi = xi - wi(i)

        IF( yr == 0.0D+00 .and. yi == 0.0D+00 ) then

           tst1 = norm
           yr = tst1
           do
             yr = 0.01D+00 * yr
             tst2 = norm + yr
             IF( tst2 <= tst1 ) then
               exit
             ENDIF
           ENDDO

        ENDIF

        call cdiv ( zzr, zzi, yr, yi, hr(i,en), hi(i,en) )
!
!  Overflow control.
!
        tr = abs ( hr(i,en) ) + abs ( hi(i,en) )

        IF( tr /= 0.0D+00 ) then

          tst1 = tr
          tst2 = tst1 + 1.0D+00 / tst1

          IF( tst2 <= tst1 ) then

            DO j = i, en
              hr(j,en) = hr(j,en)/tr
              hi(j,en) = hi(j,en)/tr
            ENDDO

          ENDIF

       ENDIF

     ENDDO

  ENDDO
!
!  End backsubstitution.
!
  enm1 = n - 1
!
!  Vectors of isolated roots.
!
  DO i = 1, n - 1

    IF( i < low .or. igh < i ) then

      DO j = i + 1, n
        zr(i,j) = hr(i,j)
        zi(i,j) = hi(i,j)
      ENDDO

    ENDIF

  ENDDO
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  DO jj = low, n - 1

     j = n + low - jj
     m = min ( j, igh )

     DO i = low, igh

        zzr = 0.0D+00
        zzi = 0.0D+00
        DO k = low, m
          zzr = zzr + zr(i,k) * hr(k,j) - zi(i,k) * hi(k,j)
          zzi = zzi + zr(i,k) * hi(k,j) + zi(i,k) * hr(k,j)
        ENDDO

        zr(i,j) = zzr
        zi(i,j) = zzi

      ENDDO

  ENDDO

  RETURN
!
!  Set error: all eigenvalues have not converged after 30*n iterations.
!
1000 CONTINUE

  ierr = en
  RETURN

END SUBROUTINE

SUBROUTINE cortb ( n, low, igh, ar, ai, ortr, orti, m, zr, zi )

!*****************************************************************************80
!
!! CORTB determines eigenvectors by undoing the CORTH transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a complex general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by CORTH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH to the order
!    of the matrix.
!
!    Input, REAL ( kind = 8 ) AR(N,IGH), AI(N,IGH), information about the 
!    unitary transformations used in the reduction by CORTH in their strict 
!    lower triangles.
!
!    Input/output, REAL ( kind = 8 ) ORTR(IGH), ORTI(IGH).  On input, further 
!    information about the transformations used in the reduction by CORTH.  On 
!    output, ORTR and ORTI have been further altered.
!
!    Input, INTEGER ( kind = 4 ) M, the number of columns of ZR and ZI to be 
!    back transformed.
!
!    Input/output, REAL ( kind = 8 ) ZR(N,M), ZI(N,M).  On input, the real and 
!    imaginary parts of the eigenvectors to be back transformed.  On output, 
!    the real and imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,igh)
  REAL ( kind = 8 ) ar(n,igh)
  REAL ( kind = 8 ) gi
  REAL ( kind = 8 ) gr
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) orti(igh)
  REAL ( kind = 8 ) ortr(igh)
  REAL ( kind = 8 ) zi(n,m)
  REAL ( kind = 8 ) zr(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  la = igh - 1

  IF( igh - 1 < low + 1 ) then
    RETURN
  ENDIF

  DO mm = low + 1, la

    mp = low + igh - mm

    IF( ar(mp,mp-1) /= 0.0D+00 .or. ai(mp,mp-1) /= 0.0D+00 ) then

      h = ar(mp,mp-1) * ortr(mp) + ai(mp,mp-1) * orti(mp)

      ortr(mp+1:igh) = ar(mp+1:igh,mp-1)
      orti(mp+1:igh) = ai(mp+1:igh,mp-1)

      DO j = 1, m

        gr = ( dot_product ( ortr(mp:igh), zr(mp:igh,j) ) &
             + dot_product ( orti(mp:igh), zi(mp:igh,j) ) ) / h

        gi = ( dot_product ( ortr(mp:igh), zi(mp:igh,j) ) &
             - dot_product ( orti(mp:igh), zr(mp:igh,j) ) ) / h

        DO i = mp, igh
          zr(i,j) = zr(i,j) + gr * ortr(i) - gi * orti(i)
          zi(i,j) = zi(i,j) + gr * orti(i) + gi * ortr(i)
        ENDDO

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE corth ( n, low, igh, ar, ai, ortr, orti )

!*****************************************************************************80
!
!! CORTH transforms a complex general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a complex general matrix, this SUBROUTINE
!    reduces a submatrix situated in rows and columns
!    LOW through IGH to upper Hessenberg form by
!    unitary similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real and
!    imaginary parts of the complex input matrix.  On output, the real and 
!    imaginary parts of the Hessenberg matrix.  Information about the unitary
!    transformations used in the reduction is stored in the remaining
!    triangles under the Hessenberg matrix.
!
!    Output, REAL ( kind = 8 ) ORTR(IGH), ORTI(IGH), further information about 
!    the transformations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) fi
  REAL ( kind = 8 ) fr
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) m,mp,low
  REAL ( kind = 8 ) orti(igh)
  REAL ( kind = 8 ) ortr(igh)
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) scale

  la = igh - 1

  IF( igh - 1 < low + 1 ) then
    RETURN
  ENDIF

  DO m = low + 1, la

    h = 0.0D+00
    ortr(m) = 0.0D+00
    orti(m) = 0.0D+00
    scale = 0.0D+00
!
!  Scale column.
!
    DO i = m, igh
      scale = scale + abs ( ar(i,m-1) ) + abs ( ai(i,m-1) )
    ENDDO

    IF( scale == 0.0D+00 ) then
      cycle
    ENDIF

    mp = m + igh

    DO ii = m, igh
      i = mp - ii
      ortr(i) = ar(i,m-1) / scale
      orti(i) = ai(i,m-1) / scale
      h = h + ortr(i) * ortr(i) + orti(i) * orti(i)
    ENDDO

    g = sqrt ( h )
    f = pythag ( ortr(m), orti(m) )

    IF( f /= 0.0D+00 ) then
      h = h + f * g
      g = g / f
      ortr(m) = ( 1.0D+00 + g ) * ortr(m)
      orti(m) = ( 1.0D+00 + g ) * orti(m)
    else
      ortr(m) = g
      ar(m,m-1) = scale
    ENDIF
!
!  Form (I-(U*Ut)/h) * A.
!
    DO j = m, n

      fr = 0.0D+00
      fi = 0.0D+00

      DO ii = m, igh
        i = mp - ii
        fr = fr + ortr(i) * ar(i,j) + orti(i) * ai(i,j)
        fi = fi + ortr(i) * ai(i,j) - orti(i) * ar(i,j)
      ENDDO

      fr = fr / h
      fi = fi / h

      ar(m:igh,j) = ar(m:igh,j) - fr * ortr(m:igh) + fi * orti(m:igh)
      ai(m:igh,j) = ai(m:igh,j) - fr * orti(m:igh) - fi * ortr(m:igh)

    ENDDO
!
!  Form (I-(U*Ut)/h) * A * (I-(U*Ut)/h)
!
    DO i = 1, igh

      fr = 0.0D+00
      fi = 0.0D+00

      DO jj = m, igh
        j = mp - jj
        fr = fr + ortr(j) * ar(i,j) - orti(j) * ai(i,j)
        fi = fi + ortr(j) * ai(i,j) + orti(j) * ar(i,j)
      ENDDO

      fr = fr / h
      fi = fi / h

      ar(i,m:igh) = ar(i,m:igh) - fr * ortr(m:igh) - fi * orti(m:igh)
      ai(i,m:igh) = ai(i,m:igh) + fr * orti(m:igh) - fi * ortr(m:igh)

    ENDDO

    ortr(m) = scale * ortr(m)
    orti(m) = scale * orti(m)
    ar(m,m-1) = - g * ar(m,m-1)
    ai(m,m-1) = - g * ai(m,m-1)

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE csroot ( xr, xi, yr, yi )

!*****************************************************************************80
!
!! CSROOT computes the complex square root of a complex quantity.
!
!  Discussion:
!
!    The branch of the square function is chosen so that
!      0.0D+00 <= YR
!    and
!      sign ( YI ) == sign ( XI )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, REAL ( kind = 8 ) XR, XI, the real and imaginary parts of the 
!    quantity whose square root is desired.
!
!    Output, REAL ( kind = 8 ) YR, YI, the real and imaginary parts of the 
!    square root.
!
  IMPLICIT NONE

  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) xi
  REAL ( kind = 8 ) xr
  REAL ( kind = 8 ) yi
  REAL ( kind = 8 ) yr

  tr = xr
  ti = xi
  s = sqrt ( 0.5D+00 * ( pythag ( tr, ti ) + abs ( tr ) ) )

  IF( 0.0D+00 <= tr ) then
    yr = s
  ENDIF

  IF( ti < 0.0D+00 ) then
    s = -s
  ENDIF

  IF( tr <= 0.0D+00 ) then
    yi = s
  ENDIF

  IF( tr < 0.0D+00 ) then
    yr = 0.5D+00 * ( ti / yi )
  else IF( 0.0D+00 < tr ) then
    yi = 0.5D+00 * ( ti / yr )
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE elmbak ( n, low, igh, a, ind, m, z )

!*****************************************************************************80
!
!! ELMBAK determines eigenvectors by undoing the ELMHES transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a real general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by ELMHES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, integers determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1 and
!    IGH equal to the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,IGH), the multipliers which were used in the
!    reduction by ELMHES in its lower triangle below the subdiagonal.
!
!    Input, INTEGER ( kind = 4 ) IND(IGH), information on the rows and columns
!    interchanged in the reduction by ELMHES.
!
!    Input, INTEGER ( kind = 4 ) M, the number of columns of Z to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the real and imaginary 
!    parts of the eigenvectors to be back transformed.  On output, the real and
!    imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,igh)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ind(igh)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) z(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  la = igh - 1

  IF( la < low + 1 ) then
    RETURN
  ENDIF

  DO mm = low + 1, la

     mp = low + igh - mm

     DO i = mp + 1, igh

       x = a(i,mp-1)
       IF( x /= 0.0D+00 ) then
         z(i,1:m) = z(i,1:m) + x * z(mp,1:m)
       ENDIF

     ENDDO

     i = ind(mp)

     IF( i /= mp ) then

       DO j = 1, m
         call r8_swap ( z(i,j), z(mp,j) )
       ENDDO

     ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE elmhes ( n, low, igh, a, ind )

!*****************************************************************************80
!
!! ELMHES transforms a real general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a real general matrix, this SUBROUTINE reduces a submatrix
!    situated in rows and columns LOW through IGH to upper Hessenberg
!    form by stabilized elementary similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, James Wilkinson,
!    ELMHES,
!    Numerische Mathematik,
!    Volume 12, pages 349-368, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing 
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, the matrix to be 
!    reduced.  On output, the Hessenberg matrix.  The multipliers
!    which were used in the reduction are stored in the
!    remaining triangle under the Hessenberg matrix.
!
!    Output, INTEGER ( kind = 4 ) IND(N), contains information on the rows
!    and columns interchanged in the reduction.  Only elements LOW through
!    IGH are used.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ind(igh)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y

  la = igh - 1

  DO m = low + 1, la

    x = 0.0D+00
    i = m

    DO j = m, igh
      IF( abs ( x ) < abs ( a(j,m-1) ) ) then
        x = a(j,m-1)
        i = j
      ENDIF
    ENDDO

    ind(m) = i
!
!  Interchange rows and columns of the matrix.
!
    IF( i /= m ) then

      DO j = m - 1, n
        call r8_swap ( a(i,j), a(m,j) )
      ENDDO

      DO j = 1, igh
        call r8_swap ( a(j,i), a(j,m) )
      ENDDO

    ENDIF

    IF( x /= 0.0D+00 ) then

      DO i = m + 1, igh

        y = a(i,m-1)

        IF( y /= 0.0D+00 ) then

          y = y / x
          a(i,m-1) = y

          a(i,m:n) = a(i,m:n) - y * a(m,m:n)

          a(1:igh,m) = a(1:igh,m) + y * a(1:igh,i)

        ENDIF

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE eltran ( n, low, igh, a, ind, z )

!*****************************************************************************80
!
!! ELTRAN accumulates similarity transformations used by ELMHES.
!
!  Discussion:
!
!    This SUBROUTINE accumulates the stabilized elementary
!    similarity transformations used in the reduction of a
!    real general matrix to upper Hessenberg form by ELMHES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Peters, James Wilkinson,
!    ELMTRANS,
!    Numerische Mathematik,
!    Volume 16, pages 181-204, 1970.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing 
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input, REAL ( kind = 8 ) A(N,IGH), the multipliers which were used in the
!    reduction by ELMHES in its lower triangle below the subdiagonal.
!
!    Input, INTEGER ( kind = 4 ) IND(IGH), information on the rows and columns
!    interchanged in the reduction by ELMHES.
!
!    Output, REAL ( kind = 8 ) Z(N,N), the transformation matrix produced in the
!    reduction by ELMHES.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,igh)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ind(igh)
  !INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) kl
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) z(n,n)
!
!  Initialize Z to the identity matrix.
!
  z(1:n,1:n) = 0.0D+00

  DO i = 1, n
    z(i,i) = 1.0D+00
  ENDDO

  kl = igh - low - 1

  IF( kl < 1 ) then
    RETURN
  ENDIF

  DO mm = 1, kl

     mp = igh - mm

     z(mp+1:igh,mp) = a(mp+1:igh,mp-1)

     i = ind(mp)

     IF( i /= mp ) then

       z(mp,mp:igh) = z(i,mp:igh)

       z(i,mp) = 1.0D+00
       z(i,mp+1:igh) = 0.0D+00

     ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE figi ( n, t, d, e, e2, ierr )

!*****************************************************************************80
!
!! FIGI transforms a real nonsymmetric tridiagonal matrix to symmetric form.
!
!  Discussion:
!
!    Given a nonsymmetric tridiagonal matrix such that the products
!    of corresponding pairs of off-diagonal elements are all
!    non-negative, this SUBROUTINE reduces it to a symmetric
!    tridiagonal matrix with the same eigenvalues.  If, further,
!    a zero product only occurs when both factors are zero,
!    the reduced matrix is similar to the original matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) T(N,3) contains the input matrix.  Its subdiagonal
!    is stored in the last N-1 positions of the first column, its diagonal in
!    the N positions of the second column, and its superdiagonal in the
!    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the symmetric 
!    matrix.
!
!    Output, REAL ( kind = 8 ) E(N), contains the subdiagonal elements of
!    the symmetric matrix in E(2:N).  E(1) is not set.
!
!    Output, REAL ( kind = 8 ) E2(N), the squares of the corresponding elements 
!    of E.  E2 may coincide with E if the squares are not needed.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    N+I, if T(I,1) * T(I-1,3) is negative,
!    -(3*N+I), if T(I,1) * T(I-1,3) is zero with one factor non-zero.  In
!      this case, the eigenvectors of the symmetric matrix are not simply
!      related to those of T and should not be sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  REAL ( kind = 8 ) t(n,3)

  ierr = 0

  DO i = 1, n

    IF( 1 < i ) then

      e2(i) = t(i,1) * t(i-1,3)

      IF( e2(i) < 0.0D+00 ) then

        ierr = n + i
        RETURN

      else IF( e2(i) == 0.0D+00 ) then

        IF( t(i,1) /= 0.0D+00 .or. t(i-1,3) /= 0.0D+00 ) then
          ierr = - 3 * n - i
          RETURN
        ENDIF

        e(i) = 0.0D+00

      else

        e(i) = sqrt ( e2(i) )

      ENDIF

    ENDIF

    d(i) = t(i,2)

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE figi2 ( n, t, d, e, z, ierr )

!*****************************************************************************80
!
!! FIGI2 transforms a real nonsymmetric tridiagonal matrix to symmetric form.
!
!  Discussion:
!
!    Given a nonsymmetric tridiagonal matrix such that the products
!    of corresponding pairs of off-diagonal elements are all
!    non-negative, and zero only when both factors are zero, this
!    SUBROUTINE reduces it to a symmetric tridiagonal matrix
!    using and accumulating diagonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) T(N,3) contains the input matrix.  Its subdiagonal
!    is stored in the last N-1 positions of the first column, its diagonal in
!    the N positions of the second column, and its superdiagonal in the
!    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the symmetric 
!    matrix.
!
!    Output, REAL ( kind = 8 ) E(N), contains the subdiagonal elements of the 
!    symmetric matrix in E(2:N).  E(1) is not set.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the transformation matrix
!    produced in the reduction.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    N+I, if T(I,1) * T(I-1,3) is negative,
!    2*N+I, if T(I,1) * T(I-1,3) is zero with one factor non-zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  !INTEGER ( kind = 4 ) j
  REAL ( kind = 8 ) t(n,3)
  REAL ( kind = 8 ) z(n,n)

  ierr = 0

  DO i = 1, n

    z(i,1:n) = 0.0D+00

    IF( i == 1 ) then

      z(i,i) = 1.0D+00

    else

      h = t(i,1) * t(i-1,3)

      IF( h < 0.0D+00 ) then

        ierr = n + i
        RETURN

      else IF( h == 0 ) then

        IF( t(i,1) /= 0.0D+00 .or. t(i-1,3) /= 0.0D+00 ) then
          ierr = 2 * n + i
          RETURN
        ENDIF

        e(i) = 0.0D+00
        z(i,i) = 1.0D+00

      else

        e(i) = sqrt ( h )
        z(i,i) = z(i-1,i-1) * e(i) / t(i-1,3)

      ENDIF

    ENDIF

    d(i) = t(i,2)

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE hqr ( n, low, igh, h, wr, wi, ierr )

!*****************************************************************************80
!
!! HQR computes all eigenvalues of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a real
!    upper Hessenberg matrix by the QR method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, Peters, James Wilkinson,
!    HQR,
!    Numerische Mathematik,
!    Volume 14, pages 219-231, 1970.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, two integers determined by 
!    BALANC.  If BALANC is not used, set LOW=1, IGH=N.
!
!    Input/output, REAL ( kind = 8 ) H(N,N), the N by N upper Hessenberg matrix.
!    Information about the transformations used in the reduction to
!    Hessenberg form by ELMHES or ORTHES, if performed, is stored
!    in the remaining triangle under the Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, no error.
!    J, the limit of 30*N iterations was reached while searching for
!      the J-th eigenvalue.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm2
  REAL ( kind = 8 ) h(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) na
  REAL ( kind = 8 ) norm
  logical notlas
  REAL ( kind = 8 ) p
  REAL ( kind = 8 ) q
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) w
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) zz

  ierr = 0
  norm = 0.0D+00
  k = 1
!
!  Store roots isolated by BALANC and compute matrix norm.
!
  DO i = 1, n

    DO j = k, n
      norm = norm + abs ( h(i,j) )
    ENDDO

    k = i
    IF( i < low .or. igh < i ) then
      wr(i) = h(i,i)
      wi(i) = 0.0D+00
    ENDIF

  ENDDO

  en = igh
  t = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 CONTINUE

  IF( en < low ) then
    RETURN
  ENDIF

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for a single small sub-diagonal element.
!
70 CONTINUE

  DO ll = low, en
    l = en + low - ll
    IF( l == low ) then
      exit
    ENDIF
    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    IF( s == 0.0D+00 ) then
      s = norm
    ENDIF
    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1) )
    IF( tst2 == tst1 ) then
      exit
    ENDIF
  ENDDO
!
!  Form shift.
!
  x = h(en,en)

  IF( l == en ) then
    go to 270
  ENDIF

  y = h(na,na)
  w = h(en,na) * h(na,en)

  IF( l == na ) then
    go to 280
  ENDIF

  IF( itn == 0 ) then
    ierr = en
    RETURN
  ENDIF
!
!  Form an exceptional shift.
!
  IF( its == 10 .or. its == 20 ) then

    t = t + x

    DO i = low, en
      h(i,i) = h(i,i) - x
    ENDDO

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75D+00 * s
    y = x
    w = -0.4375D+00 * s * s

  ENDIF

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  DO mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    IF( m == l ) then
      exit
    ENDIF

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    IF( tst2 == tst1 ) then
      exit
    ENDIF

  ENDDO

  DO i = m + 2, en
    h(i,i-2) = 0.0D+00
    IF( i /= m + 2 ) then
      h(i,i-3) = 0.0D+00
    ENDIF
  ENDDO
!
!  Double QR step involving rows l to EN and columns M to EN.
!
  DO k = m, na

    notlas = k /= na

    IF( k /= m ) then

      p = h(k,k-1)
      q = h(k+1,k-1)

      IF( notlas ) then
        r = h(k+2,k-1)
      else
        r = 0.0D+00
      ENDIF

      x = abs ( p ) + abs ( q ) + abs ( r )

      IF( x == 0.0D+00 ) then
        cycle
      ENDIF

      p = p / x
      q = q / x
      r = r / x

    ENDIF

    s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

    IF( k /= m ) then
      h(k,k-1) = - s * x
    else IF( l /= m ) then
      h(k,k-1) = - h(k,k-1)
    ENDIF

    p = p + s
    x = p / s
    y = q / s
    zz = r / s
    q = q / p
    r = r / p

    IF( .not. notlas ) then
!
!  Row modification.
!
      DO j = k, n
        p = h(k,j) + q * h(k+1,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
      ENDDO

      j = min ( en, k + 3 )
!
!  Column modification.
!
      DO i = 1, j
        p = x * h(i,k) + y * h(i,k+1)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
      ENDDO

    else
!
!  Row modification.
!
      DO j = k, n
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
        h(k+2,j) = h(k+2,j) - p * zz
      ENDDO

      j = min ( en, k + 3 )
!
!  Column modification.
!
      DO i = 1, j
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
        h(i,k+2) = h(i,k+2) - p * r
      ENDDO

    ENDIF

  ENDDO

  go to 70
!
!  One root found.
!
270 CONTINUE

  wr(en) = x + t
  wi(en) = 0.0D+00
  en = na
  go to 60
!
!  Two roots found.
!
280 CONTINUE

  p = ( y - x ) / 2.0D+00
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  x = x + t
!
!  Real root, or complex pair.
!
  IF( 0.0D+00 <= q ) then

    zz = p + sign ( zz, p )
    wr(na) = x + zz
    IF( zz == 0.0D+00 ) then
      wr(en) = wr(na)
    else
      wr(en) = x - w / zz
    ENDIF
    wi(na) = 0.0D+00
    wi(en) = 0.0D+00

  else

    wr(na) = x + p
    wr(en) = x + p
    wi(na) = zz
    wi(en) = -zz

  ENDIF

  en = enm2
  go to 60

END SUBROUTINE

SUBROUTINE hqr2 ( n, low, igh, h, wr, wi, z, ierr )

!*****************************************************************************80
!
!! HQR2 computes eigenvalues and eigenvectors of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues and eigenvectors
!    of a real upper Hessenberg matrix by the qr method.  the
!    eigenvectors of a real general matrix can also be found
!    if ELMHES and ELTRAN or ORTHES and ORTRAN have
!    been used to reduce this general matrix to Hessenberg form
!    and to accumulate the similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input/output, REAL ( kind = 8 ) H(N,N), the N by N upper Hessenberg matrix.
!    On output, the information in H has been destroyed.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts of the
!    eigenvalues.  The eigenvalues are unordered, except that complex
!    conjugate pairs of values appear consecutively, with the eigenvalue
!    having positive imaginary part listed first.  If an error exit
!    occurred, then the eigenvalues should be correct for indices
!    IERR+1 through N.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N).  On input, the transformation 
!    matrix produced by ELTRAN after the reduction by ELMHES, or by ORTRAN after
!    the reduction by ORTHES, if performed.  If the eigenvectors of the 
!    Hessenberg matrix are desired, Z must contain the identity matrix.  On 
!    output, Z contains the real and imaginary parts of the eigenvectors.
!    If the I-th eigenvalue is real, the I-th column of Z contains its
!    eigenvector.  If the I-th eigenvalue is complex with positive imaginary
!    part, the I-th and (I+1)-th columns of Z contain the real and imaginary
!    parts of its eigenvector.  The eigenvectors are unnormalized.  If an
!    error exit is made, none of the eigenvectors has been found.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!      eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm2
  REAL ( kind = 8 ) h(n,n)
  !REAL ( kind = 8 ) hnorm
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) na
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) norm
  logical notlas
  REAL ( kind = 8 ) p
  REAL ( kind = 8 ) q
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) ra
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) sa
  REAL ( kind = 8 ) t
  !REAL ( kind = 8 ) temp
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) vi
  REAL ( kind = 8 ) vr
  REAL ( kind = 8 ) w
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z(n,n)
  REAL ( kind = 8 ) zz

  ierr = 0
  norm = 0.0D+00
  k = 1
!
!  Store roots isolated by BALANC and compute the matrix norm.
!
  DO i = 1, n

    DO j = k, n
      norm = norm + abs ( h(i,j) )
    ENDDO

    k = i
    IF( i < low .or. igh < i ) then
      wr(i) = h(i,i)
      wi(i) = 0.0D+00
    ENDIF

  ENDDO

  en = igh
  t = 0.0D+00
  itn = 30 * n
!
!  Search for next eigenvalues.
!
60 CONTINUE

  IF( en < low ) then
    go to 340
  ENDIF

  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for single small sub-diagonal element.
!
70 CONTINUE

  DO ll = low, en

    l = en + low - ll

    IF( l == low ) then
      exit
    ENDIF

    s = abs ( h(l-1,l-1) ) + abs ( h(l,l) )
    IF( s == 0.0D+00 ) then
      s = norm
    ENDIF

    tst1 = s
    tst2 = tst1 + abs ( h(l,l-1) )

    IF( tst2 == tst1 ) then
      exit
    ENDIF

  ENDDO
!
!  Form shift.
!
  x = h(en,en)
  IF( l == en ) then
    go to 270
  ENDIF

  y = h(na,na)
  w = h(en,na) * h(na,en)

  IF( l == na ) then
    go to 280
  ENDIF

  IF( itn == 0 ) then
    ierr = en
    RETURN
  ENDIF
!
!  Form exceptional shift.
!
  IF( its == 10 .or. its == 20 ) then

    t = t + x

    DO i = low, en
      h(i,i) = h(i,i) - x
    ENDDO

    s = abs ( h(en,na) ) + abs ( h(na,enm2) )
    x = 0.75D+00 * s
    y = x
    w = -0.4375D+00 * s * s

  ENDIF

  its = its + 1
  itn = itn - 1
!
!  Look for two consecutive small sub-diagonal elements.
!
  DO mm = l, enm2

    m = enm2 + l - mm
    zz = h(m,m)
    r = x - zz
    s = y - zz
    p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
    q = h(m+1,m+1) - zz - r - s
    r = h(m+2,m+1)
    s = abs ( p ) + abs ( q ) + abs ( r )
    p = p / s
    q = q / s
    r = r / s

    IF( m == l ) then
      exit
    ENDIF

    tst1 = abs ( p ) * ( abs ( h(m-1,m-1) ) + abs ( zz ) + abs ( h(m+1,m+1) ) )
    tst2 = tst1 + abs ( h(m,m-1) ) * ( abs ( q ) + abs ( r ) )

    IF( tst2 == tst1 ) then
      exit
    ENDIF

  ENDDO

  DO i = m + 2, en
    h(i,i-2) = 0.0D+00
    IF( i /= m + 2 ) then
      h(i,i-3) = 0.0D+00
    ENDIF
  ENDDO
!
!  Double QR step involving rows L to EN and columns M to EN.
!
  DO k = m, na

     notlas = k /= na

     IF( k /= m ) then

       p = h(k,k-1)
       q = h(k+1,k-1)
       r = 0.0D+00
       IF( notlas ) then
         r = h(k+2,k-1)
       ENDIF

       x = abs ( p ) + abs ( q ) + abs ( r )
       IF( x == 0.0D+00 ) then
         cycle
       ENDIF

       p = p / x
       q = q / x
       r = r / x

     ENDIF

     s = sign ( sqrt ( p**2 + q**2 + r**2 ), p )

     IF( k /= m ) then
       h(k,k-1) = - s * x
     else IF( l /= m ) then
       h(k,k-1) = -h(k,k-1)
     ENDIF

     p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p

     IF( notlas ) then
       go to 225
     ENDIF
!
!  Row modification.
!
     DO j = k, n
       p = h(k,j) + q * h(k+1,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
     ENDDO

     j = min ( en, k + 3 )
!
!  Column modification.
!
     DO i = 1, j
       p = x * h(i,k) + y * h(i,k+1)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
     ENDDO
!
!  Accumulate transformations.
!
     DO i = low, igh
       p = x * z(i,k) + y * z(i,k+1)
       z(i,k) = z(i,k) - p
       z(i,k+1) = z(i,k+1) - p * q
     ENDDO

     go to 255

225  CONTINUE
!
!  Row modification.
!
     DO j = k, n
       p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
       h(k,j) = h(k,j) - p * x
       h(k+1,j) = h(k+1,j) - p * y
       h(k+2,j) = h(k+2,j) - p * zz
     ENDDO

     j = min ( en, k + 3 )
!
!  Column modification.
!
     DO i = 1, j
       p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
       h(i,k) = h(i,k) - p
       h(i,k+1) = h(i,k+1) - p * q
       h(i,k+2) = h(i,k+2) - p * r
     ENDDO
!
!  Accumulate transformations.
!
     DO i = low, igh
        p = x * z(i,k) + y * z(i,k+1) + zz * z(i,k+2)
        z(i,k) = z(i,k) - p
        z(i,k+1) = z(i,k+1) - p * q
        z(i,k+2) = z(i,k+2) - p * r
     ENDDO

255 CONTINUE

!260 CONTINUE

  ENDDO

  go to 70
!
!  One root found.
!
270 CONTINUE

  h(en,en) = x + t
  wr(en) = h(en,en)
  wi(en) = 0.0D+00
  en = na
  go to 60
!
!  Two roots found.
!
280 CONTINUE

  p = ( y - x ) / 2.0D+00
  q = p * p + w
  zz = sqrt ( abs ( q ) )
  h(en,en) = x + t
  x = h(en,en)
  h(na,na) = y + t

  IF( q < 0.0D+00 ) then
    go to 320
  ENDIF
!
!  Real pair.
!
  zz = p + sign ( zz, p )
  wr(na) = x + zz
  wr(en) = wr(na)

  IF( zz /= 0.0D+00 ) then
    wr(en) = x - w / zz
  ENDIF

  wi(na) = 0.0D+00
  wi(en) = 0.0D+00
  x = h(en,na)
  s = abs ( x ) + abs ( zz )
  p = x / s
  q = zz / s
  r = sqrt ( p**2 + q**2 )
  p = p / r
  q = q / r
!
!  Row modification.
!
  DO j = na, n
    zz = h(na,j)
    h(na,j) = q * zz + p * h(en,j)
    h(en,j) = q * h(en,j) - p * zz
  ENDDO
!
!  Column modification.
!
  DO i = 1, en
    zz = h(i,na)
    h(i,na) = q * zz + p * h(i,en)
    h(i,en) = q * h(i,en) - p * zz
  ENDDO
!
!  Accumulate transformations.
!
  DO i = low, igh
    zz = z(i,na)
    z(i,na) = q * zz + p * z(i,en)
    z(i,en) = q * z(i,en) - p * zz
  ENDDO

  go to 330
!
!  Complex pair
!
320 CONTINUE

  wr(na) = x + p
  wr(en) = x + p
  wi(na) = zz
  wi(en) = -zz

330 CONTINUE

  en = enm2
  go to 60
!
!  All roots found.
!  Backsubstitute to find vectors of upper triangular form.
!
340 CONTINUE

  IF( norm == 0.0D+00 ) then
    RETURN
  ENDIF

  DO nn = 1, n

    en = n + 1 - nn
    p = wr(en)
    q = wi(en)
    na = en - 1
!    IF( q ) 710, 600, 800
    IF( q < 0 ) THEN
      GO TO 710
    ELSE IF ( q == 0 ) THEN
      GO TO 600
    ELSE 
      GO TO 800
    ENDIF
!
!  Real vector
!
600  CONTINUE

     m = en
     h(en,en) = 1.0D+00

     IF( na == 0 ) then
       go to 800
     ENDIF

     DO ii = 1, na

        i = en - ii
        w = h(i,i) - p
        r = dot_product ( h(i,m:en), h(m:en,en) )

        IF( wi(i) < 0.0D+00 ) then
          zz = w
          s = r
          go to 700
        ENDIF

        m = i

        IF( wi(i) /= 0.0D+00 ) then
          go to 640
        ENDIF

        t = w

        IF( t == 0.0D+00 ) then

          tst1 = norm
          t = tst1

          do
            t = 0.01D+00 * t
            tst2 = norm + t
            IF( tst2 <= tst1 ) then
              exit
            ENDIF
          ENDDO

        ENDIF

        h(i,en) = -r / t
        go to 680
!
!  Solve real equations.
!
640     CONTINUE

        x = h(i,i+1)
        y = h(i+1,i)
        q = ( wr(i) - p ) * ( wr(i) - p ) + wi(i) * wi(i)
        t = ( x * s - zz * r ) / q
        h(i,en) = t

        IF( abs ( zz ) < abs ( x ) ) then
          h(i+1,en) = ( - r - w * t ) / x
        else
          h(i+1,en) = (-s - y * t ) / zz
        ENDIF
!
!  Overflow control.
!
680     CONTINUE

        t = abs ( h(i,en) )

        IF( t /= 0.0D+00 ) then

          tst1 = t
          tst2 = tst1 + 1.0D+00 / tst1

          IF( tst2 <= tst1 ) then
            h(i:en,en) = h(i:en,en) / t
          ENDIF

        ENDIF

700   CONTINUE

    ENDDO
!
!  End real vector
!
     go to 800
!
!  Complex vector
!
710  CONTINUE

     m = na
!
!  Last vector component chosen imaginary, so that the eigenvector
!  matrix is triangular.
!
     IF( abs ( h(en,na) ) > abs ( h(na,en) ) ) then

       h(na,na) = q / h(en,na)
       h(na,en) = -(h(en,en) - p) / h(en,na)

     else

       call cdiv ( 0.0D+00, -h(na,en), h(na,na)-p, q, h(na,na), h(na,en) )

     ENDIF

     h(en,na) = 0.0D+00
     h(en,en) = 1.0D+00
     enm2 = na - 1

     DO ii = 1, enm2

        i = na - ii
        w = h(i,i) - p
        ra = dot_product ( h(i,m:en), h(m:en,na) )
        sa = dot_product ( h(i,m:en), h(m:en,en) )

        IF( wi(i) < 0.0D+00 ) then
          zz = w
          r = ra
          s = sa
        ENDIF

         m = i

        IF( wi(i) == 0.0D+00 ) then
          call cdiv ( -ra, -sa, w, q, h(i,na), h(i,en) )
          go to 790
        ENDIF
!
!  Solve complex equations.
!
        x = h(i,i+1)
        y = h(i+1,i)
        vr = ( wr(i) - p ) * ( wr(i) - p ) + wi(i) * wi(i) - q * q
        vi = ( wr(i) - p ) * 2.0D+00 * q

        IF( vr == 0.0D+00 .and. vi == 0.0D+00 ) then

          tst1 = norm * ( abs ( w ) + abs ( q ) + abs ( x ) &
            + abs ( y ) + abs ( zz ) )
          vr = tst1

          do
            vr = 0.01D+00 * vr
            tst2 = tst1 + vr
            IF( tst2 <= tst1 ) then
              exit
            ENDIF
          ENDDO

        ENDIF

        call cdiv ( x*r-zz*ra+q*sa, x*s-zz*sa-q*ra, vr, vi, h(i,na), h(i,en) )

        IF( abs ( x ) > abs ( zz ) + abs ( q ) ) then
          h(i+1,na) = ( -ra - w * h(i,na) + q * h(i,en) ) / x
          h(i+1,en) = ( -sa - w * h(i,en) - q * h(i,na) ) / x
        else
          call cdiv ( -r-y*h(i,na), -s-y*h(i,en), zz, q, h(i+1,na), h(i+1,en) )
        ENDIF
!
!  Overflow control.
!
790     CONTINUE

        t = max ( abs ( h(i,na) ), abs ( h(i,en) ) )

        IF( t /= 0.0D+00 ) then
          tst1 = t
          tst2 = tst1 + 1.0D+00 / tst1
          IF( tst2 <= tst1 ) then
            h(i:en,na) = h(i:en,na) / t
            h(i:en,en) = h(i:en,en) / t
          ENDIF
        ENDIF

!795     CONTINUE

      ENDDO
!
!  End complex vector.
!
800 CONTINUE

  ENDDO
!
!  End back substitution.
!
!  Vectors of isolated roots.
!
  DO i = 1, n

    IF( i < low .or. igh < i ) then
      z(i,i:n) = h(i,i:n)
    ENDIF

  ENDDO
!
!  Multiply by transformation matrix to give vectors of original full matrix.
!
  DO jj = low, n

     j = n + low - jj
     m = min ( j, igh )

     DO i = low, igh
       z(i,j) = dot_product ( z(i,low:m), h(low:m,j) )
     ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE htrib3 ( n, a, tau, m, zr, zi )

!*****************************************************************************80
!
!! HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a complex hermitian
!    matrix by back transforming those of the corresponding
!    real symmetric tridiagonal matrix determined by HTRID3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, is the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains information about the unitary
!    transformations used in the reduction by HTRID3.
!
!    Input, REAL ( kind = 8 ) TAU(2,N), contains further information about the
!    transformations.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) ZR(N,M), ZI(N,M).  On input, ZR contains 
!    the eigenvectors to be back transformed.  On output, ZR and ZI contain
!    the real and imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) tau(2,n)
  REAL ( kind = 8 ) zi(n,m)
  REAL ( kind = 8 ) zr(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF
!
!  Transform the eigenvectors of the real symmetric tridiagonal matrix
!  to those of the hermitian tridiagonal matrix.
!
  DO k = 1, n
    DO j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    ENDDO
  ENDDO
!
!  Recover and apply the Householder matrices.
!
  DO i = 2, n

    l = i - 1
    h = a(i,i)

    IF( h /= 0.0D+00 ) then

      DO j = 1, m

        s = 0.0D+00
        si = 0.0D+00

        DO k = 1, l
          s = s + a(i,k) * zr(k,j) - a(k,i) * zi(k,j)
          si = si + a(i,k) * zi(k,j) + a(k,i) * zr(k,j)
        ENDDO

        s = ( s / h ) / h
        si = ( si / h ) / h

        zr(1:l,j) = zr(1:l,j) - s * a(i,1:l) - si * a(1:l,i)
        zi(1:l,j) = zi(1:l,j) - si * a(i,1:l) + s * a(1:l,i)

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE htribk ( n, ar, ai, tau, m, zr, zi )

!*****************************************************************************80
!
!! HTRIBK determines eigenvectors by undoing the HTRIDI transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a complex hermitian
!    matrix by back transforming those of the corresponding
!    real symmetric tridiagonal matrix determined by HTRIDI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) AR(N,N), AI(N,N), contain information about
!    the unitary transformations used in the reduction by HTRIDI in their
!    full lower triangles, except for the diagonal of AR.
!
!    Input, REAL ( kind = 8 ) TAU(2,N), contains further information about the
!    transformations.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) ZR(N,M), ZI(N,M).  On input, ZR contains 
!    the eigenvectors to be back transformed.  On output, ZR and ZI contain
!    the real and imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) tau(2,n)
  REAL ( kind = 8 ) zi(n,m)
  REAL ( kind = 8 ) zr(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF
!
!  Transform the eigenvectors of the real symmetric tridiagonal matrix to
!  those of the hermitian tridiagonal matrix.
!
  DO k = 1, n
    DO j = 1, m
      zi(k,j) = -zr(k,j) * tau(2,k)
      zr(k,j) = zr(k,j) * tau(1,k)
    ENDDO
  ENDDO
!
!  Recover and apply the Householder matrices.
!
  DO i = 2, n

    l = i - 1
    h = ai(i,i)

    IF( h /= 0.0D+00 ) then

      DO j = 1, m

        s = 0.0D+00
        si = 0.0D+00
        DO k = 1, l
          s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
          si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
        ENDDO

        s = ( s / h ) / h
        si = ( si / h ) / h

        zr(1:l,j) = zr(1:l,j) - s * ar(i,1:l) - si * ai(i,1:l)
        zi(1:l,j) = zi(1:l,j) - si * ar(i,1:l) + s * ai(i,1:l)

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE htrid3 ( n, a, d, e, e2, tau )

!*****************************************************************************80
!
!! HTRID3 tridiagonalizes a complex hermitian packed matrix.
!
!  Discussion:
!
!    This SUBROUTINE reduces a complex hermitian matrix, stored as
!    a single square array, to a real symmetric tridiagonal matrix
!    using unitary similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, the lower triangle of 
!    the complex hermitian input matrix.  The real parts of the matrix elements 
!    are stored in the full lower triangle of A, and the imaginary parts are 
!    stored in the transposed positions of the strict upper triangle of A.  No 
!    storage is required for the zero imaginary parts of the diagonal elements.
!    On output, A contains information about the unitary transformations
!    used in the reduction.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, REAL ( kind = 8 ) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) E2(N), the squares of the corresponding elements 
!    of E.  E2 may coincide with E if the squares are not needed.
!
!    Output, REAL ( kind = 8 ) TAU(2,N), contains further information about the
!    transformations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) fi
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) gi
  REAL ( kind = 8 ) h
  REAL ( kind = 8 ) hh
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) tau(2,n)

  tau(1,n) = 1.0D+00
  tau(2,n) = 0.0D+00

  DO ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
    scale = 0.0D+00

    IF( l < 1 ) then
      e(i) = 0.0D+00
      e2(i) = 0.0D+00
      go to 290
    ENDIF
!
!  Scale row.
!
     DO k = 1, l
       scale = scale + abs ( a(i,k) ) + abs ( a(k,i) )
     ENDDO

     IF( scale == 0.0D+00 ) then
       tau(1,l) = 1.0D+00
       tau(2,l) = 0.0D+00
       e(i) = 0.0D+00
       e2(i) = 0.0D+00
       go to 290
     ENDIF

      DO k = 1, l
        a(i,k) = a(i,k) / scale
        a(k,i) = a(k,i) / scale
        h = h + a(i,k) * a(i,k) + a(k,i) * a(k,i)
     ENDDO

     e2(i) = scale * scale * h
     g = sqrt ( h )
     e(i) = scale * g
     f = pythag ( a(i,l), a(l,i) )
!
!  Form next diagonal element of matrix T.
!
     IF( f /= 0.0D+00 ) then

       tau(1,l) = ( a(l,i) * tau(2,i) - a(i,l) * tau(1,i) ) / f
       si = ( a(i,l) * tau(2,i) + a(l,i) * tau(1,i) ) / f
       h = h + f * g
       g = 1.0D+00 + g / f
       a(i,l) = g * a(i,l)
       a(l,i) = g * a(l,i)

       IF( l == 1 ) then
         go to 270
       ENDIF

     else

       tau(1,l) = -tau(1,i)
       si = tau(2,i)
       a(i,l) = g

     ENDIF

     f = 0.0D+00

     DO j = 1, l

        g = 0.0D+00
        gi = 0.0D+00
!
!  Form element of A*U.
!
        DO k = 1, j - 1
          g = g + a(j,k) * a(i,k) + a(k,j) * a(k,i)
          gi = gi - a(j,k) * a(k,i) + a(k,j) * a(i,k)
        ENDDO

        g = g + a(j,j) * a(i,j)
        gi = gi - a(j,j) * a(j,i)

        DO k = j + 1, l
          g = g + a(k,j) * a(i,k) - a(j,k) * a(k,i)
          gi = gi - a(k,j) * a(k,i) - a(j,k) * a(i,k)
        ENDDO
!
!  Form element of P.
!
        e(j) = g / h
        tau(2,j) = gi / h
        f = f + e(j) * a(i,j) - tau(2,j) * a(j,i)

     ENDDO

     hh = f / ( h + h )
!
!  Form reduced A.
!
     DO j = 1, l

        f = a(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -a(j,i)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi
        a(j,j) = a(j,j) - 2.0D+00 * ( f * g + fi * gi )

        DO k = 1, j - 1
          a(j,k) = a(j,k) - f * e(k) - g * a(i,k) + fi * tau(2,k) + gi * a(k,i)
          a(k,j) = a(k,j) - f * tau(2,k) - g * a(k,i) - fi * e(k) - gi * a(i,k)
        ENDDO

     ENDDO

270  CONTINUE

     a(i,1:l) = scale * a(i,1:l)
     a(1:l,i) = scale * a(1:l,i)
     tau(2,l) = -si

290  CONTINUE

     d(i) = a(i,i)
     a(i,i) = scale * sqrt ( h )

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE htridi ( n, ar, ai, d, e, e2, tau )

!*****************************************************************************80
!
!! HTRIDI tridiagonalizes a complex hermitian matrix.
!
!  Discussion:
!
!    This SUBROUTINE reduces a complex hermitian matrix to a real symmetric
!    tridiagonal matrix using unitary similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) AR(N,N), AI(N,N).  On input, the real
!    and imaginary parts, respectively, of the complex hermitian input matrix.
!    Only the lower triangle of the matrix need be supplied.
!    On output, information about the unitary transformations used in the
!    reduction in their full lower triangles.  Their strict upper triangles
!    and the diagonal of AR are unaltered.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, REAL ( kind = 8 ) E(N), the subdiagonal elements of the tridiagonal
!    matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) E2(N), the squares of the corresponding elements 
!    of E.  E2 may coincide with E if the squares are not needed.
!
!    Output, REAL ( kind = 8 ) TAU(2,N), contains further information about the
!    transformations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) ai(n,n)
  REAL ( kind = 8 ) ar(n,n)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) fi
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) gi
  REAL ( kind = 8 ) h
  REAL ( kind = 8 ) hh
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) si
  REAL ( kind = 8 ) tau(2,n)

  tau(1,n) = 1.0D+00
  tau(2,n) = 0.0D+00

  DO i = 1, n
    d(i) = ar(i,i)
  ENDDO

  DO ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
    scale = 0.0D+00

    IF( l < 1 ) then
      e(i) = 0.0D+00
      e2(i) = 0.0D+00
      go to 290
    ENDIF
!
!  Scale row.
!
    DO k = 1, l
      scale = scale + abs ( ar(i,k) ) + abs ( ai(i,k) )
    ENDDO

    IF( scale == 0.0D+00 ) then
      tau(1,l) = 1.0D+00
      tau(2,l) = 0.0D+00
      e(i) = 0.0D+00
      e2(i) = 0.0D+00
      go to 290
    ENDIF

    ar(i,1:l) = ar(i,1:l) / scale
    ai(i,1:l) = ai(i,1:l) / scale

    DO k = 1, l
      h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
    ENDDO

    e2(i) = scale * scale * h
    g = sqrt ( h )
    e(i) = scale * g
    f = pythag ( ar(i,l), ai(i,l) )
!
!  Form next diagonal element of matrix T.
!
    IF( f /= 0.0D+00 ) then
      tau(1,l) = ( ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i) ) / f
      si = ( ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i) ) / f
      h = h + f * g
      g = 1.0D+00 + g / f
      ar(i,l) = g * ar(i,l)
      ai(i,l) = g * ai(i,l)
      IF( l == 1 ) then
        go to 270
      ENDIF
    else
      tau(1,l) = -tau(1,i)
      si = tau(2,i)
      ar(i,l) = g
    ENDIF

    f = 0.0D+00

    DO j = 1, l

      g = 0.0D+00
      gi = 0.0D+00
!
!  Form element of A*U.
!
      DO k = 1, j
        g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
        gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
      ENDDO

      DO k = j + 1, l
        g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
        gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
      ENDDO
!
!  Form element of P.
!
      e(j) = g / h
      tau(2,j) = gi / h
      f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)

    ENDDO

    hh = f / ( h + h )
!
!  Form the reduced A.
!
    DO j = 1, l

      f = ar(i,j)
      g = e(j) - hh * f
      e(j) = g
      fi = - ai(i,j)
      gi = tau(2,j) - hh * fi
      tau(2,j) = -gi

      DO k = 1, j
        ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) + fi * tau(2,k) &
          + gi * ai(i,k)
        ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) - fi * e(k) &
          - gi * ar(i,k)
      ENDDO

    ENDDO

270 CONTINUE

    ar(i,1:l) = scale * ar(i,1:l)
    ai(i,1:l) = scale * ai(i,1:l)
    tau(2,l) = -si

290 CONTINUE

    hh = d(i)
    d(i) = ar(i,i)
    ar(i,i) = hh
    ai(i,i) = scale * sqrt ( h )

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE imtql1 ( n, d, e, ierr )

!*****************************************************************************80
!
!! IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a symmetric
!    tridiagonal matrix by the implicit QL method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the matrix.  On output, the eigenvalues in ascending order.  If an error
!    exit is made, the eigenvalues are correct and ordered for indices
!    1,2,...IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, REAL ( kind = 8 ) E(N).  On input, the subdiagonal elements
!    of the matrix in its last N-1 positions.  E(1) is arbitrary.  On output,
!    E has been overwritten.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, normal RETURN,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2

  ierr = 0

  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    e(i-1) = e(i)
  ENDDO
  e(n) = 0.0D+00

  DO l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 CONTINUE

    m = l

    DO m = l, n - 1

      tst1 = abs ( d(m) ) + abs ( d(m+1) )
      tst2 = tst1 + abs ( e(m) )

      IF( tst2 == tst1 ) then
        exit
      ENDIF

    ENDDO

    p = d(l)

    IF( m == l ) then
      go to 215
    ENDIF

    IF( 30 <= j ) then
      ierr = l
      RETURN
    ENDIF

    j = j + 1
!
!  Form shift.
!
    g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
    r = pythag ( g, 1.0D+00 )
    g = d(m) - p + e(l) / ( g + sign ( r, g ) )
    s = 1.0D+00
    c = 1.0D+00
    p = 0.0D+00
    mml = m - l

    DO ii = 1, mml

      i = m - ii
      f = s * e(i)
      b = c * e(i)
      r = pythag ( f, g )
      e(i+1) = r
!
!  Recover from underflow.
!
      IF( r == 0.0D+00 ) then
        d(i+1) = d(i+1) - p
        e(m) = 0.0D+00
        go to 105
      ENDIF

      s = f / r
      c = g / r
      g = d(i+1) - p
      r = ( d(i) - g ) * s + 2.0D+00 * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b

    ENDDO

    d(l) = d(l) - p
    e(l) = g
    e(m) = 0.0D+00
    go to 105
!
!  Order the eigenvalues.
!
215 CONTINUE

    DO ii = 2, l
      i = l + 2 - ii
      IF( d(i-1) <= p ) then
        go to 270
      ENDIF
      d(i) = d(i-1)
    ENDDO

    i = 1

270 CONTINUE

    d(i) = p

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE imtql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues and eigenvectors
!    of a symmetric tridiagonal matrix by the implicit QL method.
!    The eigenvectors of a full symmetric matrix can also
!    be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the input matrix.  On output, the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct but
!    unordered for indices 1,2,...,IERR-1.
!
!    Input/output, REAL ( kind = 8 ) E(N).  On input, the subdiagonal elements
!    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
!    overwritten.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N).  On input, the transformation
!    matrix produced in the reduction by TRED2, if performed.  If the
!    eigenvectors of the tridiagonal matrix are desired, Z must contain the
!    identity matrix.  On output, Z contains orthonormal eigenvectors of the
!    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
!    contains the eigenvectors associated with the stored eigenvalues.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) t(n)
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) z(n,n)

  ierr = 0

  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    e(i-1) = e(i)
  ENDDO
  e(n) = 0.0D+00

  DO l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 CONTINUE

      m = l

      DO m = l, n - 1

        tst1 = abs ( d(m) ) + abs ( d(m+1) )
        tst2 = tst1 + abs ( e(m) )

        IF( tst2 == tst1 ) then
          exit
        ENDIF

      ENDDO

      p = d(l)

      IF( m == l ) then
        cycle
      ENDIF

      IF( 30 <= j ) then
        ierr = l
        RETURN
      ENDIF

      j = j + 1
!
!  Form shift.
!
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r = pythag ( g, 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      DO ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag ( f, g )
        e(i+1) = r
!
!  Recover from underflow.
!
        IF( r == 0.0D+00 ) then
          d(i+1) = d(i+1) - p
          e(m) = 0.0D+00
          go to 105
        ENDIF

        s = f / r
        c = g / r
        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
!
!  Form vector.
!
        DO k = 1, n
          f = z(k,i+1)
          z(k,i+1) = s * z(k,i) + c * f
          z(k,i) = c * z(k,i) - s * f
        ENDDO

      ENDDO

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    go to 105

  ENDDO
!
!  Order eigenvalues and eigenvectors.
!
  DO ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    DO j = ii, n
      IF( d(j) < p ) then
        k = j
        p = d(j)
      ENDIF
    ENDDO

    IF( k /= i ) then

      d(k) = d(i)
      d(i) = p

      t(1:n)   = z(1:n,i)
      z(1:n,i) = z(1:n,k)
      z(1:n,k) = t(1:n)

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE imtqlv ( n, d, e, e2, w, ind, ierr )

!*****************************************************************************80
!
!! IMTQLV computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a symmetric tridiagonal
!    matrix by the implicit QL method and associates with them
!    their corresponding submatrix indices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the input matrix.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E2(N).  On input, the squares of the
!    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of 
!    E2 corresponding to elements of E regarded as negligible have been
!    replaced by zero, causing the matrix to split into a direct sum of
!    submatrices.  E2(1) is also set to zero.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.  If an
!    error exit is made, the eigenvalues are correct and ordered for
!    indices 1,2,...IERR-1, but may not be the smallest eigenvalues.
!
!    Output, INTEGER ( kind = 4 ) IND(N), the submatrix indices associated with 
!    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
!    first submatrix from the top, 2 for those belonging to the second
!    submatrix, and so on.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the J-th eigenvalue has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ind(n)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) s
  INTEGER ( kind = 4 ) tag
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) w(n)

  ierr = 0
  k = 0
  tag = 0
  w(1:n) = d(1:n)
  e2(1) = 0.0D+00
  rv1(1:n-1) = e(2:n)
  rv1(n) = 0.0D+00

  DO l = 1, n

    j = 0
!
!  Look for a small sub-diagonal element.
!
105 CONTINUE

     DO m = l, n

       IF( m == n ) then
         exit
       ENDIF

       tst1 = abs ( w(m) ) + abs ( w(m+1) )
       tst2 = tst1 + abs ( rv1(m) )

       IF( tst2 == tst1 ) then
         exit
       ENDIF
!
!  Guard against underflowed element of E2.
!
       IF( e2(m+1) == 0.0D+00 ) then
         go to 125
       ENDIF

     ENDDO

!120  CONTINUE

     IF( m <= k ) then
       go to 130
     ENDIF

     IF( m /= n ) then
       e2(m+1) = 0.0D+00
     ENDIF

125  CONTINUE

     k = m
     tag = tag + 1

130  CONTINUE

     p = w(l)

     IF( m == l ) then
       go to 215
     ENDIF

     IF( 30 <= j ) then
       ierr = l
       RETURN
     ENDIF

     j = j + 1
!
!  Form shift.
!
     g = ( w(l+1) - p ) / ( 2.0D+00 * rv1(l) )
     r = pythag ( g, 1.0D+00 )
     g = w(m) - p + rv1(l) / ( g + sign ( r, g ) )
     s = 1.0D+00
     c = 1.0D+00
     p = 0.0D+00
     mml = m - l

     DO ii = 1, mml
       i = m - ii
       f = s * rv1(i)
       b = c * rv1(i)
       r = pythag ( f, g )
       rv1(i+1) = r

       IF( r == 0.0D+00 ) then
         go to 210
       ENDIF

       s = f / r
       c = g / r
       g = w(i+1) - p
       r = ( w(i) - g ) * s + 2.0D+00 * c * b
       p = s * r
       w(i+1) = g + p
       g = c * r - b
     ENDDO

     w(l) = w(l) - p
     rv1(l) = g
     rv1(m) = 0.0D+00
     go to 105
!
!  Recover from underflow.
!
210  CONTINUE

     w(i+1) = w(i+1) - p
     rv1(m) = 0.0D+00
     go to 105
!
!  Insert P = W(L) into the sorted list W(1:L-1).
!
215  CONTINUE

     DO ii = 2, l
        i = l + 2 - ii
        IF( w(i-1) <= p ) then
          go to 270
        ENDIF
        w(i) = w(i-1)
        ind(i) = ind(i-1)
     ENDDO

     i = 1

  270   CONTINUE

     w(i) = p
     ind(i) = tag

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE invit ( n, a, wr, wi, select, mm, m, z, ierr )

!*****************************************************************************80
!
!! INVIT computes eigenvectors of a real upper Hessenberg matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvectors of a real upper Hessenberg
!    matrix corresponding to specified eigenvalues, using inverse iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), the Hessenberg matrix.
!
!    Input/output, REAL ( kind = 8 ) WR(N), WI(N).  On input, the real and 
!    imaginary parts, respectively, of the eigenvalues of the matrix.  The 
!    eigenvalues must be stored in a manner identical to that of SUBROUTINE HQR,
!    which recognizes possible splitting of the matrix.  On output,
!    WR may have been altered since close eigenvalues are perturbed
!    slightly in searching for independent eigenvectors.
!
!    Input/output, logical SELECT(N).  On input, specifies the eigenvectors
!    to be found.  The eigenvector corresponding to the J-th eigenvalue is
!    specified by setting SELECT(J) to TRUE.  On output, SELECT may have been
!    altered.  If the elements corresponding to a pair of conjugate complex
!    eigenvalues were each initially set to TRUE, the program resets the
!    second of the two elements to FALSE.
!
!    Input, INTEGER ( kind = 4 ) MM, an upper bound for the number of columns
!    required to store the eigenvectors to be found.  Note that two columns are
!    required to store the eigenvector corresponding to a complex eigenvalue.
!
!    Input, INTEGER ( kind = 4 ) M, the number of columns actually used to store
!    the eigenvectors.
!
!    Output, REAL ( kind = 8 ) Z(N,MM), the real and imaginary parts of the 
!    eigenvectors.  If the next selected eigenvalue is real, the next column
!    of Z contains its eigenvector.  If the eigenvalue is complex, the next
!    two columns of Z contain the real and imaginary parts of its eigenvector.
!    The eigenvectors are normalized so that the component of largest
!    magnitude is 1.  Any vector which fails the acceptance test is set to zero.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    -(2*N+1), if more than MM columns of Z are necessary to store the
!      eigenvectors corresponding to the specified eigenvalues.
!    -K, if the iteration corresponding to the K-th value fails,
!    -(N+K), if both error situations occur.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) eps3
  REAL ( kind = 8 ) growto
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  REAL ( kind = 8 ) ilambd
  INTEGER ( kind = 4 ) ip
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) km1
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  INTEGER ( kind = 4 ) n1
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) normv
  INTEGER ( kind = 4 ) ns
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) rlambd
  REAL ( kind = 8 ) rm1(n,n)
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) rv2(n)
  INTEGER ( kind = 4 ) s
  logical select(n)
  REAL ( kind = 8 ) t
  INTEGER ( kind = 4 ) uk
  REAL ( kind = 8 ) ukroot
  REAL ( kind = 8 ) w
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z(n,mm)

  ierr = 0
  uk = 0
  s = 1
!
!  The value of IP is:
!
!   0, real eigenvalue;
!   1, first of conjugate complex pair;
!  -1, second of conjugate complex pair.
!
  ip = 0
  n1 = n - 1

  DO k = 1, n

     IF( wi(k) /= 0.0D+00 .and. 0 <= ip ) then
       ip = 1
       IF( select(k) .and. select(k+1) ) then
         select(k+1) = .false.
       ENDIF
     ENDIF

     IF( .not. select(k) ) then
       go to 960
     ENDIF

     IF( wi(k) /= 0.0D+00 ) then
       s = s + 1
     ENDIF

     IF( mm < s ) then
       go to 1000
     ENDIF

     IF( k <= uk ) then
       go to 200
     ENDIF
!
!  Check for possible splitting.
!
     DO uk = k, n
       IF( uk == n ) then
         exit
       ENDIF
       IF( a(uk+1,uk) == 0.0D+00 ) then
         exit
       ENDIF
     ENDDO
!
!  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
!
     norm = 0.0D+00
     mp = 1

     DO i = 1, uk

       x = sum ( abs ( a(i,mp:uk) ) )
       norm = max ( norm, x )
       mp = i

     ENDDO
!
!  EPS3 replaces zero pivot in decomposition and close roots are modified
!  by EPS3.
!
     IF( norm == 0.0D+00 ) then
       norm = 1.0D+00
     ENDIF

     eps3 = abs ( norm ) * epsilon ( eps3 )
!
!  GROWTO is the criterion for the growth.
!
     ukroot = uk
     ukroot = sqrt ( ukroot )
     growto = 0.1D+00 / ukroot

200  CONTINUE

     rlambd = wr(k)
     ilambd = wi(k)

     IF( k == 1 ) then
       go to 280
     ENDIF

     km1 = k - 1
     go to 240
!
!  Perturb eigenvalue if it is close to any previous eigenvalue.
!
220 CONTINUE

     rlambd = rlambd + eps3

240  CONTINUE

     DO ii = 1, km1
       i = k - ii
       IF( select(i) .and. &
            abs ( wr(i) - rlambd ) < eps3 .and. &
            abs ( wi(i) - ilambd ) < eps3 ) then
        go to 220
       ENDIF
     ENDDO

     wr(k) = rlambd
!
!  Perturb conjugate eigenvalue to match.
!
     wr(k+ip) = rlambd
!
!  Form upper Hessenberg A - rlambd*I (transposed) and initial real vector.
!
280  CONTINUE

     mp = 1

     DO i = 1, uk

        rm1(mp:uk,i) = a(i,mp:uk)

        rm1(i,i) = rm1(i,i) - rlambd
        mp = i
        rv1(i) = eps3

     ENDDO

     its = 0

     IF( ilambd /= 0.0D+00 ) then
       go to 520
     ENDIF
!
!  Real eigenvalue.
!
!  Triangular decomposition with interchanges, replacing zero pivots by eps3.
!
     DO i = 2, uk

        mp = i - 1

        IF( abs ( rm1(mp,i) ) > abs ( rm1(mp,mp) ) ) then

          DO j = mp, uk
            call r8_swap ( rm1(j,i), rm1(j,mp) )
          ENDDO

        ENDIF

        IF( rm1(mp,mp) == 0.0D+00 ) then
          rm1(mp,mp) = eps3
        ENDIF

        x = rm1(mp,i) / rm1(mp,mp)

        IF( x /= 0.0D+00 ) then
          rm1(i:uk,i) = rm1(i:uk,i) - x * rm1(i:uk,mp)
        ENDIF

      ENDDO

      IF( rm1(uk,uk) == 0.0D+00 ) then
        rm1(uk,uk) = eps3
      ENDIF
!
!  Back substitution for real vector.
!
440   CONTINUE

      DO ii = 1, uk

        i = uk + 1 - ii
        y = rv1(i)

        DO j = i + 1, uk
          y = y - rm1(j,i) * rv1(j)
        ENDDO

        rv1(i) = y / rm1(i,i)

     ENDDO

     go to 740
!
!  Complex eigenvalue.
!
!  Triangular decomposition with interchanges,
!  replacing zero pivots by EPS3.
!  Store imaginary parts in upper triangle starting at (1,3)
!
520  CONTINUE

     ns = n - s
     z(1,s-1) = - ilambd
     z(1,s) = 0.0D+00

     IF( n /= 2 ) then
       rm1(1,3) = - ilambd
       z(1,s-1) = 0.0D+00
       rm1(1,4:n) = 0.0D+00
     ENDIF

     DO i = 2, uk

        mp = i - 1
        w = rm1(mp,i)

        IF( i < n ) then
          t = rm1(mp,i+1)
        else IF( i == n ) then
          t = z(mp,s-1)
        ENDIF

        x = rm1(mp,mp) * rm1(mp,mp) + t * t

        IF( w * w <= x ) then
          go to 580
        ENDIF

        x = rm1(mp,mp) / w
        y = t / w
        rm1(mp,mp) = w

        IF( i < n ) then
          rm1(mp,i+1) = 0.0D+00
        else IF( i == n ) then
          z(mp,s-1) = 0.0D+00
        ENDIF

        DO j = i, uk

          w = rm1(j,i)
          rm1(j,i) = rm1(j,mp) - x * w
          rm1(j,mp) = w

          IF( n1 <= j ) then
            l = j - ns
            z(i,l) = z(mp,l) - y * w
            z(mp,l) = 0.0D+00
          else
            rm1(i,j+2) = rm1(mp,j+2) - y * w
            rm1(mp,j+2) = 0.0D+00
          ENDIF

        ENDDO

        rm1(i,i) = rm1(i,i) - y * ilambd

        IF( n1 <= i ) then
          l = i - ns
          z(mp,l) = -ilambd
          z(i,l) = z(i,l) + x * ilambd
        else
          rm1(mp,i+2) = -ilambd
          rm1(i,i+2) = rm1(i,i+2) + x * ilambd
        ENDIF

        go to 640

580     CONTINUE

        IF( x == 0.0D+00 ) then
          rm1(mp,mp) = eps3
          IF( i < n ) then
            rm1(mp,i+1) = 0.0D+00
          else IF( i == n ) then
            z(mp,s-1) = 0.0D+00
          ENDIF
          t = 0.0D+00
          x = eps3**2
        ENDIF

        w = w / x
        x = rm1(mp,mp) * w
        y = - t * w

        DO j = i, uk

          IF( n1 <= j ) then
            l = j - ns
            t = z(mp,l)
            z(i,l) = - x * t - y * rm1(j,mp)
          else
            t = rm1(mp,j+2)
            rm1(i,j+2) = - x * t - y * rm1(j,mp)
          ENDIF

          rm1(j,i) = rm1(j,i) - x * rm1(j,mp) + y * t

        ENDDO

        IF( n1 <= i ) then
          l = i - ns
          z(i,l) = z(i,l) - ilambd
        else
          rm1(i,i+2) = rm1(i,i+2) - ilambd
        ENDIF

640    CONTINUE

     ENDDO

     IF( n1 <= uk ) then
       l = uk - ns
       t = z(uk,l)
     else
       t = rm1(uk,uk+2)
     ENDIF

     IF( rm1(uk,uk) == 0.0D+00 .and. t == 0.0D+00 ) then
       rm1(uk,uk) = eps3
     ENDIF
!
!  Back substitution for complex vector.
!
660  CONTINUE

     DO ii = 1, uk

        i = uk + 1 - ii
        x = rv1(i)
        y = 0.0D+00

        DO j = i + 1, uk

          IF( n1 <= j ) then
            t = z(i,j-ns)
          else
            t = rm1(i,j+2)
          ENDIF

          x = x - rm1(j,i) * rv1(j) + t * rv2(j)
          y = y - rm1(j,i) * rv2(j) - t * rv1(j)

        ENDDO

        IF( i >= n1 ) then
          t = z(i,i-ns)
        else
          t = rm1(i,i+2)
        ENDIF

       call cdiv ( x, y, rm1(i,i), t, rv1(i), rv2(i) )

     ENDDO
!
!  Acceptance test for real or complex eigenvector and normalization.
!
740  CONTINUE

     its = its + 1
     norm = 0.0D+00
     normv = 0.0D+00

     DO i = 1, uk

       IF( ilambd == 0.0D+00 ) then
         x = abs ( rv1(i) )
       else
         x = pythag ( rv1(i), rv2(i) )
       ENDIF

       IF( normv < x )  then
         normv = x
         j = i
       ENDIF

       norm = norm + x

     ENDDO

     IF( norm < growto ) then
       go to 840
     ENDIF
!
!  Accept vector.
!
     x = rv1(j)
     IF( ilambd == 0.0D+00 ) then
       x = 1.0D+00 / x
     else
       y = rv2(j)
     ENDIF

     DO i = 1, uk
       IF( ilambd == 0.0D+00 ) then
         z(i,s) = rv1(i) * x
       else
         call cdiv ( rv1(i), rv2(i), x, y, z(i,s-1), z(i,s) )
       ENDIF
     ENDDO

     IF( uk == n ) then
       go to 940
     ENDIF

     j = uk + 1
     go to 900
!
!  Choose a new starting vector.
!
840  CONTINUE

     IF( uk <= its ) then
       go to 880
     ENDIF

     x = ukroot
     y = eps3 / ( x + 1.0D+00 )

     rv1(1) = eps3
     rv1(2:uk) = y

     j = uk - its + 1
     rv1(j) = rv1(j) - eps3 * x
     IF( ilambd == 0.0D+00 ) then
       go to 440
     ENDIF
     go to 660
!
!  Set error: unaccepted eigenvector.
!
880  CONTINUE

     j = 1
     ierr = - k
!
!  Set remaining vector components to zero.
!
900  CONTINUE

     DO i = j, n
       z(i,s) = 0.0D+00
       IF( ilambd /= 0.0D+00 ) then
         z(i,s-1) = 0.0D+00
       ENDIF
     ENDDO

940  CONTINUE

     s = s + 1

960  CONTINUE

     IF( ip == (-1) ) then
       ip = 0
     ENDIF

     IF( ip == 1 ) then
       ip = -1
     ENDIF

  ENDDO

  go to 1001
!
!  Set error: underestimate of eigenvector space required.
!
1000 CONTINUE

  IF( ierr /= 0 ) then
    ierr = ierr - n
  ENDIF

  IF( ierr == 0 ) then
    ierr = - ( 2 * n + 1 )
  ENDIF

1001 CONTINUE

  m = s - 1 - abs ( ip )

  RETURN

END SUBROUTINE

SUBROUTINE minfit ( nm, m, n, a, w, ip, b, ierr )

!*****************************************************************************80
!
!! MINFIT: least squares problem for a real overdetermined linear system.
!
!  Discussion:
!
!    This SUBROUTINE is part of an algorithm for solving general linear
!    systems of the form A*X=B.
!
!    It determines the singular value decomposition
!      A = U * S * V'
!    of a real M by N rectangular matrix, forming U' * B
!    rather than U.  Householder bidiagonalization and a variant of the
!    QR algorithm are used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) NM, the leading dimension of the
!    two-dimensional arrays.  NM must be at least as large as the maximum
!    of M and N.
!
!    Input, INTEGER ( kind = 4 ) M, the number of rows of A and B.
!
!    Input, INTEGER ( kind = 4 ) N, the number of columns of A, and the order
!    of V.
!
!    Input/output, REAL ( kind = 8 ) A(NM,N). On input, the rectangular
!    coefficient matrix.  On output, A has been overwritten by the orthogonal
!    matrix V of the decomposition in its first N rows and columns.  If an
!    error exit is made, the columns of V corresponding to indices of correct
!    singular values should be correct.
!
!    Output, REAL ( kind = 8 ) W(N), the singular values of A.  These are the
!    diagonal elements of S.  They are unordered.  If an error exit is made, the
!    singular values should be correct for indices IERR+1, IERR+2,...,N.
!
!    Input, INTEGER ( kind = 4 ) IP, is the number of columns of B.  IP can
!    be zero.
!
!    Input/output, REAL ( kind = 8 ) B(NM,IP).  On input, the constant column
!    matrix.  On output, B has been overwritten by U'*B.  If an error exit is
!    made, the rows of U'*B corresponding to indices of correct singular values
!    should be correct.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    K, if the K-th singular value has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) ip
  INTEGER ( kind = 4 ) n
  INTEGER ( kind = 4 ) nm

  REAL ( kind = 8 ) a(nm,n)
  REAL ( kind = 8 ) b(nm,ip)
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) i1
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) k1
  INTEGER ( kind = 4 ) kk
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) m1
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z

  ierr = 0
!
!  Householder reduction to bidiagonal form.
!
  g = 0.0D+00
  scale = 0.0D+00
  x = 0.0D+00

  DO i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = 0.0D+00
    s = 0.0D+00
    scale = 0.0D+00

    IF( i <= m ) then

      scale = sum ( abs ( a(i:m,i) ) )

      IF( scale /= 0.0D+00 ) then

        a(i:m,i) = a(i:m,i) / scale

        s = s + sum ( a(i:m,i)**2 )

        f = a(i,i)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        a(i,i) = f - g

        DO j = l, n

          s = dot_product ( a(i:m,i), a(i:m,j) )

          f = s / h
          a(i:m,j) = a(i:m,j) + f * a(i:m,i)

        ENDDO

        DO j = 1, ip

          s = dot_product ( a(i:m,i), b(i:m,j) )

          b(i:m,j) = b(i:m,j) + s * a(i:m,i) / h

        ENDDO

        a(i:m,i) = scale * a(i:m,i)

      ENDIF

    ENDIF

    w(i) = scale * g
    g = 0.0D+00
    s = 0.0D+00
    scale = 0.0D+00

    IF( i <= m .and. i /= n ) then

      scale = scale + sum ( abs ( a(i,l:n) ) )

      IF( scale /= 0.0D+00 ) then

        a(i,l:n) = a(i,l:n) / scale

        s = s + sum ( a(i,l:n)**2 )

        f = a(i,l)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        a(i,l) = f - g
        rv1(l:n) = a(i,l:n) / h

        DO j = l, m

          s = dot_product ( a(j,l:n), a(i,l:n) )

          a(j,l:n) = a(j,l:n) + s * rv1(l:n)

        ENDDO

        a(i,l:n) = scale * a(i,l:n)

      ENDIF

    ENDIF

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  ENDDO
!
!  Accumulation of right-hand transformations.
!
  DO ii = 1, n

    i = n + 1 - ii

    IF( i /= n ) then

      IF( g /= 0.0D+00 ) then

        a(l:n,i) = ( a(i,l:n) / a(i,l) ) / g

        DO j = l, n

          s = dot_product ( a(i,l:n), a(l:n,j) )

          a(l:n,j) = a(l:n,j) + s * a(l:n,i)

        ENDDO

      ENDIF

      a(i,l:n) = 0.0D+00
      a(l:n,i) = 0.0D+00

    ENDIF

    a(i,i) = 1.0D+00
    g = rv1(i)
    l = i

  ENDDO

  IF( m < n .and. ip /= 0 ) then
    m1 = m + 1
    b(m+1:n,1:ip) = 0.0D+00
  ENDIF
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  DO kk = 1, n

    k1 = n - kk
    k = k1 + 1
    its = 0
!
!  Test for splitting.
!
520 CONTINUE

    DO ll = 1, k

      l1 = k - ll
      l = l1 + 1
      tst2 = tst1 + abs ( rv1(l) )

      IF( tst2 == tst1 ) then
        go to 565
      ENDIF

      tst2 = tst1 + abs ( w(k-ll) )

      IF( tst2 == tst1 ) then
        exit
      ENDIF

    ENDDO
!
!  Cancellation of RV1(l) if l greater than 1.
!
!540 CONTINUE

    c = 0.0D+00
    s = 1.0D+00

    DO i = l, k

      f = s * rv1(i)
      rv1(i) = c * rv1(i)
      tst2 = tst1 + abs ( f)

      IF( tst2 == tst1 ) then
        exit
      ENDIF

      g = w(i)
      h = pythag ( f, g )
      w(i) = h
      c = g / h
      s = -f / h

      DO j = 1, ip
        y = b(l1,j)
        z = b(i,j)
        b(l1,j) = y * c + z * s
        b(i,j) = -y * s + z * c
      ENDDO

    ENDDO
!
!  Test for convergence.
!
565 CONTINUE

    z = w(k)

    IF( l == k ) then

      IF( z < 0.0D+00 ) then
        w(k) = - z
        a(1:n,k) = - a(1:n,k)
      ENDIF

      exit

    ENDIF
!
!  Shift from bottom 2 by 2 minor.
!
     IF( its >= 30 ) then
       ierr = k
       RETURN
     ENDIF

     its = its + 1
     x = w(l)
     y = w(k1)
     g = rv1(k1)
     h = rv1(k)
     f = 0.5D+00 * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
     g = pythag ( f, 1.0D+00 )
     f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h )
!
!  Next QR transformation.
!
     c = 1.0D+00
     s = 1.0D+00

     DO i1 = l, k1

        i = i1 + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = pythag ( f, h )
        rv1(i1) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c

        DO j = 1, n
          x = a(j,i1)
          z = a(j,i)
          a(j,i1) = x * c + z * s
          a(j,i) = -x * s + z * c
        ENDDO

        z = pythag ( f, h )
        w(i1) = z

        IF( z /= 0.0D+00 ) then
          c = f / z
          s = h / z
        ENDIF

        f = c * g + s * y
        x = -s * g + c * y

        DO j = 1, ip
          y = b(i1,j)
          z = b(i,j)
          b(i1,j) = y * c + z * s
          b(i,j) = -y * s + z * c
        ENDDO

     ENDDO

     rv1(l) = 0.0D+00
     rv1(k) = f
     w(k) = x
     go to 520

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE ortbak ( n, low, igh, a, ort, m, z )

!*****************************************************************************80
!
!! ORTBAK determines eigenvectors by undoing the ORTHES transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a real general
!    matrix by back transforming those of the corresponding
!    upper Hessenberg matrix determined by ORTHES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH equal
!    to the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,IGH), contains information about the 
!    orthogonal transformations used in the reduction by ORTHES in its strict
!    lower triangle.
!
!    Input/output, REAL ( kind = 8 ) ORT(IGH), contains further information
!    about the transformations used in the reduction by ORTHES.  On output, ORT
!    has been altered.
!
!    Input, INTEGER ( kind = 4 ) M, the number of columns of Z to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N).  On input, the real and imaginary
!    parts of the eigenvectors to be back transformed in the first M columns.
!    On output, the real and imaginary parts of the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,igh)
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) ort(igh)
  REAL ( kind = 8 ) z(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  DO mp = igh - 1, low + 1, -1

    IF( a(mp,mp-1) /= 0.0D+00 ) then

      ort(mp+1:igh) = a(mp+1:igh,mp-1)

      DO j = 1, m

        g = dot_product ( ort(mp:igh), z(mp:igh,j) )

        g = ( g / ort(mp) ) / a(mp,mp-1)

        DO i = mp, igh
          z(i,j) = z(i,j) + g * ort(i)
        ENDDO

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE orthes ( n, low, igh, a, ort )

!*****************************************************************************80
!
!! ORTHES transforms a real general matrix to upper Hessenberg form.
!
!  Discussion:
!
!    Given a real general matrix, this SUBROUTINE reduces a submatrix
!    situated in rows and columns LOW through IGH to upper Hessenberg form by
!    orthogonal similarity transformations.
!
!  Modified:
!
!    04 February 2003
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH = N.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, the matrix.  On output,
!    the Hessenberg matrix.  Information about the orthogonal transformations
!    used in the reduction is stored in the remaining triangle under the
!    Hessenberg matrix.
!
!    Output, REAL ( kind = 8 ) ORT(IGH), contains further information about the
!    transformations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) la
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) ort(igh)
  REAL ( kind = 8 ) scale

  la = igh - 1

  DO m = low + 1, la

    h = 0.0D+00
    ort(m) = 0.0D+00
    scale = 0.0D+00
!
!  Scale the column.
!
    DO i = m, igh
      scale = scale + abs ( a(i,m-1) )
    ENDDO

    IF( scale /= 0.0D+00 ) then

      mp = m + igh

      DO ii = m, igh
        i = mp - ii
        ort(i) = a(i,m-1) / scale
        h = h + ort(i) * ort(i)
      ENDDO

      g = - sign ( sqrt ( h ), ort(m) )
      h = h - ort(m) * g
      ort(m) = ort(m) - g
!
!  Form (I-(U*Ut)/h) * A.
!
      DO j = m, n

        f = 0.0D+00

        DO ii = m, igh
          i = mp - ii
          f = f + ort(i) * a(i,j)
        ENDDO

        f = f / h

        DO i = m, igh
          a(i,j) = a(i,j) - f * ort(i)
        ENDDO

      ENDDO
!
!  Form (I-(u*ut)/h) * A * (I-(u*ut)/h).
!
      DO i = 1, igh

        f = 0.0D+00
        DO jj = m, igh
          j = mp - jj
          f = f + ort(j) * a(i,j)
        ENDDO

        a(i,m:igh) = a(i,m:igh) - f * ort(m:igh) / h

      ENDDO

      ort(m) = scale * ort(m)
      a(m,m-1) = scale * g

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE ortran ( n, low, igh, a, ort, z )

!*****************************************************************************80
!
!! ORTRAN accumulates similarity transformations generated by ORTHES.
!
!  Discussion:
!
!    This SUBROUTINE accumulates the orthogonal similarity
!    transformations used in the reduction of a real general
!    matrix to upper Hessenberg form by ORTHES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) LOW, IGH, are determined by the balancing
!    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.
!
!    Input, REAL ( kind = 8 ) A(N,IGH), contains information about the 
!    orthogonal transformations used in the reduction by ORTHES in its strict 
!    lower triangle.
!
!    Input/output, REAL ( kind = 8 ) ORT(IGH), contains further information
!    about the transformations used in the reduction by ORTHES.  On output, ORT
!    has been further altered.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the transformation matrix
!    produced in the reduction by ORTHES.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) igh
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,igh)
  REAL ( kind = 8 ) g
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) kl
  INTEGER ( kind = 4 ) low
  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) mp
  REAL ( kind = 8 ) ort(igh)
  REAL ( kind = 8 ) z(n,n)
!
!  Initialize Z to the identity matrix.
!
  z(1:n,1:n) = 0.0D+00

  DO i = 1, n
    z(i,i) = 1.0D+00
  ENDDO

  kl = igh - low - 1

  IF( kl < 1 ) then
    RETURN
  ENDIF

  DO mm = 1, kl

    mp = igh - mm

    IF( a(mp,mp-1) /= 0.0D+00 ) then

      ort(mp+1:igh) = a(mp+1:igh,mp-1)

      DO j = mp, igh

        g = dot_product ( ort(mp:igh), z(mp:igh,j) )

        g = ( g / ort(mp) ) / a(mp,mp-1)

        z(mp:igh,j) = z(mp:igh,j) + g * ort(mp:igh)

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

FUNCTION pythag ( a, b ) RESULT(pyth)

!*****************************************************************************80
!
!! PYTHAG computes SQRT ( A * A + B * B ) carefully.
!
!  Discussion:
!
!    The formula
!
!      PYTHAG = sqrt ( A * A + B * B )
!
!    is reasonably accurate, but can fail if, for example, A^2 is larger
!    than the machine overflow.  The formula can lose most of its accuracy
!    if the sum of the squares is very large or very small.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Modified:
!
!    04 February 2003
!
!  Parameters:
!
!    Input, REAL ( kind = 8 ) A, B, the two legs of a right triangle.
!
!    Output, REAL ( kind = 8 ) PYTH, the length of the hypotenuse.
!
  IMPLICIT NONE

  REAL ( kind = 8 ) a
  REAL ( kind = 8 ) b
  REAL ( kind = 8 ) p
  REAL ( kind = 8 ) pyth
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) u

  p = max ( abs ( a ), abs ( b ) )

  IF( p /= 0.0D+00 ) then

    r = ( min ( abs ( a ), abs ( b ) ) / p )**2

    do

      t = 4.0D+00 + r

      IF( t == 4.0D+00 ) then
        exit
      ENDIF

      s = r / t
      u = 1.0D+00 + 2.0D+00 * s
      p = u * p
      r = ( s / u )**2 * r

    ENDDO

  ENDIF

  pyth = p

  RETURN

END FUNCTION

SUBROUTINE qzhes ( n, a, b, matz, z )

!*****************************************************************************80
!
!! QZHES carries out transformations for a generalized eigenvalue problem.
!
!  Discussion:
!
!    This SUBROUTINE is the first step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This SUBROUTINE accepts a pair of real general matrices and
!    reduces one of them to upper Hessenberg form and the other
!    to upper triangular form using orthogonal transformations.
!    it is usually followed by QZIT, QZVAL and, possibly, QZVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, the first real general
!    matrix.  On output, A has been reduced to upper Hessenberg form.  The
!    elements below the first subdiagonal have been set to zero.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, a real general matrix.
!    On output, B has been reduced to upper triangular form.  The elements
!    below the main diagonal have been set to zero.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the product of the right hand
!    transformations if MATZ is TRUE.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) lb
  logical matz
  INTEGER ( kind = 4 ) nk1
  INTEGER ( kind = 4 ) nm1
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) rho
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) u1
  REAL ( kind = 8 ) u2
  REAL ( kind = 8 ) v1
  REAL ( kind = 8 ) v2
  REAL ( kind = 8 ) z(n,n)
!
!  Set Z to the identity matrix.
!
  IF( matz ) then

    z(1:n,1:n) = 0.0D+00

    DO i = 1, n
      z(i,i) = 1.0D+00
    ENDDO

  ENDIF
!
!  Reduce B to upper triangular form.
!
  IF( n <= 1 ) then
    RETURN
  ENDIF

  nm1 = n - 1

  DO l = 1, n - 1

    l1 = l + 1

    s = sum ( abs ( b(l+1:n,l) ) )

    IF( s /= 0.0D+00 ) then

      s = s + abs ( b(l,l) )
      b(l:n,l) = b(l:n,l) / s

      r = sqrt ( sum ( b(l:n,l)**2 ) )
      r = sign ( r, b(l,l) )
      b(l,l) = b(l,l) + r
      rho = r * b(l,l)

      DO j = l + 1, n

        t = dot_product ( b(l:n,l), b(l:n,j) )

        b(l:n,j) = b(l:n,j) - t * b(l:n,l) / rho

      ENDDO

      DO j = 1, n

        t = dot_product ( b(l:n,l), a(l:n,j) )

        a(l:n,j) = a(l:n,j) - t * b(l:n,l) / rho

      ENDDO

      b(l,l) = - s * r
      b(l+1:n,l) = 0.0D+00

    ENDIF

  ENDDO
!
!  Reduce A to upper Hessenberg form, while keeping B triangular.
!
  IF( n == 2 ) then
    RETURN
  ENDIF

  DO k = 1, n - 2

     nk1 = nm1 - k

     DO lb = 1, nk1

        l = n - lb
        l1 = l + 1
!
!  Zero A(l+1,k).
!
        s = abs ( a(l,k) ) + abs ( a(l1,k) )

        IF( s /= 0.0D+00 ) then

        u1 = a(l,k) / s
        u2 = a(l1,k) / s
        r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
        v1 = - ( u1 + r) / r
        v2 = - u2 / r
        u2 = v2 / v1

        DO j = k, n
          t = a(l,j) + u2 * a(l1,j)
          a(l,j) = a(l,j) + t * v1
          a(l1,j) = a(l1,j) + t * v2
        ENDDO

        a(l1,k) = 0.0D+00

        DO j = l, n
          t = b(l,j) + u2 * b(l1,j)
          b(l,j) = b(l,j) + t * v1
          b(l1,j) = b(l1,j) + t * v2
        ENDDO
!
!  Zero B(l+1,l).
!
        s = abs ( b(l1,l1) ) + abs ( b(l1,l) )

        IF( s /= 0.0 ) then

          u1 = b(l1,l1) / s
          u2 = b(l1,l) / s
          r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
          v1 =  -( u1 + r ) / r
          v2 = -u2 / r
          u2 = v2 / v1

          DO i = 1, l1
            t = b(i,l1) + u2 * b(i,l)
            b(i,l1) = b(i,l1) + t * v1
            b(i,l) = b(i,l) + t * v2
          ENDDO

          b(l1,l) = 0.0D+00

          DO i = 1, n
            t = a(i,l1) + u2 * a(i,l)
            a(i,l1) = a(i,l1) + t * v1
            a(i,l) = a(i,l) + t * v2
          ENDDO

          IF( matz ) then

            DO i = 1, n
              t = z(i,l1) + u2 * z(i,l)
              z(i,l1) = z(i,l1) + t * v1
              z(i,l) = z(i,l) + t * v2
            ENDDO

          ENDIF

        ENDIF

      ENDIF

    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE qzit ( n, a, b, eps1, matz, z, ierr )

!*****************************************************************************80
!
!! QZIT carries out iterations to solve a generalized eigenvalue problem.
!
!  Discussion:
!
!    This SUBROUTINE is the second step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This SUBROUTINE accepts a pair of real matrices, one of them
!    in upper Hessenberg form and the other in upper triangular form.
!    It reduces the Hessenberg matrix to quasi-triangular form using
!    orthogonal transformations while maintaining the triangular form
!    of the other matrix.  It is usually preceded by QZHES and
!    followed by QZVAL and, possibly, QZVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, a real upper Hessenberg 
!    matrix.  On output, A has been reduced to quasi-triangular form.  The 
!    elements below the first subdiagonal are still zero and no two consecutive
!    subdiagonal elements are nonzero.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, a real upper triangular 
!    matrix.  On output, B is still in upper triangular form, although its 
!    elements have been altered.  The location B(N,1) is used to store EPS1 
!    times the norm of B for later use by QZVAL and QZVEC.
!
!    Input, REAL ( kind = 8 ) EPS1, a tolerance used to determine negligible
!    elements.  EPS1 = 0.0 (or negative) may be input, in which case an element
!    will be neglected only if it is less than roundoff error times the
!    norm of its matrix.  If the input EPS1 is positive, then an element
!    will be considered negligible if it is less than EPS1 times the norm
!    of its matrix.  A positive value of EPS1 may result in faster execution,
!    but less accurate results.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N).  If MATZ is FALSE, Z is not 
!    referenced.  Otherwise, on input, the transformation matrix produced in the
!    reduction by QZHES, if performed, or else the identity matrix.  On output, 
!    Z contains the product of the right hand transformations for both steps.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    J, if the limit of 30*N iterations is exhausted while the J-th
!      eigenvalue is being sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) a1
  REAL ( kind = 8 ) a11
  REAL ( kind = 8 ) a12
  REAL ( kind = 8 ) a2
  REAL ( kind = 8 ) a21
  REAL ( kind = 8 ) a22
  REAL ( kind = 8 ) a3
  REAL ( kind = 8 ) a33
  REAL ( kind = 8 ) a34
  REAL ( kind = 8 ) a43
  REAL ( kind = 8 ) a44
  REAL ( kind = 8 ) ani
  REAL ( kind = 8 ) anorm
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) b11
  REAL ( kind = 8 ) b12
  REAL ( kind = 8 ) b22
  REAL ( kind = 8 ) b33
  REAL ( kind = 8 ) b34
  REAL ( kind = 8 ) b44
  REAL ( kind = 8 ) bni
  REAL ( kind = 8 ) bnorm
  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm2
  INTEGER ( kind = 4 ) enorn
  REAL ( kind = 8 ) ep
  REAL ( kind = 8 ) eps1
  REAL ( kind = 8 ) epsa
  REAL ( kind = 8 ) epsb
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ish
  INTEGER ( kind = 4 ) itn
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) k1
  INTEGER ( kind = 4 ) k2
  INTEGER ( kind = 4 ) km1
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) ld
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) lm1
  INTEGER ( kind = 4 ) lor1
  logical matz
  INTEGER ( kind = 4 ) na
  logical notlas
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) sh
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) u1
  REAL ( kind = 8 ) u2
  REAL ( kind = 8 ) u3
  REAL ( kind = 8 ) v1
  REAL ( kind = 8 ) v2
  REAL ( kind = 8 ) v3
  REAL ( kind = 8 ) z(n,n)

  ierr = 0
!
!  Compute EPSA and EPSB.
!
  anorm = 0.0D+00
  bnorm = 0.0D+00

  DO i = 1, n

    IF( i == 1 ) then
      ani = 0.0D+00
    else
      ani = abs ( a(i,i-1) )
    ENDIF

    bni = 0.0D+00

    DO j = i, n
      ani = ani + abs ( a(i,j) )
      bni = bni + abs ( b(i,j) )
    ENDDO

    anorm = max ( anorm, ani )
    bnorm = max ( bnorm, bni )

  ENDDO

  IF( anorm == 0.0D+00 ) then
    anorm = 1.0D+00
  ENDIF

  IF( bnorm == 0.0D+00 ) then
    bnorm = 1.0D+00
  ENDIF

  ep = eps1

  IF( ep > 0.0D+00 ) then
    go to 50
  ENDIF
!
!  Use roundoff level if EPS1 is 0.
!
  ep = epsilon ( ep )

50 CONTINUE

  epsa = ep * anorm
  epsb = ep * bnorm
!
!  Reduce A to quasi-triangular form, while keeping B triangular.
!
  lor1 = 1
  enorn = n
  en = n
  itn = 30 * n
!
!  Begin QZ step.
!
60 CONTINUE

  IF( en <= 2 ) then
    IF( 1 < n ) then
      b(n,1) = epsb
    ENDIF
    RETURN
  ENDIF

  IF( .not. matz ) then
    enorn = en
  ENDIF

  its = 0
  na = en - 1
  enm2 = na - 1

70 CONTINUE

  ish = 2
!
!  Check for convergence or reducibility.
!
  DO ll = 1, en
    lm1 = en - ll
    l = lm1 + 1
    IF( l == 1 ) then
      go to 95
    ENDIF
    IF( abs ( a(l,lm1) ) <= epsa ) then
      exit
    ENDIF
  ENDDO

90 CONTINUE

  a(l,lm1) = 0.0D+00

  IF( l < na ) then
    go to 95
  ENDIF
!
!  1-by-1 or 2-by-2 block isolated.
!
  en = lm1
  go to 60
!
!  Check for small top of B.
!
95 CONTINUE

  ld = l

100 CONTINUE

  l1 = l + 1
  b11 = b(l,l)

  IF( abs ( b11 ) > epsb ) then
    go to 120
  ENDIF

  b(l,l) = 0.0D+00
  s = abs ( a(l,l) ) + abs ( a(l1,l) )
  u1 = a(l,l) / s
  u2 = a(l1,l) / s
  r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
  v1 = - ( u1 + r ) / r
  v2 = -u2 / r
  u2 = v2 / v1

  DO j = l, enorn
    t = a(l,j) + u2 * a(l1,j)
    a(l,j) = a(l,j) + t * v1
    a(l1,j) = a(l1,j) + t * v2
    t = b(l,j) + u2 * b(l1,j)
    b(l,j) = b(l,j) + t * v1
    b(l1,j) = b(l1,j) + t * v2
  ENDDO

  IF( l /= 1 ) then
    a(l,lm1) = -a(l,lm1)
  ENDIF
  lm1 = l
  l = l1
  go to 90

120 CONTINUE

  a11 = a(l,l) / b11
  a21 = a(l1,l) / b11
  IF( ish == 1 ) then
    go to 140
  ENDIF
!
!  Iteration strategy.
!
  IF( itn == 0 ) then
    go to 1000
  ENDIF

  IF( its == 10 ) then
    go to 155
  ENDIF
!
!  Determine type of shift.
!
  b22 = b(l1,l1)
  IF( abs ( b22 ) < epsb ) then
    b22 = epsb
  ENDIF

  b33 = b(na,na)
  IF( abs ( b33 ) < epsb ) then
    b33 = epsb
  ENDIF

  b44 = b(en,en)
  IF( abs ( b44 ) < epsb ) then
    b44 = epsb
  ENDIF

  a33 = a(na,na) / b33
  a34 = a(na,en) / b44
  a43 = a(en,na) / b33
  a44 = a(en,en) / b44
  b34 = b(na,en) / b44
  t = 0.5D+00 * (a43 * b34 - a33 - a44)
  r = t * t + a34 * a43 - a33 * a44

  IF( r < 0.0D+00 ) then
    go to 150
  ENDIF
!
!  Determine single shift zeroth column of A.
!
  ish = 1
  r = sqrt ( r )
  sh = -t + r
  s = -t - r
  IF( abs ( s - a44 ) < abs ( sh - a44 ) ) then
    sh = s
  ENDIF
!
!  Look for two consecutive small sub-diagonal elements of A.
!
  DO ll = ld, enm2
    l = enm2 + ld - ll
    IF( l == ld ) then
      exit
    ENDIF
    lm1 = l - 1
    l1 = l + 1
    t = a(l,l)

    IF( abs ( b(l,l) ) > epsb ) then
      t = t - sh * b(l,l)
    ENDIF

    IF( abs ( a(l,lm1) ) <= abs ( t / a(l1,l) ) * epsa ) then
      go to 100
    ENDIF

  ENDDO

140 CONTINUE

  a1 = a11 - sh
  a2 = a21

  IF( l /= ld ) then
    a(l,lm1) = -a(l,lm1)
  ENDIF

  go to 160
!
!  Determine double shift zeroth column of A.
!
150 CONTINUE

  a12 = a(l,l1) / b22
  a22 = a(l1,l1) / b22
  b12 = b(l,l1) / b22
  a1 = ( ( a33 - a11 ) * ( a44 - a11 ) - a34 * a43 + a43 * b34 * a11 ) &
    / a21 + a12 - a11 * b12
  a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34
  a3 = a(l1+1,l1) / b22
  go to 160
!
!  Ad hoc shift.
!
155 CONTINUE

  a1 = 0.0D+00
  a2 = 1.0D+00
  a3 = 1.1605D+00

  160 CONTINUE
  its = its + 1
  itn = itn - 1
  IF( .not. matz ) then
    lor1 = ld
  ENDIF
!
!  Main loop.
!
  DO k = l, na

     notlas = ( k /= na .and. ish == 2 )
     k1 = k + 1
     k2 = k + 2
     km1 = max ( k - 1, l )
     ll = min ( en, k1 + ish )

     IF( notlas ) then
       go to 190
     ENDIF
!
!  Zero A(k+1,k-1).
!
     IF( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
     ENDIF

     s = abs ( a1 ) + abs ( a2 )

     IF( s == 0.0D+00 ) then
       go to 70
     ENDIF

     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u1**2 ), u1 )
     v1 = -( u1 + r ) / r
     v2 = -u2 / r
     u2 = v2 / v1

     DO j = km1, enorn
       t = a(k,j) + u2 * a(k1,j)
       a(k,j) = a(k,j) + t * v1
       a(k1,j) = a(k1,j) + t * v2
       t = b(k,j) + u2 * b(k1,j)
       b(k,j) = b(k,j) + t * v1
       b(k1,j) = b(k1,j) + t * v2
     ENDDO

     IF( k /= l ) then
       a(k1,km1) = 0.0D+00
     ENDIF

     go to 240
!
!  Zero A(k+1,k-1) and A(k+2,k-1).
!
190  CONTINUE

     IF( k /= l ) then
       a1 = a(k,km1)
       a2 = a(k1,km1)
       a3 = a(k2,km1)
     ENDIF

     s = abs ( a1 ) + abs ( a2 ) + abs ( a3 )

     IF( s == 0.0D+00 ) then
       go to 260
     ENDIF

     u1 = a1 / s
     u2 = a2 / s
     u3 = a3 / s
     r = sign ( sqrt ( u1**2 + u2**2 + u3**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     DO j = km1, enorn
       t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
       a(k,j) = a(k,j) + t * v1
       a(k1,j) = a(k1,j) + t * v2
       a(k2,j) = a(k2,j) + t * v3
       t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
       b(k,j) = b(k,j) + t * v1
       b(k1,j) = b(k1,j) + t * v2
       b(k2,j) = b(k2,j) + t * v3
     ENDDO

     IF( k /= l ) then
       a(k1,km1) = 0.0D+00
       a(k2,km1) = 0.0D+00
     ENDIF
!
!  Zero B(k+2,k+1) and B(k+2,k).
!
     s = abs ( b(k2,k2) ) + abs ( b(k2,k1) ) + abs ( b(k2,k) )
     IF( s == 0.0D+00 ) then
       go to 240
     ENDIF

     u1 = b(k2,k2) / s
     u2 = b(k2,k1) / s
     u3 = b(k2,k) / s
     r = sign ( sqrt ( u1**2 + u2**2 + u3**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     v3 = -u3 / r
     u2 = v2 / v1
     u3 = v3 / v1

     DO i = lor1, ll
       t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
       a(i,k2) = a(i,k2) + t * v1
       a(i,k1) = a(i,k1) + t * v2
       a(i,k) = a(i,k) + t * v3
       t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
       b(i,k2) = b(i,k2) + t * v1
       b(i,k1) = b(i,k1) + t * v2
       b(i,k) = b(i,k) + t * v3
     ENDDO

     b(k2,k) = 0.0D+00
     b(k2,k1) = 0.0D+00

     IF( matz ) then

       DO i = 1, n
         t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
         z(i,k2) = z(i,k2) + t * v1
         z(i,k1) = z(i,k1) + t * v2
         z(i,k) = z(i,k) + t * v3
       ENDDO

     ENDIF
!
!  Zero B(k+1,k).
!
240  CONTINUE

     s = abs ( b(k1,k1) ) + abs ( b(k1,k) )

     IF( s /= 0.0D+00 ) then

       u1 = b(k1,k1) / s
       u2 = b(k1,k) / s
       r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
       v1 = -( u1 + r ) / r
       v2 = -u2 / r
       u2 = v2 / v1
  
       DO i = lor1, ll
         t = a(i,k1) + u2 * a(i,k)
         a(i,k1) = a(i,k1) + t * v1
         a(i,k) = a(i,k) + t * v2
         t = b(i,k1) + u2 * b(i,k)
         b(i,k1) = b(i,k1) + t * v1
         b(i,k) = b(i,k) + t * v2
       ENDDO

       b(k1,k) = 0.0D+00

       IF( matz ) then

         DO i = 1, n
           t = z(i,k1) + u2 * z(i,k)
           z(i,k1) = z(i,k1) + t * v1
           z(i,k) = z(i,k) + t * v2
         ENDDO

       ENDIF

     ENDIF

260  CONTINUE

  ENDDO

  go to 70
!
!  Set error: not all eigenvalues have converged after 30*N iterations.
!  Save EPSB for use by QZVAL and QZVEC.
!
1000 CONTINUE

  ierr = en

  IF( 1 < n ) then
    b(n,1) = epsb
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE qzval ( n, a, b, alfr, alfi, beta, matz, z )

!*****************************************************************************80
!
!! QZVAL computes eigenvalues for a generalized eigenvalue problem.
!
!  Discussion:
!
!    This SUBROUTINE is the third step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This SUBROUTINE accepts a pair of real matrices, one of them
!    in quasi-triangular form and the other in upper triangular form.
!    It reduces the quasi-triangular matrix further, so that any
!    remaining 2-by-2 blocks correspond to pairs of complex
!    eigenvalues, and RETURNs quantities whose ratios give the
!    generalized eigenvalues.  It is usually preceded by QZHES
!    and QZIT and may be followed by QZVEC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, a real upper
!    quasi-triangular matrix.  On output, A has been reduced further to a
!    quasi-triangular matrix in which all nonzero subdiagonal elements
!    correspond to pairs of complex eigenvalues.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, a real upper triangular 
!    matrix.  In addition, location B(n,1) contains the tolerance quantity EPSB
!    computed and saved in QZIT.  On output, B is still in upper triangular
!    form, although its elements have been altered.  B(N,1) is unaltered.
!
!    Output, REAL ( kind = 8 ) ALFR(N), ALFI(N), the real and imaginary parts of
!    the diagonal elements of the triangular matrix that would be obtained
!    if A were reduced completely to triangular form by unitary
!    transformations.  Non-zero values of ALFI occur in pairs, the first
!    member positive and the second negative.
!
!    Output, REAL ( kind = 8 ) BETA(N), the diagonal elements of the 
!    corresponding B, normalized to be real and non-negative.  The generalized 
!    eigenvalues are then the ratios (ALFR + I * ALFI) / BETA.
!
!    Input, logical MATZ, should be TRUE if the right hand transformations
!    are to be accumulated for later use in computing eigenvectors, and
!    to FALSE otherwise.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N), is only used if MATZ is TRUE.
!    On input, the transformation matrix produced in the reductions by QZHES
!    and QZIT, if performed, or else the identity matrix.  On output,
!    the product of the right hand transformations for all three steps.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) a1
  REAL ( kind = 8 ) a11
  REAL ( kind = 8 ) a11i
  REAL ( kind = 8 ) a11r
  REAL ( kind = 8 ) a12
  REAL ( kind = 8 ) a12i
  REAL ( kind = 8 ) a12r
  REAL ( kind = 8 ) a1i
  REAL ( kind = 8 ) a2
  REAL ( kind = 8 ) a21
  REAL ( kind = 8 ) a22
  REAL ( kind = 8 ) a22i
  REAL ( kind = 8 ) a22r
  REAL ( kind = 8 ) a2i
  REAL ( kind = 8 ) an
  REAL ( kind = 8 ) alfi(n)
  REAL ( kind = 8 ) alfr(n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) b11
  REAL ( kind = 8 ) b12
  REAL ( kind = 8 ) b22
  REAL ( kind = 8 ) beta(n)
  REAL ( kind = 8 ) bn
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) cq
  REAL ( kind = 8 ) cz
  REAL ( kind = 8 ) d
  REAL ( kind = 8 ) di
  REAL ( kind = 8 ) dr
  REAL ( kind = 8 ) e
  REAL ( kind = 8 ) ei
  INTEGER ( kind = 4 ) en
  REAL ( kind = 8 ) epsb
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) isw
  INTEGER ( kind = 4 ) j
  logical matz
  INTEGER ( kind = 4 ) na
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) sqi
  REAL ( kind = 8 ) sqr
  REAL ( kind = 8 ) ssi
  REAL ( kind = 8 ) ssr
  REAL ( kind = 8 ) szi
  REAL ( kind = 8 ) szr
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) u1
  REAL ( kind = 8 ) u2
  REAL ( kind = 8 ) v1
  REAL ( kind = 8 ) v2
  REAL ( kind = 8 ) z(n,n)

  epsb = b(n,1)
  isw = 1
!
!  Find eigenvalues of quasi-triangular matrices.
!
  DO nn = 1, n

     en = n + 1 - nn
     na = en - 1

     IF( isw == 2 ) then
       go to 505
     ENDIF

     IF( en == 1 ) then
       go to 410
     ENDIF

     IF( a(en,na) /= 0.0D+00 ) then
       go to 420
     ENDIF
!
!  1-by-1 block, one real root.
!
410  CONTINUE

     alfr(en) = a(en,en)
     IF( b(en,en) < 0.0D+00 ) then
       alfr(en) = -alfr(en)
     ENDIF

     beta(en) = abs ( b(en,en) )
     alfi(en) = 0.0D+00
     go to 510
!
!  2-by-2 block.
!
420  CONTINUE

     IF( abs ( b(na,na) ) <= epsb ) then
       a1 = a(na,na)
       a2 = a(en,na)
       go to 460
     ENDIF

     IF( abs ( b(en,en) ) <= epsb ) then
       a1 = a(en,en)
       a2 = a(en,na)
       bn = 0.0D+00
       go to 435
     ENDIF

     an = abs ( a(na,na) ) + abs ( a(na,en) ) + abs ( a(en,na) ) &
       + abs ( a(en,en) )
     bn = abs ( b(na,na) ) + abs ( b(na,en) ) + abs ( b(en,en) )
     a11 = a(na,na) / an
     a12 = a(na,en) / an
     a21 = a(en,na) / an
     a22 = a(en,en) / an
     b11 = b(na,na) / bn
     b12 = b(na,en) / bn
     b22 = b(en,en) / bn
     e = a11 / b11
     ei = a22 / b22
     s = a21 / ( b11 * b22 )
     t = ( a22 - e * b22 ) / b22

     IF( abs ( e ) > abs ( ei ) ) then
       e = ei
       t = ( a11 - e * b11 ) / b11
     ENDIF

     c = 0.5D+00 * ( t - s * b12 )
     d = c**2 + s * ( a12 - e * b12 )

     IF( d < 0.0D+00 ) then
       go to 480
     ENDIF
!
!  Two real roots.
!  Zero both A(EN,NA) and B(EN,NA).
!
     e = e + ( c + sign ( sqrt ( d ), c ) )
     a11 = a11 - e * b11
     a12 = a12 - e * b12
     a22 = a22 - e * b22

     IF( abs ( a11 ) + abs ( a12 ) >= abs ( a21 ) + abs ( a22 ) ) then
       a1 = a12
       a2 = a11
     else
       a1 = a22
       a2 = a21
     ENDIF
!
!  Choose and apply real Z.
!
435  CONTINUE

     s = abs ( a1 ) + abs ( a2 )
     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
     v1 = - ( u1 + r ) / r
     v2 = - u2 / r
     u2 = v2 / v1

     DO i = 1, en
       t = a(i,en) + u2 * a(i,na)
       a(i,en) = a(i,en) + t * v1
       a(i,na) = a(i,na) + t * v2
       t = b(i,en) + u2 * b(i,na)
       b(i,en) = b(i,en) + t * v1
       b(i,na) = b(i,na) + t * v2
     ENDDO

     IF( matz ) then

       DO i = 1, n
         t = z(i,en) + u2 * z(i,na)
         z(i,en) = z(i,en) + t * v1
         z(i,na) = z(i,na) + t * v2
       ENDDO

     ENDIF

!450  CONTINUE

     IF( bn == 0.0D+00 ) then
       go to 475
     ENDIF

     IF( abs ( e ) * bn <= an ) then
       a1 = b(na,na)
       a2 = b(en,na)
     else
       a1 = a(na,na)
       a2 = a(en,na)
     ENDIF
!
!  Choose and apply real Q.
!
460  CONTINUE

     s = abs ( a1 ) + abs ( a2 )
     IF( s == 0.0D+00 ) then
       go to 475
     ENDIF

     u1 = a1 / s
     u2 = a2 / s
     r = sign ( sqrt ( u1**2 + u2**2 ), u1 )
     v1 = -(u1 + r) / r
     v2 = -u2 / r
     u2 = v2 / v1

     DO j = na, n
       t = a(na,j) + u2 * a(en,j)
       a(na,j) = a(na,j) + t * v1
       a(en,j) = a(en,j) + t * v2
       t = b(na,j) + u2 * b(en,j)
       b(na,j) = b(na,j) + t * v1
       b(en,j) = b(en,j) + t * v2
     ENDDO

475  CONTINUE

     a(en,na) = 0.0D+00
     b(en,na) = 0.0D+00
     alfr(na) = a(na,na)
     alfr(en) = a(en,en)
     IF( b(na,na) < 0.0D+00 ) then 
       alfr(na) = -alfr(na)
     ENDIF

     IF( b(en,en) < 0.0D+00 ) then
       alfr(en) = -alfr(en)
     ENDIF

     beta(na) = abs ( b(na,na) )
     beta(en) = abs ( b(en,en) )
     alfi(en) = 0.0D+00
     alfi(na) = 0.0D+00
     go to 505
!
!  Two complex roots.
!
480  CONTINUE

     e = e + c
     ei = sqrt ( -d )
     a11r = a11 - e * b11
     a11i = ei * b11
     a12r = a12 - e * b12
     a12i = ei * b12
     a22r = a22 - e * b22
     a22i = ei * b22

     IF( abs ( a11r ) + abs ( a11i ) + abs ( a12r ) + abs ( a12i ) >= &
            abs ( a21 ) + abs ( a22r ) + abs ( a22i ) ) then
       a1 = a12r
       a1i = a12i
       a2 = -a11r
       a2i = -a11i
     else
       a1 = a22r
       a1i = a22i
       a2 = -a21
       a2i = 0.0D+00
     ENDIF
!
!  Choose complex Z.
!
     cz = sqrt ( a1**2 + a1i**2 )

     IF( cz /= 0.0D+00 ) then
       szr = ( a1 * a2 + a1i * a2i) / cz
       szi = ( a1 * a2i - a1i * a2) / cz
       r = sqrt ( cz**2 + szr**2 + szi**2 )
       cz = cz / r
       szr = szr / r
       szi = szi / r
     else
       szr = 1.0D+00
       szi = 0.0D+00
     ENDIF

     IF( an >= ( abs ( e ) + ei ) * bn ) then
       a1 = cz * b11 + szr * b12
       a1i = szi * b12
       a2 = szr * b22
       a2i = szi * b22
     else
       a1 = cz * a11 + szr * a12
       a1i = szi * a12
       a2 = cz * a21 + szr * a22
       a2i = szi * a22
     ENDIF
!
!  Choose complex Q.
!
     cq = sqrt ( a1**2 + a1i**2 )

     IF( cq /= 0.0D+00 ) then
       sqr = ( a1 * a2 + a1i * a2i ) / cq
       sqi = ( a1 * a2i - a1i * a2 ) / cq
       r = sqrt ( cq**2 + sqr**2 + sqi**2 )
       cq = cq / r
       sqr = sqr / r
       sqi = sqi / r
     else
       sqr = 1.0D+00
       sqi = 0.0D+00
     ENDIF
!
!  Compute diagonal elements that would result if transformations were applied.
!
     ssr = sqr * szr + sqi * szi
     ssi = sqr * szi - sqi * szr
     i = 1
     tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22
     ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
     dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
     di = cq * szi * b12 + ssi * b22
     go to 503

502  CONTINUE

     i = 2
     tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22
     ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
     dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
     di = -ssi * b11 - sqi * cz * b12

503  CONTINUE

     t = ti * dr - tr * di

     IF( t < 0.0D+00 ) then
       j = en
     else
       j = na
     ENDIF

     r = sqrt ( dr**2 + di**2 )
     beta(j) = bn * r
     alfr(j) = an * (tr * dr + ti * di) / r
     alfi(j) = an * t / r

     IF( i == 1 ) then
       go to 502
     ENDIF

505  CONTINUE

     isw = 3 - isw

510  CONTINUE

  ENDDO

  b(n,1) = epsb

  RETURN

END SUBROUTINE

SUBROUTINE qzvec ( n, a, b, alfr, alfi, beta, z )

!*****************************************************************************80
!
!! QZVEC computes eigenvectors for a generalized eigenvalue problem.
!
!  Discussion:
!
!    This SUBROUTINE is the optional fourth step of the QZ algorithm
!    for solving generalized matrix eigenvalue problems.
!
!    This SUBROUTINE accepts a pair of real matrices, one of them in
!    quasi-triangular form (in which each 2-by-2 block corresponds to
!    a pair of complex eigenvalues) and the other in upper triangular
!    form.  It computes the eigenvectors of the triangular problem and
!    transforms the results back to the original coordinate system.
!    it is usually preceded by QZHES, QZIT, and QZVAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains a real upper quasi-triangular 
!    matrix.  Its subdiagonal elements provide information about the storage of
!    the complex eigenvectors.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, a real upper triangular
!    matrix.  In addition, location B(N,1) contains the tolerance quantity EPSB
!    computed and saved in QZIT.  On output, B has been destroyed.
!
!    Input, REAL ( kind = 8 ) ALFR(N), ALFI(N), BETA(N), vectors whose ratios
!      ( ALFR + I * ALFI ) / BETA
!    are the generalized eigenvalues.  They are usually obtained from QZVAL.
!
!    Input/output, REAL ( kind = 8 ) Z(N,N).  On input, the transformation
!    matrix produced in the reductions by QZHES, QZIT, and QZVAL, if performed.
!    If the eigenvectors of the triangular problem are desired, Z must contain
!    the identity matrix.  On output, Z contains the real and imaginary parts of
!    the eigenvectors:
!    If ALFI(I) == 0.0, the I-th eigenvalue is real and the I-th column of Z
!    contains its eigenvector.
!    If ALFI(I) > 0.0, the eigenvalue is the first of a complex pair and the
!    I-th and (I+1)-th columns of Z contain its eigenvector.
!    If ALFI(I) < 0.0, the eigenvalue is the second of a complex pair and the
!    (I-1)-th and I-th columns of Z contain the conjugate of its eigenvector.
!    Each eigenvector is normalized so that the modulus of its largest
!    component is 1.0D+00 .
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) alfi(n)
  REAL ( kind = 8 ) alfm
  REAL ( kind = 8 ) alfr(n)
  REAL ( kind = 8 ) almi
  REAL ( kind = 8 ) almr
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) beta(n)
  REAL ( kind = 8 ) betm
  REAL ( kind = 8 ) d
  REAL ( kind = 8 ) di
  REAL ( kind = 8 ) dr
  INTEGER ( kind = 4 ) en
  INTEGER ( kind = 4 ) enm2
  REAL ( kind = 8 ) epsb
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) isw
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) na
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) q
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) ra
  REAL ( kind = 8 ) rr
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) sa
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) t1
  REAL ( kind = 8 ) t2
  REAL ( kind = 8 ) ti
  REAL ( kind = 8 ) tr
  REAL ( kind = 8 ) w
  REAL ( kind = 8 ) w1
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z(n,n)
  REAL ( kind = 8 ) z1
  REAL ( kind = 8 ) zz

  epsb = b(n,1)
  isw = 1

  DO nn = 1, n

     en = n + 1 - nn
     na = en - 1

     IF( isw == 2 ) then
       go to 795
     ENDIF

     IF( alfi(en) /= 0.0D+00 ) then
       go to 710
     ENDIF
!
!  Real vector.
!
     m = en
     b(en,en) = 1.0D+00

     IF( na == 0 ) then
       go to 800
     ENDIF

     alfm = alfr(m)
     betm = beta(m)

     DO ii = 1, na

        i = en - ii
        w = betm * a(i,i) - alfm * b(i,i)
        r = 0.0D+00

        DO j = m, en
          r = r + ( betm * a(i,j) - alfm * b(i,j) ) * b(j,en)
        ENDDO

        IF( i == 1 .or. isw == 2 ) then
          go to 630
        ENDIF

        IF( betm * a(i,i-1) == 0.0D+00 ) then
          go to 630
        ENDIF

        zz = w
        s = r
        go to 690

630     CONTINUE

        m = i

        IF( isw == 2 ) then
          go to 640
        ENDIF
!
!  Real 1-by-1 block.
!
        IF( w == 0.0D+00 ) then
          t = epsb
        else
          t = w
        ENDIF

        b(i,en) = - r / t
        go to 700
!
!  Real 2-by-2 block.
!
640     CONTINUE

        x = betm * a(i,i+1) - alfm * b(i,i+1)
        y = betm * a(i+1,i)
        q = w * zz - x * y
        t = ( x * s - zz * r ) / q
        b(i,en) = t

        IF( abs ( x ) <= abs ( zz ) ) then
          go to 650
        ENDIF

        b(i+1,en) = (-r - w * t) / x

        go to 690

650     CONTINUE

        b(i+1,en) = (-s - y * t) / zz

690     CONTINUE

        isw = 3 - isw

700     CONTINUE

     ENDDO
!
!  End real vector.
!
     go to 800
!
!  Complex vector.
!
710  CONTINUE

     m = na
     almr = alfr(m)
     almi = alfi(m)
     betm = beta(m)
!
!  Last vector component chosen imaginary so eigenvector matrix is triangular.
!
     y = betm * a(en,na)
     b(na,na) = -almi * b(en,en) / y
     b(na,en) = ( almr * b(en,en) - betm * a(en,en) ) / y
     b(en,na) = 0.0D+00
     b(en,en) = 1.0D+00
     enm2 = na - 1

     DO ii = 1, enm2

        i = na - ii
        w = betm * a(i,i) - almr * b(i,i)
        w1 = -almi * b(i,i)
        ra = 0.0D+00
        sa = 0.0D+00

        DO j = m, en
          x = betm * a(i,j) - almr * b(i,j)
          x1 = -almi * b(i,j)
          ra = ra + x * b(j,na) - x1 * b(j,en)
          sa = sa + x * b(j,en) + x1 * b(j,na)
        ENDDO

        IF( i == 1 .or. isw == 2 ) then
          go to 770
        ENDIF

        IF( betm * a(i,i-1) == 0.0D+00 ) then
          go to 770
        ENDIF

        zz = w
        z1 = w1
        r = ra
        s = sa
        isw = 2
        go to 790
770     CONTINUE

        m = i
        IF( isw == 2 ) then
          go to 780
        ENDIF
!
!  Complex 1-by-1 block.
!
        tr = -ra
        ti = -sa

773     CONTINUE

        dr = w
        di = w1
!
!  Complex divide (t1,t2) = (tr,ti) / (dr,di),
!
775     CONTINUE

        IF( abs ( di ) > abs ( dr ) ) then
          go to 777
        ENDIF

        rr = di / dr
        d = dr + di * rr
        t1 = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        go to ( 787, 782 ), isw

777     CONTINUE

        rr = dr / di
        d = dr * rr + di
        t1 = ( tr * rr + ti ) / d
        t2 = ( ti * rr - tr ) / d
        go to ( 787, 782 ), isw
!
!  Complex 2-by-2 block.
!
780     CONTINUE

        x = betm * a(i,i+1) - almr * b(i,i+1)
        x1 = -almi * b(i,i+1)
        y = betm * a(i+1,i)
        tr = y * ra - w * r + w1 * s
        ti = y * sa - w * s - w1 * r
        dr = w * zz - w1 * z1 - x * y
        di = w * z1 + w1 * zz - x1 * y
        IF( dr == 0.0D+00 .and. di == 0.0D+00 ) then
          dr = epsb
        ENDIF
        go to 775

782     CONTINUE

        b(i+1,na) = t1
        b(i+1,en) = t2
        isw = 1

        IF( abs ( y ) > abs ( w ) + abs ( w1 ) ) then
          go to 785
        ENDIF

        tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
        ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
        go to 773

785     CONTINUE

        t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en) ) / y
        t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na) ) / y

787     CONTINUE

        b(i,na) = t1
        b(i,en) = t2

790     CONTINUE

     ENDDO
!
!  End complex vector.
!
795   CONTINUE

      isw = 3 - isw

800   CONTINUE

  ENDDO
!
!  End back substitution.
!  Transform to original coordinate system.
!
  DO jj = 1, n

     j = n + 1 - jj

     DO i = 1, n

        zz = 0.0D+00

        DO k = 1, j
          zz = zz + z(i,k) * b(k,j)
        ENDDO

        z(i,j) = zz

      ENDDO

  ENDDO
!
!  Normalize so that modulus of largest component of each vector is 1.
!  (ISW is 1 initially from before).
!
  DO j = 1, n

     d = 0.0D+00
     IF( isw == 2 ) then
       go to 920
     ENDIF

     IF( alfi(j) /= 0.0D+00 ) then
       go to 945
     ENDIF

     DO i = 1, n
       d = max ( d, abs ( z(i,j) ) )
     ENDDO

     z(1:n,j) = z(1:n,j) / d

     go to 950

920  CONTINUE

     DO i = 1, n
       r = abs ( z(i,j-1) ) + abs ( z(i,j) )
       IF( r /= 0.0D+00 ) then
         r = r * sqrt ( ( z(i,j-1) / r )**2 + ( z(i,j) / r )**2 )
       ENDIF
       d = max ( d, r )
     ENDDO

     z(1:n,j-1) = z(1:n,j-1) / d
     z(1:n,j) = z(1:n,j) / d

945  CONTINUE

     isw = 3 - isw

950  CONTINUE

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, REAL ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  IMPLICIT NONE

  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z

  z = x
  x = y
  y = z

  RETURN

END SUBROUTINE

SUBROUTINE r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) M, the number of rows in A.
!
!    Input, INTEGER ( kind = 4 ) N, the number of columns in A.
!
!    Input, REAL ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  RETURN

END SUBROUTINE

SUBROUTINE r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, REAL ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, INTEGER ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, INTEGER ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ), parameter :: incx = 5
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) i2hi
  INTEGER ( kind = 4 ) i2lo
  INTEGER ( kind = 4 ) ihi
  INTEGER ( kind = 4 ) ilo
  INTEGER ( kind = 4 ) inc
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) j2
  INTEGER ( kind = 4 ) j2hi
  INTEGER ( kind = 4 ) j2lo
  INTEGER ( kind = 4 ) jhi
  INTEGER ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  DO j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    DO j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    ENDDO

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    DO i = i2lo, i2hi

      DO j2 = 1, inc

        j = j2lo - 1 + j2

        IF( a(i,j) == REAL ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        ENDIF

      ENDDO

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, INTEGER ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, REAL ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ), parameter :: i4_huge = 2147483647
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) seed
  REAL ( kind = 8 ) r(m,n)

  DO j = 1, n

    DO i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      IF( seed < 0 ) then
        seed = seed + i4_huge
      ENDIF

      r(i,j) = REAL ( seed, kind = 8 ) * 4.656612875D-10

    ENDDO
  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the number of components of the vector.
!
!    Input, REAL ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n)
  INTEGER ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  DO i = 1, n
    write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the number of components of the vector.
!
!    Input, REAL ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a1(n)
  REAL ( kind = 8 ) a2(n)
  INTEGER ( kind = 4 ) i
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  IF( all ( a1(1:n) == aint ( a1(1:n) ) ) .and. &
       all ( a2(1:n) == aint ( a2(1:n) ) ) ) then
    DO i = 1, n
      write ( *, '(i8,2i8)' ) i, int ( a1(i) ), int ( a2(i) )
    ENDDO
  else IF( all ( abs ( a1(1:n) ) < 1000000.0D+00 ) .and. &
            all ( abs ( a2(1:n) ) < 1000000.0D+00 ) ) then
    DO i = 1, n
      write ( *, '(i8,2f14.6)' ) i, a1(i), a2(i)
    ENDDO
  else
    DO i = 1, n
      write ( *, '(i8,2g14.6)' ) i, a1(i), a2(i)
    ENDDO
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE ratqr ( n, eps1, d, e, e2, m, w, ind, bd, type, idef, ierr )

!*****************************************************************************80
!
!! RATQR computes selected eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the algebraically smallest or largest
!    eigenvalues of a symmetric tridiagonal matrix by the
!    rational QR method with Newton corrections.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) EPS1.  On input, a theoretical absolute
!    error tolerance for the computed eigenvalues.  If the input EPS1 is
!    non-positive, or indeed smaller than its default value, it is reset at
!    each iteration to the respective default value, namely, the product of
!    the relative machine precision and the magnitude of the current eigenvalue
!    iterate.  The theoretical absolute error in the K-th eigenvalue is usually
!    not greater than K times EPS1.  On output, EPS1 is unaltered unless it has
!    been reset to its (last) default value.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the input matrix.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E2(N).  On input, E2(2:N-1) contains the
!    squares of the corresponding elements of E, and E2(1) is arbitrary.  On
!    output, elements of E2 corresponding to elements of E regarded as
!    negligible have been replaced by zero, causing the matrix to split into
!    a direct sum of submatrices.  E2(1) is set to 0.0D+00 if the smallest
!    eigenvalues have been found, and to 2.0D+00 if the largest eigenvalues
!    have been found.  E2 is otherwise unaltered (unless overwritten by BD).
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvalues to be found.
!
!    Output, REAL ( kind = 8 ) W(M), the M algebraically smallest eigenvalues in
!    ascending order, or the M largest eigenvalues in descending order.
!    If an error exit is made because of an incorrect specification of IDEF,
!    no eigenvalues are found.  If the Newton iterates for a particular
!    eigenvalue are not monotone, the best estimate obtained is RETURNed
!    and IERR is set.  W may coincide with D.
!
!    Outpt, integer IND(N), contains in its first M positions the submatrix
!    indices associated with the corresponding eigenvalues in W:
!    1 for eigenvalues belonging to the first submatrix from the top, 2 for
!    those belonging to the second submatrix, and so on.
!
!    Output, REAL ( kind = 8 ) BD(N), contains refined bounds for the
!    theoretical errors of the corresponding eigenvalues in W.  These bounds
!    are usually within the tolerance specified by EPS1.  BD may coincide
!    with E2.
!
!    Input, INTEGER ( kind = 4 ) IDEF, should be set to 1 if the input matrix
!    is known to be positive definite, to -1 if the input matrix is known to
!    be negative  definite, and to 0 otherwise.
!
!    Input, logical TYPE, should be set to TRUE if the smallest eigenvalues
!    are to be found, and to FALSE if the largest eigenvalues are to be found.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    6*N+1, if IDEF is set to 1 and TYPE to .true. when the matrix is not
!      positive definite, or if IDEF is set to -1 and TYPE to .false.
!      when the matrix is not negative definite,
!    5*N+K, if successive iterates to the K-th eigenvalue are not monotone
!      increasing, where K refers to the last such occurrence.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) bd(n)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) delta
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) ep
  REAL ( kind = 8 ) eps1
  REAL ( kind = 8 ) err
  REAL ( kind = 8 ) f
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) idef
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ind(n)
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jdef
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  !INTEGER ( kind = 4 ) k1
  INTEGER ( kind = 4 ) m
  REAL ( kind = 8 ) p
  REAL ( kind = 8 ) q
  REAL ( kind = 8 ) qp
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) tot
  logical type
  REAL ( kind = 8 ) w(n)

  ierr = 0
  jdef = idef
  w(1:n) = d(1:n)

  IF( .not. type ) then
    j = 1
    go to 400
  ENDIF

40 CONTINUE

  err = 0.0D+00
  s = 0.0D+00
!
!  Look for small sub-diagonal entries and define initial shift
!  from lower Gerschgorin bound.
!
!  Copy E2 array into BD.
!
  tot = w(1)
  q = 0.0D+00
  j = 0

  DO i = 1, n

     p = q

     IF( i == 1 ) then
       go to 60
     ENDIF

     IF( p > ( abs ( d(i) ) + abs (  d(i-1) ) ) * epsilon ( p ) ) then
       go to 80
     ENDIF

60   CONTINUE

     e2(i) = 0.0D+00

80   CONTINUE

     bd(i) = e2(i)
!
!  Count also if element of E2 has underflowed.
!
     IF( e2(i) == 0.0D+00 ) then
       j = j + 1
     ENDIF

     ind(i) = j
     q = 0.0D+00
     IF( i /= n ) then
       q = abs ( e(i+1) )
     ENDIF

     tot = min ( w(i) - p - q, tot )

  ENDDO

  IF( jdef == 1 .and. tot < 0.0D+00 ) then
    go to 140
  ENDIF

  w(1:n) = w(1:n) - tot

  go to 160

140 CONTINUE

  tot = 0.0D+00

160 CONTINUE

  DO k = 1, m
!
!  Next QR transformation.
!
180  CONTINUE

     tot = tot + s
     delta = w(n) - s
     i = n
     f = abs ( tot ) * epsilon ( f )
     IF( eps1 < f ) then
       eps1 = f
     ENDIF

     IF( delta > eps1 ) then
       go to 190
     ENDIF

     IF( delta < (-eps1) ) then
       ierr = 6 * n + 1
       RETURN
     ENDIF

     go to 300
!
!  Replace small sub-diagonal squares by zero to reduce the incidence of
!  underflows.
!
190  CONTINUE

     DO j = k + 1, n
       IF( bd(j) <= ( abs (  w(j) + w(j-1) ) * epsilon ( bd(j) ) ) ** 2 ) then
         bd(j) = 0.0D+00
       ENDIF
     ENDDO

     f = bd(n) / delta
     qp = delta + f
     p = 1.0D+00

     DO ii = 1, n - k

       i = n - ii
       q = w(i) - s - f
       r = q / qp
       p = p * r + 1.0D+00
       ep = f * r
       w(i+1) = qp + ep
       delta = q - ep

       IF( delta > eps1 ) then
         go to 220
       ENDIF

       IF( delta < (-eps1) ) then
         ierr = 6 * n + 1
         RETURN
       ENDIF

       go to 300

220    CONTINUE

       f = bd(i) / q
       qp = delta + f
       bd(i+1) = qp * ep

     ENDDO

     w(k) = qp
     s = qp / p

     IF( tot + s > tot ) then
       go to 180
     ENDIF
!
!  Set error: irregular end of iteration.
!  Deflate minimum diagonal element.
!
     ierr = 5 * n + k
     s = 0.0D+00
     delta = qp

     DO j = k, n
       IF( w(j) <= delta ) then
         i = j
         delta = w(j)
       ENDIF
     ENDDO
!
!  Convergence.
!
300  CONTINUE

     IF( i < n ) then
       bd(i+1) = bd(i) * f / qp
     ENDIF

     ii = ind(i)

     DO jj = 1, i - k
       j = i - jj
       w(j+1) = w(j) - s
       bd(j+1) = bd(j)
       ind(j+1) = ind(j)
     ENDDO

     w(k) = tot
     err = err + abs ( delta)
     bd(k) = err
     ind(k) = ii

  ENDDO

  IF( type ) then
    RETURN
  ENDIF

  f = bd(1)
  e2(1) = 2.0D+00
  bd(1) = f
  j = 2
!
!  Negate elements of W for largest values.
!
400 CONTINUE

  w(1:n) = - w(1:n)
  jdef = -jdef

  IF( j == 1 ) then
    go to 40
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rebak ( n, b, dl, m, z )

!*****************************************************************************80
!
!! REBAK determines eigenvectors by undoing the REDUC transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a generalized
!    symmetric eigensystem by back transforming those of the
!    derived symmetric matrix determined by REDUC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) B(N,N), contains information about the similarity
!    transformation (Cholesky decomposition) used in the reduction by REDUC
!    in its strict lower triangle.
!
!    Input, REAL ( kind = 8 ) DL(N), further information about the 
!    transformation.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the eigenvectors to be
!    back transformed in its first M columns.  On output, the transformed
!    eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) dl(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  REAL ( kind = 8 ) z(n,m)

  DO j = 1, m
    DO i = n, 1, -1
      z(i,j) = ( z(i,j) - dot_product ( b(i+1:n,i), z(i+1:n,j) ) ) / dl(i)
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE rebakb ( n, b, dl, m, z )

!*****************************************************************************80
!
!! REBAKB determines eigenvectors by undoing the REDUC2 transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a generalized
!    symmetric eigensystem by back transforming those of the
!    derived symmetric matrix determined by REDUC2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) B(N,N), contains information about the similarity
!    transformation (Cholesky decomposition) used in the reduction by REDUC2
!    in its strict lower triangle.
!
!    Input, REAL ( kind = 8 ) DL(N), further information about the 
!    transformation.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the eigenvectors to be 
!    back transformed in its first M columns.  On output, the transformed
!    eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) dl(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  REAL ( kind = 8 ) z(n,m)

  DO j = 1, m

    DO i = n, 1, -1

      z(i,j) = dl(i) * z(i,j) + dot_product ( b(i,1:i-1), z(1:i-1,j) )

    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE reduc ( n, a, b, dl, ierr )

!*****************************************************************************80
!
!! REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.
!
!  Discussion:
!
!    This SUBROUTINE reduces the generalized symmetric eigenproblem
!    ax=(lambda)bx, where B is positive definite, to the standard
!    symmetric eigenproblem using the Cholesky factorization of B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.  If the
!    Cholesky factor L of B is already available, N should be prefixed with a
!    minus sign.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, A contains a real
!    symmetric matrix.  Only the full upper triangle of the matrix need be
!    supplied.  On output, A contains in its full lower triangle the full lower
!    triangle of the symmetric matrix derived from the reduction to the
!    standard form.  The strict upper triangle of a is unaltered.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, the real symmetric
!    input matrix.  Only the full upper triangle of the matrix need be supplied.
!    If N is negative, the strict lower triangle of B contains, instead, the
!    strict lower triangle of its Cholesky factor L.  In any case, on output,
!    B contains in its strict lower triangle the strict lower triangle of
!    its Cholesky factor L.  The full upper triangle of B is unaltered.
!
!    Input/output, REAL ( kind = 8 ) DL(N).  If N is negative, then the DL
!    contains the diagonal elements of L on input.  In any case, DL will contain
!    the diagonal elements of L on output,
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    7*N+1, if B is not positive definite.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) dl(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y

  ierr = 0
  nn = abs ( n )
!
!  Form L in the arrays B and DL.
!
  DO i = 1, n

     DO j = i, n

        x = b(i,j)

        DO k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        ENDDO

        IF( j == i ) then

          IF( x <= 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'REDUC - Fatal error!'
            write ( *, '(a)' ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            RETURN
          ENDIF

          y = sqrt ( x )
          dl(i) = y
        else
          b(j,i) = x / y
        ENDIF

    ENDDO

  ENDDO
!
!  Form the transpose of the upper triangle of INV(L)*A
!  in the lower triangle of the array A.
!
  DO i = 1, nn

     y = dl(i)

     DO j = i, nn

        x = a(i,j)

        DO k = 1, i - 1
          x = x - b(i,k) * a(j,k)
        ENDDO

        a(j,i) = x / y

      ENDDO

  ENDDO
!
!  Pre-multiply by INV(L) and overwrite.
!
  DO j = 1, nn

     DO i = j, nn

        x = a(i,j)

        DO k = j, i - 1
          x = x - a(k,j) * b(i,k)
        ENDDO

        DO k = 1, j - 1
          x = x - a(j,k) * b(i,k)
        ENDDO

        a(i,j) = x / dl(i)

    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE reduc2 ( n, a, b, dl, ierr )

!*****************************************************************************80
!
!! REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.
!
!  Discussion:
!
!    This SUBROUTINE reduces the generalized symmetric eigenproblems
!    abx=(lambda)x or bay=(lambda)y, where B is positive definite,
!    to the standard symmetric eigenproblem using the Cholesky
!    factorization of B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.  If the
!    Cholesky factor L of B is already available, N should be prefixed with a
!    minus sign.
!
!    Input/output, REAL ( kind = 8 ) A(N,N).  On input, A contains a real
!    symmetric matrix.  Only the full upper triangle of the matrix need be
!    supplied.  On output, A contains in its full lower triangle the full lower
!    triangle of the symmetric matrix derived from the reduction to the
!    standard form.  The strict upper triangle of a is unaltered.
!
!    Input/output, REAL ( kind = 8 ) B(N,N).  On input, the real symmetric
!    input matrix.  Only the full upper triangle of the matrix need be supplied.
!    If N is negative, the strict lower triangle of B contains, instead, the
!    strict lower triangle of its Cholesky factor L.  In any case, on output,
!    B contains in its strict lower triangle the strict lower triangle of
!    its Cholesky factor L.  The full upper triangle of B is unaltered.
!
!    Input/output, REAL ( kind = 8 ) DL(N).  If N is negative, then the DL
!    contains the diagonal elements of L on input.  In any case, DL will contain
!    the diagonal elements of L on output,
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    7*N+1, if B is not positive definite.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) dl(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) nn
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y

  ierr = 0
  nn = abs ( n )
!
!  Form L in the arrays B and DL.
!
  DO i = 1, n

     DO j = i, n

        x = b(i,j)

        DO k = 1, i - 1
          x = x - b(i,k) * b(j,k)
        ENDDO

        IF( j == i ) then

          IF( x <= 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'REDUC2 - Fatal error!'
            write ( *, '(a)' ) '  The matrix is not positive definite.'
            ierr = 7 * n + 1
            RETURN
          ENDIF

          y = sqrt ( x )
          dl(i) = y

        else

          b(j,i) = x / y

        ENDIF

    ENDDO

  ENDDO
!
!  Form the lower triangle of A*L in the lower triangle of A.
!
  DO i = 1, nn

     DO j = 1, i

        x = a(j,i) * dl(j)

        DO k = j + 1, i
          x = x + a(k,i) * b(k,j)
        ENDDO

        DO k = i + 1, nn
          x = x + a(i,k) * b(k,j)
        ENDDO

        a(i,j) = x

     ENDDO

  ENDDO
!
!  Pre-multiply by L' and overwrite.
!
  DO i = 1, nn

    y = dl(i)

    DO j = 1, i

      x = y * a(i,j)

      DO k = i + 1, nn
        x = x + a(k,j) * b(k,i)
      ENDDO

      a(i,j) = x

    ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE rg ( n, a, wr, wi, matz, z, ierr )

!*****************************************************************************80
!
!! RG computes eigenvalues and eigenvectors of a real general matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) A(N,N), the real general matrix.  On 
!    output, A has been overwritten.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) WR(N), WI(N), the real and imaginary parts,
!    respectively, of the eigenvalues.  Complex conjugate pairs of eigenvalues
!    appear consecutively with the eigenvalue having the positive imaginary
!    part first.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the real and imaginary parts of 
!    the eigenvectors if MATZ is not zero.  If the J-th eigenvalue is real, the
!    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
!    complex with positive imaginary part, the J-th and (J+1)-th columns of
!    Z contain the real and imaginary parts of its eigenvector.  The
!    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
!
!    Output, INTEGER ( kind = 4 ) IERR, an error completion code described in 
!    the documentation for HQR and HQR2.  The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) fv1(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) is1
  INTEGER ( kind = 4 ) is2
  INTEGER ( kind = 4 ) iv1(n)
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) wi(n)
  REAL ( kind = 8 ) wr(n)
  REAL ( kind = 8 ) z(n,n)

  call balanc ( n, a, is1, is2, fv1 )

  call elmhes ( n, is1, is2, a, iv1 )

  IF( matz == 0 ) then

    call hqr ( n, is1, is2, a, wr, wi, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

  else

    call eltran ( n, is1, is2, a, iv1, z )

    call hqr2 ( n, is1, is2, a, wr, wi, z, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

    call balbak ( n, is1, is2, fv1, n, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rgg ( n, a, b, alfr, alfi, beta, matz, z, ierr )

!*****************************************************************************80
!
!! RGG: eigenvalues/vectors for the generalized problem A*x = lambda*B*x.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real general generalized eigenproblem
!
!      A * x = lambda * B * x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.
!
!    Input/output, REAL ( kind = 8 ) A(N,N), B(N,N), the two real general 
!    matrices.  On output, A and B have been overwritten.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) ALFR(N), ALFI(N), the real and imaginary parts,
!    respectively, of the numerators of the eigenvalues.
!
!    Output, REAL ( kind = 8 ) BETA(N), the denominators of the eigenvalues,
!    which are thus given by the ratios (ALFR + I * ALFI ) / BETA.
!    Complex conjugate pairs of eigenvalues appear consecutively
!    with the eigenvalue having the positive imaginary part first.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the real and imaginary parts of 
!    the eigenvectors if MATZ is not zero.  If the J-th eigenvalue is real, the
!    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
!    complex with positive imaginary part, the J-th and (J+1)-th columns of
!    Z contain the real and imaginary parts of its eigenvector.  The
!    conjugate of this vector is the eigenvector for the conjugate eigenvalue.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set equal to an error completion code
!    described in the documentation for QZIT.  The normal completion
!    code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) alfi(n)
  REAL ( kind = 8 ) alfr(n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) beta(n)
  REAL ( kind = 8 ) eps1
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  logical tf
  REAL ( kind = 8 ) z(n,n)

  eps1 = 0.0D+00

  IF( matz == 0 ) then
    tf = .false.
  else
    tf = .true.
  ENDIF

  call qzhes ( n, a, b, tf, z )

  call qzit ( n, a, b, eps1, tf, z, ierr )

  IF( ierr /= 0 ) then
    RETURN
  ENDIF

  call qzval ( n, a, b, alfr, alfi, beta, tf, z )

  IF( matz /= 0 ) then
    call qzvec ( n, a, b, alfr, alfi, beta, z )
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rs ( n, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RS computes eigenvalues and eigenvectors of real symmetric matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), the real symmetric matrix.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set equal to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  IF( matz == 0 ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsb ( n, mb, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric band matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) MB, the half band width of the matrix,
!    defined as the number of adjacent diagonals, including the principal
!    diagonal, required to specify the non-zero portion of the lower triangle
!    of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,MB), contains the lower triangle of the real
!    symmetric band matrix.  Its lowest subdiagonal is stored in the last N+1-MB
!    positions of the first column, its next subdiagonal in the last
!    N+2-MB positions of the second column, further subdiagonals similarly,
!    and finally its principal diagonal in the N positions of the last
!    column.  Contents of storages not part of the matrix are arbitrary.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired,
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mb
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,mb)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  logical tf
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  IF( mb <= 0 ) then
    ierr = 12 * n
    RETURN
  ENDIF

  IF( n < mb ) then
    ierr = 12 * n
    RETURN
  ENDIF

  IF( matz == 0 ) then

    tf = .false.

    call bandr ( n, mb, a, w, fv1, fv2, tf, z )

    call tqlrat ( n, w, fv2, ierr )

  else

    tf = .true.

    call bandr ( n, mb, a, w, fv1, fv1, tf, z )

    call tql2 ( n, w, fv1, z, ierr )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsg ( n, a, b, w, matz, z, ierr )

!*****************************************************************************80
!
!! RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real symmetric generalized eigenproblem  ax = (lambda)bx.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Modified:
!
!    04 February 2003
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains a real symmetric matrix.
!
!    Input, REAL ( kind = 8 ) B(N,N), contains a positive definite real
!    symmetric matrix.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  call reduc ( n, a, b, fv2, ierr )

  IF( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSG - Fatal error!'
    write ( *, '(a)' ) '  Error RETURN from REDUC.'
    RETURN
  ENDIF

  IF( matz == 0 ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSG - Warning!'
      write ( *, '(a)' ) '  Error RETURN from TQLRAT!'
      RETURN
    ENDIF

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSG - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from TQL2!'
      RETURN
    ENDIF

    call rebak ( n, b, fv2, n, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsgab ( n, a, b, w, matz, z, ierr )

!*****************************************************************************80
!
!! RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real symmetric generalized eigenproblem  abx = (lambda)x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains a real symmetric matrix.
!
!    Input, REAL ( kind = 8 ) B(N,N), contains a positive definite real
!    symmetric matrix.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ 
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  call reduc2 ( n, a, b, fv2, ierr )

  IF( ierr /= 0 ) then
    RETURN
  ENDIF

  IF( matz == 0 ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

    call rebak ( n, b, fv2, n, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsgba ( n, a, b, w, matz, z, ierr )

!*****************************************************************************80
!
!! RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    for the real symmetric generalized eigenproblem:
!
!      B * A * x = lambda * x
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrices A and B.
!
!    Input, REAL ( kind = 8 ) A(N,N), a real symmetric matrix.
!
!    Input, REAL ( kind = 8 ) B(N,N), a positive definite symmetric matrix.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired, 
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) b(n,n)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  call reduc2 ( n, a, b, fv2, ierr )

  IF( ierr /= 0 ) then
    RETURN
  ENDIF

  IF( matz == 0 ) then

    call tred1 ( n, a, w, fv1, fv2 )

    call tqlrat ( n, w, fv2, ierr )

  else

    call tred2 ( n, a, w, fv1, z )

    call tql2 ( n, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      RETURN
    ENDIF

    call rebakb ( n, b, fv2, n, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsm ( n, a, w, m, z, ierr )

!*****************************************************************************80
!
!! RSM computes eigenvalues, some eigenvectors, real symmetric matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find all of the eigenvalues and some of the eigenvectors
!    of a real symmetric matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), the symmetric matrix.
!
!    Input, INTEGER ( kind = 4 ) M, specifies the number of eigenvectors to 
!    compute.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,M), contains the orthonormal eigenvectors
!    associated with the first M eigenvalues.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT, IMTQLV and
!    TINVIT.  The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) fwork1(n)
  REAL ( kind = 8 ) fwork2(n)
  REAL ( kind = 8 ) fwork3(n)
  !REAL ( kind = 8 ) fwork4(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) iwork(n)
  INTEGER ( kind = 4 ) k1
  INTEGER ( kind = 4 ) k2
  INTEGER ( kind = 4 ) k3
  INTEGER ( kind = 4 ) k4
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,m)

  k1 = 1
  k2 = k1 + n
  k3 = k2 + n
  k4 = k3 + n

  IF( m <= 0 ) then

    call tred1 ( n, a, w, fwork1, fwork2 )

    call tqlrat ( n, w, fwork2, ierr )

  else

    call tred1 ( n, a, fwork1, fwork2, fwork3 )

    call imtqlv ( n, fwork1, fwork2, fwork3, w, iwork, ierr )

    call tinvit ( n, fwork1, fwork2, fwork3, m, w, iwork, z, ierr )

    call trbak1 ( n, a, fwork2, m, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rsp ( n, nv, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RSP computes eigenvalues and eigenvectors of real symmetric packed matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of
!    SUBROUTINEs from the eigensystem SUBROUTINE package (eispack)
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric packed matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) NV, the dimension of the array A, which
!    must be at least (N*(N+1))/2.
!
!    Input, REAL ( kind = 8 ) A(NV), contains the lower triangle of the
!    real symmetric packed matrix stored row-wise.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired,
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ is 
!    nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for TQLRAT and TQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n
  INTEGER ( kind = 4 ) nv

  REAL ( kind = 8 ) a(nv)
  REAL ( kind = 8 ) fv1(n)
  REAL ( kind = 8 ) fv2(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  IF( ( n * ( n + 1 ) ) / 2 > nv ) then
    ierr = 20 * n
    RETURN
  ENDIF

  call tred3 ( n, nv, a, w, fv1, fv2 )

  IF( matz == 0 ) then

    call tqlrat ( n, w, fv2, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSP - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from TQLRAT.'
      RETURN
    ENDIF

  else

    z(1:n,1:n) = 0.0D+00

    DO i = 1, n
      z(i,i) = 1.0D+00
    ENDDO

    call tql2 ( n, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSP - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from TQL2.'
      RETURN
    ENDIF

    call trbak3 ( n, nv, a, n, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rspp ( n, nv, a, w, matz, z, ierr, m, type )

!*****************************************************************************80
!
!! RSPP computes some eigenvalues/vectors, real symmetric packed matrix.
!
!  Discussion:
!
!    This routine calls the appropriate routines for the following problem:
!
!    Given a symmetric matrix A, which is stored in a packed mode, find
!    the M smallest or largest eigenvalues, and corresponding eigenvectors.
!
!    The routine RSP RETURNs all eigenvalues and eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of A, the number of rows and
!    columns in the original matrix.
!
!    Input, INTEGER ( kind = 4 ) NV, is the of the array A as specified in the
!    calling program.  NV must not be less than N*(N+1)/2.
!
!    Input, REAL ( kind = 8 ) A((N*(N+1))/2), on input the lower triangle of the
!    real symmetric matrix, stored row-wise in the vector,
!    in the order A(1,1), / A(2,1), A(2,2), / A(3,1), A(3,2), A(3,3)/
!    and so on.
!
!    Output, REAL ( kind = 8 ) W(M), the eigenvalues requested.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is set to 0 if only eigenvalues are
!    desired.  Otherwise it is set to any non-zero integer for both eigenvalues
!    and eigenvectors.
!
!    Output, REAL ( kind = 8 ) Z(N,M), the eigenvectors.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag from RATQR.  IERR=0 on
!    normal RETURN.  IERR nonzero, in this case, means that the algorithm broke
!    down while computing an eigenvalue.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvalues to be found.
!
!    Input, logical TYPE, set to .true. if the smallest eigenvalues
!    are to be found, or .false. if the largest ones are sought.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n
  INTEGER ( kind = 4 ) nv

  REAL ( kind = 8 ) a(nv)
  REAL ( kind = 8 ) bd(n)
  REAL ( kind = 8 ) eps1
  INTEGER ( kind = 4 ) idef
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) iwork(n)
  INTEGER ( kind = 4 ) matz
  logical type
  REAL ( kind = 8 ) w(m)
  REAL ( kind = 8 ) work1(n)
  REAL ( kind = 8 ) work2(n)
  REAL ( kind = 8 ) work3(n)
  REAL ( kind = 8 ) z(n,m)
!
!  IDEF =
!    -1 if the matrix is known to be negative definite,
!    +1 if the matrix is known to be positive definite, or
!    0 otherwise.
!
  idef = 0
!
!  Reduce to symmetric tridiagonal form.
!
  call tred3 ( n, nv, a, work1, work2, work3 )
!
!  Find the eigenvalues.
!
  eps1 = 0.0D+00

  call ratqr ( n, eps1, work1, work2, work3, m, w, iwork, &
    bd, type, idef, ierr )

  IF( ierr /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSPP - Fatal error!'
    write ( *, '(a)' ) '  Error RETURN from RATQR.'
    RETURN
  ENDIF
!
!  Find eigenvectors for the first M eigenvalues.
!
  IF( matz /= 0 ) then

    call tinvit ( n, work1, work2, work3, m, w, iwork, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RSPP - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from TINVIT.'
      RETURN
    ENDIF
!
!  Reverse the transformation.
!
    call trbak3 ( n, nv, a, m, z )

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rst ( n, w, e, matz, z, ierr )

!*****************************************************************************80
!
!! RST computes eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of SUBROUTINEs
!    to find the eigenvalues and eigenvectors (if desired)
!    of a real symmetric tridiagonal matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) W(N).  On input, the diagonal elements
!    of the real symmetric tridiagonal matrix.  On output, the eigenvalues in
!    ascending order.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the matrix in
!    E(2:N).  E(1) is arbitrary.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is zero if only eigenvalues are desired,
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for IMTQL1 and IMTQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) e(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  IF( matz == 0 ) then

    call imtql1 ( n, w, e, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RST - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from IMTQL1.'
      RETURN
    ENDIF

  else

    z(1:n,1:n) = 0.0D+00

    DO i = 1, n
      z(i,i) = 1.0D+00
    ENDDO

    call imtql2 ( n, w, e, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RST - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from IMTQL2.'
      RETURN
    ENDIF

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE rt ( n, a, w, matz, z, ierr )

!*****************************************************************************80
!
!! RT computes eigenvalues/vectors, real sign-symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE calls the recommended sequence of SUBROUTINEs
!    to find the eigenvalues and eigenvectors (if desired)
!    of a special real tridiagonal matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains the special real tridiagonal
!    matrix in its first three columns.  The subdiagonal elements are stored
!    in the last N-1 positions of the first column, the diagonal elements
!    in the second column, and the superdiagonal elements in the first N-1
!    positions of the third column.  Elements A(1,1) and A(N,3) are arbitrary.
!
!    Input, INTEGER ( kind = 4 ) MATZ, is 0 if only eigenvalues are desired,
!    and nonzero if both eigenvalues and eigenvectors are desired.
!
!    Output, REAL ( kind = 8 ) W(N), the eigenvalues in ascending order.
!
!    Output, REAL ( kind = 8 ) Z(N,N), contains the eigenvectors, if MATZ
!    is nonzero.
!
!    Output, INTEGER ( kind = 4 ) IERR, is set to an error
!    completion code described in the documentation for IMTQL1 and IMTQL2.
!    The normal completion code is zero.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) fv1(n)
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) matz
  REAL ( kind = 8 ) a(n,3)
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) z(n,n)

  IF( matz == 0 ) then

    call figi ( n, a, w, fv1, fv1, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from FIGI.'
      RETURN
    ENDIF

    call imtql1 ( n, w, fv1, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from IMTQL1.'
      RETURN
    ENDIF

  else

    call figi2 ( n, a, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from FIGI2.'
      RETURN
    ENDIF

    call imtql2 ( n, w, fv1, z, ierr )

    IF( ierr /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RT - Fatal error!'
      write ( *, '(a)' ) '  Error RETURN from IMTQL2.'
      RETURN
    ENDIF

  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE svd ( m, n, a, w, matu, u, matv, v, ierr )

!*****************************************************************************80
!
!! SVD computes the singular value decomposition for a real matrix.
!
!  Discussion:
!
!    This SUBROUTINE determines the singular value decomposition
!
!      A = U * S * V'
!
!    of a real M by N rectangular matrix.  Householder bidiagonalization
!    and a variant of the QR algorithm are used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Golub, Christian Reinsch,
!    Numerische Mathematik,
!    Volume 14, 1970, pages 403-420.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) M, the number of rows of A and U.
!
!    Input, INTEGER ( kind = 4 ) N, the number of columns of A and U, and
!    the order of V.
!
!    Input, REAL ( kind = 8 ) A(M,N), the M by N matrix to be decomposed.
!
!    Output, REAL ( kind = 8 ) W(N), the singular values of A.  These are the
!    diagonal elements of S.  They are unordered.  If an error exit is
!    made, the singular values should be correct for indices
!    IERR+1, IERR+2,..., N.
!
!    Input, logical MATU, should be set to TRUE if the U matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, REAL ( kind = 8 ) U(M,N), contains the matrix U, with orthogonal
!    columns, of the decomposition, if MATU has been set to TRUE.  Otherwise
!    U is used as a temporary array.  U may coincide with A.
!    If an error exit is made, the columns of U corresponding
!    to indices of correct singular values should be correct.
!
!    Input, logical MATV, should be set to TRUE if the V matrix in the
!    decomposition is desired, and to FALSE otherwise.
!
!    Output, REAL ( kind = 8 ) V(N,N), the orthogonal matrix V of the
!    decomposition if MATV has been set to TRUE.  Otherwise V is not referenced.
!    V may also coincide with A if U is not needed.  If an error
!    exit is made, the columns of V corresponding to indices of
!    correct singular values should be correct.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    K, if the K-th singular value has not been determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(m,n)
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) i1
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) kk
  INTEGER ( kind = 4 ) k1
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) ll
  INTEGER ( kind = 4 ) l1
  logical matu
  logical matv
  INTEGER ( kind = 4 ) mn
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) u(m,n)
  REAL ( kind = 8 ) v(n,n)
  REAL ( kind = 8 ) w(n)
  REAL ( kind = 8 ) x
  REAL ( kind = 8 ) y
  REAL ( kind = 8 ) z

  ierr = 0
  u(1:m,1:n) = a(1:m,1:n)
!
!  Householder reduction to bidiagonal form.
!
  g = 0.0D+00
  scale = 0.0D+00
  x = 0.0D+00

  DO i = 1, n

    l = i + 1
    rv1(i) = scale * g
    g = 0.0D+00
    s = 0.0D+00
    scale = 0.0D+00

    IF( i <= m ) then

      scale = sum ( abs ( u(i:m,i) ) )

      IF( scale /= 0.0D+00 ) then

        u(i:m,i) = u(i:m,i) / scale

        s = sum ( u(i:m,i)**2 )

        f = u(i,i)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        u(i,i) = f - g

        IF( i /= n ) then

          DO j = l, n
            s = dot_product ( u(i:m,i), u(i:m,j) )
            u(i:m,j) = u(i:m,j) + s * u(i:m,i) / h
          ENDDO

        ENDIF

        u(i:m,i) = scale * u(i:m,i)

      ENDIF

    ENDIF

    w(i) = scale * g
    g = 0.0D+00
    s = 0.0D+00
    scale = 0.0D+00

    IF( i <= m .and. i /= n ) then

      scale = sum ( abs ( u(i,l:n) ) )

      IF( scale /= 0.0D+00 ) then

        u(i,l:n) = u(i,l:n) / scale
        s = sum ( u(i,l:n)**2 )
        f = u(i,l)
        g = - sign ( sqrt ( s ), f )
        h = f * g - s
        u(i,l) = f - g
        rv1(l:n) = u(i,l:n) / h

        IF( i /= m ) then

          DO j = l, m

            s = dot_product ( u(j,l:n), u(i,l:n) )

            u(j,l:n) = u(j,l:n) + s * rv1(l:n)

          ENDDO

        ENDIF

        u(i,l:n) = scale * u(i,l:n)

      ENDIF

    ENDIF

    x = max ( x, abs ( w(i) ) + abs ( rv1(i) ) )

  ENDDO
!
!  Accumulation of right-hand transformations.
!
  IF( matv ) then

    DO i = n, 1, -1

      IF( i /= n ) then

         IF( g /= 0.0D+00 ) then

          v(l:n,i) = ( u(i,l:n) / u(i,l) ) / g

          DO j = l, n

            s = dot_product ( u(i,l:n), v(l:n,j) )

            v(l:n,j) = v(l:n,j) + s * v(l:n,i)

          ENDDO

        ENDIF

        v(i,l:n) = 0.0D+00
        v(l:n,i) = 0.0D+00

      ENDIF

      v(i,i) = 1.0D+00
      g = rv1(i)
      l = i

    ENDDO

  ENDIF
!
!  Accumulation of left-hand transformations.
!
  IF( matu ) then

    mn = min ( m, n )

    DO i = min ( m, n ), 1, -1

      l = i + 1
      g = w(i)

      IF( i /= n ) then
        u(i,l:n) = 0.0D+00
      ENDIF

      IF( g /= 0.0D+00 ) then

        IF( i /= mn ) then

          DO j = l, n
            s = dot_product ( u(l:m,i), u(l:m,j) )
            f = ( s / u(i,i) ) / g
            u(i:m,j) = u(i:m,j) + f * u(i:m,i)
          ENDDO

        ENDIF

        u(i:m,i) = u(i:m,i) / g

      else

        u(i:m,i) = 0.0D+00

      ENDIF

      u(i,i) = u(i,i) + 1.0D+00

    ENDDO

  ENDIF
!
!  Diagonalization of the bidiagonal form.
!
  tst1 = x

  DO kk = 1, n

     k1 = n - kk
     k = k1 + 1
     its = 0
!
!  Test for splitting.
!
520  CONTINUE

     DO ll = 1, k

       l1 = k - ll
       l = l1 + 1
       tst2 = tst1 + abs ( rv1(l) )

       IF( tst2 == tst1 ) then
         go to 565
       ENDIF

       tst2 = tst1 + abs ( w(l1) )

       IF( tst2 == tst1 ) then
         exit
       ENDIF

     ENDDO
!
!  Cancellation of rv1(l) if L greater than 1.
!
     c = 0.0D+00
     s = 1.0D+00

     DO i = l, k

       f = s * rv1(i)
       rv1(i) = c * rv1(i)
       tst2 = tst1 + abs ( f )

       IF( tst2 == tst1 ) then
         go to 565
       ENDIF

       g = w(i)
       h = pythag ( f, g )
       w(i) = h
       c = g / h
       s = -f / h

       IF( matu ) then

         DO j = 1, m
           y = u(j,l1)
           z = u(j,i)
           u(j,l1) = y * c + z * s
           u(j,i) = -y * s + z * c
         ENDDO

       ENDIF

    ENDDO
!
!  Test for convergence.
!
565 CONTINUE

    z = w(k)

    IF( l == k ) then
      go to 650
    ENDIF
!
!  Shift from bottom 2 by 2 minor.
!
    IF( its >= 30 ) then
      ierr = k
      RETURN
    ENDIF

    its = its + 1
    x = w(l)
    y = w(k1)
    g = rv1(k1)
    h = rv1(k)
    f = 0.5D+00 * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y )
    g = pythag ( f, 1.0D+00 )
    f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + sign ( g, f ) ) - h )
!
!  Next QR transformation.
!
    c = 1.0D+00
    s = 1.0D+00

    DO i1 = l, k1

      i = i1 + 1
      g = rv1(i)
      y = w(i)
      h = s * g
      g = c * g
      z = pythag ( f, h )
      rv1(i1) = z
      c = f / z
      s = h / z
      f = x * c + g * s
      g = -x * s + g * c
      h = y * s
      y = y * c

      IF( matv ) then

        DO j = 1, n
          x = v(j,i1)
          z = v(j,i)
          v(j,i1) = x * c + z * s
          v(j,i) = -x * s + z * c
        ENDDO

      ENDIF

      z = pythag ( f, h )
      w(i1) = z
!
!  Rotation can be arbitrary if Z is zero.
!
      IF( z /= 0.0D+00 ) then
        c = f / z
        s = h / z
      ENDIF

      f = c * g + s * y
      x = - s * g + c * y

      IF( matu ) then

        DO j = 1, m
          y = u(j,i1)
          z = u(j,i)
          u(j,i1) = y * c + z * s
          u(j,i) = - y * s + z * c
        ENDDO

      ENDIF

    ENDDO

    rv1(l) = 0.0D+00
    rv1(k) = f
    w(k) = x
    go to 520
!
!  Convergence.
!
650 CONTINUE

    IF( z <= 0.0D+00 ) then

      w(k) = - z

      IF( matv ) then
        v(1:n,k) = - v(1:n,k)
      ENDIF

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE


! SUBROUTINE timestamp ( )
! 
! !*****************************************************************************80
! !
! !! TIMESTAMP prints the current YMDHMS date as a time stamp.
! !
! !  Example:
! !
! !    31 May 2001   9:45:54.872 AM
! !
! !  Licensing:
! !
! !    This code is distributed under the GNU LGPL license.
! !
! !  Modified:
! !
! !    18 May 2013
! !
! !  Author:
! !
! !    John Burkardt
! !
! !  Parameters:
! !
! !    None
! !
!   IMPLICIT NONE
! 
!   character ( len = 8 ) ampm
!   INTEGER ( kind = 4 ) d
!   INTEGER ( kind = 4 ) h
!   INTEGER ( kind = 4 ) m
!   INTEGER ( kind = 4 ) mm
!   character ( len = 9 ), parameter, dimension(12) :: month = (/ &
!     'January  ', 'February ', 'March    ', 'April    ', &
!     'May      ', 'June     ', 'July     ', 'August   ', &
!     'September', 'October  ', 'November ', 'December ' /)
!   INTEGER ( kind = 4 ) n
!   INTEGER ( kind = 4 ) s
!   INTEGER ( kind = 4 ) values(8)
!   INTEGER ( kind = 4 ) y
! 
!   call date_and_time ( values = values )
! 
!   y = values(1)
!   m = values(2)
!   d = values(3)
!   h = values(5)
!   n = values(6)
!   s = values(7)
!   mm = values(8)
! 
!   IF( h < 12 ) then
!     ampm = 'AM'
!   else IF( h == 12 ) then
!     IF( n == 0 .and. s == 0 ) then
!       ampm = 'Noon'
!     else
!       ampm = 'PM'
!     ENDIF
!   else
!     h = h - 12
!     IF( h < 12 ) then
!       ampm = 'PM'
!     else IF( h == 12 ) then
!       IF( n == 0 .and. s == 0 ) then
!         ampm = 'Midnight'
!       else
!         ampm = 'AM'
!       ENDIF
!     ENDIF
!   ENDIF
! 
!   write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
!     d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
! 
!   RETURN
! end

SUBROUTINE tinvit ( n, d, e, e2, m, w, ind, z, ierr )

!*****************************************************************************80
!
!! TINVIT computes eigenvectors from eigenvalues, real tridiagonal symmetric.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvectors of a tridiagonal
!    symmetric matrix corresponding to specified eigenvalues,
!    using inverse iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    B Smith, J Boyle, J Dongarra, B Garbow, Y Ikebe, V Klema, C Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the matrix.
!
!    Input, REAL ( kind = 8 ) E(N), contains the subdiagonal elements of
!    the input matrix in E(2:N).  E(1) is arbitrary.
!
!    Input, REAL ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E, with zeros corresponding to negligible elements of E.
!    E(I) is considered negligible if it is not larger than the product of
!    the relative machine precision and the sum of the magnitudes of D(I)
!    and D(I-1).  E2(1) must contain 0.0D+00 if the eigenvalues are in
!    ascending order, or 2.0D+00 if the eigenvalues are in descending order.
!    If BISECT, TRIDIB, or IMTQLV has been used to find the eigenvalues,
!    their output E2 array is exactly what is expected here.
!
!    Input, INTEGER ( kind = 4 ) M, the number of specified eigenvalues.
!
!    Input, REAL ( kind = 8 ) W(M), the eigenvalues.
!
!    Input, INTEGER ( kind = 4 ) IND(M), the submatrix indices associated with 
!    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
!    first submatrix from the top, 2 for those belonging to the second
!    submatrix, and so on.
!
!    Output, REAL ( kind = 8 ) Z(N,M), the associated set of orthonormal
!    eigenvectors.  Any vector which fails to converge is set to zero.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
!      converge in 5 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) eps2
  REAL ( kind = 8 ) eps3
  REAL ( kind = 8 ) eps4
  INTEGER ( kind = 4 ) group
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ind(m)
  INTEGER ( kind = 4 ) ip
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  REAL ( kind = 8 ) norm
  REAL ( kind = 8 ) order
  INTEGER ( kind = 4 ) p
  !REAL ( kind = 8 ) pythag
  INTEGER ( kind = 4 ) q
  INTEGER ( kind = 4 ) r
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) rv2(n)
  REAL ( kind = 8 ) rv3(n)
  REAL ( kind = 8 ) rv4(n)
  REAL ( kind = 8 ) rv6(n)
  INTEGER ( kind = 4 ) s
  INTEGER ( kind = 4 ) tag
  REAL ( kind = 8 ) u
  REAL ( kind = 8 ) uk
  REAL ( kind = 8 ) v
  REAL ( kind = 8 ) w(m)
  REAL ( kind = 8 ) x0
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) xu
  REAL ( kind = 8 ) z(n,m)

  ierr = 0

  IF( m == 0 ) then
    RETURN
  ENDIF

  u = 0.0D+00
  x0 = 0.0D+00

  tag = 0
  order = 1.0D+00 - e2(1)
  q = 0
!
!  Establish and process next submatrix.
!
100 CONTINUE

  p = q + 1

  DO q = p, n
    IF( q == n ) then
      exit
    ENDIF
    IF( e2(q+1) == 0.0D+00 ) then
      exit
    ENDIF
  ENDDO
!
!  Find vectors by inverse iteration.
!
!140 CONTINUE

  tag = tag + 1
  s = 0

  DO r = 1, m

     IF( ind(r) /= tag ) then
       go to 920
     ENDIF

     its = 1
     x1 = w(r)

     IF( s /= 0 ) then
       go to 510
     ENDIF
!
!  Check for isolated root.
!
     xu = 1.0D+00

     IF( p == q ) then
       rv6(p) = 1.0D+00
       go to 870
     ENDIF

     norm = abs ( d(p) )
     ip = p + 1

     DO i = p + 1, q
       norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
     ENDDO
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by EPS3,
!  EPS4 is taken very small to avoid overflow.
!
     eps2 = 0.001D+00 * norm
     eps3 = abs ( norm ) * epsilon ( eps3 )
     uk = q - p + 1
     eps4 = uk * eps3
     uk = eps4 / sqrt ( uk )
     s = p

505 CONTINUE

     group = 0
     go to 520
!
!  Look for close or coincident roots.
!
510  CONTINUE

     IF( abs ( x1 - x0 ) >= eps2 ) then
       go to 505
     ENDIF

     group = group + 1

     IF( order * (x1 - x0) <= 0.0D+00 ) then
       x1 = x0 + order * eps3
     ENDIF
!
!  Elimination with interchanges and initialization of vector.
!
520  CONTINUE

     v = 0.0D+00

     DO i = p, q

        rv6(i) = uk

        IF( i == p ) then
          go to 560
        ENDIF

        IF( abs ( e(i) ) < abs ( u ) ) then
          go to 540
        ENDIF

        xu = u / e(i)
        rv4(i) = xu
        rv1(i-1) = e(i)
        rv2(i-1) = d(i) - x1
        rv3(i-1) = 0.0D+00
        IF( i /= q ) then
          rv3(i-1) = e(i+1)
        ENDIF
        u = v - xu * rv2(i-1)
        v = - xu * rv3(i-1)
        cycle

540     CONTINUE

        xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = 0.0D+00

560     CONTINUE

        u = d(i) - x1 - xu * v
        IF( i /= q ) then
          v = e(i+1)
        ENDIF

     ENDDO

     IF( u == 0.0D+00 ) then
       u = eps3
     ENDIF

     rv1(q) = u
     rv2(q) = 0.0D+00
     rv3(q) = 0.0D+00
!
!  Back substitution.
!
600   CONTINUE

  DO ii = p, q
    i = p + q - ii
    rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
    v = u
    u = rv6(i)
  ENDDO
!
!  Orthogonalize with respect to previous members of group.
!
     j = r

     DO jj = 1, group

       do

         j = j - 1

         IF( ind(j) == tag ) then
           exit
         ENDIF

       ENDDO

       xu = dot_product ( rv6(p:q), z(p:q,j) )

       rv6(p:q) = rv6(p:q) - xu * z(p:q,j)

     ENDDO

     norm = sum ( abs ( rv6(p:q) ) )

     IF( norm >= 1.0D+00 ) then
       go to 840
     ENDIF
!
!  Forward substitution.
!
     IF( its == 5 ) then
       go to 830
     ENDIF

     IF( norm == 0.0D+00 ) then
       rv6(s) = eps4
       s = s + 1
       IF( q < s ) then
         s = p
       ENDIF
       go to 780
     ENDIF

!740  CONTINUE

     xu = eps4 / norm
     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate.
!
780  CONTINUE
!
!  If RV1(I-1) == E(I), a row interchange was performed earlier in the
!  triangularization process.
!
     DO i = ip, q

       u = rv6(i)

       IF( rv1(i-1) == e(i) ) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       ENDIF

       rv6(i) = u - rv4(i) * rv6(i-1)

     ENDDO

     its = its + 1
     go to 600
!
!  Set error: non-converged eigenvector.
!
830  CONTINUE

     ierr = -r
     xu = 0.0D+00
     go to 870
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
840  CONTINUE

     u = 0.0D+00
     DO i = p, q
       u = pythag ( u, rv6(i) )
     ENDDO

     xu = 1.0D+00 / u

870  CONTINUE

     z(1:n,r) = 0.0D+00
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1

920  CONTINUE

  ENDDO

  IF( q < n ) then
    go to 100
  ENDIF

  RETURN

END SUBROUTINE

SUBROUTINE tql1 ( n, d, e, ierr )

!*****************************************************************************80
!
!! TQL1 computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a symmetric tridiagonal
!    matrix by the QL method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  References:
!
!    Bowdler, Martin, Reinsch, James Wilkinson,
!    Numerische Mathematik,
!    Volume 11, 1968, pages 293-306.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, is the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) D(N).
!    On input, the diagonal elements of the matrix.
!    On output, the eigenvalues in ascending order.
!    If an error exit is made, the eigenvalues are correct and
!    ordered for indices 1, 2,... IERR-1, but may not be
!    the smallest eigenvalues.
!
!    Input/output, REAL ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, normal RETURN,
!    J, if the J-th eigenvalue has not been determined after
!    30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) c2
  REAL ( kind = 8 ) c3
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) dl1
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) el1
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) l2
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) s2
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2

  ierr = 0
  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    e(i-1) = e(i)
  ENDDO

  f = 0.0D+00
  tst1 = 0.0D+00
  e(n) = 0.0D+00

  DO l = 1, n

    j = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    DO m = l, n

      tst2 = tst1 + abs ( e(m) )

      IF( tst2 == tst1 ) then
        exit
      ENDIF

    ENDDO

    IF( m == l ) then
      go to 210
    ENDIF

130 CONTINUE

    IF( j >= 30 ) then
      ierr = l
      RETURN
    ENDIF

    j = j + 1
!
!  Form the shift.
!
    l1 = l + 1
    l2 = l1 + 1
    g = d(l)
    p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
    r = pythag ( p, 1.0D+00 )
    d(l) = e(l) / ( p + sign ( r, p ) )
    d(l1) = e(l) * ( p + sign ( r, p ) )
    dl1 = d(l1)
    h = g - d(l)

    d(l2:n) = d(l2:n) - h

    f = f + h
!
!  QL transformation.
!
    p = d(m)
    c = 1.0D+00
    c2 = c
    el1 = e(l1)
    s = 0.0D+00
    mml = m - l

    DO ii = 1, mml
      c3 = c2
      c2 = c
      s2 = s
      i = m - ii
      g = c * e(i)
      h = c * p
      r = pythag ( p, e(i) )
      e(i+1) = s * r
      s = e(i) / r
      c = p / r
      p = c * d(i) - s * g
      d(i+1) = h + s * ( c * g + s * d(i) )
    ENDDO

    p = - s * s2 * c3 * el1 * e(l) / dl1
    e(l) = s * p
    d(l) = c * p
    tst2 = tst1 + abs ( e(l) )
    IF( tst2 > tst1 ) then
      go to 130
    ENDIF

210 CONTINUE

    p = d(l) + f
!
!  Order the eigenvalues.
!
    DO ii = 2, l
      i = l + 2 - ii
      IF( p >= d(i-1) ) then
        go to 270
      ENDIF
      d(i) = d(i-1)
    ENDDO

    i = 1

270 CONTINUE

    d(i) = p

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE tql2 ( n, d, e, z, ierr )

!*****************************************************************************80
!
!! TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues and eigenvectors of a symmetric
!    tridiagonal matrix by the QL method.  The eigenvectors of a full
!    symmetric matrix can also be found if TRED2 has been used to reduce this
!    full matrix to tridiagonal form.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Bowdler, Martin, Reinsch, James Wilkinson,
!    TQL2,
!    Numerische Mathematik,
!    Volume 11, pages 293-306, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) D(N).  On input, the diagonal elements of
!    the matrix.  On output, the eigenvalues in ascending order.  If an error
!    exit is made, the eigenvalues are correct but unordered for indices
!    1,2,...,IERR-1.
!
!    Input/output, REAL ( kind = 8 ) E(N).  On input, E(2:N) contains the
!    subdiagonal elements of the input matrix, and E(1) is arbitrary.
!    On output, E has been destroyed.
!
!    Input, REAL ( kind = 8 ) Z(N,N).  On input, the transformation matrix
!    produced in the reduction by TRED2, if performed.  If the eigenvectors of
!    the tridiagonal matrix are desired, Z must contain the identity matrix.
!    On output, Z contains the orthonormal eigenvectors of the symmetric
!    tridiagonal (or full) matrix.  If an error exit is made, Z contains
!    the eigenvectors associated with the stored eigenvalues.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, normal RETURN,
!    J, if the J-th eigenvalue has not been determined after
!    30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) c2
  REAL ( kind = 8 ) c3
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) dl1
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) el1
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) l2
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) s2
  REAL ( kind = 8 ) t
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) z(n,n)

  ierr = 0

  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    e(i-1) = e(i)
  ENDDO

  f = 0.0D+00
  tst1 = 0.0D+00
  e(n) = 0.0D+00

  DO l = 1, n

    j = 0
    h = abs ( d(l) ) + abs ( e(l) )
    tst1 = max ( tst1, h )
!
!  Look for a small sub-diagonal element.
!
    DO m = l, n
      tst2 = tst1 + abs ( e(m) )
      IF( tst2 == tst1 ) then
        exit
      ENDIF
    ENDDO

    IF( m /= l ) then

      do

        IF( 30 <= j ) then
          ierr = l
          RETURN
        ENDIF

        j = j + 1
!
!  Form shift.
!
        l1 = l + 1
        l2 = l1 + 1
        g = d(l)
        p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
        r = pythag ( p, 1.0D+00 )
        d(l) = e(l) / ( p + sign ( r, p ) )
        d(l1) = e(l) * ( p + sign ( r, p ) )
        dl1 = d(l1)
        h = g - d(l)
        d(l2:n) = d(l2:n) - h
        f = f + h
!
!  QL transformation.
!
        p = d(m)
        c = 1.0D+00
        c2 = c
        el1 = e(l1)
        s = 0.0D+00
        mml = m - l

        DO ii = 1, mml

          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
!
!  Form vector.
!
          DO k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * h
            z(k,i) = c * z(k,i) - s * h
          ENDDO

        ENDDO

        p = - s * s2 * c3 * el1 * e(l) / dl1
        e(l) = s * p
        d(l) = c * p
        tst2 = tst1 + abs ( e(l) )

        IF( tst2 <= tst1 ) then
          exit
        ENDIF

      ENDDO

    ENDIF

    d(l) = d(l) + f

  ENDDO
!
!  Order eigenvalues and eigenvectors.
!
  DO ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    DO j = ii, n

      IF( d(j) < p ) then
        k = j
        p = d(j)
      ENDIF

    ENDDO

    IF( k /= i ) then

      d(k) = d(i)
      d(i) = p

      DO j = 1, n
        t      = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = t
      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE tqlrat ( n, d, e2, ierr )

!*****************************************************************************80
!
!! TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds the eigenvalues of a symmetric
!    tridiagonal matrix by the rational QL method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    C Reinsch,
!    Algorithm 464, TQLRAT,
!    Communications of the ACM,
!    Volume 16, page 689, 1973.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) D(N).  On input, D contains the diagonal
!    elements of the matrix.  On output, D contains the eigenvalues in ascending
!    order.  If an error exit was made, then the eigenvalues are correct
!    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
!
!    Input/output, REAL ( kind = 8 ) E2(N), contains in positions 2 through N 
!    the squares of the subdiagonal elements of the matrix.  E2(1) is
!    arbitrary.  On output, E2 has been overwritten by workspace
!    information.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for no error,
!    J, if the J-th eigenvalue could not be determined after 30 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) b
  REAL ( kind = 8 ) c
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) l1
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) mml
  REAL ( kind = 8 ) p
  !REAL ( kind = 8 ) pythag
  REAL ( kind = 8 ) r
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) t

  ierr = 0

  IF( n == 1 ) then
    RETURN
  ENDIF

  DO i = 2, n
    e2(i-1) = e2(i)
  ENDDO

  f = 0.0D+00
  t = 0.0D+00
  e2(n) = 0.0D+00

  DO l = 1, n

     j = 0
     h = abs ( d(l) ) + sqrt ( e2(l) )

     IF( t <= h ) then

       t = h
       b = abs ( t ) * epsilon ( b )
       c = b * b

     ENDIF
!
!  Look for small squared sub-diagonal element.
!
     DO m = l, n
       IF( e2(m) <= c ) then
         exit
       ENDIF
     ENDDO

     IF( m /= l ) then

       do

         IF( 30 <= j ) then
           ierr = l
           RETURN
         ENDIF

         j = j + 1
!
!  Form shift.
!
         l1 = l + 1
         s = sqrt ( e2(l) )
         g = d(l)
         p = ( d(l1) - g ) / ( 2.0D+00 * s )
         r = pythag ( p, 1.0D+00 )
         d(l) = s / ( p + sign ( r, p ) )
         h = g - d(l)
         d(l1:n) = d(l1:n) - h
         f = f + h
!
!  Rational QL transformation.
!
         g = d(m)
         IF( g == 0.0D+00 ) then
           g = b
         ENDIF

         h = g
         s = 0.0D+00
         mml = m - l

         DO ii = 1, mml
           i = m - ii
           p = g * h
           r = p + e2(i)
           e2(i+1) = s * r
           s = e2(i) / r
           d(i+1) = h + s * ( h + d(i) )
           g = d(i) - e2(i) / g
           IF( g == 0.0D+00 ) then
             g = b
           ENDIF
           h = g * p / r
         ENDDO

         e2(l) = s * g
         d(l) = h
!
!  Guard against underflow in convergence test.
!
         IF( h == 0.0D+00 ) then
           exit
         ENDIF

         IF( abs ( e2(l) ) <= abs ( c / h ) ) then
           exit
         ENDIF

         e2(l) = h * e2(l)

         IF( e2(l) == 0.0D+00 ) then
            exit
         ENDIF

       ENDDO

     ENDIF

     p = d(l) + f
!
!  Order the eigenvalues.
!
     DO i = l, 1, -1
       IF( i == 1 ) then
         d(i) = p
         exit
       else IF( d(i-1) <= p ) then
         d(i) = p
         exit
       ENDIF
       d(i) = d(i-1)
     ENDDO

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE trbak1 ( n, a, e, m, z )

!*****************************************************************************80
!
!! TRBAK1 determines eigenvectors by undoing the TRED1 transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a real symmetric
!    matrix by back transforming those of the corresponding
!    symmetric tridiagonal matrix determined by TRED1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), contains information about the orthogonal
!    transformations used in the reduction by TRED1 in its strict lower
!    triangle.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is arbitrary.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the eigenvectors to be
!    back transformed.  On output, the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) e(n)
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) j
  !INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) z(n,m)

  IF( m <= 0 ) then
    RETURN
  ENDIF

  IF( n <= 1 ) then
    RETURN
  ENDIF

  DO i = 2, n

    l = i - 1

    IF( e(i) /= 0.0D+00 ) then

      DO j = 1, m

        s = dot_product ( a(i,1:l), z(1:l,j) )

        s = ( s / a(i,l) ) / e(i)

        z(1:l,j) = z(1:l,j) + s * a(i,1:l)

      ENDDO

    ENDIF

  ENDDO

  CONTINUE

  RETURN

END SUBROUTINE

SUBROUTINE trbak3 ( n, nv, a, m, z )

!*****************************************************************************80
!
!! TRBAK3 determines eigenvectors by undoing the TRED3 transformation.
!
!  Discussion:
!
!    This SUBROUTINE forms the eigenvectors of a real symmetric
!    matrix by back transforming those of the corresponding
!    symmetric tridiagonal matrix determined by TRED3.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) NV, the dimension of the array paramater A,
!    which must be at least N*(N+1)/2.
!
!    Input, REAL ( kind = 8 ) A(NV), information about the orthogonal
!    transformations used in the reduction by TRED3.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvectors to be back
!    transformed.
!
!    Input/output, REAL ( kind = 8 ) Z(N,M).  On input, the eigenvectors to be 
!    back transformed.  On output, the transformed eigenvectors.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) nv

  REAL ( kind = 8 ) a(nv)
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ik
  INTEGER ( kind = 4 ) iz
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  INTEGER ( kind = 4 ) n
  REAL ( kind = 8 ) s
  REAL ( kind = 8 ) z(n,m)

  IF( m == 0 ) then
    RETURN
  ENDIF

  DO i = 2, n

    l = i - 1
    iz = ( i * l ) / 2
    ik = iz + i
    h = a(ik)

    IF( h /= 0.0D+00 ) then

      DO j = 1, m

        s = 0.0D+00
        ik = iz

        DO k = 1, l
          ik = ik + 1
          s = s + a(ik) * z(k,j)
        ENDDO

        s = ( s / h ) / h
        ik = iz

        DO k = 1, l
          ik = ik + 1
          z(k,j) = z(k,j) - s * a(ik)
        ENDDO

      ENDDO

    ENDIF

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE tred1 ( n, a, d, e, e2 )

!*****************************************************************************80
!
!! TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    The routine reduces a real symmetric matrix to a symmetric
!    tridiagonal matrix using orthogonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED1,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix A.
!
!    Input/output, REAL ( kind = 8 ) A(N,N), on input, contains the real
!    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
!    On output, A contains information about the orthogonal transformations
!    used in the reduction in its strict lower triangle.
!    The full upper triangle of A is unaltered.
!
!    Output, REAL ( kind = 8 ) D(N), contains the diagonal elements of the
!    tridiagonal matrix.
!
!    Output, REAL ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) E2(N), contains the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) scale

  d(1:n) = a(n,1:n)

  DO i = 1, n
    a(n,i) = a(i,i)
  ENDDO

  DO ii = 1, n

    i = n + 1 - ii
    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    IF( scale == 0.0D+00 ) then

      DO j = 1, l
        d(j) = a(l,j)
        a(l,j) = a(i,j)
        a(i,j) = 0.0D+00
      ENDDO

      e(i) = 0.0D+00
      e2(i) = 0.0D+00

      cycle

    ENDIF

    d(1:l) = d(1:l) / scale

    DO k = 1, l
      h = h + d(k)**2
    ENDDO

    e2(i) = h * scale**2
    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g

    IF( 1 <= l ) then
!
!  Form A * U.
!
      e(1:l) = 0.0D+00

      DO j = 1, l

        f = d(j)
        g = e(j) + a(j,j) * f

        DO k = j + 1, l
          g = g + a(k,j) * d(k)
          e(k) = e(k) + a(k,j) * f
        ENDDO

        e(j) = g

      ENDDO
!
!  Form P.
!
      f = 0.0D+00

      DO j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
      ENDDO

      h = f / ( h + h )
!
!  Form Q.
!
      e(1:l) = e(1:l) - h * d(1:l)
!
!  Form reduced A.
!
      DO j = 1, l

        f = d(j)
        g = e(j)

        a(j:l,j) = a(j:l,j) - f * e(j:l) - g * d(j:l)

      ENDDO

    ENDIF

    DO j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
    ENDDO


  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE tred2 ( n, a, d, e, z )

!*****************************************************************************80
!
!! TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    This SUBROUTINE reduces a real symmetric matrix to a
!    symmetric tridiagonal matrix using and accumulating
!    orthogonal similarity transformations.
!
!    A and Z may coincide, in which case a single storage area is used
!    for the input of A and the output of Z.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2012
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED2,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, REAL ( kind = 8 ) A(N,N), the real symmetric input matrix.  Only the
!    lower triangle of the matrix need be supplied.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, REAL ( kind = 8 ) E(N), contains the subdiagonal elements of the
!    tridiagonal matrix in E(2:N).  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) Z(N,N), the orthogonal transformation matrix
!    produced in the reduction.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) a(n,n)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  REAL ( kind = 8 ) hh
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) scale
  REAL ( kind = 8 ) z(n,n)

  DO i = 1, n
    z(i:n,i) = a(i:n,i)
  ENDDO

  d(1:n) = a(n,1:n)

  DO ii = 2, n

    i = n + 2 - ii
    l = i - 1
    h = 0.0D+00
!
!  Scale row.
!
    scale = sum ( abs ( d(1:l) ) )

    IF( scale == 0.0D+00 ) then

      e(i) = d(l)

      DO j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0D+00
        z(j,i) = 0.0D+00
      ENDDO

      d(i) = 0.0D+00

      cycle

    ENDIF

    d(1:l) = d(1:l) / scale

    h = h + dot_product ( d(1:l), d(1:l) )

    f = d(l)
    g = - sign ( sqrt ( h ), f )
    e(i) = scale * g
    h = h - f * g
    d(l) = f - g
!
!  Form A*U.
!
    e(1:l) = 0.0D+00

    DO j = 1, l

      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f

      DO k = j + 1, l
        g = g + z(k,j) * d(k)
        e(k) = e(k) + z(k,j) * f
      ENDDO

      e(j) = g

    ENDDO
!
!  Form P.
!
    e(1:l) = e(1:l) / h

    f = dot_product ( e(1:l), d(1:l) )

    hh = 0.5D+00 * f / h
!
!  Form Q.
!
    e(1:l) = e(1:l) - hh * d(1:l)
!
!  Form reduced A.
!
    DO j = 1, l

      f = d(j)
      g = e(j)

      z(j:l,j) = z(j:l,j) - f * e(j:l) - g * d(j:l)

      d(j) = z(l,j)
      z(i,j) = 0.0D+00

    ENDDO

    d(i) = h

  ENDDO
!
!  Accumulation of transformation matrices.
!
  DO i = 2, n

    l = i - 1
    z(n,l) = z(l,l)
    z(l,l) = 1.0D+00
    h = d(i)

    IF( h /= 0.0D+00 ) then

      d(1:l) = z(1:l,i) / h

      DO j = 1, l

        g = dot_product ( z(1:l,i), z(1:l,j) )

        DO k = 1, l
          z(k,j) = z(k,j) - g * d(k)
        ENDDO

      ENDDO

    ENDIF

    z(1:l,i) = 0.0D+00

  ENDDO

  d(1:n) = z(n,1:n)

  z(n,1:n-1) = 0.0D+00
  z(n,n) = 1.0D+00

  e(1) = 0.0D+00

  RETURN

END SUBROUTINE

SUBROUTINE tred3 ( n, nv, a, d, e, e2 )

!*****************************************************************************80
!
!! TRED3: transform real symmetric packed matrix to symmetric tridiagonal form.
!
!  Discussion:
!
!    This SUBROUTINE reduces a real symmetric matrix, stored as
!    a one-dimensional array, to a symmetric tridiagonal matrix
!    using orthogonal similarity transformations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Martin, Reinsch, James Wilkinson,
!    TRED3,
!    Numerische Mathematik,
!    Volume 11, pages 181-195, 1968.
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input, INTEGER ( kind = 4 ) NV, the dimension of A, which must be at least
!    (N*(N+1))/2.
!
!    Input/output, REAL ( kind = 8 ) A(NV).  On input, the lower triangle of
!    the real symmetric matrix, stored row-wise.  On output, information about
!    the orthogonal transformations used in the reduction.
!
!    Output, REAL ( kind = 8 ) D(N), the diagonal elements of the tridiagonal
!    matrix.
!
!    Output, REAL ( kind = 8 ) E(N), the subdiagonal elements of the tridiagonal
!    matrix in E(2:N).  E(1) is set to zero.
!
!    Output, REAL ( kind = 8 ) E2(N),  the squares of the corresponding
!    elements of E.  E2 may coincide with E if the squares are not needed.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) n
  INTEGER ( kind = 4 ) nv

  REAL ( kind = 8 ) a(nv)
  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) f
  REAL ( kind = 8 ) g
  REAL ( kind = 8 ) h
  REAL ( kind = 8 ) hh
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) iz
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jk
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) scale

  DO ii = 1, n

     i = n + 1 - ii
     l = i - 1
     iz = ( i * l ) / 2
     h = 0.0D+00
     scale = 0.0D+00
!
!  Scale row.
!
     DO k = 1, l
       iz = iz + 1
       d(k) = a(iz)
       scale = scale + abs ( d(k) )
     ENDDO

     IF( scale == 0.0D+00 ) then
       e(i) = 0.0D+00
       e2(i) = 0.0D+00
       go to 290
     ENDIF

     DO k = 1, l
       d(k) = d(k) / scale
       h = h + d(k)**2
     ENDDO

     e2(i) = scale * scale * h
     f = d(l)
     g = - sign ( sqrt ( h ), f )
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
     a(iz) = scale * d(l)

     IF( l == 1 ) then
       go to 290
     ENDIF

     jk = 1

     DO j = 1, l

        f = d(j)
        g = 0.0D+00

        DO k = 1, j - 1
          g = g + a(jk) * d(k)
          e(k) = e(k) + a(jk) * f
          jk = jk + 1
        ENDDO

        e(j) = g + a(jk) * f
        jk = jk + 1

     ENDDO
!
!  Form P.
!
     e(1:l) = e(1:l) / h
     f = dot_product ( e(1:l), d(1:l) )
     hh = f / ( h + h )
!
!  Form Q.
!
     e(1:l) = e(1:l) - hh * d(1:l)
     jk = 1
!
!  Form reduced A.
!
     DO j = 1, l
       f = d(j)
       g = e(j)
       DO k = 1, j
         a(jk) = a(jk) - f * e(k) - g * d(k)
         jk = jk + 1
       ENDDO
     ENDDO

290  CONTINUE

     d(i) = a(iz+1)
     a(iz+1) = scale * sqrt ( h )

!300  CONTINUE

  ENDDO

  RETURN

END SUBROUTINE

SUBROUTINE tridib ( n, eps1, d, e, e2, lb, ub, m11, m, w, ind, ierr )

!*****************************************************************************80
!
!! TRIDIB computes some eigenvalues of a real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvalues of a tridiagonal
!    symmetric matrix between specified boundary indices,
!    using bisection.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 April 2011
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) EPS1.  On input, an absolute error
!    tolerance for the computed eigenvalues.  It should be chosen commensurate
!    with relative perturbations in the matrix elements of the order of the
!    relative machine precision.  If the input EPS1 is non-positive, it
!    is reset for each submatrix to a default value, namely, minus the
!    product of the relative machine precision and the 1-norm of the submatrix.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the input matrix.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E2(N).  On input, the squares of the
!    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
!    E2 corresponding to elements of E regarded as negligible, have been
!    replaced by zero, causing the matrix to split into a direct sum of
!    submatrices.  E2(1) is also set to zero.
!
!    Input, INTEGER ( kind = 4 ) M11, the lower boundary index for the desired
!    eigenvalues.
!
!    Input, INTEGER ( kind = 4 ) M, the number of eigenvalues desired.  The
!    upper boundary index M22 is then obtained as M22 = M11 + M - 1.
!
!    Output, REAL ( kind = 8 ) LB, UB, define an interval containing exactly
!    the desired eigenvalues.
!
!    Output, REAL ( kind = 8 ) W(M), the eigenvalues between indices M11 and M22
!    in ascending order.
!
!    Output, INTEGER ( kind = 4 ) IND(M), the submatrix indices associated with
!    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
!    first submatrix from the top, 2 for those belonging to the second
!    submatrix, and so on.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, for normal RETURN,
!    3*N+1, if multiple eigenvalues at index M11 make unique selection
!      impossible,
!    3*N+2, if multiple eigenvalues at index M22 make unique selection
!      impossible.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) eps1
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ind(m)
  INTEGER ( kind = 4 ) isturm
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) k
  INTEGER ( kind = 4 ) l
  REAL ( kind = 8 ) lb
  INTEGER ( kind = 4 ) m1
  INTEGER ( kind = 4 ) m11
  INTEGER ( kind = 4 ) m2
  INTEGER ( kind = 4 ) m22
  INTEGER ( kind = 4 ) p
  INTEGER ( kind = 4 ) q
  INTEGER ( kind = 4 ) r
  REAL ( kind = 8 ) rv4(n)
  REAL ( kind = 8 ) rv5(n)
  INTEGER ( kind = 4 ) s
  REAL ( kind = 8 ) t1
  REAL ( kind = 8 ) t2
  INTEGER ( kind = 4 ) tag
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) u
  REAL ( kind = 8 ) ub
  REAL ( kind = 8 ) v
  REAL ( kind = 8 ) w(m)
  REAL ( kind = 8 ) x0
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) xu

  ierr = 0
  tag = 0
  xu = d(1)
  x0 = d(1)
  s = 0
  u = 0.0D+00
!
!  Look for small sub-diagonal entries and determine an
!  interval containing all the eigenvalues.
!
  DO i = 1, n

     x1 = u

     IF( i == n ) then
       u = 0.0D+00
     else
       u = abs ( e(i+1) )
     ENDIF

     xu = min ( xu, d(i) - ( x1 + u ) )
     x0 = max ( x0, d(i) + ( x1 + u ) )

     IF( 1 < i ) then
       tst1 = abs ( d(i) ) + abs ( d(i-1) )
       tst2 = tst1 + abs ( e(i) )
       IF( tst2 <= tst1 ) then
         e2(i) = 0.0D+00
       ENDIF
     else
       e2(i) = 0.0D+00
     ENDIF

  ENDDO

  x1 = n
  x1 = x1 * max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )
  xu = xu - x1
  t1 = xu
  x0 = x0 + x1
  t2 = x0
!
!  Determine an interval containing exactly the desired eigenvalues.
!
  p = 1
  q = n
  m1 = m11 - 1

  IF( m1 == 0 ) then
    go to 75
  ENDIF

  isturm = 1

50 CONTINUE

  v = x1
  x1 = xu + ( x0 - xu ) * 0.5D+00

  IF( x1 == v ) then
    go to 980
  ENDIF

  go to 320

60 CONTINUE

  IF( s < m1 ) then
    xu = x1
    go to 50
  else IF( m1 < s ) then
    x0 = x1
    go to 50
  ENDIF

  xu = x1
  t1 = x1

75 CONTINUE

  m22 = m1 + m

  IF( m22 /= n ) then
    x0 = t2
    isturm = 2
    go to 50
  ENDIF

  go to 90

80 CONTINUE

  IF( s < m22 ) then
    xu = x1
    go to 50
  else IF( m22 < s ) then
    x0 = x1
    go to 50
  ENDIF

   t2 = x1

90 CONTINUE

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the
!  Gerschgorin bounds.
!
100 CONTINUE

  IF( r == m ) then
    go to 1001
  ENDIF

  tag = tag + 1
  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0D+00

  DO q = p, n

    x1 = u
    u = 0.0D+00
    v = 0.0D+00

    IF( q < n ) then
      u = abs ( e(q+1) )
      v = e2(q+1)
    ENDIF

    xu = min ( d(q) - ( x1 + u ), xu )
    x0 = max ( d(q) + ( x1 + u ), x0 )

    IF( v == 0.0D+00 ) then
      exit
    ENDIF

  ENDDO

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )

  IF( eps1 <= 0.0D+00 ) then
    eps1 = -x1
  ENDIF

  IF( p /= q ) then
    go to 180
  ENDIF
!
!  Check for isolated root within interval.
!
  IF( d(p) < t1 .or. t2 <= d(p) ) then
    go to 940
  ENDIF

  m1 = p
  m2 = p
  rv5(p) = d(p)
  go to 900

180 CONTINUE

  x1 = x1 * ( q - p + 1 )
  lb = max ( t1, xu - x1 )
  ub = min ( t2, x0 + x1 )
  x1 = lb
  isturm = 3
  go to 320

200 CONTINUE

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

220 CONTINUE

  m2 = s
  IF( m2 < m1 ) then
    go to 940
  ENDIF
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for the K-th eigenvalue.
!
  k = m2

250 CONTINUE

  xu = lb

  DO ii = m1, k

    i = m1 + k - ii
    IF( xu < rv4(i) ) then
      xu = rv4(i)
      exit
    ENDIF

  ENDDO

  x0 = min ( x0, rv5(k) )
!
!  Next bisection step.
!
300  CONTINUE

     x1 = ( xu + x0 ) * 0.5D+00

     IF( ( x0 - xu ) <= abs ( eps1) ) then
       go to 420
     ENDIF

     tst1 = 2.0D+00 * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + (x0 - xu)

     IF( tst2 == tst1 ) then
       go to 420
     ENDIF
!
!  Sturm sequence.
!
320  CONTINUE

     s = p - 1
     u = 1.0D+00

     DO i = p, q

       IF( u == 0.0D+00 ) then
         v = abs ( e(i) ) / epsilon ( v )
         IF( e2(i) == 0.0D+00 ) then
           v = 0.0D+00
         ENDIF
       else
         v = e2(i) / u
       ENDIF

       u = d(i) - x1 - v

       IF( u < 0.0D+00 ) then
         s = s + 1
       ENDIF

     ENDDO

     go to (60,80,200,220,360), isturm
!
!  Refine intervals.
!
360  CONTINUE

     IF( k <= s ) then
       go to 400
     ENDIF

     xu = x1

     IF( m1 <= s ) then
       go to 380
     ENDIF

     rv4(m1) = x1
     go to 300

380  CONTINUE

     rv4(s+1) = x1
     rv5(s) = min ( rv5(s), x1 )
     go to 300

400  CONTINUE

     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420  CONTINUE

  rv5(k) = x1
  k = k - 1
  IF( m1 <= k ) then
    go to 250
  ENDIF
!
!  Order eigenvalues tagged with their submatrix associations.
!
900 CONTINUE

  s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1

  DO l = 1, r

     IF( s < j ) then
       go to 910
     ENDIF

     IF( m2 < k ) then
       go to 940
     ENDIF

     IF( w(l) <= rv5(k) ) then
       go to 915
     ENDIF

     DO ii = j, s
       i = l + s - ii
       w(i+1) = w(i)
       ind(i+1) = ind(i)
     ENDDO

910  CONTINUE

     w(l) = rv5(k)
     ind(l) = tag
     k = k + 1
     go to 920

915  CONTINUE

     j = j + 1

920  CONTINUE

  ENDDO

940 CONTINUE

  IF( q < n ) then
    go to 100
  ENDIF

  go to 1001
!
!  Set error: interval cannot be found containing exactly the
!  desired eigenvalues.
!
980 CONTINUE

  ierr = 3 * n + isturm

1001 CONTINUE

  lb = t1
  ub = t2

  RETURN

END SUBROUTINE

SUBROUTINE tsturm ( n, eps1, d, e, e2, lb, ub, mm, m, w, z, ierr )

!*****************************************************************************80
!
!! TSTURM computes some eigenvalues/vectors, real symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This SUBROUTINE finds those eigenvalues of a tridiagonal
!    symmetric matrix which lie in a specified interval and their
!    associated eigenvectors, using bisection and inverse iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 October 2009
!
!  Author:
!
!    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
!    Klema, Moler.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    James Wilkinson, Christian Reinsch,
!    Handbook for Automatic Computation,
!    Volume II, Linear Algebra, Part 2,
!    Springer, 1971,
!    ISBN: 0387054146,
!    LC: QA251.W67.
!
!    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
!    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
!    Matrix Eigensystem Routines, EISPACK Guide,
!    Lecture Notes in Computer Science, Volume 6,
!    Springer Verlag, 1976,
!    ISBN13: 978-3540075462,
!    LC: QA193.M37.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, REAL ( kind = 8 ) EPS1.  On input, an absolute error
!    tolerance for the computed eigenvalues.  It should be chosen commensurate
!    with relative perturbations in the matrix elements of the order of the
!    relative machine precision.  If the input EPS1 is non-positive, it
!    is reset for each submatrix to a default value, namely, minus the
!    product of the relative machine precision and the 1-norm of the submatrix.
!
!    Input, REAL ( kind = 8 ) D(N), the diagonal elements of the input matrix.
!
!    Input, REAL ( kind = 8 ) E(N), the subdiagonal elements of the input matrix
!    in E(2:N).  E(1) is arbitrary.
!
!    Input/output, REAL ( kind = 8 ) E2(N).  On input, the squares of the
!    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
!    E2 corresponding to elements of E regarded as negligible have been
!    replaced by zero, causing the matrix to split into a direct sum of
!    submatrices.  E2(1) is also set to zero.
!
!    Input, REAL ( kind = 8 ) LB, UB, define the interval to be searched for
!    eigenvalues.  If LB is not less than UB, no eigenvalues will be found.
!
!    Input, INTEGER ( kind = 4 ) MM, an upper bound for the number of
!    eigenvalues in the interval.  If more than MM eigenvalues are determined
!    to lie in the interval, an error RETURN is made with no values or vectors
!    found.
!
!    Output, INTEGER ( kind = 4 ) M, the number of eigenvalues determined to lie
!    in (LB, UB).
!
!    Output, REAL ( kind = 8 ) W(M), the eigenvalues in ascending order if the
!    matrix does not split.  If the matrix splits, the eigenvalues are in
!    ascending order for each submatrix.  If a vector error exit is made, W
!    contains those values already found.
!
!    Output, REAL ( kind = 8 ) Z(N,MM), the associated set of orthonormal
!    eigenvectors.  If an error exit is made, Z contains those vectors already
!    found.
!
!    Output, INTEGER ( kind = 4 ) IERR, error flag.
!    0, normal RETURN.
!    3*N+1, if M exceeds MM.
!    4*N+R, if the eigenvector corresponding to the R-th
!      eigenvalue fails to converge in 5 iterations.
!
  IMPLICIT NONE

  INTEGER ( kind = 4 ) mm
  INTEGER ( kind = 4 ) n

  REAL ( kind = 8 ) d(n)
  REAL ( kind = 8 ) e(n)
  REAL ( kind = 8 ) e2(n)
  REAL ( kind = 8 ) eps1
  REAL ( kind = 8 ) eps2
  REAL ( kind = 8 ) eps3
  REAL ( kind = 8 ) eps4
  INTEGER ( kind = 4 ) group
  INTEGER ( kind = 4 ) i
  INTEGER ( kind = 4 ) ierr
  INTEGER ( kind = 4 ) ii
  INTEGER ( kind = 4 ) ip
  INTEGER ( kind = 4 ) isturm
  INTEGER ( kind = 4 ) its
  INTEGER ( kind = 4 ) j
  INTEGER ( kind = 4 ) jj
  INTEGER ( kind = 4 ) k
  REAL ( kind = 8 ) lb
  INTEGER ( kind = 4 ) m
  INTEGER ( kind = 4 ) m1
  INTEGER ( kind = 4 ) m2
  REAL ( kind = 8 ) norm
  INTEGER ( kind = 4 ) p
  !REAL ( kind = 8 ) pythag
  INTEGER ( kind = 4 ) q
  INTEGER ( kind = 4 ) r
  REAL ( kind = 8 ) rv1(n)
  REAL ( kind = 8 ) rv2(n)
  REAL ( kind = 8 ) rv3(n)
  REAL ( kind = 8 ) rv4(n)
  REAL ( kind = 8 ) rv5(n)
  REAL ( kind = 8 ) rv6(n)
  INTEGER ( kind = 4 ) s
  REAL ( kind = 8 ) t1
  REAL ( kind = 8 ) t2
  REAL ( kind = 8 ) tst1
  REAL ( kind = 8 ) tst2
  REAL ( kind = 8 ) u
  REAL ( kind = 8 ) ub
  REAL ( kind = 8 ) uk
  REAL ( kind = 8 ) v
  REAL ( kind = 8 ) w(mm)
  REAL ( kind = 8 ) x0
  REAL ( kind = 8 ) x1
  REAL ( kind = 8 ) xu
  REAL ( kind = 8 ) z(n,mm)

  ierr = 0
  s = 0
  t1 = lb
  t2 = ub
!
!  Look for small sub-diagonal entries.
!
  e2(1) = 0.0D+00

  DO i = 2, n

    tst1 = abs ( d(i) ) + abs ( d(i-1) )
    tst2 = tst1 + abs ( e(i) )

    IF( tst2 <= tst1 ) then
      e2(i) = 0.0D+00
    ENDIF

  ENDDO
!
!  Determine the number of eigenvalues in the interval.
!
  p = 1
  q = n
  x1 = ub
  isturm = 1
  go to 320

60 CONTINUE

  m = s
  x1 = lb
  isturm = 2
  go to 320

80 CONTINUE

  m = m - s

  IF( mm < m ) then
    go to 980
  ENDIF

  q = 0
  r = 0
!
!  Establish and process next submatrix, refining interval by the
!  Gerschgorin bounds.
!
100 CONTINUE

  IF( r == m ) then
    go to 1001
  ENDIF

  p = q + 1
  xu = d(p)
  x0 = d(p)
  u = 0.0D+00

  DO q = p, n

     x1 = u
     u = 0.0D+00
     v = 0.0D+00

     IF( q /= n ) then
       u = abs ( e(q+1) )
       v = e2(q+1)
     ENDIF

     xu = min ( d(q) - ( x1 + u ), xu )
     x0 = max ( d(q) + ( x1 + u ), x0 )

     IF( v == 0.0D+00 ) then
       exit
     ENDIF

  ENDDO

  x1 = max ( abs ( xu ), abs ( x0 ) ) * epsilon ( x1 )

  IF( eps1 <= 0.0D+00 ) then
    eps1 = -x1
  ENDIF

  IF( p /= q ) then
    go to 180
  ENDIF
!
!  Check for isolated root within interval.
!
  IF( d(p) < t1 .or. t2 <= d(p) ) then
    go to 940
  ENDIF

  r = r + 1

  z(1:n,r) = 0.0D+00

  w(r) = d(p)
  z(p,r) = 1.0D+00
  go to 940

180 CONTINUE

  u = q - p + 1
  x1 = u * x1
  lb = max ( t1, xu - x1 )
  ub = min ( t2, x0 + x1 )
  x1 = lb
  isturm = 3
  go to 320

200 CONTINUE

  m1 = s + 1
  x1 = ub
  isturm = 4
  go to 320

220 CONTINUE

  m2 = s
  IF( m2 < m1 ) then
    go to 940
  ENDIF
!
!  Find roots by bisection.
!
  x0 = ub
  isturm = 5

  rv5(m1:m2) = ub
  rv4(m1:m2) = lb
!
!  Loop for K-th eigenvalue.
!
  k = m2

250 CONTINUE

  xu = lb

  DO ii = m1, k

    i = m1 + k - ii

    IF( xu < rv4(i) ) then
      xu = rv4(i)
      exit
    ENDIF

  ENDDO

!280 CONTINUE

  x0 = min ( x0, rv5(k) )
!
!  Next bisection step.
!
300 CONTINUE

     x1 = ( xu + x0 ) * 0.5D+00

     IF( ( x0 - xu ) <= abs ( eps1 ) ) then
       go to 420
     ENDIF

     tst1 = 2.0D+00 * ( abs ( xu ) + abs ( x0 ) )
     tst2 = tst1 + (x0 - xu)

     IF( tst2 == tst1 ) then
       go to 420
     ENDIF
!
!  Sturm sequence.
!
320  CONTINUE

     s = p - 1
     u = 1.0D+00

     DO i = p, q

        IF( u /= 0.0D+00 ) then
          go to 325
        ENDIF

        v = abs ( e(i) ) / epsilon ( v )
        IF( e2(i) == 0.0D+00 ) then
          v = 0.0D+00
        ENDIF

        go to 330

325     CONTINUE

        v = e2(i) / u
330     CONTINUE

        u = d(i) - x1 - v
        IF( u < 0.0D+00 ) then
          s = s + 1
        ENDIF

     ENDDO

     go to ( 60,80,200,220,360 ), isturm
!
!  Refine intervals.
!
360  CONTINUE

     IF( k <= s ) then
       go to 400
     ENDIF

     xu = x1

     IF( m1 <= s ) then
       go to 380
     ENDIF

     rv4(m1) = x1
     go to 300

380  CONTINUE

     rv4(s+1) = x1
     IF( x1 < rv5(s) ) then
       rv5(s) = x1
     ENDIF
     go to 300

400  CONTINUE

     x0 = x1
     go to 300
!
!  K-th eigenvalue found.
!
420  CONTINUE

  rv5(k) = x1
  k = k - 1

  IF( m1 <= k ) then
    go to 250
  ENDIF
!
!  Find vectors by inverse iteration.
!
  norm = abs ( d(p) )
  ip = p + 1

  DO i = ip, q
    norm = max ( norm, abs ( d(i) ) + abs ( e(i) ) )
  ENDDO
!
!  EPS2 is the criterion for grouping,
!  EPS3 replaces zero pivots and equal roots are modified by eps3,
!  EPS4 is taken very small to avoid overflow.
!
  eps2 = 0.001D+00 * norm
  eps3 = abs ( norm ) * epsilon ( eps3 )
  uk = q - p + 1
  eps4 = uk * eps3
  uk = eps4 / sqrt ( uk )
  group = 0
  s = p

  DO k = m1, m2

     r = r + 1
     its = 1
     w(r) = rv5(k)
     x1 = rv5(k)
!
!  Look for close or coincident roots.
!
     IF( k /= m1 ) then
       IF( eps2 <= x1 - x0 ) then
         group = -1
       ENDIF
       group = group + 1
       IF( x1 <= x0 ) then
         x1 = x0 + eps3
       ENDIF
     ENDIF
!
!  Elimination with interchanges and initialization of vector.
!
!520  CONTINUE

     v = 0.0D+00

     DO i = p, q

        rv6(i) = uk

        IF( i == p ) then
          go to 560
        ENDIF

        IF( abs ( u ) <= abs ( e(i) ) ) then
          xu = u / e(i)
          rv4(i) = xu
          rv1(i-1) = e(i)
          rv2(i-1) = d(i) - x1
          rv3(i-1) = 0.0D+00
          IF( i /= q ) then
            rv3(i-1) = e(i+1)
          ENDIF
          u = v - xu * rv2(i-1)
          v = -xu * rv3(i-1)
          cycle
        ENDIF

!540     CONTINUE

        xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = 0.0D+00

560     CONTINUE

        u = d(i) - x1 - xu * v

        IF( i /= q ) then
          v = e(i+1)
        ENDIF

     ENDDO

     IF( u == 0.0D+00 ) then
       u = eps3
     ENDIF

     rv1(q) = u
     rv2(q) = 0.0D+00
     rv3(q) = 0.0D+00
!
!  Back substitution.
!
600  CONTINUE

     DO ii = p, q
        i = p + q - ii
        rv6(i) = ( rv6(i) - u * rv2(i) - v * rv3(i) ) / rv1(i)
        v = u
        u = rv6(i)
     ENDDO
!
!  Orthogonalize with respect to previous members of group.
!
     DO jj = 1, group
        j = r - group - 1 + jj
        xu = dot_product ( rv6(p:q), z(p:q,j) )
        rv6(p:q) = rv6(p:q) - xu * z(p:q,j)
     ENDDO

!700  CONTINUE

     norm = sum ( abs ( rv6(p:q) ) )

     IF( 1.0D+00 <= norm ) then
       go to 840
     ENDIF
!
!  Forward substitution.
!
     IF( its == 5 ) then
       ierr = 4 * n + r
       go to 1001
     ENDIF

     IF( norm == 0.0D+00 ) then
       rv6(s) = eps4
       s = s + 1
       IF( q < s ) then
         s = p
       ENDIF
       go to 780
     ENDIF

!740  CONTINUE

    xu = eps4 / norm

     rv6(p:q) = rv6(p:q) * xu
!
!  Elimination operations on next vector iterate.
!
780    CONTINUE
!
!  If rv1(i-1) == e(i), a row interchange was performed earlier in the
!  triangularization process.
!
     DO i = p, q

       u = rv6(i)

       IF( rv1(i-1) == e(i) ) then
         u = rv6(i-1)
         rv6(i-1) = rv6(i)
       ENDIF

       rv6(i) = u - rv4(i) * rv6(i-1)

     ENDDO

     its = its + 1
     go to 600
!
!  Normalize so that sum of squares is 1 and expand to full order.
!
840  CONTINUE

     u = 0.0D+00

     DO i = p, q
       u = pythag ( u, rv6(i) )
     ENDDO

     xu = 1.0D+00 / u

     z(1:n,r) = 0.0D+00
     z(p:q,r) = rv6(p:q) * xu

     x0 = x1

  ENDDO

940 CONTINUE

  IF( q < n ) then
    go to 100
  ENDIF

  go to 1001
!
!  Set error: underestimate of number of eigenvalues in interval.
!
980 CONTINUE

  ierr = 3 * n + 1

1001 CONTINUE

  lb = t1
  ub = t2

  RETURN

END SUBROUTINE


END MODULE EISPACK