MODULE EPI_MATH

  !USE cdf_normal_mod              !! Part of CDFLIB (DCDFLIB IN F90)
  !USE cdf_nc_chisq_mod            !! Part of CDFLIB (DCDFLIB IN F90)
  !USE cdf_chisq_mod               !! Part of CDFLIB (DCDFLIB IN F90)
  USE EISPACK                     !! EISPACK LIBRARY
  USE EPI_UTILS                   !! MODULE WITH VARIOUS ROUTINES
  USE LOGISTICREG                 !! LOGISTIC REGRESSION MODULE

#ifdef _OPENMP
  USE PAR_ZIG_MOD 
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PARAMETERS AND INTERFACES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERFACES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE IsEven
     MODULE PROCEDURE IsEven_ikn, IsEven_iks, IsEven_ikb !, IsEvenL
  END INTERFACE

  INTERFACE IsOdd
     MODULE PROCEDURE IsOdd_ikn, IsOdd_iks, IsOdd_ikb !, IsOddL
  END INTERFACE

  INTERFACE Mean
     MODULE PROCEDURE Mean_r, Mean_ikn, Mean_iks, Mean_ikb
  END INTERFACE

  INTERFACE Divide
     MODULE PROCEDURE Divide_ikn, Divide_iks, Divide_ikb
  END INTERFACE

  INTERFACE Invert
     MODULE PROCEDURE Invert1, InvertN, InvertNM
  END INTERFACE
  
  INTERFACE RandNumber
     MODULE PROCEDURE RandNumberN, RandNumber1
  END INTERFACE
  
  INTERFACE CountTrue
     MODULE PROCEDURE CountTrue1, CountTrue2
  END INTERFACE

  INTERFACE TRAN
     MODULE PROCEDURE TRAN_R, TRAN_I
  END INTERFACE
  
  INTERFACE IsMultiple
     MODULE PROCEDURE IsMultiple_ikn, IsMultiple_iks, IsMultiple_ikb
  END INTERFACE
  
  INTERFACE QuickSort
     MODULE PROCEDURE QuickSortReal, QuickSortInt
  END INTERFACE

  INTERFACE MultiQSort
     MODULE PROCEDURE MultiQSortReal, MultiQSortInt
  END INTERFACE

  INTERFACE 
    FUNCTION FMIN(AX,BX,F,TOL) RESULT(FMIN_0)
      REAL(KIND=8) :: AX
      REAL(KIND=8) :: BX
      REAL(KIND=8) :: F
      EXTERNAL F
      REAL(KIND=8) :: TOL
      REAL(KIND=8) :: FMIN_0
    END FUNCTION FMIN
  END INTERFACE 

!   INTERFACE 
!     FUNCTION FindMin(A,B,FUN,TOL) RESULT(MIN)
!       EXTERNAL FUN
!       REAL(dpp), INTENT(IN) :: A
!       REAL(dpp), INTENT(IN) :: B
!       REAL(dpp), INTENT(IN) :: FUN
!       REAL(dpp), INTENT(IN) :: TOL
!       REAL(dpp)             :: MIN
!       REAL(dpd)             :: A1, B1, FUN1, TOL1, MIN1
!     END FUNCTION FindMin
!   END INTERFACE 

  
 CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SUBROUTINES AND FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION FindMin(fun, a, b) RESULT(min)
  IMPLICIT NONE
  EXTERNAL fun, fmin
  REAL(dpp), INTENT(IN) :: a, b
  REAL(dpd)             :: fun, fmin, ax, bx, tol
  REAL(dpp)             :: min
  
  tol = ten_d*EPSILON(zero_d)  
  ax = DBLE(a) + tol  
  bx = DBLE(b) - tol
  
  print *,"FindMin"
  
  !! Check for invalid bounds
  IF(ax >= bx) &
    CALL Prnt0("Interval must have positive length! (FindMin)", Q=.TRUE.)  
  
  print *,"XXX2"
  !! Find the minimum
  min = REAL(fmin(ax, bx, fun, tol), dpp)
  print *,"XXX3"
  
  RETURN
  
END FUNCTION FindMin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION dpsqrt(x)
  IMPLICIT NONE
  INTRINSIC SQRT
  REAL(dpp), INTENT(IN) :: x
  REAL(dpp)             :: dpsqrt
  
  dpsqrt = REAL(SQRT(DBLE(x)), dpp)
  
  RETURN
  
END FUNCTION dpsqrt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Mean_r(x) RESULT(m)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: x(:)
  REAL(dpp)             :: m
  
  m = SUM(x) / SIZE(x)
  
END FUNCTION Mean_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Mean_iks(x) RESULT(m)
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: x(:)
  REAL(dpp)                :: m
  
  m = SUM(REAL(x, dpp)) / SIZE(x)
  
END FUNCTION Mean_iks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Mean_ikn(x) RESULT(m)
  IMPLICIT NONE
  INTEGER(ikn), INTENT(IN) :: x(:)
  REAL(dpp)                :: m
  
  m = SUM(REAL(x, dpp)) / SIZE(x)
  
END FUNCTION Mean_ikn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Mean_ikb(x) RESULT(m)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: x(:)
  REAL(dpp)                :: m
  
  m = SUM(REAL(x, dpp)) / SIZE(x)
  
END FUNCTION Mean_ikb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Sd(x) RESULT(stdev)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: x(:)
  REAL(dpp)             :: stdev, n
  
  n = SIZE(x)
  IF(n<2) THEN
    stdev = -one
  ELSE
    stdev = SUM((x - SUM(x)/n)**2) / (n-1)
  ENDIF
  
END FUNCTION Sd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION TRAN_R(A)
!! Tranposes a real matrix A
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: A(:,:)
  REAL(dpp)             :: TRAN_R(SIZE(A,2),SIZE(A,1))
  
  TRAN_R = TRANSPOSE(A)
  
END FUNCTION TRAN_R

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION TRAN_I(A)
!! Tranposes an integer matrix A
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: A(:,:)
  INTEGER             :: TRAN_I(SIZE(A,2),SIZE(A,1))
  
  TRAN_I = TRANSPOSE(A)
  
END FUNCTION TRAN_I

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Diag(v)
!! Returns diagonal matrix with diagonal equal to values in 'v'
!! (Right now works only with vectors)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: v(:)
  REAL(dpp)             :: Diag(SIZE(v),SIZE(v))
  INTEGER               :: j
  
    Diag = 0
    DO j=1,SIZE(v,1)
        Diag(j,j) = v(j)
    ENDDO  

END FUNCTION Diag
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Id(v)
!! Returns the unit matrix of dim v if v is natural, throws error otherwise
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: v
  REAL(dpp)           :: Id(v,v)
  INTEGER             :: j
  
  IF(v<=0) CALL Prnt0("Id cannot take a non-positive argument.",Q=.TRUE.)  

  Id = zero
  DO j=1,v
    Id(j,j) = one
  ENDDO  
  
END FUNCTION Id
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsMultiple_iks(i,j) RESULT(mult)
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: i, j
  LOGICAL                  :: mult
    
  mult = .FALSE.
  IF(INT(i/j, iks)*j == i) mult = .TRUE.
  RETURN
    
END FUNCTION IsMultiple_iks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsMultiple_ikn(i,j) RESULT(mult)
  IMPLICIT NONE
  INTEGER(ikn), INTENT(IN) :: i, j
  LOGICAL                  :: mult
    
  mult = .FALSE.
  IF(INT(i/j, ikn)*j == i) mult = .TRUE.
  RETURN
    
END FUNCTION IsMultiple_ikn
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsMultiple_ikb(i,j) RESULT(mult)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: i
  INTEGER, INTENT(IN)      :: j
  LOGICAL                  :: mult
    
  mult = .FALSE.
  IF(INT(i/j, ikb)*j == i) mult = .TRUE.
  RETURN
    
END FUNCTION IsMultiple_ikb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsEven_iks(i) RESULT(even)
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: i
  LOGICAL                  :: even
    
  even = INT(i/2, iks)*2 == i
  RETURN  

END FUNCTION IsEven_iks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsEven_ikn(i) RESULT(even)
  IMPLICIT NONE
  INTEGER(ikn), INTENT(IN) :: i
  LOGICAL                  :: even
    
  even = INT(i/2, ikn)*2 == i
  RETURN  

END FUNCTION IsEven_ikn
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsEven_ikb(i) RESULT(even)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: i
  LOGICAL                  :: even
    
  even = INT(i/2, ikb)*2 == i
  RETURN  

END FUNCTION IsEven_ikb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! FUNCTION IsEvenL(i) RESULT(even)
!   IMPLICIT NONE
!   INTEGER(ikl), INTENT(IN) :: i
!   LOGICAL                  :: even
!     
!   even = INT(i/2, ikl)*2 == i
!   RETURN  
! 
! END FUNCTION IsEvenL
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsOdd_iks(i) RESULT(odd)
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: i
  LOGICAL                  :: odd
    
  odd = INT((i+1)/2, iks)*2 == i+1
  RETURN
    
END FUNCTION IsOdd_iks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsOdd_ikn(i) RESULT(odd)
  IMPLICIT NONE
  INTEGER(ikn), INTENT(IN) :: i
  LOGICAL             :: odd
    
  odd = INT((i+1)/2, ikn)*2 == i+1
  RETURN
    
END FUNCTION IsOdd_ikn
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsOdd_ikb(i) RESULT(odd)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: i
  LOGICAL                  :: odd
    
  odd = INT((i+1)/2, ikb)*2 == i+1
  RETURN
    
END FUNCTION IsOdd_ikb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! FUNCTION IsOddL(i) RESULT(odd)
!   IMPLICIT NONE
!   INTEGER(ikl), INTENT(IN) :: i
!   LOGICAL                  :: odd
!     
!   odd = INT((i+1)/2, ikl)*2 == i
!   RETURN
!     
! END FUNCTION IsOddL
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION CountOccurrence(A, v) RESULT(n)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: A(:), v
  INTEGER               :: n, i

  n = 0
  DO i=1,SIZE(A)
    IF(A(i) /= v) CYCLE
    n = n + 1
  ENDDO

END FUNCTION CountOccurrence
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION CountTrue1(A) RESULT(n)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: A(:)
  INTEGER             :: n, i

  n = 0
  DO i=1,SIZE(A)
    IF(.NOT.A(i)) CYCLE
    n = n + 1
  ENDDO

END FUNCTION CountTrue1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION CountTrue2(A) RESULT(n)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: A(:,:)
  INTEGER             :: n, i, j

  n = 0
  DO i=1,SIZE(A,1)
    DO j=1,SIZE(A,2)
      IF(.NOT.A(i,j)) CYCLE
      n = n + 1
    ENDDO
  ENDDO

END FUNCTION CountTrue2
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Multiply(a, b) RESULT(prod)
!! Multiply two reals in an underflow-safe way
  IMPLICIT NONE
  REAL(dpp), PARAMETER  :: eps = SQRT(tiny(one))
  REAL(dpp), INTENT(IN) :: a, b
  REAL(dpp)             :: prod
  
  prod = a * b
  RETURN 

  prod = zero 
  IF(a==zero .OR. b==zero) RETURN
  prod = eps 
  IF(ABS(a) < eps .AND. ABS(b)<eps) RETURN
  IF(LOG10(ABS(a))+LOG10(ABS(b)) < 100*one) RETURN

  prod = a * b 
  
END FUNCTION Multiply
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Divide_ikn(a, b) RESULT(frac)
!! Divides two integers
  IMPLICIT NONE
  INTEGER(ikn), INTENT(IN) :: a, b
  REAL(dpp)                :: frac
  
  frac = zero
  IF(a /= 0) frac = REAL(REAL(a, dpp) * Invert(REAL(b, dpp)), dpp) 
  
END FUNCTION Divide_ikn
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Divide_iks(a, b) RESULT(frac)
!! Divides two integers
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: a, b
  REAL(dpp)                :: frac
  
  frac = zero
  IF(a /= 0) frac = REAL(REAL(a, dpp) * Invert(REAL(b, dpp)), dpp) 
  
END FUNCTION Divide_iks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Divide_ikb(a, b) RESULT(frac)
!! Divides two integers
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: a, b
  REAL(dpp)                :: frac
  
  frac = zero
  IF(a /= 0) frac = REAL(REAL(a, dpp) * Invert(REAL(b, dpp)), dpp) 
  
END FUNCTION Divide_ikb
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Invert1(v) RESULT(inv)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: v
  REAL(dpp)             :: inv
  
  inv = zero
  IF(ABS(v) > epstol) inv = one / v
  
END FUNCTION Invert1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION InvertN(v) RESULT(inv)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: v(:)
  REAL(dpp)             :: inv(SIZE(v))
  INTEGER               :: j
  
  IF(SIZE(v)<=0) &
    CALL Prnt0("Argument must have positive length (InvertN).", Q=.TRUE.)  

  inv = zero
  DO j=1,SIZE(v)
    IF(ABS(v(j)) > epstol) inv(j) = one / v(j)
  ENDDO  
  
END FUNCTION InvertN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION InvertNM(v) RESULT(inv)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: v(:,:)
  REAL(dpp)             :: inv(SIZE(v,1),SIZE(v,2))
  INTEGER               :: i,j
  
  IF(SIZE(v)<=0) CALL PrntE("Argument must have positive length.", Q=.TRUE.)  

  inv = zero
  DO i=1,SIZE(v,1)
    DO j=1,SIZE(v,2)
      IF(ABS(v(i,j)) > epstol) inv(i,j) = one / v(i,j)
    ENDDO
  ENDDO  
  
END FUNCTION InvertNM
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION pnorm(x, m, s) RESULT(p)
!! Returns the p-value of distribution N(m,s^2)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)           :: x
  REAL(dpp), INTENT(IN), OPTIONAL :: m, s
  REAL(dpp)                       :: p, y
  
  !! Work with local variable
  y = x
  
  !! Center and/or standardize
  IF(PRESENT(m)) y = y - m
  IF(PRESENT(s)) y = y / s

  !! Return standard normal distribution p-value 
  p = REAL(ERFC(DBLE(- y / SQRT(two_d))) / two, dpp)
  
  RETURN
  
END FUNCTION pnorm
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION qnorm(p, m, s, status) RESULT(x)
!! Returns the p-value of distribution N(m,s^2)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: p
  REAL(dpp), INTENT(IN), OPTIONAL  :: m, s
  REAL(dpp)                        :: x, s1, m1
  !REAL(dpcdf)                      :: p1, s1, m1
  INTEGER, INTENT(INOUT), OPTIONAL :: status
  INTEGER                          :: status1
  
  !! Return standard normal distribution quantile 
  CALL ppnd16(p, x, status1)
  
  !! Transform the variance and mean
  m1 = zero 
  s1 = one 
  IF(PRESENT(m)) m1 = m 
  IF(PRESENT(s)) s1 = s 
  x = x*s1 + m1

  !! Check for errors
  IF(status1 /= 0) THEN
    CALL PrntE("Non-zero status in qnorm (code "//i2cp(status1)//&
                    ")! See source code of module cdf_normal_mod for help.")
    IF(PRESENT(status)) THEN
      status = status1
      RETURN
    ELSE
      CALL Prnt0(Q=.TRUE.)
    ENDIF
  ENDIF

  !!! Return normal distribution quantile 
  !m1 = REAL(zero, dpcdf) 
  !s1 = REAL(one, dpcdf) 
  !IF(PRESENT(m)) m1 = REAL(m, dpcdf) 
  !IF(PRESENT(s)) s1 = REAL(s, dpcdf) 
  !p1 = REAL(p, dpcdf) 
  !x = REAL(INV_NORMAL(p1, mean=m1, sd=s1), dpp)
  
  RETURN
  
END FUNCTION qnorm
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION norm_pval(x, m, s) RESULT(p)
!! Returns the p-value of distribution N(m,s^2)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)           :: x
  REAL(dpp), INTENT(IN), OPTIONAL :: m, s
  REAL(dpp)                       :: p, y
  
  y = x
  
  !! Center and/or standardize
  IF(PRESENT(m)) y = y - m
  IF(PRESENT(s)) y = y / s

  !! Return standard normal distribution p-value 
  p = MAX(REAL(ERFC( DBLE( ABS(y) / SQRT(two_d) ) ), dpp), ten**(-99))
  
  RETURN
  
END FUNCTION norm_pval
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION pchisq(x, df, ncp, status) RESULT(p)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: x
  INTEGER, INTENT(IN)              :: df
  REAL(dpp), INTENT(IN), OPTIONAL  :: ncp
  INTEGER, INTENT(INOUT), OPTIONAL :: status
  REAL(dpp)                        :: p, dfr
  INTEGER                          :: ifault

  ifault = 0
  !! Call the (non-central) chisquare distribution function
  IF(PRESENT(ncp)) THEN
    dfr = REAL(df, dpp) 
    p = chi2nc(x, dfr, ncp, ifault)
  ELSE
    p = chi_squared(x, df)
  ENDIF
  
  !! Print error or return status number
  IF(ifault/=0) THEN
    IF(PRESENT(status)) THEN
      status = 3000 + ifault
    ELSE
      CALL PrntE("Non-zero return code in pchisq (error code "//&
                      i2cp(ifault)//")!", Q=.TRUE.)
    ENDIF
  ELSE
    IF(PRESENT(status)) status = 0
  ENDIF 

  RETURN
  
END FUNCTION pchisq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION qchisq(p, df, status, silent) RESULT(x)
!! Returns the p-value of distribution N(m,s^2)
  IMPLICIT NONE
  REAL(dpp), PARAMETER             :: min_p = 0.000000000002_dpp, &
                                      max_p = 0.999999999998_dpp
  !REAL(dpp), PARAMETER             :: min_p = 0.00000002_dpp, &
  !                                    max_p = 0.99999998_dpp
  REAL(dpp), INTENT(IN)            :: p
  INTEGER, INTENT(IN)              :: df
  INTEGER, INTENT(INOUT), OPTIONAL :: status
  LOGICAL, INTENT(IN), OPTIONAL    :: silent
  REAL(dpp)                        :: x, p1, dfr, g, a, b 
  INTEGER                          :: ifault
  LOGICAL                          :: warn
  !REAL(dpcdf)                      :: p1, df1, dummy
  
  warn = .TRUE.
  IF(PRESENT(silent)) THEN
    warn = .NOT.silent
  ENDIF
  
  p1 = p
  !! Check for too big or small probability
  IF(p1<=min_p .OR. p1>=max_p) THEN
    IF(warn) &
      CALL PrntE("Argument 'p' out of range (qchisq). The value of 'p'"//&
             " must be between "//TRIM(r2c(min_p))//" and "//TRIM(r2c(max_p))//&
             " while the current value is "//TRIM(r2c(p))//"! Approximating"//&
             " the p-quantile by the nearest allowed value.")
    p1 = MIN(MAX(p1, min_p), max_p)
  ENDIF
              
  !! Check for negative degrees of freedom
  IF(df<=0) &
    CALL PrntE("Argument 'df' must be non-negative while the current"//&
                    " value of 'df' is "//i2cp(df)//".", Q=.TRUE.)

  !! Call the chisquare quantile function
  dfr = REAL(df, dpp)
  g = lngamma(dfr / two)
  x = ppchi2(p1, dfr, g, ifault)

  !! If maximum number of iterations error occurred approximate the chi-square 
  !! distribution function on (0,0.15) by (x^2)/8.4
  IF(ifault == 3 .AND. df==4 .AND. p<three*tenth/two) THEN
    a = one / 11
    b = 1.8_dpp
    x = (p1/a)**(one/b)
    ifault = 0
  ENDIF

  !! Check the ifault and print errors
  IF(ifault /= 0) THEN
    IF(PRESENT(status)) THEN
      IF(ifault == 1) THEN
        status = 2001
      ELSEIF(ifault == 2) THEN
        status = 2002
      ELSEIF(ifault == 3) THEN
        status = 2003
      ELSE
        status = 2099
      ENDIF
    ELSE
      IF(ifault == 1) THEN
        CALL PrntE("p must be between 0.000002 and 0.999998 (qchisq).", &
                        Q=.TRUE.)
      ELSEIF(ifault == 2) THEN
        CALL PrntE("Negative degrees of freedom (qchisq).", Q=.TRUE.)
      ELSEIF(ifault == 3) THEN
        CALL PrntE("Maximum number of iterations exceeded! (p="//&
                        TRIM(r2c(p))//", v="//i2cp(df)//", g="//&
                        TRIM(r2c(g))//")", Q=.TRUE.)
      ELSE
        CALL PrntE("Non-zero ifault in qchisq (code "//i2cp(ifault)//&
                        ")! See source code of module cdf_chisq_mod for help.",&
                        Q=.TRUE.)
      ENDIF
    ENDIF
  ELSE
    IF(PRESENT(status)) status = 0
  ENDIF

  RETURN
  
END FUNCTION qchisq
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE RandNumber1(rnumber)
!! Generates n random numbers either by a standard RANDOM_NUMBER intrinsic
!! or by an OPENMP optimized 
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT) :: rnumber 
#ifdef _OPENMP
  INTEGER                :: thread, OMP_GET_THREAD_NUM
#endif
  
#ifdef _OPENMP
  
  !! Use the openmp routine par_uni
  thread = OMP_GET_THREAD_NUM()
  rnumber = par_uni(thread)
  
#else
  
  !! Use fortran intrinsic function RANDOM_NUMBER
  CALL RANDOM_NUMBER(rnumber)
  
#endif

  RETURN   

END SUBROUTINE RandNumber1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
SUBROUTINE RandNumberN(rnumbers)
!! Generates n random numbers either by a standard RANDOM_NUMBER intrinsic
!! or by an OPENMP optimized 
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT) :: rnumbers(:)
  INTEGER                :: i 
  
  DO i=1,SIZE(rnumbers)
    CALL RandNumber1(rnumbers(i))
  ENDDO

  RETURN   

END SUBROUTINE RandNumberN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
SUBROUTINE RandSample(S,n,low,high)
!! Generates a random integer sequence: S(1), S(2), ... , S(N)
!! such that each element is in the closed interval (low,high) and
!! sampled with replacement.
  IMPLICIT NONE
  INTEGER   :: n, S(n), low, high, IX, i
  REAL(dpp) :: U, X
  
  !! Sample N times uniform and rescale it and round it  
  DO i=1,N
    CALL RandNumber(U)
    X = REAL((high+1)-low, dpp)*U + REAL(low, dpp)
    IX = INT(X)
    IF(X<0 .AND. IX/=X) IX = INT(X - one)
    S(i) = IX
  ENDDO
  
  RETURN

END SUBROUTINE RandSample
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
SUBROUTINE Permutate(P)
!! Generates a random permutation P of length n of 0's and 1's 
!! with delta the ratio of 1's.
  IMPLICIT NONE
  INTEGER, PARAMETER     :: nr = 100
  INTEGER, INTENT(INOUT) :: P(:)
  REAL(dpp)              :: U(nr)
  INTEGER                :: i, j, ipj, k, temp, M, n
       
  n = SIZE(P)
  !! Permutate while generating up to nr uniform(0,1) numbers at a time
  DO i=1,n,nr
    M = MIN(n-i+1,nr)
    CALL RandNumber(U(1:M))
    DO j=1,M
      ipj = i+j-1
      k = INT(U(j)*(n-ipj+1))+ipj
      temp = P(ipj)
      P(ipj) = P(k)
      P(k) = temp
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE Permutate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
SUBROUTINE GetRandPerm(P,n)
!! Generates a random permutation of 1, 2, ..., n 
  IMPLICIT NONE
  INTEGER, PARAMETER   :: nr = 100
  INTEGER, INTENT(IN)  :: n
  INTEGER, INTENT(OUT) :: P(n)
  REAL(dpp)            :: U(nr)
  INTEGER              :: i, j, ipj, k, temp, M
       
    DO i=1,n
      P(i)=i
    END DO

    !! Permutate while generating up to nr U(0,1) numbers at a time.
    DO i=1,n,nr
       M = MIN(n-i+1,nr)
       CALL RandNumber(U(1:M))
       DO j=1,M
         ipj = i+j-1
         k = INT(U(j)*(n-ipj+1))+ipj
         temp = P(ipj)
         P(ipj) = P(k)
         P(k) = temp
       ENDDO
    ENDDO

    RETURN

END SUBROUTINE GetRandPerm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetOwnSeed(seed, seed_size)
! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
  IMPLICIT NONE
#ifdef __INTEL_COMPILER
  INTEGER, EXTERNAL :: GETPID
#endif
  INTEGER, INTENT(OUT) :: seed(:)
  INTEGER, INTENT(IN)  :: seed_size
  INTEGER              :: i, dt(8), pid, t(2), s
  INTEGER(8)           :: count, tms
            
  CALL system_clock(count)
  IF(count /= 0) THEN
    t = transfer(count, t)
  ELSE
    CALL date_and_time(values=dt)
    tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
         + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
         + dt(3) * 24 * 60 * 60 * 60 * 1000 &
         + dt(5) * 60 * 60 * 1000 &
         + dt(6) * 60 * 1000 + dt(7) * 1000 &
         + dt(8)
    t = transfer(tms, t)
  ENDIF
  s = ieor(t(1), t(2))
  pid = GETPID() + 1099279 ! Add a prime
  s = ieor(s, pid)
  IF(seed_size >= 3) THEN
    seed(1) = t(1) + 36269
    seed(2) = t(2) + 72551
    seed(3) = pid
    IF(seed_size > 3) THEN
       seed(4:) = s + 37 * (/ (i, i = 0, seed_size - 4) /)
    ENDIF
  ELSE
    seed = s + 37 * (/ (i, i = 0, seed_size - 1 ) /)
  ENDIF
               
END SUBROUTINE GetOwnSeed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SetSeed(seed, init, get, put, useed_int, useed_par, start, const, &
                   nthreads, announce)
!! Sets or stores the current seed for random number generation for both the
!! internal and the par_zig (OPENMP enabled) random number generators.
!! If OPENMP available and the number of threads is larger than 1 use
!! module PAR_ZIG_MOD initialization. Otherwise, use ordinary way of
!! initializing the random seed.
  IMPLICIT NONE
#ifdef __INTEL_COMPILER
  INTEGER, EXTERNAL                             :: GETPID
#endif
  INTEGER, PARAMETER                            :: c0 = 37 
#ifdef _OPENMP
  INTEGER, PARAMETER                            :: grainsize = 32
#endif
  INTEGER, INTENT(INOUT), ALLOCATABLE, OPTIONAL :: seed(:), useed_int(:), &
                                                   useed_par(:)
  LOGICAL, INTENT(IN), OPTIONAL                 :: init, get, put, announce
  !INTEGER(ikb), INTENT(IN), OPTIONAL            :: const
  INTEGER, INTENT(IN), OPTIONAL                 :: const, start, nthreads
  INTEGER                                       :: i, j, a, pid, c, nthreads1, &
                                                   nseed_int, nseed_par
  INTEGER, ALLOCATABLE                          :: lseed(:), lseed_int(:), &
                                                   lseed_par(:)  
  CHARACTER(mntl)                               :: txt, txt2
  
  nthreads1 = 1
  IF(PRESENT(nthreads)) nthreads1 = MAX(nthreads,1)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                       INITIALIZE NEW SEED                              !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IF(PRESENT(init)) THEN
    IF(init) THEN

      !! Set seed size to the number of threads for (par_zig)
      nseed_par = nthreads1 

      !! Get the seed size for internal random number generator
      CALL RANDOM_SEED(SIZE=nseed_int)

      !! Set the local seed sizes
      CALL ResizeVar(lseed_int, nseed_int, -1)
      CALL ResizeVar(lseed_par, nseed_par, -1)
      
      !! If user specified a positive value in 'useed_int', use it as the seed
      IF(PRESENT(useed_int)) THEN
        IF(ALLOCATED(useed_int)) THEN
          IF(ANY(useed_int>0)) THEN
            IF(SIZE(useed_int)==1) THEN
              lseed_int = useed_int(1) + 7 * (/ (i, i=0,nseed_int-1) /)
            ELSEIF(SIZE(useed_int) < nseed_int) THEN
              CALL PrntE("User specified seed is too short. Length "//TRIM(i2cp(nseed_int))//" required.")
            ELSE
              lseed_int = useed_int(1:nseed_int)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      
      !! If user specified a positive value in 'useed_par', use it as the seed
      IF(PRESENT(useed_par)) THEN
        IF(ALLOCATED(useed_par)) THEN
          IF(ANY(useed_par>0)) THEN
            IF(SIZE(useed_par)==1) THEN
              lseed_par = useed_par(1) + 7 * (/ (i, i=0,nseed_par-1) /)
            ELSEIF(SIZE(useed_par) < nseed_par) THEN
              CALL PrntE("User specified seed is too short. Length "//TRIM(i2cp(nseed_par))//" required.")
            ELSE
              lseed_par = useed_par(1:nseed_par)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
     
      !! Generate the local seeds for the two random number generators
      DO j=1,2
      
        !! Otherwise generate seed based on system clock
        IF((j==1 .AND. ANY(lseed_int<=0)) .OR. (j==2 .AND. ANY(lseed_par<=0))) THEN
        
          !! Set the base constant
          c = c0
          IF(PRESENT(const)) THEN
            IF(const/=0) c = const
          ENDIF
          
          !! Set the initial constant: IF start present and non-zero, the value 
          !! will be equal to start, otherwise use system clock to set it
          CALL SYSTEM_CLOCK(a)
          pid = GETPID()
          a = INT(abs( mod((a*181)*((pid-83)*359), 104729) ))
        
          IF(PRESENT(start)) THEN
            IF(start/=0) a = start
          ENDIF
        
          !! Set lseed to a time based value
          IF(j==1) lseed_int = a + c * (/ (i+5 , i=1,nseed_int) /)
          IF(j==2) lseed_par = a + c * (/ (i+5 , i=1,nseed_par) /)
          
        ENDIF
        
      ENDDO

      !! Initialize the random seed based on lseed_int value
#ifdef _OPENMP

      !! This uses external seed initialization
      CALL par_zigset(nthreads1, lseed_par, grainsize)
    
#endif

!#else
    
      !! This will initialize seed in G95 and ifort
      !CALL RANDOM_SEED()
      !CALL RANDOM_SEED(GET = lseed_int)
      !!! THIS IS FOR GFORTRAN
      !CALL GetOwnSeed(lseed_int, nseed_int)
      CALL RANDOM_SEED(PUT=lseed_int)
      
!#endif
    
    !! Announce what the seed is
    IF(PRESENT(announce)) THEN
      IF(announce) THEN

#ifdef _OPENMP
        txt = " (OPENMP)"
        txt2 = "--seed-omp"
        CALL ResizeVar(lseed, nseed_par, -1)
        lseed = lseed_par
#else
        txt = ""
        txt2 = "--seed"
        CALL ResizeVar(lseed, nseed_int, -1)
        lseed = lseed_int 
#endif

        !CALL Prnt("Random seed 'start' and 'step': "//TRIM(i2cp(a))//" "//TRIM(i2cp(c)))
        txt = "Random seed"//TRIM(txt)//" initialized to:"
        DO i=1,SIZE(lseed)
          txt = TRIM(txt)//" "//TRIM(i2cp(lseed(i)))
        ENDDO
        CALL Prnt(TRIM(txt))
        
        txt = "Note: The seed can be set via '"//TRIM(txt2)//" "//TRIM(i2cp(lseed(1)))
        DO i=2,SIZE(lseed)
          txt = TRIM(txt)//":"//TRIM(i2cp(lseed(i)))
        ENDDO
        txt = TRIM(txt)//"'"
        CALL Prnt(txt)
        
        DEALLOCATE(lseed)
        
      ENDIF
    ENDIF      
  
    DEALLOCATE(lseed_int, lseed_par)
  
    ENDIF    
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                       RETRIEVE THE CURRENT SEED                        !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Store the current seed
  IF(PRESENT(get)) THEN
    IF(get) THEN
    
      !! Get the current seed
      IF(PRESENT(seed)) THEN

#ifdef _OPENMP

        CALL ResizeVar(seed, SIZE(par_seed))
        seed = par_seed
        
#else
    
        !! Read the size of the seed
        CALL RANDOM_SEED(SIZE=nseed_int)
      
        CALL ResizeVar(seed, nseed_int)
      
        !! Read the value
        CALL RANDOM_SEED(GET=seed)
        
#endif
      
      ENDIF
      
    ENDIF
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                          SET THE GIVEN SEED                            !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Set the random seed to SEED and destroy SEED
  IF(PRESENT(put)) THEN
    IF(put) THEN
    
      IF(.NOT.PRESENT(seed)) &
        CALL PrntE("Missing 'seed' (SetSeed)", Q=.TRUE.)
      IF(.NOT.ALLOCATED(seed)) &
        CALL PrntE("'seed' not allocated (SetSeed)", Q=.TRUE.)
        
#ifdef _OPENMP

      IF(SIZE(seed)<nthreads1) &
        CALL PrntE("'seed' too small (SetSeed)", Q=.TRUE.)

      CALL par_zigset(SIZE(seed), seed, grainsize)
      txt = " (OPENMP)"

#else

      IF(SIZE(seed)/=12) &
        CALL PrntE("'seed' must have length 12 (SetSeed)", Q=.TRUE.)
      
      CALL RANDOM_SEED(PUT=seed)
      txt = ""
      
#endif
            
      !! Announce what the seed is
      IF(PRESENT(announce)) THEN
        IF(announce) THEN
          txt = "Random seed"//TRIM(txt)//" set to:"
          DO i=1,SIZE(seed)
            txt = TRIM(txt)//" "//TRIM(i2cp(seed(i)))
          ENDDO
          CALL Prnt(txt, skip1=1, skip2=1)
        ENDIF
      ENDIF
      
    ENDIF   
  ENDIF

  RETURN

END SUBROUTINE SetSeed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PseudoInv(M, I, Sq, SqI, tol, ierr)
  !! Returns a pseudoinverse matrix of a square matrix M and possibly also
  !! a square root and square root of the pseudoinverse
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: M(:,:)
  REAL(dpp), INTENT(OUT), OPTIONAL :: I(:,:), Sq(:,:), SqI(:,:)
  REAL(dpp), INTENT(IN), OPTIONAL  :: tol
  INTEGER, INTENT(INOUT), OPTIONAL :: ierr
  REAL(dpp)                        :: M_e(SIZE(M,1)), M_v(SIZE(M,1),SIZE(M,2))
  REAL(dpp)                        :: E(SIZE(M,1)), eps!, trace, sume
  INTEGER                          :: k, l!, ierr

  IF(SIZE(M,1)/=SIZE(M,2)) THEN
    CALL PrntE("Argument is not a square matrix (PseudoInv).")
    IF(PRESENT(ierr)) THEN
      ierr = 20000
      RETURN
    ELSE
      CALL Prnt0(Q=.TRUE.)
    ENDIF
  ENDIF
  
  !! Set the tol for zero values
  eps = epstol
  IF(PRESENT(tol)) THEN
    IF(tol>zero) eps = tol
  ENDIF
  
  !print *,"           X1"
  
  !! Get spectral decomposition of V
  !CALL GetEigen(M, M_e, M_v, ierr)
  CALL GetEigRS(M, M_e, M_v, ierr)
  IF(ierr/=0) RETURN

  !print *,"           X2"
  !! Get spectral decomposition of V
  !M_v = M
  !CALL kaiser(M_v, SIZE(M_v,1), SIZE(M_v,2), M_e, trace, sume, ierr)

  !! Check for (truly) negative eigenvalues which defy positive semidefiniteness
  IF(ANY(M_e < -eps)) THEN
    CALL PrntE("Cannot compute square-root of M unless it is positive semidefinite (PseudoInverse).")
    WRITE(*,'(A)') "Eigenvalues:"
    WRITE(*,'(ES11.3)') M_e
    WRITE(10,*) M_e
    WRITE(*,'(A)') "Tolerance:"
    WRITE(*,'(ES11.3)') eps
    WRITE(*,'(A)') "Matrix M:"
    DO k=1,SIZE(M,1)
      DO l=1,SIZE(M,2)
        WRITE(*,'(ES11.3)', ADVANCE='NO') M(k,l)
        WRITE(11,*) M(k,l)
      ENDDO
      WRITE(*,'(A)') ""
    ENDDO
    WRITE(*,'(A)') "Matrix M_v:"
    DO k=1,SIZE(M_v,1)
      DO l=1,SIZE(M_v,2)
        WRITE(*,'(ES11.3)', ADVANCE='NO') M_v(k,l)
        WRITE(12,*) M_v(k,l)
      ENDDO
      WRITE(*,'(A)') ""
    ENDDO

    CALL GetEigRS(M, M_e, M_v, ierr)
    WRITE(*,'(A)') "Eigenvalues:"
    WRITE(*,'(ES11.3)') M_e
    WRITE(20,*) M_e
    WRITE(*,'(A)') "Tolerance:"
    WRITE(*,'(ES11.3)') eps
    WRITE(*,'(A)') "Matrix M:"
    DO k=1,SIZE(M,1)
      DO l=1,SIZE(M,2)
        WRITE(*,'(ES11.3)', ADVANCE='NO') M(k,l)
        WRITE(21,*) M(k,l)
      ENDDO
      WRITE(*,'(A)') ""
    ENDDO
    WRITE(*,'(A)') "Matrix M_v:"
    DO k=1,SIZE(M_v,1)
      DO l=1,SIZE(M_v,2)
        WRITE(*,'(ES11.3)', ADVANCE='NO') M_v(k,l)
        WRITE(22,*) M_v(k,l)
      ENDDO
      WRITE(*,'(A)') ""
    ENDDO

    CALL Prnt0(Q=.TRUE.)
    IF(PRESENT(ierr)) THEN
      ierr = 20001
      RETURN
    ELSE
      CALL Prnt0(Q=.TRUE.)
    ENDIF
  ENDIF
  
  !! Set the zero tolerance
  eps = MIN(MAX(eps, MAXVAL(ABS(M_e)) * eps), 1e-4)

  if(.false.) then
    print *,""
    print *,""
    print *,"-----------------------------------------------"
    print *,"eps=", eps
    print *,"M_e=", M_e
  endif

  !! Nulify tiny values in eigenvectors and eigenvalues
  DO k=1,SIZE(M_v,1)
    IF(ABS(M_e(k)) <= eps) M_e(k) = zero
    DO l=1,SIZE(M_v,2)
      IF(ABS(M_v(k,l)) <= eps) M_v(k,l) = zero
    ENDDO
  ENDDO

  if(.false.) then
    print *,"M_e=", M_e
    print *,""
    if(present(I) .and. present(Sq) .and. present(SqI)) then
      print *,""
      print *,"eps=", eps
      stop
    endif
  endif
  
  !! Compute pseudoinverse and/or squareroot and/or squareroot pseudoinverse
  IF(PRESENT(I) .OR. PRESENT(Sq) .OR. PRESENT(SqI)) THEN
    
    !! Compute the pseudoinverse of M (spectral decomposition)
    IF(PRESENT(I)) THEN
      E = Invert(ABS(M_e))
      I = MATMUL(MATMUL(M_v, Diag(E)), TRANSPOSE(M_v))
      if(.false.) then
        print *,"M_v=", MAXVAL(M_v)
        print *,"E=", E
        print *,"I=", I(SIZE(I,1),:)
      endif
    ENDIF

    !! Compute squareroot matrix of M (spectral decomposition)
    IF(PRESENT(Sq)) THEN
      E = REAL(SQRT(DBLE(ABS(M_e))), dpp)
      Sq = MATMUL(MATMUL(M_v, Diag(E)), TRANSPOSE(M_v))
    ENDIF

    !! Compute square-root of the pseudoinverse of M (spectral decomposition)
    IF(PRESENT(SqI)) THEN
      E = REAL(SQRT(DBLE(Invert(ABS(M_e)))), dpp)
      SqI = MATMUL(MATMUL(M_v, Diag(E)), TRANSPOSE(M_v))
    ENDIF
    
  ENDIF

  !print *,"           X3"
  RETURN
    
END SUBROUTINE PseudoInv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SquareRoot(M, M_Sq, tol, ierr)
  !! Returns a pseudoinverse matrix of a square matrix M
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: M(:,:)
  REAL(dpp), OPTIONAL, INTENT(IN)  :: tol
  REAL(dpp), INTENT(OUT)           :: M_Sq(:,:)
  INTEGER, INTENT(INOUT), OPTIONAL :: ierr
  REAL(dpp)                        :: M_e(SIZE(M,1)), M_v(SIZE(M,1),SIZE(M,2)),&
                                      E(SIZE(M,1)), eps!, trace, sume
  INTEGER                          :: k

  IF(SIZE(M,1)/=SIZE(M,2)) THEN
    CALL PrntE("Argument is not a square matrix (SquareRoot).")
    IF(PRESENT(ierr)) THEN
      ierr = 10001
      RETURN
    ELSE
      CALL Prnt0(Q=.TRUE.)
    ENDIF
  ENDIF
  
  eps = epstol
  IF(PRESENT(tol)) THEN
    IF(tol>zero) eps = tol
  ENDIF
  
  !! Get spectral decomposition of V
  !CALL GetEigen(M, M_e, M_v, ierr)
  CALL GetEigRS(M, M_e, M_v, ierr)
  IF(ierr/=0) RETURN

  !! Get spectral decomposition of V
  !M_v = M
  !CALL kaiser(M_v, SIZE(M_v,1), SIZE(M_v,2), M_e, trace, sume, ierr)

  
  !! Check for (truly) negative eigenvalues which defy positive 
  !! semidefiniteness
  IF(ANY(M_e < -eps)) THEN
    CALL PrntE("Argument is not a positive semidefinite matrix"//&
                   " (SquareRoot).")
    PRINT *,"Eigenvalues:", M_e
    PRINT *,"Tolerance:", eps
    IF(PRESENT(ierr)) THEN
      ierr = 10002
      RETURN
    ELSE
      CALL Prnt0(Q=.TRUE.)
    ENDIF
  ENDIF

  !! Nulify tiny eigenvalues
  DO k=1,SIZE(M_e)
    IF(M_e(k) < eps) M_e(k) = zero
  ENDDO

  !! Take square root of only non-negative eigenvalues
  E = REAL(SQRT(DBLE(ABS(M_e))), dpp)
  
  !! Compute the pseudoinverse from spectral decomposition
  M_Sq = MATMUL(MATMUL(M_v,Diag(E)),TRANSPOSE(M_v))
  
  RETURN
    
END SUBROUTINE SquareRoot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION FindDet(matrix, n) RESULT(det)
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)       :: n
  REAL(dpp), DIMENSION(n,n) :: matrix
  REAL(dpp)                 :: det, m, temp
  INTEGER                   :: i, j, k, l
  LOGICAL                   :: DetExists = .TRUE.

    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
      IF(matrix(k,k) == 0) THEN
        DetExists = .FALSE.
        DO i = k+1, n
          IF(matrix(i,k) /= 0) THEN
            DO j = 1, n
              temp = matrix(i,j)
              matrix(i,j)= matrix(k,j)
              matrix(k,j) = temp
            END DO
            DetExists = .TRUE.
            l=-l
            EXIT
          ENDIF
        END DO
        IF(DetExists .EQV. .FALSE.) THEN
          det = 0
          return
        ENDIF
      ENDIF
      DO j = k+1, n
        m = matrix(j,k)/matrix(k,k)
        DO i = k+1, n
          matrix(j,i) = matrix(j,i) - m*matrix(k,i)
        END DO
      END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    det = l
    DO i = 1, n
      det = det * matrix(i,i)
    END DO
        
END FUNCTION FindDet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!! INCLUDE THE EISPACK LIBRARY ROUTINES USED BY 'GetEigRS' !!!

!#include "epi_inc_eispack.F90"

SUBROUTINE GetEigRS(M, D, V, ierr)
!! This subroutine calculates the eigenvalues and eigenvectors of a real 
!! symmetric matrix M using EISPACK's subroutine 'rs'.
  IMPLICIT NONE

  REAL(dpp), INTENT(IN)            :: M(:,:)
  REAL(dpp), INTENT(OUT)           :: D(:), V(:,:)
  !INTEGER                          :: NROT
  INTEGER, INTENT(INOUT), OPTIONAL :: ierr
  REAL(dpd)                        :: A(SIZE(M,1),SIZE(M,2))
  INTEGER                          :: N, get_ev, ierror

  !print *,"   E1"
  !! Check symmetry
  N = SIZE(M,1)
  IF(N /= SIZE(M,2)) CALL PrntE("M is not square (GetEigRS)",Q=.TRUE.)
  
  !! Make a local copy of M
  A(1:N,1:N) = M(1:N,1:N)
  
  !! Tell rs to also get the eigenvectors
  get_ev = 1
  
  !! Calculate the eigenvalues and eigenvectors
  CALL rs(N, A, D, get_ev, V, ierror)
  
  !! Check for error code
  IF(ierror /= 0) THEN
    IF(PRESENT(ierr)) THEN
      ierr = 30002
    ELSE
      CALL PrntE("Non-zero code '"//i2cp(ierror)//" returned by 'rs' (GetEigRS).", Q=.TRUE.)
    ENDIF
  ENDIF
  
  !print *,"   E2"
  RETURN
  
END SUBROUTINE GetEigRS   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetEigen(M, D, V, ierr, niter)
  !*************************************************************
  !* This subroutine computes all eigenvalues and eigenvectors *
  !* of a real symmetric square matrix A(N,N). On output, D(N) *
  !* returns the eigenvalues of matrix A. V(N,N) contains      *
  !* the eigenvectors of A by columns. The normalization to    *
  !* unity is made by main program before printing results.    *
  !* NROT could return the number of Jacobi matrix             *
  !* rotations which were required.                            *
  !* --------------------------------------------------------- *
  !* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
  !*       chap. 11, pages 346-348" [BIBLI 08].                *
  !* --------------------------------------------------------- *
  !* Edited by Jakub Pecanka for the purposes of EpiDetector   *
  !*************************************************************
  IMPLICIT NONE

  INTEGER, PARAMETER               :: maxitr = 300
  REAL(dpd), PARAMETER             :: & !small_l = 1.0E-12, 
                                      zero_l = 0.0_dpd, &
                                      half_l = 0.5_dpd, &
                                      one_l = 1.0_dpd, &
                                      ten_l = 10.0_dpd, &
                                      eps = SQRT(SQRT(tiny(one_l)))
  REAL(dpp), INTENT(IN)            :: M(:,:)
  REAL(dpp), INTENT(OUT)           :: D(:), V(:,:)
  INTEGER                          :: NROT
  INTEGER, INTENT(OUT), OPTIONAL   :: niter
  INTEGER, INTENT(INOUT), OPTIONAL :: ierr
  REAL(dpd)                        :: A(SIZE(M,1),SIZE(M,2)), &
                                      B(SIZE(M,1)), Z(SIZE(M,1)), &
                                      c,g,h,s,sm,t,tau,theta,tresh
  INTEGER                          :: N, i, j, ip, iq

  !! Work with local variables
  A = REAL(M, dpd)           
  N = SIZE(M,1)
  
  theta = zero_l
  
  IF(N /= SIZE(M,2)) CALL PrntE("M is not square (GetEigen)",Q=.TRUE.)
  
  !! initialize V to an identity matrix
  V = zero_l
  DO j=1,N
    V(j,j) = one_l
  ENDDO  
  
  !! initialize vector Z to zero
  Z = zero_l
      
  !! initialize vector B to the diagonal of A
  DO ip=1,N
    B(ip) = A(ip,ip)
  ENDDO
  
  !! initialize vector D to the diagonal of A
  D = REAL(B, dpp)
  
  NROT=0
  DO i=1,maxitr
  
    sm=zero_l
    DO ip=1,N-1     !sum off-diagonal elements
      DO iq=ip+1,N
        sm = sm + ABS(A(ip,iq))
      ENDDO
    ENDDO

    IF(sm==zero_l) RETURN  !normal return

    IF(i > 4) THEN
      tresh = (2*one_l/ten_l) * sm**2
    ELSE
      tresh = zero_l
    ENDIF

    DO ip=1,N-1
      DO iq=ip+1,N
        g = 100 * one_l * ABS(A(ip,iq))
        
        ! after 4 sweeps, skip the rotation IF the off-diagonal element is small
        IF((i > 4) .AND. (ABS(D(ip))+g == ABS(D(ip))) .AND. &
        (ABS(D(iq))+g == ABS(D(iq)))) THEN
          A(ip,iq) = zero_l
        ELSEIF(ABS(A(ip,iq)) > tresh) THEN

          h = REAL(D(iq)-D(ip), dpd)
          
          IF(ABS(h)+g == ABS(h)) THEN
            t = A(ip,iq) / h
          ELSE
            theta = half_l*h/A(ip,iq)  
            t = one_l / (ABS(theta)+SQRT(one_l+theta**2))
            IF(theta > zero_l) t = -t
          ENDIF

          c = one_l / SQRT(one_l+t**2)
          s = t*c
          tau = s/(one_l+c)
          h = t * A(ip,iq)
          Z(ip) = Z(ip) - h
          Z(iq) = Z(iq) + h
          D(ip) = D(ip) - REAL(h, dpp)
          D(iq) = D(iq) + REAL(h, dpp)
          A(ip,iq) = zero_l

          DO j=1,ip-1
            g = A(j,ip)
            h = A(j,iq)
            !! Update A in a underflow-safe way
            A(j,ip) = g - Multiply(s, h + g*tau)
            A(j,iq) = h + Multiply(s, g - h*tau)
            
          ENDDO

          DO j=ip+1, iq-1
            g = A(ip,j)
            h = A(j,iq)
            !! Update A in a underflow-safe way
            A(ip,j) = g - Multiply(s, h + g*tau)
            A(j,iq) = h + Multiply(s, g - h*tau)
          ENDDO                     

          DO j=iq+1,N
            g = A(ip,j)
            h = A(iq,j)
            !! Update A in a underflow-safe way
            A(ip,j) = g - Multiply(s, h + g*tau)
            A(iq,j) = h + Multiply(s, g - h*tau)
          ENDDO                 

          DO j=1,N
            g = REAL(V(j,ip), dpd)
            h = REAL(V(j,iq), dpd)
            V(j,ip) = REAL(g - Multiply(s, h + g*tau), dpp)
            V(j,iq) = REAL(h + Multiply(s, g - h*tau), dpp)
          ENDDO                 

          NROT = NROT+1
          
        ENDIF !IF((i > 4)...
        
      ENDDO ! DO iq=ip+1,N
    ENDDO ! DO ip=1,N-1
    
    DO ip=1,N
      B(ip) = B(ip) + Z(ip)
      D(ip) = B(ip)
      Z(ip) = zero_l
    ENDDO

  ENDDO !main i loop
  
  IF(PRESENT(niter)) niter = NROT

  !! Should have returned by now, print error
  CALL PrntE("Maximum iterations ("//i2cp(maxitr)//") reached (GetEigen).")
  IF(PRESENT(ierr)) THEN
    ierr = 30001
  ELSE
    CALL Prnt0(Q=.TRUE.)
  ENDIF

  RETURN

END SUBROUTINE GetEigen   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE kaiser(a, nrows, n, eigenv, trace, sume, ier)
  !*************************************************************************
  !  EIGENVALUES AND VECTORS OF A SYMMETRIC +VE DEFINITE MATRIX,
  !  USING KAISER'S METHOD.
  !  REFERENCE: KAISER,H.F. 'THE JK METHOD: A PROCEDURE FOR FINDING THE
  !  EIGENVALUES OF A REAL SYMMETRIC MATRIX', COMPUT.J., VOL.15, 271-273, 1972.
  
  !  ARGUMENTS:-
  !  A       = INPUT, AN ARRAY CONTAINING THE MATRIX
  !            OUTPUT, THE COLUMNS OF A CONTAIN THE NORMALIZED EIGENVECTORS
  !            OF A.   N.B. A IS OVERWRITTEN !
  !  NROWS   = INPUT, THE FIRST DIMENSION OF A IN THE CALLING PROGRAM.
  !  N       = INPUT, THE ORDER OF A, I.E. NO. OF COLUMNS.
  !            N MUST BE <= NROWS.
  !  EIGENV()= OUTPUT, A VECTOR CONTAINING THE ORDERED EIGENVALUES.
  !  TRACE   = OUTPUT, THE TRACE OF THE INPUT MATRIX.
  !  SUME    = OUTPUT, THE SUM OF THE EIGENVALUES COMPUTED.
  !            N.B. ANY SYMMETRIC MATRIX MAY BE INPUT, BUT IF IT IS NOT +VE
  !            DEFINITE, THE ABSOLUTE VALUES OF THE EIGENVALUES WILL BE FOUND.
  !            IF TRACE = SUME, THEN ALL OF THE EIGENVALUES ARE POSITIVE
  !            OR ZERO.   IF SUME > TRACE, THE DIFFERENCE IS TWICE THE SUM OF
  !            THE EIGENVALUES WHICH HAVE BEEN GIVEN THE WRONG SIGNS !
  !  IER     = OUTPUT, ERROR INDICATOR
  !             = 0 NO ERROR
  !             = 1 N > NROWS OR N < 1
  !             = 2 FAILED TO CONVERGE IN 10 ITERATIONS
  
  !  LATEST REVISION - 6 September 1990
  !  Fortran 90 version - 20 November 1998
  !*************************************************************************
  
  IMPLICIT NONE
  
  REAL(dpp), INTENT(IN OUT) :: a(:,:)
  INTEGER, INTENT(IN)        :: nrows
  INTEGER, INTENT(IN)        :: n
  REAL(dpp), INTENT(OUT)    :: eigenv(:)
  REAL(dpp), INTENT(OUT)    :: trace
  REAL(dpp), INTENT(OUT)    :: sume
  INTEGER, INTENT(OUT)       :: ier
  
  ! Local variables
  
  REAL(dpp), PARAMETER :: small_l = 1.0E-12, zero_l = 0.0D0, half_l = 0.5D0,&
                          one_l = 1.0D0
  INTEGER              :: i, iter, j, k, ncount, nn
  REAL(dpp)            :: absp, absq, COS, ctn, eps, halfp, p, q, SIN, ss, &
                          TAN, temp, xj, xk

    !   CALCULATE CONVERGENCE TOLERANCE, EPS.
    !   CALCULATE TRACE.   INITIAL SETTINGS.
    
    ier = 1
    IF(n < 1 .OR. n > nrows) RETURN
    ier = 0
    iter = 0
    trace = zero_l
    ss = zero_l
    DO j = 1,n
      trace = trace + a(j,j)
      DO i = 1,n
        ss = ss + a(i,j)**2
      END DO
    END DO
    sume = zero_l
    eps = small_l*ss/n
    nn = n*(n-1)/2
    ncount = nn
    
    !   ORTHOGONALIZE PAIRS OF COLUMNS J & K, K > J.
    
    20 DO j = 1,n-1
      DO k = j+1,n
        
    !   CALCULATE PLANAR ROTATION REQUIRED
        
        halfp = zero_l
        q = zero_l
        DO i = 1,n
          xj = a(i,j)
          xk = a(i,k)
          halfp = halfp + xj*xk
          q = q + (xj+xk) * (xj-xk)
        END DO
        p = halfp + halfp
        absp = ABS(p)
        
    !   IF P is very small, the vectors are almost orthogonal.
    !   Skip the rotation if Q >= 0 (correct ordering).
        
        IF(absp < eps .AND. q >= zero_l) THEN
          ncount = ncount - 1
          IF(ncount <= 0) GO TO 160
          CYCLE
        ENDIF
        
    !   Rotation needed.
        
        absq = ABS(q)
        IF(absp <= absq) THEN
          TAN = absp/absq
          COS = one_l/SQRT(one_l + TAN*TAN)
          SIN = TAN*COS
        ELSE
          ctn = absq/absp
          SIN = one_l/SQRT(one_l + ctn*ctn)
          COS = ctn*SIN
        ENDIF
        COS = SQRT((one_l + COS)*half_l)
        SIN = SIN/(COS + COS)
        IF(q < zero_l) THEN
          temp = COS
          COS = SIN
          SIN = temp
        ENDIF
        IF(p < zero_l) SIN = -SIN
        
    !   PERFORM ROTATION
        
        DO i = 1,n
          temp = a(i,j)
          a(i,j) = temp*COS + a(i,k)*SIN
          a(i,k) = -temp*SIN + a(i,k)*COS
        END DO
      END DO
    END DO
    ncount = nn
    iter = iter + 1
    IF(iter < 10) GO TO 20
    ier = 2
    
    !   CONVERGED, OR GAVE UP AFTER 10 ITERATIONS
    
    160 DO j = 1,n
      temp = SUM( a(1:n,j)**2 )
      eigenv(j) = SQRT(temp)
      sume = sume + eigenv(j)
    END DO
    
    !   SCALE COLUMNS TO HAVE UNIT LENGTH
    
    DO j = 1,n
      IF(eigenv(j) > zero_l) THEN
        temp = one_l/eigenv(j)
      ELSE
        temp = zero_l
      ENDIF
      a(1:n,j) = a(1:n,j)*temp
    END DO
    
    RETURN

END SUBROUTINE kaiser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!
! Source : http://en.wikibooks.org/wiki/Algorithm_Implementation/
!                 Sorting/Quicksort#FORTRAN_90.2F95
!
! ***********************************
! *
  SUBROUTINE QuickSortReal(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * IF present Ipt, a pointer with the 
! * changes is returned in Ipt.
! ***********************************
 
  TYPE Limits
     INTEGER :: Ileft, Iright
  END TYPE Limits

  ! For a list with Isw number of elements or
  ! less use Insrt
  INTEGER, PARAMETER :: Isw = 10

  REAL(kind=dpp), INTENT(INOUT) :: X(:)
  INTEGER, INTENT(OUT), OPTIONAL :: Ipt(:)

  INTEGER :: I, Ipvn, Ileft, Iright, ISpos, ISmax
  INTEGER, ALLOCATABLE :: IIpt(:)
  TYPE (Limits), ALLOCATABLE :: Stack(:)

  ALLOCATE(Stack(SIZE(X)))

  Stack(:)%Ileft = 0
  IF(PRESENT(Ipt)) THEN
     
    FORALL(I=1:SIZE(Ipt)) Ipt(I) = I
    
    ! Initialize the stack
    Ispos = 1
    Ismax = 1
    Stack(ISpos)%Ileft  = 1
    Stack(ISpos)%Iright = SIZE(X)
    
    DO WHILE (Stack(ISpos)%Ileft /= 0)

      Ileft = Stack(ISPos)%Ileft
      Iright = Stack(ISPos)%Iright
      
      IF(Iright-Ileft <= Isw) THEN
        CALL InsrtLC(X, Ipt, Ileft,Iright)
      ELSE
        Ipvn = ChoosePiv(X, Ileft, Iright)
        Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
        
        Stack(ISmax+1)%Ileft = Ileft
        Stack(ISmax+1)%Iright = Ipvn-1
        Stack(ISmax+2)%Ileft = Ipvn + 1
        Stack(ISmax+2)%Iright = Iright
        ISmax = ISmax + 2
      ENDIF
  
      ISpos = ISpos + 1
      IF(ISpos>SIZE(X)) EXIT

    ENDDO

  ELSE

    ! Initialize the stack
    Ispos = 1
    Ismax = 1
    Stack(ISpos)%Ileft  = 1
    Stack(ISpos)%Iright = SIZE(X)
    
    ALLOCATE(IIpt(SIZE(X)))
    DO WHILE (Stack(ISpos)%Ileft /= 0)

      Ileft = Stack(ISPos)%Ileft
      Iright = Stack(ISPos)%Iright
      
      IF(Iright-Ileft <= Isw) THEN
         CALL InsrtLC(X, IIpt, Ileft, Iright)
      ELSE
         Ipvn = ChoosePiv(X, Ileft, Iright)
         Ipvn = Partition(X, Ileft, Iright, Ipvn)
  
         Stack(ISmax+1)%Ileft = Ileft
         Stack(ISmax+1)%Iright = Ipvn-1
         Stack(ISmax+2)%Ileft = Ipvn + 1
         Stack(ISmax+2)%Iright = Iright
         ISmax = ISmax + 2
      ENDIF
      
      ISpos = ISpos + 1
      IF(ISpos>SIZE(X)) EXIT
        
     ENDDO
     DEALLOCATE(IIpt)

  ENDIF

  DEALLOCATE(Stack)

  RETURN
 
  CONTAINS
 
    ! ***********************************
    INTEGER FUNCTION ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright) for QuickSortReal.
    ! * This routine chooses the median of the first, last and mid 
    ! * element of the list.
    ! ***********************************
 
      REAL(kind=dpp), INTENT(IN) :: XX(:)
      INTEGER, INTENT(IN) :: IIleft, IIright
 
      REAL(kind=dpp) :: XXcp(3)
      INTEGER :: IIpt(3), IImd
 
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)
 
      CALL InsrtLC(XXcp, IIpt, 1, 3)
 
      SELECT CASE (IIpt(2))
      CASE (1)
         IIpv = IIleft
      CASE (2)
         IIpv = IImd
      CASE DEFAULT ! Same as 'CASE (3)'
         IIpv = IIright
      END SELECT
 
      RETURN
    END FUNCTION ChoosePiv
 
    ! ***********************************
    SUBROUTINE InsrtLC(XX, IIpt, II_left, II_right)
    ! ***********************************
    ! * Perform an insertion sort of the list 
    ! * XX(:) between index values II_left and II_right.
    ! * IIpt(:) returns the permutations made to sort.
    ! ***********************************
 
      REAL(kind=dpp), INTENT(INOUT) :: XX(:)
      INTEGER, INTENT(INOUT) :: IIpt(:)
      INTEGER, INTENT(IN) :: II_left, II_right
 
      REAL(kind=dpp) :: RRtmp
      INTEGER :: II, JJ
 
      DO II = II_left+1, II_right
        RRtmp = XX(II)
        DO JJ = II-1, 1, -1
          IF(RRtmp < XX(JJ)) THEN
            XX(JJ+1) = XX(JJ)
            CALL Swap_IN(IIpt, JJ, JJ+1)
          ELSE
            EXIT
          ENDIF
        END DO
        XX(JJ+1) = RRtmp
      ENDDO
 
      RETURN
    END SUBROUTINE InsrtLC
 
END SUBROUTINE QuickSortReal

! ***********************************

SUBROUTINE QuickSortInt(X, Ipt)
  INTEGER, INTENT(INOUT) :: X(:)
  INTEGER, INTENT(OUT), OPTIONAL :: Ipt(:)
  REAL(dpp) :: Y(SIZE(X))
  INTEGER :: Ord(SIZE(X)), i
  
  Y = REAL(X, dpp)
  Ord = (/ (i, i=1,SIZE(X)) /)
  CALL QuickSortReal(Y, Ord)
  X = X(Ord)
  IF(PRESENT(Ipt)) THEN
    Ipt = Ord
  ENDIF
  RETURN

END SUBROUTINE QuickSortInt

! ***********************************
! *
INTEGER FUNCTION Partition(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by QuickSortReal.
! ***********************************

  REAL(kind=dpp), INTENT(INOUT) :: X(:)
  INTEGER, INTENT(IN) :: Ileft, Iright, Ipv
  INTEGER, INTENT(INOUT), OPTIONAL :: Ipt(:)

  REAL(kind=dpp) :: Rpv
  INTEGER :: I

  Rpv = X(Ipv)
  CALL Swap(X, Ipv, Iright)
  IF(PRESENT(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
  Ipvfn = Ileft

  IF(PRESENT(Ipt))  THEN
     DO I = Ileft, Iright-1
        IF(X(I) <= Rpv) THEN
           CALL Swap(X, I, Ipvfn)
           CALL Swap_IN(Ipt, I, Ipvfn)
           Ipvfn = Ipvfn + 1
        ENDIF
     END DO
  ELSE
     DO I = Ileft, Iright-1
        IF(X(I) <= Rpv) THEN
           CALL Swap(X, I, Ipvfn)
           Ipvfn = Ipvfn + 1
        ENDIF
     END DO
  ENDIF

  CALL Swap(X, Ipvfn, Iright)
  IF(PRESENT(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)

  RETURN
END FUNCTION Partition
 
! ***********************************
! *
SUBROUTINE Swap(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************

  REAL(kind=dpp), INTENT(INOUT) :: X(:)
  INTEGER, INTENT(IN) :: I, J

  REAL(kind=dpp) :: Itmp

  Itmp = X(I)
  X(I) = X(J)
  X(J) = Itmp

  RETURN
END SUBROUTINE Swap
 
! ***********************************
! *
SUBROUTINE Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
  INTEGER, INTENT(INOUT) :: X(:)
  INTEGER, INTENT(IN) :: I, J

  INTEGER :: Itmp

  Itmp = X(I)
  X(I) = X(J)
  X(J) = Itmp

  RETURN
END SUBROUTINE Swap_IN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE MultiQSortReal(X)

  IMPLICIT NONE
  REAL(dpp), INTENT(INOUT) :: X(:,:)
  INTEGER                  :: order(SIZE(X,1)), i, n, start
  
  n = SIZE(X,1)
  IF(n==1) RETURN
  
  !! Sort the first column of the array
  order = (/ (i, i=1,n) /)
  CALL qsort(X(:,1), n, order)
  !CALL QuickSort(X(:,1), order)
  
  !! Reorder the second column to match the newly ordered first column
  X(:,2) = X(order,2)
  
  !! Sort the second column within runs of same values in the first column
  start = 1
  DO i = 2,n

    !! Run continues -> go to the next one
    IF(i<n .AND. X(i,1)==X(start,1)) CYCLE

    !! Run ended due to a different value in 1st column -> sort
    IF(X(i,1)/=X(start,1)) THEN
      CALL qsort(X(start:i-1,2), i-start, order(start:i-1))
      !CALL QuickSort(X(start:i-1,2))
    !! Run ended because were at the last element
    ELSE
      CALL qsort(X(start:i,2), i-start+1, order(start:i))
      !CALL QuickSort(X(start:i,2))
    ENDIF

    start = i

  ENDDO
  
  RETURN
  
END SUBROUTINE MultiQSortReal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE MultiQSortInt(X)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: X(:,:)
  REAL(dpp)              :: Y(SIZE(X,1),SIZE(X,2))
  
  Y = X
  CALL MultiQSortReal(Y)
  X = INT(Y+half)
  
  RETURN
  
END SUBROUTINE MultiQSortInt

END MODULE EPI_MATH
