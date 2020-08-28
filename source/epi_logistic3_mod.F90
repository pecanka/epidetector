MODULE LogisticReg
! A package for fitting linear logistic models by iteratively re-weighted
! least squares.

! The model fitted is that the probability of a `success' is:

!      p = 1 - 1/[1 + exp(b0 + b1.X1 + b2.X2 + ... + bk.Xk)]

! where X1, X2, ... , Xk are the predictor variables, and the coefficients
! b0, b1, b2, ..., bk are to be determined.

! N.B. The residual variance used here in estimating  standard errors is the
!      larger of  1.0 and that from the weighted least squares calculations;
!      it is not the theoretical residual variance (1.0) assuming a binomial
!      distribution about the logistic curve.     If the fit of the logistic
!      is poor, the standard errors of the coefficients in the logistic will
!      be much larger.

! The calculation of chi-squared was corrected on 9 August 2003.
! My thanks to Zhang Guan

! By Alan Miller
! amiller @ bigpond.net.au
! users.bigpond.net.au/amiller/

! Latest revision - 9 August 2003

!!!! USE lsq

IMPLICIT NONE

  !!! MOVED FROM LSQ !!!

  ! Note. dp is being set to give at least 12 decimal digit representation of
  !       floating point numbers.   This should be adequate for most problems 
  !       except the fitting of polynomials. dp is being set here so that the 
  !       same code can be run on PCs and Unix systems,    which will usually 
  !       represent floating-point numbers in   'double precision', and other
  !       other systems with larger word lengths,     which will give similar
  !       accuracy in `single precision'.

  INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(12,60)

  !INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(16,60)
  !INTEGER, PARAMETER           :: dp = SELECTED_REAL_KIND(13)     ! returns 8
  !INTEGER, PARAMETER           :: dp = KIND(1.0D0)
  
  !!! MOVED FROM LSQ !!!
  
  LOGICAL, PARAMETER           :: report_logistic_errors = .FALSE.

PUBLIC  :: logistic

CONTAINS


SUBROUTINE logistic(ngroups, x, k, s, n, chisq, devnce, ndf, beta, se_beta,  &
                    ier, cov_beta, fit, stdres, deveps)

! Input arguments:

! ngroups     The number of groups of observations.

! x(:,:)      x(i,j) contains the value of the j-th predictor for group i.
!             The X-values are the same for all the n(i) trials in group i.

! k           The number of predictors.   N.B. A constant will be fitted in
!             the model; do not include it in the k predictor variables.

! s(:), n(:)  s(i) is the number of `successes' out of the n(i) `trials' in
!             the i-th group.   In many applications, each n(i) = 1, that is
!             each case will be treated individually, and then s(i) = 0 or 1
!             for each of the two possible outcomes.

  INTEGER, INTENT(IN)    :: ngroups, k, s(:), n(:)
  REAL (dp), INTENT(IN)  :: x(:,:)

! Output arguments:

! chisq       The value of chi-squared on exit, when a model has been fitted.

! devnce      The deviance on exit, when a model has been fitted.

! ndf         Number of degrees of freedom.

! beta(0:)    The fitted coefficients in the logistic model.

! se_beta(0:) Approximate standard errors of the beta coefficients.

! ier         Error indicator
!             = 0 successful termination
!             = 1 if ngroups < 2 or ndf < 0
!             = 2 if any n(i) < 0
!             = 3 if any s(i) < 0
!             = 4 if any s(i) > n(i)
!             = 5 if any X-variable is constant
!             = 6 if a singularity is detected
!             = 7 if any beta(i) is tending to +/- infinity
!             = 8 failure to converge in maxiter iterations

  REAL (dp), INTENT(OUT)  :: chisq, devnce, beta(0:), se_beta(0:)
  INTEGER, INTENT(OUT)    :: ndf, ier

! Optional output arguments:

! cov_beta(0:,0:)     Approximate covariance matrix of the fitted coefficients.

! fit(:)              The fitted probabilities of a success for each group.

! stdres(:)           Vector of standardized residuals.

! deveps              Determines precision of comparisons with zero

  REAL (dp), INTENT(OUT), OPTIONAL  :: cov_beta(0:,0:), fit(:), stdres(:)
  REAL (dp), INTENT(IN), OPTIONAL   :: deveps

! ************************************ !
! Variables moved here from module LSQ !
! ************************************ !

  INTEGER                      :: nobs, ncol, r_dim
  INTEGER, ALLOCATABLE         :: vorder(:), row_ptr(:)
  LOGICAL                      :: tol_set, rss_set
  
  REAL (dp), ALLOCATABLE       :: d(:), rhs(:), r(:), tol(:), rss(:)
  REAL (dp)                    :: zero = 0.0_dp, one = 1.0_dp, vsmall
  REAL (dp)                    :: sserr

! ************************************ !
! Variables moved here from module LSQ !
! ************************************ !

! Local variables

! maxiter     Number of iterations after which it produces error 8

  INTEGER    :: i, iter, j, ncov, pos, maxiter
  REAL (dp)  :: propn(ngroups), p(ngroups), wt(ngroups), xrow(0:k), db(0:k), &
                bnew(0:k), dev_new, xb, pnew(ngroups), wnew(ngroups), a, b,  &
                range(k), var, e(ngroups), hii, eps
  LOGICAL    :: lindep(0:k)
  REAL (dp), ALLOCATABLE  :: covmat(:)

    ! Set up values every time logistic is called
    tol_set = .FALSE.
    rss_set = .FALSE.
    
    ! Initial checks
    
    ier = 0
    ndf = ngroups - k - 1
    IF (ngroups < 2 .OR. ndf < 0) THEN
      ier = 1
      RETURN
    END IF
    IF (ANY(n(1:ngroups) < 0)) THEN
      ier = 2
      RETURN
    END IF
    IF (ANY(s(1:ngroups) < 0)) THEN
      ier = 3
      RETURN
    END IF
    IF (ANY(s(1:ngroups) > n(1:ngroups))) THEN
      ier = 4
      RETURN
    END IF
    
    ! Calculate ranges of the X-variables to use in testing whether any beta
    ! is tending to +/- infinity.  Also test that no variable is constant.
    
    DO i = 1, k
      a = MAXVAL(x(1:ngroups,i))
      b = MINVAL(x(1:ngroups,i))
      range(i) = a - b
      IF (range(i) < EPSILON(0.0_dp) * (ABS(a) + ABS(b))) THEN
        ier = 5
        RETURN
      END IF
    END DO
    
    ! Start with all beta's = 0 and weights = 1.
    
    beta(0:k) = 0.0_dp
    wt(1:ngroups) = 1.0_dp
    p(1:ngroups) = 0.5_dp
    
    ! propn stores the sample proportions, i.e. s(i) / n(i)
    propn(1:ngroups) = REAL(s(1:ngroups), KIND=dp) / n(1:ngroups)
    iter = 1
    
    maxiter = 20
    IF(.NOT.PRESENT(deveps)) THEN
      eps = 0.0001_dp
    ELSE
      eps = deveps
    ENDIF
    
    ! Start of iterative cycle
    
    DO
      CALL startup(k, .TRUE.)
      DO i = 1, ngroups
        IF (iter == 1) THEN
          xrow(0) = 0.25_dp
          xrow(1:k) = 0.25_dp*x(i, 1:k)
        ELSE
          xrow(0) = p(i)*(1.0_dp - p(i))
          xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
        END IF
        CALL includ(wt(i), xrow, propn(i)-p(i))
      END DO

    ! Test for a singularity
    
      IF (.NOT. tol_set) CALL tolset()
      CALL sing(lindep, ier)
    
      IF (ier /= 0) THEN
        DO i = 1, k
          IF (lindep(i) .AND. report_logistic_errors) THEN
            WRITE(*, '(a, i6, a)') ' Variable number ', i,  &
                                   ' is linearly dependent upon earlier variables'
          END IF
        END DO
        ier = 6
        CALL cleanup()
        RETURN
      END IF
    
      CALL regcf(db, k+1, ier)
      10 bnew = beta(0:k) + db
    
    ! Calculate new p(i)'s, weights & deviance
    
      dev_new = 0.0_dp
      DO i = 1, ngroups
        xb = DOT_PRODUCT( x(i,1:k), bnew(1:k) ) + bnew(0)
        xb = EXP(xb)
        pnew(i) = xb / (1.0_dp + xb)
        wnew(i) = REAL(n(i), KIND=dp) / (pnew(i)*(1.0_dp - pnew(i)))
        IF (iter == 1) wnew(i) = SQRT(wnew(i))
        IF (s(i) > 0) dev_new = dev_new + s(i)*LOG(propn(i)/pnew(i))
        IF (s(i) < n(i)) dev_new = dev_new +   &
                               (n(i)-s(i))*LOG((1.0_dp-propn(i))/(1.0_dp-pnew(i)))
      END DO
      dev_new = 2 * dev_new
    
    ! If deviance has increased, reduce the step size.
    
      IF (iter > 2) THEN
        IF (dev_new > devnce*(1+eps)) THEN
          db = 0.5_dp * db
          GO TO 10
        END IF
      END IF
    
    ! Replace betas, weights & p's with new values
    
      beta(0:k) = bnew(0:k)
      wt = wnew
      p(1:ngroups) = pnew
    
    ! Test for convergence
    
      IF (iter > 2 .AND. devnce - dev_new < eps) EXIT
      devnce = dev_new
      iter = iter + 1
      IF (iter > maxiter) THEN
        ier = 8
        CALL cleanup()
        RETURN
      END IF
    
    ! Test for a very large beta
    
      DO i = 1, k
        IF (ABS(beta(i))*range(i) > 30.0_dp .AND. report_logistic_errors) THEN
          WRITE(*, '(a, i4, a)') ' Coefficient for variable no.', i,  &
                                 ' tending to infinity'
          ier = 7
          CALL cleanup()
          RETURN
        END IF
      END DO
    
    END DO
    
    e = n(1:ngroups)*p(1:ngroups)
    chisq = SUM( (s(1:ngroups) - e)**2 / (e * (1.0_dp - p(1:ngroups))) )
    devnce = dev_new
    
    ! Calculate the approximate covariance matrix for the beta's, if ndf > 0.
    
    IF (ndf > 0) THEN
      ncov = (k+1)*(k+2)/2
      ALLOCATE( covmat(ncov) )
      var = 0.0_dp
      CALL cov(k+1, var, covmat, ncov, se_beta, ier)
      IF (var < 1.0_dp) THEN
        covmat = covmat / var
        se_beta = se_beta / SQRT(var)
      END IF
    
      IF(PRESENT(cov_beta)) THEN
        pos = 1
        DO i = 0, k
          cov_beta(i,i) = covmat(pos)
          pos = pos + 1
          DO j = i+1, k
            cov_beta(i,j) = covmat(pos)
            cov_beta(j,i) = covmat(pos)
            pos = pos + 1
          END DO
        END DO
      END IF
    END IF
    
    IF(PRESENT(fit)) fit(1:ngroups) = p
    
    IF(PRESENT(stdres)) THEN
      DO i = 1, ngroups
        xrow(0) = p(i)*(1.0_dp - p(i))
        xrow(1:k) = p(i)*(1.0_dp - p(i))*x(i, 1:k)
        hii = 0.0_dp
        CALL hdiag(xrow, k+1, hii, ier)
        stdres(i) = (s(i) - n(i)*p(i)) /   &
                    SQRT(n(i)*p(i)*(1.0_dp - p(i))*(1.0_dp - hii))
      END DO
    END IF
    
    IF (ALLOCATED(covmat)) DEALLOCATE( covmat )

    CALL cleanup()
    RETURN

!!! INCLUDE THE CONTENT OF MODIFIED MODULE LSQ INSIDE THIS SUBROUTINE
!!! TO MAKE ITS VARIABLES LOCAL TO THIS SUBROUTINE FOR THE PURPOSES
!!! OF PARALLEL RUNNING OF MULTIPLE COMPUTATION THREADS (OPENMP)

#include "epi_inc_lsq.F90"

!!! END OF INCLUSION

END SUBROUTINE logistic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Chisquare distribution function functions follow  !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION chi_squared(chi2, ndf) RESULT(prob)
! Previously called chi_squared
! Calculate the chi-squared distribution function
! ndf  = number of degrees of freedom
! chi2 = chi-squared value
! prob = probability of a chi-squared value <= chi2 (i.e. the left-hand
!        tail area)

  INTEGER, INTENT(IN)    :: ndf
  REAL (dp), INTENT(IN)  :: chi2
  REAL (dp)              :: prob
  
  ! Local variables
  REAL (dp) :: half = 0.5_dp, x, p
  
    x = half * chi2
    p = half * REAL(ndf, KIND=dp)
    prob = gammad(x, p)
    RETURN

END FUNCTION chi_squared

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION chi2nc(x, f, theta, status) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-06-16  Time: 00:34:40

! ALGORITHM AS 275 APPL.STATIST. (1992), VOL.41, NO.2

! Computes the noncentral chi-square distribution function with positive
! real degrees of freedom f and nonnegative noncentrality parameter theta

IMPLICIT NONE
!INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp), INTENT(IN)  :: f
REAL (dp), INTENT(IN)  :: theta
INTEGER, INTENT(OUT), OPTIONAL :: status        !! Added by Jakub
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: ifault                   
LOGICAL    :: flag
REAL (dp)  :: lam, n, u, v, x2, f2, t, term, bound

! EXTERNAL alogam

INTEGER, PARAMETER    :: itrmax = 150 !itrmax = 50
REAL (dp), PARAMETER  :: errmax = 1.0E-12, zero = 0.0_dp, one = 1.0_dp,  &
                         two = 2.0_dp

IF(PRESENT(status)) status = 0

fn_val = x
ifault = 2
IF (f <= zero .OR. theta < zero) GO TO 999
ifault = 3
IF (x < zero) GO TO 999
ifault = 0
IF (x == zero) RETURN
lam = theta / two

!       Evaluate the first term

n = one
u = EXP(-lam)
v = u
x2 = x / two
f2 = f / two
t = x2 ** f2 * EXP(-x2) / EXP(lngamma(f2 + one))
term = v * t
fn_val = term

!       Check if (f+2n) is greater than x

flag = .false.
10 IF (f + two * n - x <= zero) GO TO 30

!       Find the error bound and check for convergence

flag = .true.
20 bound = t * x / (f + two * n - x)
IF (bound > errmax .AND. INT(n) <= itrmax) GO TO 30
IF (bound > errmax) ifault = 1
GO TO 999

!       Evaluate the next term of the expansion and then the partial sum

30 u = u * lam / n
v = v + u
t = t * x / (f + two * n)
term = v * t
fn_val = fn_val + term
n = n + one
IF (flag) GO TO 20
GO TO 10

999 CONTINUE

IF(ifault /= 0) THEN
  IF(PRESENT(status)) THEN
    status = ifault
  ELSE
    WRITE(*, *) 'Error in CHI2NC, IFAULT = ', ifault, '. Arguments: x =', x, &
                ', f =',f,', theta=',theta,'. Computed return value: fn_val =',&
                fn_val
  ENDIF
ENDIF

RETURN

END FUNCTION chi2nc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION ppchi2(p, v, g, IFAULT) RESULT(fn_val)

! N.B. Argument IFAULT has been removed.

! This version by Alan Miller
! amiller @ bigpond.net.au
! Latest revision - 27 October 2000

!  Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35

!  To evaluate the percentage points of the chi-squared
!  probability distribution function.

!  p must lie in the range 0.000002 to 0.999998,
!  v must be positive,
!  g must be supplied and should be equal to ln(gamma(v/2.0))

!  Incorporates the suggested changes in AS R85 (vol.40(1), pp.233-5, 1991)
!  which should eliminate the need for the limited range for p above,
!  though these limits have not been removed from the routine.

!  If IFAULT = 4 is returned, the result is probably as accurate as
!  the machine will allow.

!  Auxiliary routines required: PPND = AS 111 (or AS 241) and GAMMAD = AS 239.

IMPLICIT NONE
!INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: p
REAL (dp), INTENT(IN)  :: v
REAL (dp), INTENT(IN)  :: g
INTEGER, INTENT(OUT)   :: IFAULT
REAL (dp)              :: fn_val

!INTERFACE
!  FUNCTION gammad(x, p) RESULT(fn_val)
!    IMPLICIT NONE
!    INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(12, 60)
!    REAL (dp), INTENT(IN) :: x, p
!    REAL (dp)             :: fn_val
!  END FUNCTION gammad
!
!  SUBROUTINE ppnd16 (p, normal_dev, ifault)
!    IMPLICIT NONE
!    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
!    REAL (dp), INTENT(IN)   :: p
!    INTEGER, INTENT(OUT)    :: ifault
!    REAL (dp), INTENT(OUT)  :: normal_dev
!  END SUBROUTINE ppnd16
!END INTERFACE

! Local variables

REAL (dp)  :: a, b, c, p1, p2, q, s1, s2, s3, s4, s5, s6, t, x, xx
INTEGER    :: i, if1

INTEGER, PARAMETER    :: maxit = 150 !maxit = 20
REAL (dp), PARAMETER  :: aa = 0.6931471806_dp, e = 0.5e-06_dp,         &
                         pmin = 0.000000000002_dp, pmax = 0.999999999998_dp, &
                         !pmin = 0.0000000002_dp, pmax = 0.9999999998_dp, &
                         !pmin = 0.000002_dp, pmax = 0.999998_dp,       &
                         zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,   &
                         two = 2.0_dp, three = 3.0_dp, six = 6.0_dp,   &
                         c1 = 0.01_dp, c2 = 0.222222_dp, c3 = 0.32_dp, &
                         c4 = 0.4_dp, c5 = 1.24_dp, c6 = 2.2_dp,       &
                         c7 = 4.67_dp, c8 = 6.66_dp, c9 = 6.73_dp,     &
                         c10 = 13.32_dp, c11 = 60.0_dp, c12 = 70.0_dp, &
                         c13 = 84.0_dp, c14 = 105.0_dp, c15 = 120.0_dp, &
                         c16 = 127.0_dp, c17 = 140.0_dp, c18 = 175.0_dp, &
                         c19 = 210.0_dp, c20 = 252.0_dp, c21 = 264.0_dp, &
                         c22 = 294.0_dp, c23 = 346.0_dp, c24 = 420.0_dp, &
                         c25 = 462.0_dp, c26 = 606.0_dp, c27 = 672.0_dp, &
                         c28 = 707.0_dp, c29 = 735.0_dp, c30 = 889.0_dp, &
                         c31 = 932.0_dp, c32 = 966.0_dp, c33 = 1141.0_dp, &
                         c34 = 1182.0_dp, c35 = 1278.0_dp, c36 = 1740.0_dp, &
                         c37 = 2520.0_dp, c38 = 5040.0_dp

!       Test arguments and initialise

IFAULT = 0

fn_val = -one
IF (p < pmin .OR. p > pmax) THEN
  IFAULT = 1
  !WRITE(*, *) 'Error in PPCHI2: p must be between 0.000002 & 0.999998'
  RETURN
END IF
IF (v <= zero) THEN
  IFAULT = 2
  !WRITE(*, *) 'Error in PPCHI2: Number of deg. of freedom <= 0'
  RETURN
END IF

xx = half * v
c = xx - one

!       Starting approximation for small chi-squared

IF (v < -c5 * LOG(p)) THEN
  fn_val = (p * xx * EXP(g + xx * aa)) ** (one/xx)
  IF (fn_val < e) GO TO 6
  GO TO 4
END IF

!       Starting approximation for v less than or equal to 0.32

IF (v > c3) GO TO 3
fn_val = c4
a = LOG(one-p)

2 q = fn_val
p1 = one + fn_val * (c7+fn_val)
p2 = fn_val * (c9 + fn_val * (c8 + fn_val))
t = -half + (c7 + two * fn_val) / p1 - (c9 + fn_val * (c10 + three * fn_val)) / p2
fn_val = fn_val - (one - EXP(a + g + half * fn_val + c * aa) * p2 / p1) / t
IF (ABS(q / fn_val - one) > c1) GO TO 2
GO TO 4

!       Call to algorithm AS 241 - note that p has been tested above.

3 CALL ppnd16(p, x, if1)

!       Starting approximation using Wilson and Hilferty estimate

p1 = c2 / v
fn_val = v * (x * SQRT(p1) + one - p1) ** 3

!       Starting approximation for p tending to 1

IF (fn_val > c6 * v + six) fn_val = -two * (LOG(one-p) - c * LOG(half * fn_val) + g)

!       Call to algorithm AS 239 and calculation of seven term Taylor series

4 DO i = 1, maxit
  q = fn_val
  p1 = half * fn_val
  p2 = p - gammad(p1, xx)

  t = p2 * EXP(xx * aa + g + p1 - c * LOG(fn_val))
  b = t / fn_val
  a = half * t - b * c
  s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 + c11 * a))))) / c24
  s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 * a)))) / c37
  s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37
  s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 + c36 * a))) / c38
  s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37
  s6 = (c15 + c * (c23 + c16 * c)) / c38
  fn_val = fn_val + t * (one + half * t * s1 - b * c * (s1 - b *   &
           (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))))
  IF (ABS(q / fn_val - one) > e) RETURN
END DO

  IFAULT = 3
  !WRITE(*, *) 'Error in PPCHI2: Max. number of iterations ',maxit,'exceeded'//&
  !            ' for arguments: p =',p,', v =',v,', g =',g

6 RETURN
END FUNCTION ppchi2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


FUNCTION gammad(x, p) RESULT(gamma_prob)

!  ALGORITHM AS239  APPL. STATIST. (1988) VOL. 37, NO. 3

!  Computation of the Incomplete Gamma Integral

!  Auxiliary functions required: ALNORM = algorithm AS66 (included) & LNGAMMA

!  Converted to be compatible with ELF90 by Alan Miller
!  N.B. The return parameter IFAULT has been removed as ELF90 allows only
!  one output parameter from functions.   An error message is issued instead.

! This revision - 15 October 1996

  REAL (dp), INTENT(IN) :: x, p
  REAL (dp)             :: gamma_prob
  
  !     Local variables
  
  REAL (dp) :: pn1, pn2, pn3, pn4, pn5, pn6, tol = 1.d-14, oflo = 1.d+37,  &
               xbig = 1.d+8, arg, c, rn, a, b, one = 1._dp, zero = 0._dp, an, &
               two = 2._dp, elimit = -88._dp, plimit = 1000._dp, three = 3._dp, &
               nine = 9._dp
  
    gamma_prob = zero
    
    !      Check that we have valid values for X and P
    
    IF (p <= zero .OR. x < zero) THEN
      WRITE(*, *)'Error: Function gammad. 1st argument < 0 or 2nd argument <= 0'
      RETURN
    END IF
    IF (x == zero) RETURN
    
    !      Use a normal approximation if P > PLIMIT
    
    IF (p > plimit) THEN
      pn1 = three * SQRT(p) * ((x / p) ** (one / three) + one / (nine * p) - one)
      gamma_prob = alnorm(pn1, .false.)
      RETURN
    END IF
    
    !      If X is extremely large compared to P then set gamma_prob = 1
    
    IF (x > xbig) THEN
      gamma_prob = one
      RETURN
    END IF
    
    IF (x <= one .OR. x < p) THEN
    
    !      Use Pearson's series expansion.
    !      (Note that P is not large enough to force overflow in LNGAMMA)
    
      arg = p * LOG(x) - x - lngamma(p + one)
      c = one
      gamma_prob = one
      a = p
      DO
        a = a + one
        c = c * x / a
        gamma_prob = gamma_prob + c
        IF (c < tol) EXIT
      END DO
      arg = arg + LOG(gamma_prob)
      gamma_prob = zero
      IF (arg >= elimit) gamma_prob = EXP(arg)
    
    ELSE
    
    !      Use a continued fraction expansion
    
      arg = p * LOG(x) - x - lngamma(p)
      a = one - p
      b = a + x + one
      c = zero
      pn1 = one
      pn2 = x
      pn3 = x + one
      pn4 = x * b
      gamma_prob = pn3 / pn4
      DO
        a = a + one
        b = b + two
        c = c + one
        an = a * c
        pn5 = b * pn3 - an * pn1
        pn6 = b * pn4 - an * pn2
        IF (ABS(pn6) > zero) THEN
          rn = pn5 / pn6
          IF (ABS(gamma_prob - rn) <= MIN(tol, tol * rn)) EXIT
          gamma_prob = rn
        END IF
    
        pn1 = pn3
        pn2 = pn4
        pn3 = pn5
        pn4 = pn6
        IF (ABS(pn5) >= oflo) THEN
    
      !      Re-scale terms in continued fraction if terms are large
    
          pn1 = pn1 / oflo
          pn2 = pn2 / oflo
          pn3 = pn3 / oflo
          pn4 = pn4 / oflo
        END IF
      END DO
      arg = arg + LOG(gamma_prob)
      gamma_prob = one
      IF (arg >= elimit) gamma_prob = one - EXP(arg)
    END IF
    
    RETURN
END FUNCTION gammad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION alnorm(x, upper) RESULT(norm_prob)

!  Algorithm AS66 Applied Statistics (1973) vol.22, no.3

!  Evaluates the tail area of the standardised normal curve
!  from x to infinity if upper is .true. or
!  from minus infinity to x if upper is .false.

  REAL (dp), INTENT(IN)  :: x
  LOGICAL, INTENT(IN)    :: upper
  REAL (dp)              :: norm_prob
  
  
  ! Local variables
  REAL (dp) :: zero = 0.0_dp, one = 1.0_dp, half = 0.5_dp
  REAL (dp) :: con = 1.28_dp, z, y, ltone = 7.0_dp, utzero = 18.66_dp
  REAL (dp) :: p = 0.398942280444_dp, q = 0.39990348504_dp,   &
               r = 0.398942280385_dp, a1 = 5.75885480458_dp,  &
               a2 = 2.62433121679_dp, a3 = 5.92885724438_dp,  &
               b1 = -29.8213557807_dp, b2 = 48.6959930692_dp, &
               c1 = -3.8052D-8, c2 = 3.98064794D-4,         &
               c3 = -0.151679116635_dp, c4 = 4.8385912808_dp, &
               c5 = 0.742380924027_dp, c6 = 3.99019417011_dp, &
               d1 = 1.00000615302_dp, d2 = 1.98615381364_dp,  &
               d3 = 5.29330324926_dp, d4 = -15.1508972451_dp, &
               d5 = 30.789933034_dp
  LOGICAL   :: up

    up = upper
    z = x
    IF(z >=  zero) GO TO 10
    up = .NOT. up
    z = -z
    10 IF(z <= ltone .OR. up .AND. z <= utzero) GO TO 20
    norm_prob = zero
    GO TO 40
    20 y = half*z*z
    IF(z > con) GO TO 30
    
    norm_prob = half - z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
    GO TO 40
    30 norm_prob = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
    40 IF(.NOT. up) norm_prob = one - norm_prob
    RETURN

END FUNCTION alnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

  REAL(dp), INTENT(IN) :: z
  REAL(dp)             :: lanczos
  
  ! Local variables
  
  REAL(dp)  :: a(9) = (/ 0.9999999999995183_dp, 676.5203681218835_dp, &
                        -1259.139216722289_dp, 771.3234287757674_dp, &
                        -176.6150291498386_dp, 12.50734324009056_dp, &
                        !-0.1385710331296526_dp, 0.9934937113930748D-05, &
                        ! 0.1659470187408462D-06 /), zero = 0._dp,   &
                        -0.1385710331296526_dp, 0.00009934937113930748_dp, &
                         0.000001659470187408462_dp /), zero = 0._dp,   &
                         one = 1._dp, lnsqrt2pi = 0.9189385332046727_dp, &
                         half = 0.5_dp, sixpt5 = 6.5_dp, seven = 7._dp, tmp
  INTEGER   :: j

    lanczos = zero
    IF (z <= zero) THEN
      WRITE(*, *)'Error: zero or -ve argument for lngamma'
      RETURN
    ENDIF
    
    tmp = z + seven
    DO j = 9, 2, -1
      lanczos = lanczos + a(j)/tmp
      tmp = tmp - one
    END DO
    lanczos = lanczos + a(1)
    lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
    RETURN

END FUNCTION lngamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ppnd16 (p, normal_dev, ifault)

! ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

! Produces the normal deviate Z corresponding to a given lower
! tail area of P; Z is accurate to about 1 part in 10**16.

! The hash sums below are the sums of the mantissas of the
! coefficients.   They are included for use in checking
! transcription.

! This ELF90-compatible version by Alan Miller - 20 August 1996
! N.B. The original algorithm is as a function; this is a subroutine

IMPLICIT NONE

!INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)   :: p
INTEGER, INTENT(OUT)    :: ifault
REAL (dp), INTENT(OUT)  :: normal_dev

! Local variables

REAL (dp) :: zero = 0.d0, one = 1.d0, half = 0.5d0,  &
             split1 = 0.425d0, split2 = 5.d0, const1 = 0.180625d0, &
             const2 = 1.6d0, q, r

! Coefficients for P close to 0.5

REAL (dp) :: a0 = 3.3871328727963666080D0, &
             a1 = 1.3314166789178437745D+2, &
             a2 = 1.9715909503065514427D+3, &
             a3 = 1.3731693765509461125D+4, &
             a4 = 4.5921953931549871457D+4, &
             a5 = 6.7265770927008700853D+4, &
             a6 = 3.3430575583588128105D+4, &
             a7 = 2.5090809287301226727D+3, &
             b1 = 4.2313330701600911252D+1, &
             b2 = 6.8718700749205790830D+2, &
             b3 = 5.3941960214247511077D+3, &
             b4 = 2.1213794301586595867D+4, &
             b5 = 3.9307895800092710610D+4, &
             b6 = 2.8729085735721942674D+4, &
             b7 = 5.2264952788528545610D+3
! HASH SUM AB    55.8831928806149014439

! Coefficients for P not close to 0, 0.5 or 1.

REAL (dp) :: c0 = 1.42343711074968357734D0, &
             c1 = 4.63033784615654529590D0, &
             c2 = 5.76949722146069140550D0, &
             c3 = 3.64784832476320460504D0, &
             c4 = 1.27045825245236838258D0, &
             c5 = 2.41780725177450611770D-1, &
             c6 = 2.27238449892691845833D-2, &
             c7 = 7.74545014278341407640D-4, &
             d1 = 2.05319162663775882187D0, &
             d2 = 1.67638483018380384940D0, &
             d3 = 6.89767334985100004550D-1, &
             d4 = 1.48103976427480074590D-1, &
             d5 = 1.51986665636164571966D-2, &
             d6 = 5.47593808499534494600D-4, &
             d7 = 1.05075007164441684324D-9
! HASH SUM CD    49.33206503301610289036

! Coefficients for P near 0 or 1.

REAL (dp) :: e0 = 6.65790464350110377720D0, &
             e1 = 5.46378491116411436990D0, &
             e2 = 1.78482653991729133580D0, &
             e3 = 2.96560571828504891230D-1, &
             e4 = 2.65321895265761230930D-2, &
             e5 = 1.24266094738807843860D-3, &
             e6 = 2.71155556874348757815D-5, &
             e7 = 2.01033439929228813265D-7, &
             f1 = 5.99832206555887937690D-1, &
             f2 = 1.36929880922735805310D-1, &
             f3 = 1.48753612908506148525D-2, &
             f4 = 7.86869131145613259100D-4, &
             f5 = 1.84631831751005468180D-5, &
             f6 = 1.42151175831644588870D-7, &
             f7 = 2.04426310338993978564D-15
! HASH SUM EF    47.52583317549289671629

ifault = 0
q = p - half
IF (ABS(q) <= split1) THEN
  r = const1 - q * q
  normal_dev = q * (((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0) / &
           (((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + one)
  RETURN
ELSE
  IF (q < zero) THEN
    r = p
  ELSE
    r = one - p
  END IF
  IF (r <= zero) THEN
    ifault = 1
    normal_dev = zero
    RETURN
  END IF
  r = SQRT(-LOG(r))
  IF (r <= split2) THEN
    r = r - const2
    normal_dev = (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0) / &
             (((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + one)
  ELSE
    r = r - split2
    normal_dev = (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0) / &
             (((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + one)
  END IF
  IF (q < zero) normal_dev = - normal_dev
  RETURN
END IF
END SUBROUTINE ppnd16

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE qsort(a, n, t)

!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

!     DP PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.

IMPLICIT NONE

INTEGER, INTENT(IN)     :: n
REAL(dp), INTENT(INOUT) :: a(n)
INTEGER, INTENT(INOUT)  :: t(n)

!     Local Variables

INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL(dp)               :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE qsortInt(a, n, t)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: n
  INTEGER, INTENT(INOUT) :: a(n), t(n)
  REAL(dp)               :: b(n)
  
  b = a
  CALL qsort(b, n, t)
  a = INT(b + 0.5_dp)
  
  RETURN
  
END SUBROUTINE qsortInt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE LogisticReg
