! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random 
! sequence. (Gib Bogle)
!--------------------------------------------------------------------------

MODULE PAR_ZIG_MOD

#ifdef _OPENMP
  USE OMP_LIB
  
  IMPLICIT NONE
  
  PRIVATE
  
  INTEGER,  PARAMETER  ::  DP = SELECTED_REAL_KIND( 12, 60 )
  REAL(DP), PARAMETER  ::  m1 = 2147483648.0_DP, &
                           m2 = 2147483648.0_DP, &
                           half = 0.5_DP
  REAL(DP)             ::  dn0 = 3.442619855899_DP, & 
                           tn0 = 3.442619855899_DP, &
                           vn = 0.00991256303526217_DP, &
                           q, &
                           de0 = 7.697117470131487_DP, &
                           te0 = 7.697117470131487_DP, &
                           ve = 0.003949659822581572_DP
!  INTEGER,  SAVE       ::  iz, jz, jsr = 123456789, kn(0:127),              &
!                           ke(0:255), hz
!  REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
!  LOGICAL,  SAVE       ::  initialized = .FALSE.

  INTEGER, SAVE               :: par_n = 0, par_step
  INTEGER, ALLOCATABLE, SAVE  :: par_jsr(:), par_kn(:,:), par_ke(:,:), par_seed(:) 
  REAL(DP), ALLOCATABLE, SAVE :: par_wn(:,:), par_fn(:,:), par_we(:,:), par_fe(:,:)

  PUBLIC  :: par_zigset, par_shr3, par_uni, par_rnor, par_rexp, par_seed


  CONTAINS


SUBROUTINE par_zigset(npar, par_jsrseed, grainsize)
!! None: Indices of par_jsrseed need to start from 0:
  INTEGER, INTENT(IN) :: npar, par_jsrseed(0:), grainsize
  INTEGER  :: i, kpar
  REAL(DP) :: dn, tn, de, te

  par_n = npar
  par_step = grainsize

  ! First we need to allocate all the non-volatile arrays with the size npar
  IF(ALLOCATED(par_jsr))  DEALLOCATE(par_jsr)
  IF(ALLOCATED(par_kn))   DEALLOCATE(par_kn)
  IF(ALLOCATED(par_ke))   DEALLOCATE(par_ke)
  IF(ALLOCATED(par_wn))   DEALLOCATE(par_wn)
  IF(ALLOCATED(par_fn))   DEALLOCATE(par_fn)
  IF(ALLOCATED(par_we))   DEALLOCATE(par_we)
  IF(ALLOCATED(par_fe))   DEALLOCATE(par_fe)
  IF(ALLOCATED(par_seed)) DEALLOCATE(par_seed)
  ALLOCATE(par_jsr(0:npar*par_step))
  ALLOCATE(par_kn(0:127,0:npar-1))
  ALLOCATE(par_ke(0:255,0:npar-1))
  ALLOCATE(par_wn(0:127,0:npar-1))
  ALLOCATE(par_fn(0:127,0:npar-1))
  ALLOCATE(par_we(0:255,0:npar-1))
  ALLOCATE(par_fe(0:255,0:npar-1))
  ALLOCATE(par_seed(npar))
  
  !! Store the seed globally
  par_seed = par_jsrseed

  ! Now treat each instance separately
  DO kpar = 0,npar-1
  
    !  Set the seed
    par_jsr(kpar*par_step) = par_jsrseed(kpar)

    !  Tables for RNOR
    dn = dn0
    tn = tn0
    q = vn*EXP(half*dn*dn)
    par_kn(0,kpar) = INT((dn/q)*m1)
    par_kn(1,kpar) = 0
    par_wn(0,kpar) = q/m1
    par_wn(127,kpar) = dn/m1
    par_fn(0,kpar) = 1.0_DP
    par_fn(127,kpar) = EXP( -half*dn*dn )
    DO  i = 126, 1, -1
      dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )    ! dn
      par_kn(i+1,kpar) = INT((dn/tn)*m1)
      tn = dn                            ! tn
      par_fn(i,kpar) = EXP(-half*dn*dn)
      par_wn(i,kpar) = dn/m1
    END DO

    !  Tables for REXP
    de = de0
    te = te0
    q = ve*EXP( de )
    par_ke(0,kpar) = INT((de/q)*m2)
    par_ke(1,kpar) = 0
    par_we(0,kpar) = q/m2
    par_we(255,kpar) = de/m2
    par_fe(0,kpar) = 1.0_DP
    par_fe(255,kpar) = EXP( -de )
    DO  i = 254, 1, -1
      de = -LOG( ve/de + EXP( -de ) )                ! de
      par_ke(i+1,kpar) = INT(m2 * (de/te))
      te = de                            ! te
      par_fe(i,kpar) = EXP( -de )
      par_we(i,kpar) = de/m2
    END DO
    
  ENDDO
  
  RETURN
  
END SUBROUTINE par_zigset

!  Generate random 32-bit integers
FUNCTION par_shr3(kpar) RESULT( ival )
!! Note by Jakub Pecanka: This function relies on integer overlow when shifting
!! the bits by ISHFT, which can cause runtime error when compiled with checks
!! for integer overflow and for example if compiled by Gfortran 4.5.0 with
!! --ffpe-trap=overflow it can be very hard to locate a problem. Therefore, 
!! I modified this function by replacing the line
!!        ival = jz + jsr
!! with
!!        ival = INT(INT(jz, 8) + INT(jsr, 8), 4)
!! which achieves the same thing as the integer overflow.
 
  INTEGER :: ival, kpar
  INTEGER :: jz, jsr
  
  jsr = par_jsr(kpar*par_step)
  jz = jsr
  jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
  jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
  jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
  par_jsr(kpar*par_step) = jsr

  ival = INT(INT(jz, 8) + INT(jsr, 8), 4)
  
  RETURN
END FUNCTION par_shr3



!  Generate uniformly distributed random numbers, sequence kpar
FUNCTION par_uni(kpar) RESULT( fn_val )
  INTEGER  :: kpar
  REAL(DP) :: fn_val

  IF(kpar >= par_n) THEN
    WRITE(*,*) 'Thread number exceeds initialized max: ',kpar+1,">",par_n
    STOP
  ENDIF
  fn_val = half + 0.2328306e-9_DP * par_shr3(kpar)
  RETURN
END FUNCTION par_uni



!  Generate random normals, sequence kpar
FUNCTION par_rnor(kpar) RESULT( fn_val )
  REAL(DP)             ::  fn_val
  INTEGER :: kpar

  REAL(DP), PARAMETER  ::  r = 3.442620_DP
  REAL(DP)             ::  x, y
  INTEGER :: iz, hz

!   IF( .NOT. initialized ) CALL zigset( jsr )
  if (kpar >= par_n) then
    write(*,*) 'thread number exceeds initialized max: ',kpar,par_n
    stop
  endif

  hz = par_shr3(kpar)
  iz = IAND( hz, 127 )
  IF( ABS( hz ) < par_kn(iz,kpar) ) THEN
    fn_val = hz * par_wn(iz,kpar)
  ELSE
    DO
       IF( iz == 0 ) THEN
          DO
             x = -0.2904764_DP* LOG( par_uni(kpar) )
             y = -LOG( par_uni(kpar) )
             IF( y+y >= x*x ) EXIT
          END DO
          fn_val = r+x
          IF( hz <= 0 ) fn_val = -fn_val
          RETURN
       END IF
       x = hz * par_wn(iz,kpar)
       IF( par_fn(iz,kpar) + par_uni(kpar)*(par_fn(iz-1,kpar)-par_fn(iz,kpar)) < EXP(-half*x*x) ) THEN
          fn_val = x
          RETURN
       END IF
       hz = par_shr3(kpar)
       iz = IAND( hz, 127 )
       IF( ABS( hz ) < par_kn(iz,kpar) ) THEN
          fn_val = hz * par_wn(iz,kpar)
          RETURN
       END IF
    END DO
  END IF
  RETURN
END FUNCTION par_rnor



!  Generate random exponentials, sequence kpar
FUNCTION par_rexp(kpar) RESULT( fn_val )
   REAL(DP) ::  fn_val
   INTEGER  :: kpar

   REAL(DP) ::  x
   INTEGER  :: iz, jz

!   IF( .NOT. initialized ) CALL Zigset( jsr )
  if (kpar >= par_n) then
    write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
    stop
  endif
   jz = par_shr3(kpar)
   iz = IAND( jz, 255 )
   IF( ABS( jz ) < par_ke(iz,kpar) ) THEN
      fn_val = ABS(jz) * par_we(iz,kpar)
      RETURN
   END IF
   DO
      IF( iz == 0 ) THEN
         fn_val = 7.69711 - LOG(par_uni(kpar) )
         RETURN
      END IF
      x = ABS( jz ) * par_we(iz,kpar)
      IF( par_fe(iz,kpar) + par_uni(kpar)*(par_fe(iz-1,kpar) - par_fe(iz,kpar)) < EXP( -x ) ) THEN
         fn_val = x
         RETURN
      END IF
      jz = par_shr3(kpar)
      iz = IAND( jz, 255 )
      IF( ABS( jz ) < par_ke(iz,kpar) ) THEN
         fn_val = ABS( jz ) * par_we(iz,kpar)
         RETURN
      END IF
   END DO
   RETURN
END FUNCTION par_rexp

#endif

END MODULE PAR_ZIG_MOD
