MODULE EPI_SIMUL

  !USE EPI_UTILS
  USE EPI_MATH

  IMPLICIT NONE

  CHARACTER(mstl)             :: model_names(25) = &
                                  (/"A", "B", "C", "D", "E", "F", "G", "H", &
                                    "I", "J", "K", "L", "M", "N", "O", "P", &
                                    "Q", "R" , "S" , "T", "U", "V", "W", "X", "Z" /)   
                                 ! Defines names of interaction models, now
                                 ! these are simply letters but they could
                                 ! be given more meaningful values (up to
                                 ! length 'mstl', which is currently 50)

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GetModelName(model) RESULT(name)
!! This function returns the name of a model specified by a number. 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: model
  CHARACTER(mstl)     :: name
  
  IF(model>0 .AND. model <= SIZE(model_names)) THEN
    name = model_names(model)
  ELSE  
    name = "unknown model"
  ENDIF
  
  RETURN
    
END FUNCTION GetModelName

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION GetModelNumber(model) RESULT(number)
!! This function return the number of a model specified by a number. 
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: model
  INTEGER                  :: number, i
  
  !! Find the model name among defined model names
  number = -1
  DO i=1,SIZE(model_names)
    IF(model == model_names(i)) number = i
  ENDDO
  
  !! Throw an error if model name not found
  IF(number == -1) &
    CALL PrntE("Unknown model '"//TRIM(model)//"'! Execution stopped.", &
                    Q=.TRUE., premature=.FALSE., skip1=1)
  
  RETURN
    
END FUNCTION GetModelNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GenerateCaCo(nca, nco, frA, frB, frC, fixed_freq, OR, OR1, &
                        OR2, LD, coef, prev, nmarkers, nLDpairs, HWE, binary, &
                        X, GENca, GENco, prev_actual, nindiv_total)   
!! Generate genotypes for the given number of GENca and GENco. 
!! OnlyMinorAlleleAtRisk determines whether it is always the minor allele that
!! is at risk or whether it is randomly selected between minor and major
!! Number of rows of 'GENca' and 'GENco' must be 'nca' and 'nco' and
!! the number of columns must be 4*(2+nmarkers+nLDpairs) for both and the  
!! 'length of character' for both should be 1. 
  IMPLICIT NONE
  LOGICAL, PARAMETER                  :: OnlyMinorAlleleAtRisk = .TRUE., &
                                         def_HWE = .TRUE.
  LOGICAL, INTENT(IN)                 :: fixed_freq
  REAL(dpp), INTENT(IN)               :: frA, frB, frC,&
                                         prev, OR, OR1, OR2, LD, coef(:,:)
  INTEGER, INTENT(IN)                 :: nca, nco, nmarkers, nLDpairs
  LOGICAL, INTENT(IN), OPTIONAL       :: binary, HWE

  INTEGER(iks), INTENT(OUT), OPTIONAL :: X(:,:)
  CHARACTER(*), INTENT(OUT), OPTIONAL :: GENca(:,:), GENco(:,:)
  REAL(dpp), INTENT(OUT), OPTIONAL    :: prev_actual
  INTEGER, INTENT(OUT), OPTIONAL      :: nindiv_total
  
  INTEGER(iks)                        :: ng1, ng2, Y(nca+nco,2)
  INTEGER                             :: nca_kept, nco_kept, i, j, m, n, &
                                         nindiv, nca_total
  CHARACTER(2)                        :: lg1, lg2
  REAL(dpp)                           :: randval, affprob, minfr1, &
                                         minfr2, minfr3, & 
                                         fr1, fr2, fr3
  LOGICAL                             :: binary1, HWE1

  !! Determine whether the HWE assumption is to be used or not
  HWE1 = def_HWE
  IF(PRESENT(HWE)) HWE1 = HWE
  
  IF(PRESENT(prev_actual)) prev_actual = zero
  IF(PRESENT(nindiv_total)) nindiv_total = 0
  
  !! If no output variables present, no output wanted
  IF(.NOT.(PRESENT(GENca) .OR. PRESENT(GENco) .OR. PRESENT(X))) RETURN
  
  !! Based on presence and value of binary determine whether the output
  !! will be binary or plain
  binary1 = .FALSE.
  IF(PRESENT(binary)) binary1 = binary

  !! Read in constants
  minfr1 = MAX(frA, zero)
  minfr2 = MAX(frB, zero)
  minfr3 = MAX(frC, zero)

  !! Initialize counters
  nindiv = 0        
  nca_kept = 0
  nco_kept = 0
  nca_total = 0
  
  !! You can use either fixed allele frequencies or a lower cutoff to a uniform
  !! minor allele frequency distribution 
  IF(fixed_freq) THEN
    !! Use fixed allele freqs
    fr1 = minfr1
    fr2 = minfr2
    fr3 = minfr3
  ELSE
    !! Use uniform allele frequency distribution with a lower cutoff...
    CALL RandNumber(fr1)
    CALL RandNumber(fr2)
    fr1 = fr1 * (half-minfr1) + minfr1
    fr2 = fr2 * (half-minfr2) + minfr2
  ENDIF

  !! Repeat the process of generating random genotypes and deciding if an 
  !! individual is case or control until the required sample size is reached
  DO WHILE (nca_kept < nca .OR. nco_kept < nco)
  
    !! Increase the counter of total number of "observed" individuals
    nindiv = nindiv + 1

    !! If OnlyMinorAlleleAtRisk is .FALSE., then randomly decide whether MINOR 
    !! or MAJOR alleles are at risk for the current pair
    IF(.NOT.OnlyMinorAlleleAtRisk) THEN
      CALL RandNumber(randval)
      IF(randval >= half) fr1 = one-fr1
      CALL RandNumber(randval) 
      IF(randval >= half) fr2 = one-fr2
    ENDIF
    
    !! Generate the genotypes at the first TWO loci assuming Hardy-Weinberg 
    !! equilibrium but allowing for a user specified amount of LD
    CALL GetGenotypes(ng1, ng2, LD, fr1, fr2, HWE1, lg1, lg2)
                           
    !! Get the probability of being a case for given genotypes (ng1, ng2)
    CALL GetProbs(affprob, coef, ng1, ng2, prev, OR, OR1, OR2, fr1, fr2)

    !! Decide if an individual is case/control based on the two-locus genotype
    !! If randval is SMALLER than affprob, it is a CASE
    CALL RandNumber(randval)
    IF(randval <= affprob) THEN 

      !! Increase the counter of total number of "observed" GENca
      nca_total = nca_total + 1
    
      !! If all required CASES already simulated, skip to the next cycle
      IF(nca_kept >= nca) CYCLE
      
      !! Otherwise, save the genotype as a CASE
      nca_kept = nca_kept+1
      
      !! If GENca present, save the genotype
      IF(PRESENT(GENca)) THEN
        GENca(nca_kept,1) = lg1(1:1)
        GENca(nca_kept,2) = lg1(2:2)
        GENca(nca_kept,3) = lg2(1:1)
        GENca(nca_kept,4) = lg2(2:2)
      ENDIF
      !! If X present, save the numerical representation of the genotype 
      IF(PRESENT(X)) THEN
        Y(nca_kept,1) = ng1
        Y(nca_kept,2) = ng2
      ENDIF 

    !! Otherwise, if randval is BIGGER than affprob, it is a CONTROL
    ELSE

      !! If all required GENca CONTROLS simulated, skip to the next cycle
      IF(nco_kept >= nco) CYCLE
      !! Otherwise, save the genotype as a CONTROL
      nco_kept = nco_kept + 1
      !! If GENco present, save the genotype
      IF(PRESENT(GENco)) THEN
        GENco(nco_kept,1) = lg1(1:1)
        GENco(nco_kept,2) = lg1(2:2)
        GENco(nco_kept,3) = lg2(1:1)
        GENco(nco_kept,4) = lg2(2:2)
      ENDIF
      !! If X present, save the numerical representation of the genotype 
      IF(PRESENT(X)) THEN
        Y(nco_kept+nca,1) = ng1
        Y(nco_kept+nca,2) = ng2
      ENDIF 

    ENDIF

  ENDDO
  
  !! Save/convert Y into output variable X
  IF(PRESENT(X)) THEN
    IF(binary1) THEN
      CALL Ped2Bed(Y, X(:,1:2), .TRUE., .TRUE.)
    ELSE
      X(:,1:2) = Y
    ENDIF
  ENDIF
  
  !! Number of pairs
  n = nLDpairs + 1    !! +1 is for the "causal" pair 

!!***********************!!
!! START PARALLEL REGION !!
!!***********************!!
!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP FIRSTPRIVATE(nca, nco, nmarkers, fixed_freq, LD, n) &
!!$OMP FIRSTPRIVATE(minfr1, minfr2, fr1, fr2) & 
!!$OMP FIRSTPRIVATE(fr3) &
!!$OMP PRIVATE(i, j, m, lg1, lg2, ng1, ng2) &
!!$OMP PRIVATE(randval) &
!!$OMP SHARED(GENca, GENco, X)

!!$OMP DO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Generate neutral marker genotypes for the SNP pairs that should be in LD 
  !! both for GENco and GENca
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO m = 2, n
    DO i = 1, nca+nco
      !! Get allele frequency if NOT FIXED. Use unif. distr. with a lower cutoff
      IF(.NOT.fixed_freq) THEN
        CALL RandNumber(fr1)
        CALL RandNumber(fr2)
        fr1 = fr1*(half-minfr1) + minfr1
        fr2 = fr2*(half-minfr2) + minfr2
      ENDIF

      !! Get the genotypes
      CALL GetGenotypes(ng1, ng2, LD, fr1, fr2, HWE1, lg1, lg2)

      !! Append the genotypes, first cases, then controls
      IF(i <= nca) THEN
!!$OMP CRITICAL
        !! If GENca present, save the genotype
        IF(PRESENT(GENca)) THEN
          GENca(i,4*(m-1)+1) = lg1(1:1)
          GENca(i,4*(m-1)+2) = lg1(2:2)
          GENca(i,4*(m-1)+3) = lg2(1:1)
          GENca(i,4*(m-1)+4) = lg2(2:2)
        ENDIF
        !! If X present, save the numerical representation of the genotype 
        IF(PRESENT(X)) THEN
          Y(i,1) = ng1
          Y(i,2) = ng2
        ENDIF 
!!$OMP END CRITICAL
      ENDIF

      IF(i > nca) THEN
         j = i-nca
!!$OMP CRITICAL
        !! If GENco present, save the genotype
        IF(PRESENT(GENco)) THEN
          GENco(j,4*(m-1)+1) = lg1(1:1)
          GENco(j,4*(m-1)+2) = lg1(2:2)
          GENco(j,4*(m-1)+3) = lg2(1:1)
          GENco(j,4*(m-1)+4) = lg2(2:2)
        ENDIF
        !! If X present, save the numerical representation of the genotype 
        IF(PRESENT(X)) THEN
          Y(i,1) = ng1
          Y(i,2) = ng2
        ENDIF 
!!$OMP END CRITICAL
      ENDIF

    ENDDO !! i: loop over individuals
    
    !! Save Y into output variable X
    IF(PRESENT(X)) THEN
      IF(binary1) THEN
        CALL Ped2Bed(Y, X(:,2*(m-1)+1:2*(m-1)+2), .TRUE., .TRUE.)
      ELSE
        X(:,2*(m-1)+1:2*(m-1)+2) = Y
      ENDIF
    ENDIF

  ENDDO !! m: loop over LD pairs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Generate neutral marker genotypes for the SNP pairs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO m = 1, nmarkers
    DO i = 1, nca+nco
      !! Generate neutral marker genotypes for the same case and control 
      !! individuals generated in the previous step
      IF(.NOT.fixed_freq) THEN
        CALL RandNumber(fr3)
        fr3 = fr3*(half-minfr3) + minfr3
      ENDIF
      
      !! Get random number to determine what the genotypes are
      CALL RandNumber(randval)

      !! First cases
       IF(i <= nca) THEN
!!$OMP CRITICAL
          !! If GENca present, save the genotype
        IF(PRESENT(GENca)) THEN
          IF(randval <= fr3**2) THEN
            GENca(i,4*n+2*(m-1)+1) = "G"
            GENca(i,4*n+2*(m-1)+2) = "G"
          ELSEIF(randval > one-(one-fr3)**2) THEN
            GENca(i,4*n+2*(m-1)+1) = "A"
            GENca(i,4*n+2*(m-1)+2) = "A"
          ELSE
            GENca(i,4*n+2*(m-1)+1) = "A"
            GENca(i,4*n+2*(m-1)+2) = "G"
          ENDIF
        ENDIF
        !! If X present, save the numerical representation of the genotype 
        IF(PRESENT(X)) THEN
          IF(randval <= fr3**2) THEN
            Y(i,1) = 2_iks
          ELSEIF(randval > one-(one-fr3)**2) THEN
            Y(i,1) = 0_iks
          ELSE
            Y(i,1) = 1_iks
          ENDIF
        ENDIF
!!$OMP END CRITICAL
      ENDIF
      
      !! Then controls      
       IF(i > nca) THEN
         j = i-nca
!!$OMP CRITICAL
        !! If GENco present, save the genotype
        IF(PRESENT(GENco)) THEN
          IF(randval <= fr3**2) THEN
            GENco(j,4*n+2*(m-1)+1) = "G"
            GENco(j,4*n+2*(m-1)+2) = "G"
          ELSEIF(randval > one-(one-fr3)**2) THEN
            GENco(j,4*n+2*(m-1)+1) = "A"
            GENco(j,4*n+2*(m-1)+2) = "A"
          ELSE
            GENco(j,4*n+2*(m-1)+1) = "A"
            GENco(j,4*n+2*(m-1)+2) = "G"
          ENDIF
        ENDIF
        !! If X present, save the numerical representation of the genotype 
        IF(PRESENT(X)) THEN
          IF(randval <= fr3**2) THEN
            Y(i,1) = 2_iks
          ELSEIF(randval > one-(one-fr3)**2) THEN
            Y(i,1) = 0_iks
          ELSE
            Y(i,1) = 1_iks
          ENDIF
        ENDIF
!!$OMP END CRITICAL
      ENDIF

    ENDDO !! i: loop over individuals
  
    !! Save Y into output variable X
    IF(PRESENT(X)) THEN
      IF(binary1) THEN
        CALL Ped2Bed(Y(:,1:1), X(:,2*n+m:2*n+m), .TRUE., .TRUE.)
      ELSE
        X(:,2*n+m) = Y(:,1)
      ENDIF
    ENDIF

  ENDDO !! m: loop over neutral pairs
    
!!$OMP END DO
!!$OMP END PARALLEL
!!***********************!!
!!  END PARALLEL REGION  !!
!!***********************!!

  IF(PRESENT(prev_actual)) prev_actual = REAL(nca_total, dpp) / nindiv
  IF(PRESENT(nindiv_total)) nindiv_total = nindiv

  RETURN
    
END SUBROUTINE GenerateCaCo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetGenotypes(ng1, ng2, LD, fr1, fr2, HWE, lg1, lg2)
!! Returns genotypes for two loci generated with or without the assumption of
!! Hardy-Weinberg equilibrium (HWE). Returns the numerical genotypes in ng1 and 
!! ng2 and the "letter" representation in lg1 and lg2. The numerical
!! representation is obtained by counting the minor allele at each locus, which
!! in our notation is the allele G.
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)               :: LD, fr1, fr2
  INTEGER(iks), INTENT(OUT)           :: ng1, ng2
  LOGICAL, INTENT(IN)                 :: HWE
  CHARACTER(2), INTENT(OUT), OPTIONAL :: lg1, lg2
  REAL(dpp)                           :: P0, P11, P12, P21, T1, T2, T3, v1, v2

  !! Get random number that decide genotypes
  CALL RandNumber(v1)
  CALL RandNumber(v2)
  
  !! Generate random genotypes WHILE assuming HARDY-WEINBERG EQUILIBRIUM
  IF(HWE) THEN

    ng1 = 1_iks
    IF(v1<=fr1**2)          ng1 = 2_iks
    IF(v1>one-(one-fr1)**2) ng1 = 0_iks

    ng2 = 1_iks
    IF(v2<=fr2**2)          ng2 = 2_iks
    IF(v2>one-(one-fr2)**2) ng2 = 0_iks

  !! Generate random genotypes WITHOUT assuming HARDY-WEINBERG EQUILIBRIUM
  ELSE
  
    !! Calculate haplotype frequencies based on the user-specified squared  
    !! correlation coefficient r^2 (=|LD|) and the allele frequencies...
    P0  = SQRT( ABS(LD) * fr1 * (one-fr1) * fr2 * (one-fr2) )
    P11 =  P0 + fr1 * fr2
    P12 = -P0 + fr1 * (one-fr2)
    P21 = -P0 + (one-fr1) * fr2
    
    T1 = P11
    T2 = T1 + P12
    T3 = T2 + P21
    
    !! This looks complicated but you need to draw two random haplotypes 
    !! (basically the parental haplotypes) that determine the two-locus genotype
    !! Default genotype is 1-1, which would be when 
    !!      v1<=T1 .AND. T3<v2                                 or
    !!      T1 < v1 .AND. v1<=T2 .AND. T2 < v2 .AND. v2<=T3    or
    !!      T2 < v1 .AND. v1<=T3 .AND. T1 < v2 .AND. v2<=T2    or
    !!      T3 < v1 .AND. v2<=T1
    ng1 = 1_iks; ng2 = 1_iks                                  ! AGAG

    !! Otherwise the genotypes change to: 
    IF(v1<=T1 .AND. v2<=T1) THEN
      ng1 = 2_iks; ng2 = 2_iks;                                 ! GGGG
    ENDIF
    IF(v1<=T1 .AND. T1<v2 .AND. v2<=T2) ng1 = 2_iks              ! GGAG
    IF(v1<=T1 .AND. T2<v2 .AND. v2<=T3) ng2 = 2_iks              ! AGGG
    IF(T1 < v1 .AND. v1<=T2 .AND. v2<=T1) ng1 = 2_iks            ! GGAG
    IF(T1 < v1 .AND. v1<=T2 .AND. T1 < v2 .AND. v2<=T2) THEN
      ng1 = 2_iks; ng2 = 0_iks                                 ! GGAA
    ENDIF
    IF(T1 < v1 .AND. v1<=T2 .AND. T3 < v2) ng2 = 0_iks          ! AGAA
    IF(T2 < v1 .AND. v1<=T3 .AND. v2<=T1)  ng2 = 2_iks            ! AGGG
    IF(T2 < v1 .AND. v1<=T3 .AND. T2 < v2 .AND. v2<=T3) THEN
      ng1 = 0_iks; ng2 = 2_iks                                 ! AAGG
    ENDIF
    IF(T2 < v1 .AND. v1<=T3 .AND. T3 < v2) ng1 = 0_iks          ! AAAG
    IF(T3 < v1 .AND. T1 < v2 .AND. v2<=T2) ng2 = 0_iks          ! AGAA
    IF(T3 < v1 .AND. T2 < v2 .AND. v2<=T3) ng1 = 0_iks          ! AAAG
    IF(T3 < v1 .AND. T3 < v2) THEN
      ng1 = 0_iks; ng2 = 0_iks                                ! AAAA
    ENDIF
    
    !! Reverse one of the genotypes if correlation should be negative
    IF(LD<zero) THEN
      ng2 = 2_iks - ng2
    ENDIF

  ENDIF
  
  !! "Letter" representation of genotypes
  IF(PRESENT(lg1)) THEN
    SELECT CASE(ng1)
      CASE(0_iks); lg1 = "AA" 
      CASE(1_iks);  lg1 = "AG" 
      CASE(2_iks);  lg1 = "GG"
    END SELECT
  ENDIF 
  IF(PRESENT(lg2)) THEN
    SELECT CASE(ng2)
      CASE(0_iks); lg2 = "AA" 
      CASE(1_iks);  lg2 = "AG" 
      CASE(2_iks);  lg2 = "GG"
    END SELECT 
  ENDIF 
  
  RETURN
    
END SUBROUTINE GetGenotypes

!!! ************************************************************************ !!!
!!! OLDER VERSION OF GetGenotypes follows ...                                !!!
!!! ************************************************************************ !!!
! SUBROUTINE GetGenotypes(g1, g2, LD, fr1, fr2, ng1, &
!                         ng2, HWE)
! !! Returns genotypes for two loci generated with or without the assumption of
! !! Hardy-Weinberg equilibrium. Returns the genotypes in g1 and g2 and
! !! the numerical representation in ng1 and ng2. The numerical
! !! representation is obtained by counting the minor allele at each locus, which
! !! in our notation is the allele G.
!   IMPLICIT NONE
!   CHARACTER(2), INTENT(OUT)          :: g1, g2
!   REAL(dpp), INTENT(IN)              :: LD, fr1, fr2
!   INTEGER(iks), INTENT(OUT), OPTIONAL :: ng1, ng2
!   LOGICAL, INTENT(IN)                :: HWE
!   REAL(dpp)                          :: P0, P11, P12, P21, T1, T2, T3, val1, &
!                                         val2
!   CHARACTER(4)                       :: both
! 
!   !! Get random number that decide genotypes
!   CALL RandNumber(val1)
!   CALL RandNumber(val2)
! 
!   !! Generate random genotypes WHILE assuming HARDY-WEINBERG EQUILIBRIUM
!   IF(HWE) THEN
! 
!     IF(val1<=fr1**2)                                    g1 = "GG" 
!     IF(val1>one-(one-fr1)**2)                           g1 = "AA"   
!     IF(val1>fr1**2 .AND. val1<one-(one-fr1)**2) g1 = "AG"
! 
!     IF(val2<=fr2**2)                                    g2 = "GG" 
!     IF(val2>one-(one-fr2)**2)                           g2 = "AA"   
!     IF(val2>fr2**2 .AND. val2<one-(one-fr2)**2) g2 = "AG"
! 
!   !! Generate random genotypes WITHOUT assuming HARDY-WEINBERG EQUILIBRIUM
!   ELSE
!     !! Calculate haplotype frequencies based on the user-specified r-squared and 
!     !! the allele frequencies...
!     P0  = SQRT( LD*fr1*(one-fr1)*fr2*(one-fr2) )
!     P11 =  P0 + fr1 * fr2
!     P12 = -P0 + fr1 * (one-fr2)
!     P21 = -P0 + (one-fr1) * fr2
!     
!     T1 = P11
!     T2 = T1 + P12
!     T3 = T2 + P21
!     
!     !! This looks complicated but you need to draw two random haplotypes 
!     !! (basically the parental haplotypes) that determine the two-locus 
!     !! genotype
!     IF(val1<=T1 .AND. val2<=T1)                                 both = "GGGG"
!     IF(val1<=T1 .AND. T1 < val2 .AND. val2<=T2)                 both = "GGAG"
!     IF(val1<=T1 .AND. T2 < val2 .AND. val2<=T3)                 both = "AGGG"
!     IF(val1<=T1 .AND. T3 < val2)                                both = "AGAG"
!     
!     IF(T1 < val1 .AND. val1<=T2 .AND. val2<=T1)                 both = "GGAG"
!     IF(T1 < val1 .AND. val1<=T2 .AND. T1 < val2 .AND. val2<=T2) both = "GGAA"
!     IF(T1 < val1 .AND. val1<=T2 .AND. T2 < val2 .AND. val2<=T3) both = "AGAG"
!     IF(T1 < val1 .AND. val1<=T2 .AND. T3 < val2)                both = "AGAA"
!     
!     IF(T2 < val1 .AND. val1<=T3 .AND. val2<=T1)                 both = "AGGG"
!     IF(T2 < val1 .AND. val1<=T3 .AND. T1 < val2 .AND. val2<=T2) both = "AGAG"
!     IF(T2 < val1 .AND. val1<=T3 .AND. T2 < val2 .AND. val2<=T3) both = "AAGG"
!     IF(T2 < val1 .AND. val1<=T3 .AND. T3 < val2)                both = "AAAG"
!   
!     IF(T3 < val1 .AND. val2<=T1)                                both = "AGAG"
!     IF(T3 < val1 .AND. T1 < val2 .AND. val2<=T2)                both = "AGAA"
!     IF(T3 < val1 .AND. T2 < val2 .AND. val2<=T3)                both = "AAAG"
!     IF(T3 < val1 .AND. T3 < val2)                               both = "AAAA"
! 
!     !! Genetic representation of genotypes
!     g1 = both(1:2)
!     g2 = both(3:4)
!   ENDIF
!   
!   !! Numeric representation of genotypes
!   IF(PRESENT(ng1)) THEN
!     ng1 = 0_iks 
!     IF(g1(1:1) == "G") ng1 = ng1 + 1_iks 
!     IF(g1(2:2) == "G") ng1 = ng1 + 1_iks
!   ENDIF 
!   IF(PRESENT(ng2)) THEN
!     ng2 = 0_iks 
!     IF(g2(1:1) == "G") ng2 = ng2 + 1_iks 
!     IF(g2(2:2) == "G") ng2 = ng2 + 1_iks
!   ENDIF
!   
!   RETURN
!     
! END SUBROUTINE GetGenotypes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetInteractionTable(coef, model)
!! Returns in affprob the probability of being a case for given parameters
!! NOTE: A new model is added by adding an aditional case below and changing 
!!       'model_names' at the top of this module  
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT) :: coef(3,3)
  INTEGER, INTENT(IN)    :: model

  SELECT CASE (model)
  
  !! Multiplicative risk at both loci
  CASE(1) !! Model "A"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, one, two, &
                                 zero, two, four /), (/3,3/) ) )
  !! Only AG-GG is at risk
  CASE(2) !! Model "B"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, zero, two, &
                                 zero, zero, zero /), (/3,3/) ) )
  !! Only the single heterozygous are at risk
  CASE(3) !! Model "C"
    coef = TRANSPOSE(RESHAPE( (/ zero, one, zero, &
                                 one, zero, one, &
                                 zero, one, zero /), (/3,3/) ) )
  !! Only AG-GG is at high risk
  CASE(4) !! Model "D"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, one, &
                                 four, zero, two, &
                                 zero, zero, one /), (/3,3/) ) )
  !! Only AA-AA, AA-GG, GG-AA, GG-GG are at risk
  CASE(5) !! Model "E"
    coef = TRANSPOSE(RESHAPE( (/ one, zero, one, &
                                 zero, zero, zero, &
                                 one, zero, one /), (/3,3/) ) )
  !! Only AG-?? are at risk
  CASE(6) !! Model "F"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 one, one, one, &
                                 zero, zero, zero /), (/3,3/) ) )
  !! AA-AA, AA-GG, GG-AG are at risk, AA-AG, GG-AA, GG-GG are at high risk
  CASE(7) !! Model "G"
    coef = TRANSPOSE(RESHAPE( (/ one, two, one, &
                                 zero, zero, zero, &
                                 two, one, two /), (/3,3/) ) )
  !! Only AA-AG, AG-AG, AG-GG, GG-AG are at risk
  CASE(8) !! Model "H"
    coef = TRANSPOSE(RESHAPE( (/ zero, one, zero, &
                                 zero, one, three, &
                                 zero, one, zero /), (/3,3/) ) )
  !! Dominant risk at both loci
  CASE(9) !! Model "I"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, one, one, &
                                 zero, one, one /), (/3,3/) ) )
  !! Only the double heterozygous are at (strong) risk
  CASE(10) !! Model "J"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, two, zero, &
                                 zero, zero, three /), (/3,3/) ) )
  CASE(11) !! Model "K"
    coef = TRANSPOSE(RESHAPE( (/ one, zero, zero, &
                                 zero, two, zero, &
                                 zero, zero, two /), (/3,3/) ) )
  CASE(12) !! Model "L"
    coef = TRANSPOSE(RESHAPE( (/ four, zero, zero, &
                                 zero, one, zero, &
                                 zero, zero, three /), (/3,3/) ) )
  CASE(13) !! Model "M"
    coef = TRANSPOSE(RESHAPE( (/ four, zero, one, &
                                 zero, one, zero, &
                                 zero, two, three /), (/3,3/) ) )
  CASE(14) !! Model "N"
    coef = TRANSPOSE(RESHAPE( (/ one, zero, one, &
                                 one, two, one, &
                                 zero, five, three /), (/3,3/) ) )
  CASE(15) !! Model "O"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, one, two, &
                                 zero, two, three /), (/3,3/) ) )
  CASE(16) !! Model "P"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, zero, one, &
                                 zero, one, one /), (/3,3/) ) )
  CASE(17) !! Model "Q"
    coef = TRANSPOSE(RESHAPE( (/ one, zero, three, &
                                 zero, zero, zero, &
                                 two, zero, four /), (/3,3/) ) )
  !! Only AG-?? are at risk
  CASE(18) !! Model "R"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 one, three, one, &
                                 zero, zero, zero /), (/3,3/) ) )

  CASE(19) !! Model "S"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, zero, three, &
                                 zero, three, four /), (/3,3/) ) )

  CASE(20) !! Model "T"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, zero, zero, &
                                 zero, zero, four /), (/3,3/) ) )

  CASE(21) !! Model "U"
    coef = TRANSPOSE(RESHAPE( (/ four, two, zero, &
                                 two, one, zero, &
                                 zero, zero, zero /), (/3,3/) ) )

  CASE(22) !! Model "V"
    coef = TRANSPOSE(RESHAPE( (/ one, two, four, &
                                 two, five, six, &
                                 four, six, ten /), (/3,3/) ) )
  CASE(23) !! Model "W"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, two, one, &
                                 zero, one, four /), (/3,3/) ) )
  CASE(24) !! Model "X"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, two, three, &
                                 zero, one, five /), (/3,3/) ) )
                                 
  CASE(25) !! Model "Z"
    coef = TRANSPOSE(RESHAPE( (/ zero, zero, zero, &
                                 zero, zero, 1.2_dpp, &
                                 zero, 1.2_dpp, 2.5_dpp /), (/3,3/) ) )
   CASE DEFAULT
    CALL PrntE("Invalid model selection '"//TRIM(GetModelName(model))//"'!", &
                 Q=.TRUE., skip1=1)
  END SELECT
  
  RETURN
  
END SUBROUTINE GetInteractionTable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetProbs(affprob, coef, g1, g2, prev, OR, OR1, OR2, m1, m2)
!! Returns in affprob the probability of being a case for the given parameters
!! the given values of numerical genotypes g1, g2 
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)  :: affprob
  INTEGER(iks), INTENT(IN) :: g1, g2
  REAL(dpp), INTENT(IN)   :: coef(:,:), prev, OR, OR1, OR2, m1, m2
  REAL(dpp)               :: beta0, beta1, beta2, beta3

  !! Assign LRM coefficients
  beta1 = LOG(OR1)
  beta2 = LOG(OR2)
  beta3 = LOG(OR)
  
  !! Define the logistic model intercept so that approximately the population
  !! prevalence of the case status is about 'prev'  
  !beta0 = LOG(prev) - m1 * LOG(OR1) - m2 * LOG(OR2) - m1 * m2 * LOG(OR)
  beta0 = prev - (m1 - m1) - (m2 - m2)
  
  !! Calculate the logistic probability based on the given model (via coef)
  affprob = one/(one+EXP(-beta0 - g1*beta1 - g2*beta2 - coef(g1+1,g2+1)*beta3))

  RETURN
  
END SUBROUTINE GetProbs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SimulateMap(MAPA, nloci, ncols, chcol, rscol, dscol, bpcol, &
                       a1col, a2col, binary, ch_start, rs_start, sglchr)
  IMPLICIT NONE
  CHARACTER(*), INTENT(OUT), ALLOCATABLE :: MAPA(:,:)
  !INTEGER, INTENT(IN)           :: nloci(:)
  INTEGER, INTENT(IN)           :: nloci, ncols
  INTEGER, INTENT(IN), OPTIONAL :: chcol, rscol, dscol, bpcol, a1col, a2col, &
                                   ch_start, rs_start
  LOGICAL, INTENT(IN), OPTIONAL :: binary, sglchr
  INTEGER                       :: chr, rs, k, chcol1, rscol1, dscol1, &
                                   bpcol1, a1col1, a2col1, & !, ii, jj
                                   ch_start1, rs_start1, chr1
  LOGICAL                       :: binary1, sglchr1
                
  binary1 = .FALSE.
  sglchr1 = .TRUE.
  chcol1 = def_map_chcol
  rscol1 = def_map_rscol
  dscol1 = def_map_dscol
  bpcol1 = def_map_bpcol
  a1col1 = def_map_a1col 
  a2col1 = def_map_a2col
  ch_start1 = 0
  rs_start1 = 1
  IF(PRESENT(binary)) binary1 = binary
  IF(PRESENT(sglchr)) sglchr1 = sglchr
  IF(PRESENT(chcol)) chcol1 = chcol
  IF(PRESENT(rscol)) rscol1 = rscol
  IF(PRESENT(dscol)) dscol1 = dscol
  IF(PRESENT(bpcol)) bpcol1 = bpcol
  IF(PRESENT(a1col)) a1col1 = a1col
  IF(PRESENT(a2col)) a2col1 = a2col
  IF(PRESENT(ch_start)) ch_start1 = ch_start
  IF(PRESENT(rs_start)) rs_start1 = rs_start
  
  !! Check for proper size of MAPA
  IF(.NOT.ALLOCATED(MAPA)) &
    CALL ResizeVar(MAPA, nloci, ncols, "")
  IF(SIZE(MAPA,1)/=nloci .OR. SIZE(MAPA,2)/=ncols) &
    CALL ResizeVar(MAPA, nloci, ncols, "")
    
  !! Assign default values to MAPA
  MAPA(:,chcol1) = "0"
  MAPA(:,rscol1) = "rs0"
  IF(dscol1<=SIZE(MAPA,2)) MAPA(:,dscol1) = "0"
  IF(bpcol1<=SIZE(MAPA,2)) MAPA(:,bpcol1) = "0"   
  
  !! First assign "fake" allele names
  IF(MAX(a1col1, a2col1)<=SIZE(MAPA,2)) THEN
    WHERE(MAPA(:,a1col1)=="") MAPA(:,a1col1) = "G"
    WHERE(MAPA(:,a2col1)=="") MAPA(:,a2col1) = "A"
  ENDIF
    
  !! Then the "fake" chromosome numbers and rs numbers
  chr = 0
  rs = 0
  DO k=1,nloci
    !! Assign "fake" chromosome and (unique) rs numbers
    chr1 = ch_start1 + chr
    MAPA(k,chcol1) = i2cp(chr1)
    MAPA(k,rscol1) = "rs"//TRIM(i2cp(rs_start1 + rs))
    rs = rs + 1
    IF(.NOT.sglchr) chr = chr + 1
    IF(chr1 == 22) chr = 0
  ENDDO

  RETURN
  
END SUBROUTINE SimulateMap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE WriteSimulOutput(CAs, COs, nca, nco, pedfile, mapfile, nomapfile, &
                            onlygenes, chck_exist, nocausalpair)
!! Write the data to PLINK format peding and map files. Note that the two 
!! "causal" markers are now the first two on the first chromosome!
  IMPLICIT NONE
  CHARACTER, INTENT(IN)         :: CAs(:,:), COs(:,:)
  CHARACTER(*), INTENT(INOUT)   :: pedfile, mapfile
  INTEGER, INTENT(IN)           :: nca, nco
  LOGICAL, INTENT(IN)           :: nomapfile, onlygenes
  LOGICAL, INTENT(IN), OPTIONAL :: chck_exist, nocausalpair
  LOGICAL                       :: chck_exist1
  CHARACTER                     :: sts, position1
  INTEGER                       :: i, m, nloci, chr, loci1
  
  chck_exist1 = .FALSE.
  IF(PRESENT(chck_exist)) chck_exist1 = chck_exist
  
  position1 = "A"
  IF(chck_exist1) position1 = "R"
   
  CALL OpenFile(UNIT=uped, FILE=pedfile, ACTION='W', STATUS='U', FORM='F',&
                POSITION=position1, ACCESS='S', chck_exist=chck_exist1)
  
  !! Write all loci unless nocausalpair present and true
  loci1 = 1
  IF(PRESENT(nocausalpair)) THEN
    IF(nocausalpair) loci1 = 5
  ENDIF

  !! Get the number of loci (CAs has even number of columns)
  nloci = SIZE(CAs,2)/2
  !! First write cases
  sts = "2"
  DO i = 1, nca + nco 
    IF(.NOT.onlygenes) WRITE(uped,'(I0,9A)', ADVANCE='NO') &
                                  i,tab,"1",tab,"0",tab,"0",tab,"1",tab
    !! Now write controls
    IF(i>nca) sts = "1"
    WRITE(uped,'(A)', ADVANCE='NO') sts
    DO m = loci1, 2*nloci
      IF(i<=nca) WRITE(uped,'(4A)', ADVANCE='NO') tab, CAs(i,m)
      IF(i>nca)  WRITE(uped,'(4A)', ADVANCE='NO') tab, COs(i-nca,m)
    ENDDO 
    WRITE(uped,'(A)') ""
  ENDDO
                               
  CLOSE(uped)
  
  IF(2*nloci<loci1) CALL PrntW("No genetical data saved into ped file!", skip1=1)
    
  IF(nomapfile) RETURN

  !! Write the mapping info
  CALL OpenFile(UNIT=umap, FILE=mapfile, ACTION='W', STATUS='U', &
                FORM='F', POSITION="A", ACCESS='S')
  chr=0
  DO m = 1, nloci
    chr = chr+1
    IF(chr>23) chr=1
    WRITE(umap,'(I0,2A,I0,3A,I0)') &
            chr,tab,"rs",m,tab,"0",tab,m
            !INT(m/(nneutralloci+2)*21.4)+1,tab,"rs",m,tab,"0",tab,m
            
  ENDDO
  CLOSE(umap)
  
  RETURN
  
END SUBROUTINE WriteSimulOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CheckAlleleCount(X, count_minor, count_major)
!! Checks whether the numerical representation of genotypes in X counts the
!! minor/major allele  
  IMPLICIT NONE
  INTEGER(iks), INTENT(INOUT)   :: X(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: count_minor, count_major
  LOGICAL                       :: count_minor1
  INTEGER                       :: i, nalleles, nrow, ncol

  !! If both count_minor and count_major are present with the same value 
  !! it means there is an error in the code
  IF(PRESENT(count_minor) .AND. PRESENT(count_major)) THEN
    IF(count_minor .EQV. count_major) & 
      CALL PrntE("Contradictory allele type (CheckAlleleCount).", Q=.TRUE.)
  ENDIF

  !! Minor allele is the default one
  count_minor1 = .TRUE.
  !! Choose type base on input
  IF(PRESENT(count_minor)) count_minor1 = count_minor
  IF(PRESENT(count_major)) count_minor1 = .NOT.count_major
  
  nrow = SIZE(X,1)
  ncol = SIZE(X,2)
  !! Check for alleles count and if the wrong allele counted, revert X
  DO i=1,ncol
    nalleles = SUM(INT(X(:,i))) - nrow
    IF((count_minor1.AND.nalleles>0) .OR. (.NOT.count_minor1.AND.nalleles<0)) &
      X(:,i) = 2_iks - X(:,i)
  ENDDO

  RETURN
  
END SUBROUTINE CheckAlleleCount

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Bed2Ped(V, U, bed_major, ped_minor, chr, sex, maf, NM, maf_keep)
!! Converts plink BED representation genotypes stored in V to a 0-1-2 format
!! returned in U. In BED bytes are read backwards, which means that 10 is stored 
!! as 01 vice-versa. Thus, NA is coded as 10=1, 0 alleles as 00=0, 1 allele as
!! 01=2 and 2 as 11=3.
!! Summary : 00  Homozygote "1"/"1"
!!           01  Heterozygote
!!           11  Homozygote "2"/"2"
!!           10  Missing genotype
!! Example : 01101100
!!           HGFEDCBA
!! 
!!                 AB   00  -- homozygote (first)
!!               CD     11  -- other homozygote (second)
!!             EF       01  -- heterozygote (third)
!!           GH         10  -- missing genotype (fourth)
!!
!! Parameter chr, sex, maf, NM are optional and specify the chromosome of the
!! current marker (chr), the vector of sex indicators (sex), request for MAF
!! to be computed (maf), number of alleles based on which MAF was computed (NM)
!! and the possibility to indicated which individuals should be used for the
!! calculation of MAF (maf_keep, which should be a boolean vector of the same 
!! length as U) 

  IMPLICIT NONE
  INTEGER(1), INTENT(IN)             :: V(:)
  INTEGER(iks), INTENT(OUT)          :: U(:)
  LOGICAL, INTENT(IN), OPTIONAL      :: bed_major, ped_minor
  INTEGER, INTENT(IN), OPTIONAL      :: chr
  INTEGER(iks), INTENT(IN), OPTIONAL :: sex(:)
  REAL(dpp), INTENT(OUT), OPTIONAL   :: maf
  INTEGER, INTENT(OUT), OPTIONAL     :: NM
  LOGICAL, INTENT(IN), OPTIONAL      :: maf_keep(:)
  INTEGER(2)                         :: VV(SIZE(V)) !, UU
  INTEGER(iks)                       :: U0(4*SIZE(V)) !, U1, U2, U3, U4
  INTEGER                            :: maf_nallele_relev, i, j, sU, sV
  REAL(dpp)                          :: maf_nallele_valid
  LOGICAL                            :: bed_major1, ped_minor1
  
  !! Select major/minor format
  bed_major1 = .TRUE.
  ped_minor1 = .TRUE.
  IF(PRESENT(bed_major)) bed_major1 = bed_major
  IF(PRESENT(ped_minor)) ped_minor1 = ped_minor
  
  !! Store size of U (i.e., sample size) in sU 
  sU = SIZE(U)
  sV = SIZE(V)
  
  !! Make sure integer(1) V is properly stored in integer(2) VV
  VV = V
  WHERE(V<0) VV = INT(256, 2) + V
  
  !! Convert from binary represenation to usual integer representation by 
  !! reading the bits from V in the order 0:1, 2:3, 4:5, 6:7 which mean reading
  !! from the back (smallest power position) of V(1), then going to 8:9, 10:11, 
  !! 12:13, 14:15, which are the bits of V(2), and so on 
  !U0 = INT((/ (( IBITS(VV(i), j+1, 1) + 2*IBITS(VV(i), j, 1), j=0,6,2), i=1,sV) /), iks)
  U0 = INT((/ ((IBITS(VV(i), j, 2), j=0,6,2), i=1,sV) /), iks)
  U = U0(1:sU)
  
!   !! Assign values to U based on V
!   DO i=1,SIZE(V)
!     !! Store V(i) into U with ...
!     IF(V(i)>=0) THEN
!       UU = V(i)
!     ELSE
!       UU = INT(256, 2) + V(i)
!     ENDIF
!     !! Extract the individual pairs of bits
!     U4 = INT( UU / 4**3 , 1)
!     U3 = INT( ( UU - U4*4**3) / 4**2 , 1)
!     U2 = INT( ( UU - U4*4**3 - U3*4**2 ) / 4 , 1)
!     U1 = INT( ( UU - U4*4**3 - U3*4**2 - U2*4 ) , 1)
!     !! Store it to the next (up to) 4 elements in U
!     j = 4*(i-1)
!     IF(j+1<=sU) U(j+1) = U1 
!     IF(j+2<=sU) U(j+2) = U2 
!     IF(j+3<=sU) U(j+3) = U3 
!     IF(j+4<=sU) U(j+4) = U4
!   ENDDO
!
!   print *,""
!   DO i=1,sU
!     print *, U(i), U0(i), U(i) - U0(i)
!   ENDDO
!   print *,"U-U0=", SUM(INT(ABS(U-U0(1:sU))))
!   STOP

  !! Change representation and revert 01 and 10
  !DO i=1,sU
  !  !! Change values in U: 1 -> NumNA, 2 -> 1, 3 -> 2 (and 0 -> 0) 
  !  IF(U(i)==1_iks) THEN
  !    U(i) = NumNA
  !  ELSEIF(U(i)==2_iks) THEN
  !    U(i) = 1_iks
  !  ELSEIF(U(i)==3_iks) THEN 
  !    U(i) = 2_iks
  !  ENDIF
  !
  !  !IF(reverse_counts) THEN
  !  !  IF(U(i)==0_iks) THEN
  !  !    U(i) = 2_iks 
  !  !  ELSEIF(U(i)==2_iks) THEN
  !  !    U(i) = 0_iks
  !  !  ENDIF
  !  !ENDIF
  !  
  !ENDDO
  
  !! Change values in U: 1 -> -1 (=NumNA), 2 -> 1, 3 -> 2 (and 0 -> 0) 
  U = U - 1_iks
  WHERE(U <= 0_iks) U = - U - 1_iks
  
  !! Decide whether to reverse U or not. It needs to be reversed iF bed_major1
  !! and ped_minor1 are logically equivalent, meaning both true or both false
  IF(bed_major1 .EQV. ped_minor1) THEN
    U = 2_iks - U
    WHERE(U == 3_iks) U = NumNA
  ENDIF

  !! Count the number of non-NA alleles (2 per locus). For "special chromosomes"
  !! this is corrected below
  IF(PRESENT(maf) .OR. PRESENT(NM)) THEN
    IF(PRESENT(maf_keep)) THEN
      IF(SIZE(maf_keep)/=SIZE(U)) &
        CALL PrntE("Incompatible maf_keep! Size "//i2c(SIZE(U))//" expected while "//&
                     i2cp(SIZE(maf_keep))//" found.")
      maf_nallele_relev = 2 * COUNT(maf_keep .AND. U/=NumNA)
    ELSE 
      !maf_nallele_relev = 2 * (sU - COUNT(U==NumNA))
      maf_nallele_relev = 2 * COUNT(U/=NumNA)
    ENDIF
  ENDIF 
  !write(*,'I0') U

  !! Look at special chromosomes (only if chr or sex are not missing)
  IF(.NOT.PRESENT(chr) .OR. .NOT.PRESENT(sex)) GOTO 10
  IF(ANY(chr /= special_chr)) GOTO 10
  
  !! Do a "special chromosomes check"
  IF(ANY(chr==(/sex_chr, mitochondrial_chr/))) THEN
    
    DO i=1,sU
          
      !! For chromosomes 23, 24 and mitochondrial chrom. the males must be homo-  
      !! zygous, otherwise invalid, and the "two" alleles are counted only once      
      IF(sex(i) == NumMa .AND. U(i)/=NumNA) THEN
        IF(U(i)==1_iks) U(i) = NumNA
        IF(U(i)==2_iks) U(i) = 1_iks
        IF(PRESENT(maf) .OR. PRESENT(NM)) THEN
          maf_nallele_relev = maf_nallele_relev - 1
          IF(U(i)==NumNA) maf_nallele_relev = maf_nallele_relev - 1
        ENDIF
      ENDIF
  
      !! For chromosome 24 females must be ignored altogether
      IF(sex(i) == NumFe .AND. chr==sex_chr(2)) THEN
        IF((PRESENT(maf) .OR. PRESENT(NM)) .AND. U(i)/=NumNA) &
          maf_nallele_relev = maf_nallele_relev - 2
        U(i) = NumNA
      ENDIF
  
      !! Do a "special chromosome check"
      !! For chromosomes 23, 24 and mitochondrial chrom. the males must be  
      !! homozygous, otherwise invalid, and the "two" alleles are counted only 
      !! once      
      !IF(sex(i) == NumMa .AND. ANY(chr==sex_chr)) THEN
      !  IF(U(i)/=NumNA) maf_nallele_relev = maf_nallele_relev - 1
      !  IF(U(i)==1_iks) THEN
      !    maf_nallele_relev = maf_nallele_relev - 1
      !    U(i) = NumNA
      !  ENDIF
      !  IF(U(i)==2_iks) U(i) = 1_iks
      !ENDIF
      !! Then check mitochondrial chromosomes, where only homozygous individuals
      !! are considered and the "two" alleles are counted only once      
      !IF(chr == mitochondrial_chr) THEN
      !  IF(U(i)/=NumNA) maf_nallele_relev = maf_nallele_relev - 1
      !  IF(U(i)==1_iks) THEN
      !    maf_nallele_relev = maf_nallele_relev - 1
      !    U(i) = NumNA
      !  ENDIF
      !  IF(U(i)==2_iks) U(i) = 1_iks
      !ENDIF
       
    ENDDO
    
  ENDIF

  10 CONTINUE
    
  !! Return also the allele frequency
  IF(PRESENT(maf)) THEN
  
    IF(PRESENT(maf_keep)) THEN
      maf_nallele_valid = SUM(MAX(zero, REAL(U(PACK((/(i, i=1,SIZE(U))/), maf_keep)), dpp)))
    ELSE
      maf_nallele_valid = SUM(MAX(zero, REAL(U, dpp)))
    ENDIF
    
    maf = zero
    IF(maf_nallele_relev>0) maf = maf_nallele_valid / maf_nallele_relev
    IF(.NOT.ped_minor1) maf = one - maf
  
  ENDIF
  
  !! And the number of relevant alleles
  IF(PRESENT(NM)) NM = maf_nallele_relev
  
  RETURN

END SUBROUTINE Bed2Ped

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Ped2Bed(Y, X, bed_major, ped_minor)
!! Coverts the 0-1-2 format in Y to plink BED representation in X
!! In the BED format the bytes are read backwards, which means that 10 (=2) is  
!! 01 (=1) and vice-versa. Thus, NA is coded as 01=1, 0 as 00=0, 1 as 10=2 and 
!! 2 as 11=3.
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN)  :: Y(:,:)
  INTEGER(iks), INTENT(OUT) :: X(:,:)
  LOGICAL, INTENT(IN)       :: bed_major, ped_minor
  INTEGER(iks)              :: Y1(SIZE(Y,1))
  INTEGER                   :: i, j, iloci, sY, z
  LOGICAL                   :: reverse_counts
  
  !! Check for dimension match
  IF(SIZE(Y,2)/=SIZE(X,2)) &
    CALL PrntE("Array dimensions mismatch (Ped2Bed)!", Q=.TRUE.)
    
  !! Get first array dimension
  sY = SIZE(Y,1)
  
  !! Decide whether counts of alleles should be reversed; By default
  !! (bed_major .AND. ped_minor) is true, thus reverse_counts is true
  reverse_counts = .FALSE.
  IF((bed_major .AND. ped_minor) .OR. .NOT.(bed_major .OR. ped_minor)) &
    reverse_counts = .TRUE.
  
  !! Covert data for each locus to a binary format
  DO iloci=1,SIZE(Y,2)
  
    Y1 = Y(:,iloci)
    
    !! If binary format should contain the counts of MAJOR ALLELES (default),
    !! revert the values in Y1 so that it contains the counts of MAJOR ALLELES
    !! Then change values in Y1: 2 -> 3, 1 -> 2, NumNA -> 1 so that they can
    !! be transformed to binary format which stores values in reverse order,
    !! i.e. 01 -> 10, 10 -> 01
    !! Then assign values to X based on Y, where the binary format orders the loci
    !! inside one byte as (4th-3rd-2nd-1st)
    z = 0
    j = 0
    DO i=1,sY

      !! First revert if necessary 2's and 0's
      !IF(Y1(i)>=0 .AND. reverse_counts) Y1(i) = 2_iks - Y1(i)
      IF(reverse_counts) THEN
        IF(Y1(i)==0_iks) THEN
          Y1(i) = 2_iks 
        ELSEIF(Y1(i)==2_iks) THEN
          Y1(i) = 0_iks
        ENDIF
      ENDIF 
        
      !! Make sure NumNA is coded as 01=1, 1 as 10=2 and 2 as 11=3 (ab=a*2+b)
      IF(Y1(i)==NumNA) THEN
        Y1(i) = 1_iks
      ELSEIF(Y1(i)==1_iks) THEN 
        Y1(i) = 2_iks 
      ELSEIF(Y1(i)==2_iks) THEN
        Y1(i) = 3_iks
      ENDIF 

      !! Then store the counts in X
      z = z + Y1(i)*4**j
      j = j + 1
      IF(j==4 .OR. i==sY) THEN
        X((i+3)/4, iloci) = INT(z,1)
        z = 0
        j = 0
      ENDIF
    ENDDO

  ENDDO
  
  RETURN

END SUBROUTINE Ped2Bed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_SIMUL
