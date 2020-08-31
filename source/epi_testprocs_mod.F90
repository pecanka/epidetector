MODULE EPI_TESTPROC

  USE EPI_MATH         !! VARIOUS MATHEMATICAL ROUTINES
  USE EPI_SIMUL        !! NEEDED TO GET THE INTERACTION TABLE (FUNCTION GetInteractionTable)
  
  IMPLICIT NONE

  !! Default Accessibility
  PRIVATE
  
  !! Public statements
  PUBLIC :: ScanData, & 
            GetSubsample, &
            GetAutoLevel, &
            TestS1, &
            TestS2, &
            GetBeta
  
  LOGICAL, PARAMETER :: PF_two_sided = .TRUE.

  REAL(dpp)          :: PF_beta, PF_level, PF_K, PF_lambda, PF_slope
  LOGICAL            :: PF_level_unreliable

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!                    PUBLIC FUNCTIONS AND SUBROUTINES                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ScanData(n, s, x, y, nv, ierr, & 
             !! OPTIONAL parameters from here (except for delta really) !!
             d_co, d_ca, WT1, s1, x1, y1, nv1, s2, x2, y2, nv2, cco1, cco2, &
             cca1, cca2, cor_all, c0, min_cnt, min_ca, min_co, subs, &
             fix_subs, sex, excl_male, excl_fema)
  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: n
  INTEGER(iks), INTENT(INOUT)              :: x(:), y(:), s(:)
  INTEGER, INTENT(OUT)                     :: nv, ierr 
  !! OPTIONAL parameters from here (except for d_co and d_ca really) !!
  INTEGER(iks), INTENT(OUT), OPTIONAL      :: x1(n), y1(n), s1(n), &
                                              x2(n), y2(n), s2(n)
  INTEGER(iks), INTENT(INOUT), OPTIONAL    :: subs(n), sex(n)
  INTEGER, INTENT(IN), OPTIONAL            :: WT1
  INTEGER, INTENT(OUT), OPTIONAL           :: nv1, nv2
  REAL(dpp), INTENT(OUT), TARGET, OPTIONAL :: cco1(dmS,dmS), cco2(dmS,dmS),&
                                              cca1(dmS,dmS), cca2(dmS,dmS)
  REAL(dpp), INTENT(IN), OPTIONAL          :: d_co, d_ca, c0, &
                                              min_cnt, min_ca, min_co
  LOGICAL, INTENT(IN), OPTIONAL            :: fix_subs, excl_male, &
                                              excl_fema, cor_all

  INTEGER(iks)       :: subs1(n)
  REAL(dpp)          :: c01, min_cnt1, min_ca1, min_co1
  INTEGER            :: n1, n2, g1, g2, j, N_co, N_ca, N_co1, N_ca1, N_co2, &
                        N_ca2
  LOGICAL            :: low_cell_count, fix_subs1, cor_all1
  REAL(dpp), POINTER :: cc(:,:)
  
  !!#################!!
  !INTEGER(iks) :: x9(n), y9(n), s9(n), sex9(n), subs9(n)
  !!#################!!
  
  !! Reset error indicator
  ierr = 0
  
  !! Resolve input concerning whether fixed subsample will be used in all tests
  fix_subs1 = .FALSE.
  IF(PRESENT(fix_subs)) fix_subs1 = fix_subs
                   
  !! Default ss is zeros (no individuals used)
  subs1 = 0_iks

  !! If fixed sample is to be used subsample must be present. If it is, use it.
  IF(fix_subs1) THEN
    
    IF(.NOT.PRESENT(subs)) CALL PrntE("Variable subs missing! (ScanData)", Q=.TRUE.)
    
    subs1 = subs
  
  !! If subsample is not fixed, generate it based on d_co and d_ca
  ELSEIF(PRESENT(d_co) .AND. PRESENT(d_ca)) THEN
    
    CALL GetSubsample(subs1, s, n, d_co=d_co, d_ca=d_ca)
  
  !! If d_co or d_ca is missing, then throw an error
  ELSE
    CALL PrntE("S1 subsample is not fixed but d_co and/or"//&
               " d_ca are/is missing! (ScanData).", Q=.TRUE.)
  ENDIF
  
  !!#############!!
  !x9 = INT(x, iks)
  !y9 = INT(y, iks)
  !s9 = INT(s, iks)
  !sex9 = INT(sex, iks)
  !subs9 = INT(subs1, iks)
  
  !CALL CleanUp(n, s9, x9, y9, subs9, nv, n1, n2, ierr, sex9, excl_male, excl_fema)
  
  !x = INT(x9, iks)
  !y = INT(s9, iks)
  !s = INT(s9, iks)
  !sex = INT(sex9, iks)
  !subs1 = INT(subs9, iks)
  !!#############!!

  !! Remove the NA values from the data
  CALL CleanUp(n, s, x, y, subs1, nv, n1, n2, ierr, sex, excl_male, excl_fema)
  
  !!! Return modified optional variables
  IF(PRESENT(subs)) subs = subs1
  IF(PRESENT(nv1))  nv1 = n1
  IF(PRESENT(nv2))  nv2 = n2

  !! Check for no data
  IF(nv==0) THEN
    ierr = err_no_cc
    !! Assign error indicating values in case early return occurred
    IF(PRESENT(x1))   x1 = NumEmpty
    IF(PRESENT(y1))   y1 = NumEmpty
    IF(PRESENT(s1))   s1 = NumEmpty
    IF(PRESENT(x2))   x2 = NumEmpty
    IF(PRESENT(y2))   y2 = NumEmpty
    IF(PRESENT(s2))   s2 = NumEmpty
    IF(PRESENT(cco1)) cco1 = NumEmpty
    IF(PRESENT(cco2)) cco2 = NumEmpty
    IF(PRESENT(cca1)) cca1 = NumEmpty
    IF(PRESENT(cca2)) cca2 = NumEmpty
    RETURN
  ENDIF
  
  !! Check for missing variables and if any are, just return
  IF(.NOT.PRESENT(WT1))  RETURN
  IF(.NOT.PRESENT(x1))   RETURN
  IF(.NOT.PRESENT(y1))   RETURN
  IF(.NOT.PRESENT(s1))   RETURN
  IF(.NOT.PRESENT(x2))   RETURN
  IF(.NOT.PRESENT(y2))   RETURN
  IF(.NOT.PRESENT(s2))   RETURN
  IF(.NOT.PRESENT(nv1))  RETURN
  IF(.NOT.PRESENT(nv2))  RETURN
  IF(.NOT.PRESENT(cco1)) RETURN
  IF(.NOT.PRESENT(cco2)) RETURN
  IF(.NOT.PRESENT(cca1)) RETURN
  IF(.NOT.PRESENT(cca2)) RETURN
  
  !! Get the genotype contingency tables
  CALL GetCounts(nv, s(1:nv), x(1:nv), y(1:nv), subs1(1:nv), n1, s1(1:n1), &
                 x1(1:n1), y1(1:n1), n2, s2(1:n2), x2(1:n2), y2(1:n2), cco1, &
                 cco2, cca1, cca2)
                 
  if(.false.) then
  
  print *,""              
  print *,""              
  print *,"x==0", ANY(x(1:nv)==0_iks)
  print *,"x==1", ANY(x(1:nv)==1_iks)
  print *,"x==2", ANY(x(1:nv)==2_iks)
  print *,""              
  print *,"x1==0", ANY(x1(1:n1)==0_iks)
  print *,"x1==1", ANY(x1(1:n1)==1_iks)
  print *,"x1==2", ANY(x1(1:n1)==2_iks)
  print *,""              
  print *,"x2==0", ANY(x2(1:n2)==0_iks)
  print *,"x2==1", ANY(x2(1:n2)==1_iks)
  print *,"x2==2", ANY(x2(1:n2)==2_iks)
  print *,""              
  print *,"y==0", ANY(y(1:nv)==0_iks)
  print *,"y==1", ANY(y(1:nv)==1_iks)
  print *,"y==2", ANY(y(1:nv)==2_iks)
  print *,""              
  print *,"y1==0", ANY(y1(1:n1)==0_iks)
  print *,"y1==1", ANY(y1(1:n1)==1_iks)
  print *,"y1==2", ANY(y1(1:n1)==2_iks)
  print *,""              
  print *,"y2==0", ANY(y2(1:n2)==0_iks)
  print *,"y2==1", ANY(y2(1:n2)==1_iks)
  print *,"y2==2", ANY(y2(1:n2)==2_iks)
  print *,""              
    
  
  print *,""              
  print *,"cco=", cco1(1,:)+cco2(1,:)
  print *,"cco=", cco1(2,:)+cco2(2,:)
  print *,"cco=", cco1(3,:)+cco2(3,:)
  print *,""              
  print *,"cco1=", cco1(1,:)
  print *,"cco1=", cco1(2,:)
  print *,"cco1=", cco1(3,:)
  print *,""              
  print *,"cco2=", cco2(1,:)
  print *,"cco2=", cco2(2,:)
  print *,"cco2=", cco2(3,:)
  print *,""              
  print *,"cca=", cca1(1,:)+cca2(1,:)
  print *,"cca=", cca1(2,:)+cca2(2,:)
  print *,"cca=", cca1(3,:)+cca2(3,:)
  print *,""              
  print *,"cca1=", cca1(1,:)
  print *,"cca1=", cca1(2,:)
  print *,"cca1=", cca1(3,:)
  print *,""              
  print *,"cca2=", cca2(1,:)
  print *,"cca2=", cca2(2,:)
  print *,"cca2=", cca2(3,:)
  print *,""              
  
  endif
  
  min_cnt1 = zero; IF(PRESENT(min_cnt)) min_cnt1 = min_cnt 

  !! Check for too low counts in any of the groups
  low_cell_count = .FALSE.
  IF(min_cnt1 > 0) THEN
    IF(ANY(cco1 + cco2 < min_cnt1))          low_cell_count = .TRUE.
    IF(ANY(cca1 + cca2 < min_cnt1))          low_cell_count = .TRUE.
    IF(ANY(cco1 < min_cnt1) .AND. WT1/=T1ca) low_cell_count = .TRUE.
    IF(ANY(cca1 < min_cnt1) .AND. WT1/=T1co) low_cell_count = .TRUE.
  ENDIF
  
  !! Return low cell count error number
  IF(low_cell_count) ierr = err_low_cell

  c01 = zero; IF(PRESENT(c0)) c01 = c0 

  !! Yates type correction for zero cell counts and low cell count error removal
  IF(min_cnt1>zero .OR. c01>zero) THEN
    
    cor_all1 = .FALSE.; IF(PRESENT(cor_all)) cor_all1 = cor_all
     
    DO j=1,4
      
      !! Associate pointers to what type of counts will be corrected
      IF(j==1) cc => cco1
      IF(j==2) cc => cco2
      IF(j==3) cc => cca1
      IF(j==4) cc => cca2
  
      !! Check if any correction needed. If so, do it.
      IF(MINVAL(cc) >= min_cnt1) CYCLE
      
      !! First type of correction: add cell count correction to all cells
      IF(cor_all1) THEN
        IF(c01 > zero) THEN
          cc = cc + c01
        ELSE
          cc = MAX(cc, min_cnt1)
        ENDIF
      !! Second type of correction: make all cells at least min_cnt1
      ELSE
        
        DO g1=1,SIZE(cc,1)
          DO g2=1,SIZE(cc,2)
            IF(cc(g1,g2) >= min_cnt1) CYCLE
            IF(c01 > zero) THEN
              cc(g1,g2) = cc(g1,g2) + c01
            ELSE
              cc(g1,g2) = MAX(cc(g1,g2), min_cnt1)
            ENDIF
          ENDDO 
        ENDDO
        
      ENDIF
      ierr = err_low_cell_cor
  
    ENDDO
    
  ENDIF
  
  !! Check for too low sample sizes
  N_co1 = INT(SUM(cco1))
  N_ca1 = INT(SUM(cca1))
  N_co2 = INT(SUM(cco2))
  N_ca2 = INT(SUM(cca2))
  
  min_co1 = zero 
  min_ca1 = zero
  IF(PRESENT(min_co)) min_co1 = min_co 
  IF(PRESENT(min_ca)) min_ca1 = min_ca
   
  !! Sample size checks (CURRENTLY DISABLED)
  !IF(.FALSE.) THEN
  !  do_DS1 = .FALSE.
  !  IF(PRESENT(do_DS))   do_DS1 = do_DS
  ! 
  !  IF(WT1==T1cc .AND. MIN(N_co1,N_ca1) < MIN(min_co1,min_ca1)) &
  !    ierr = err_low_ss_S1
  !  IF(WT1==T1co .AND. N_co1 < min_co1) &
  !    ierr = err_low_ss_S1
  !  IF(WT1==T1ca .AND. N_ca1 < min_ca1) &
  !    ierr = err_low_ss_S1
  !  IF(WT1==T1cc .AND. N_co1 < min_co1 .AND. N_ca1 < min_ca1) &
  !    ierr = err_low_ss_S1
  !
  !  IF(WT1==T1cc .AND. (cntrgrp1>0 .OR. do_DS1)) THEN
  !    IF(N_co2 < min_co1) ierr = err_low_ss_S2
  !    IF(N_ca2 < min_ca1) ierr = err_low_ss_S2
  !  ENDIF
  !
  !  both1 = (cntrgrp1==1 .OR. cntrgrp1==3)
  !  both2 = (cntrgrp1==2 .OR. cntrgrp1==3)
  !
  !  IF(WT1==T1co .AND. cntrgrp1>0) THEN
  !    IF(both1 .AND. N_co2 < min_co1) ierr = err_low_ss_S2
  !    IF(both2 .AND. N_ca2 < min_ca1) ierr = err_low_ss_S2
  !  ENDIF
  !
  !  IF(WT1==T1ca .AND. cntrgrp1>0) THEN
  !    IF(both1 .AND. N_ca2 < min_ca1) ierr = err_low_ss_S2
  !    IF(both2 .AND. N_co2 < min_co1) ierr = err_low_ss_S2
  !  ENDIF
  !ENDIF
    
  !! Count the controls and cases
  N_co = N_co1 + N_co2
  N_ca = N_ca1 + N_ca2
  
  !! Check for less than minimum count
  IF(N_co < min_co1) ierr = err_low_co
  IF(N_ca < min_ca1) ierr = err_low_ca
  IF(N_co < min_co1 .AND. N_ca < min_ca1) ierr = err_low_cc

  !! If no controls or no cases, return error 
  IF(N_co==0) ierr = err_no_co
  IF(N_ca==0) ierr = err_no_ca
  IF(MAX(N_co,N_ca)==0) ierr = err_no_cc
  
  RETURN

  CONTAINS
!END SUBROUTINE ScanData
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE CleanUp(n, s, x, y, ss, nv, n1, n2, ierr, sex, excl_male, excl_fema)
  ! Removes records with invalid genotypes or status
  ! s ... status, x,y ... genotype 1 and 2, nv ... number of valid records,
  ! ierr ... error code, ss ... subsample, sex ... sex indicator
    IMPLICIT NONE
    !!#################!!
    !INTEGER, PARAMETER :: iks = 4
    !!#################!!
    INTEGER, INTENT(IN)                   :: n 
    !INTEGER(iks), INTENT(INOUT)           :: s(n), x(n), y(n), ss(n)
    INTEGER(iks), INTENT(INOUT)           :: s(:), x(:), y(:), ss(:)
    INTEGER, INTENT(OUT)                  :: nv, n1, n2, ierr 
    INTEGER(iks), INTENT(INOUT), OPTIONAL :: sex(n)
    LOGICAL, INTENT(IN), OPTIONAL         :: excl_male, excl_fema
  
    INTEGER(iks)                          :: sex1(n), tx(n), ty(n), ts(n), &
                                             tss(n), tsex1(n)
    INTEGER                               :: j, N_co, N_ca, code
    LOGICAL                               :: excl_male1, excl_fema1, isv(n)
                                           
    !! Reset error indicator
    ierr = 0
    
    !! Define vector of sex
    sex1 = 0_iks
    excl_male1 = .FALSE.
    excl_fema1 = .FALSE.
    IF(PRESENT(sex)) sex1 = sex
    IF(PRESENT(excl_male)) excl_male1 = excl_male
    IF(PRESENT(excl_fema)) excl_fema1 = excl_fema
    
    !! Select how the data is cleaned 
    !! The suspicion is that the most elegant code=2 is the slowest, code=1 should 
    !! be the fastest and is short too, while code=3 is in fact the fastest
    code = 3
    IF(code==1) THEN
    
      !! Make sure all non-valid (or excluded) individuals have the right NA value
      N_co = 0; N_ca = 0; nv = 0; n1 = 0
      DO j=1,n
      
        !! Identify invalid genotypes or status or exclude males/females 
        !IF(x(j)<0_iks .OR. x(j)>2_iks) THEN;      s(j) = NumNA; CYCLE; ENDIF
        !IF(y(j)<0_iks .OR. y(j)>2_iks) THEN;      s(j) = NumNA; CYCLE; ENDIF
        !IF(excl_male1 .AND. sex1(j)==NumMa) THEN; s(j) = NumNA; CYCLE; ENDIF
        !IF(excl_fema1 .AND. sex1(j)==NumFe) THEN; s(j) = NumNA; CYCLE; ENDIF
        !IF(s(j)/=NumCo .AND. s(j)/=NumCa) THEN;   s(j) = NumNA; CYCLE; ENDIF
        IF(x(j)<0_iks .OR. x(j)>2_iks)      CYCLE
        IF(y(j)<0_iks .OR. y(j)>2_iks)      CYCLE
        IF(excl_male1 .AND. sex1(j)==NumMa) CYCLE
        IF(excl_fema1 .AND. sex1(j)==NumFe) CYCLE
        IF(s(j)/=NumCo .AND. s(j)/=NumCa)   CYCLE
  
        !! Count controls and cases
        IF(s(j)==NumCo) N_co = N_co + 1
        IF(s(j)==NumCa) N_ca = N_ca + 1
        
        !! Count S1 individuals
        IF(ss(j)==1_iks) n1 = n1 + 1
        
        !! Shift values closer to the beginning of the vectors
        nv = nv + 1
        IF(nv==j) CYCLE
        x(nv)    = x(j)
        y(nv)    = y(j)
        s(nv)    = s(j)
        ss(nv)   = ss(j)
        sex1(nv) = sex1(j)
        
      ENDDO
    
      !nv = N_co + N_ca
      n2 = nv - n1
      
    ELSEIF(code==2) THEN
    
      !! Make sure all non-valid (or excluded) individuals have the right NA value
      !WHERE(x<0_iks .OR. x>2_iks .OR. y<0_iks .OR. y>2_iks .OR. s/=NumCo .AND. s/=NumCa) s = NumNA
      WHERE(x<0_iks .OR. x>2_iks .OR. y<0_iks .OR. y>2_iks) s = NumNA
      WHERE(s/=NumCo .AND. s/=NumCa) s = NumNA
      
      !! If males/females should be excluded, set their genotypes/phenotypes to NA 
      IF(excl_male1) WHERE(sex1==NumMa) s = NumNA
      IF(excl_fema1) WHERE(sex1==NumFe) s = NumNA
    
      !! Count controls and cases
      N_co = COUNT(s==NumCo)
      N_ca = COUNT(s==NumCa)
      nv = N_co + N_ca
      
      !! Extract the valid values and put them at the beginning of each vector
      IF(nv>0 .AND. nv<SIZE(x)) THEN
        isv = s/=NumNA
        x(1:nv)    = PACK(x, isv)
        y(1:nv)    = PACK(y, isv)
        s(1:nv)    = PACK(s, isv)
        ss(1:nv)   = PACK(ss, isv)
        sex1(1:nv) = PACK(sex1, isv)
      ENDIF
    
      !! Count the valid individuals both stages
      n1 = COUNT(ss(1:nv)==1_iks)
      n2 = nv - n1
      
    ELSEIF(code==3) THEN
    
      !! Reset controls and cases counter 
      N_co = 0; N_ca = 0; nv = 0; n1 = 0
      
      !! Eliminate missing values and store input in tx, ty, ts
      DO j=1,SIZE(x)
    
        !! Check for missing values in x, y or s
        IF(x(j) < 0_iks .OR. x(j) > 2_iks) CYCLE
        IF(y(j) < 0_iks .OR. y(j) > 2_iks) CYCLE
        IF(ALL(s(j) /= (/NumCo, NumCa/))) CYCLE
        
        !! Check for excluding males/females
        IF(excl_male1 .AND. sex1(j) == NumMa) CYCLE
        IF(excl_fema1 .AND. sex1(j) == NumFe) CYCLE
        
        !! Count controls and cases
        IF(s(j) == NumCo) N_co = N_co + 1
        IF(s(j) == NumCa) N_ca = N_ca + 1
        
        !! Count S1 individuals
        IF(ss(j)==1_iks) n1 = n1 + 1
        
        !! Assign valid observations into tx, ty, ts and ts
        nv = nv + 1
        tx(nv) = x(j)
        ty(nv) = y(j)
        ts(nv) = s(j)
        tss(nv) = ss(j)
        tsex1(nv) = sex1(j)
      
      ENDDO
      
      !! Nulify the original variables
      !CALL SetVal(NumNA, x, y, s, sex1, ss) 
      IF(nv<n) THEN
        x(nv+1:n)  = NumNA
        y(nv+1:n)  = NumNA
        s(nv+1:n)  = NumNA
        sex1(nv+1:n) = NumNA
        ss(nv+1:n) = NumNA
      ENDIF
      !x  = NumNA
      !y  = NumNA
      !s  = NumNA
      !sex1 = NumNA
      !ss = NumNA
      
      !! Reassign the data into original variables
      IF(nv>0) THEN
        x(1:nv)  = tx(1:nv)
        y(1:nv)  = ty(1:nv)
        s(1:nv)  = ts(1:nv)
        ss(1:nv) = tss(1:nv)
        sex1(1:nv) = tsex1(1:nv)
      ENDIF
  
      !! Calculate the number of S2 valid observations
      n2 = nv - n1
      
    ELSE
      STOP "Wrong 'code' in CleanUp."
    ENDIF
    
    !! Return modified sex and subs
    IF(PRESENT(sex)) sex = sex1
    
    RETURN
    
  END SUBROUTINE CleanUp
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  SUBROUTINE GetCounts(n, s, x, y, ss, n1, s1, x1, y1, n2, s2, x2, y2, cco1, &
                       cco2, cca1, cca2)
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: n, n1, n2 
    INTEGER(iks), INTENT(IN)  :: x(n), y(n), s(n), ss(n)
    INTEGER(iks), INTENT(OUT)   :: x1(n1), y1(n1), s1(n1), x2(n2), y2(n2), s2(n2)
    REAL(dpp), INTENT(OUT)      :: cco1(dmS,dmS), cco2(dmS,dmS), cca1(dmS,dmS), & 
                                   cca2(dmS,dmS)
  
    INTEGER                     :: g1, g2, j, m1, m2, code
    LOGICAL                     :: use_phase1, isS1(n), isCo1(n1), isCo2(n2), &
                                   x1g1(n1), x2g1(n2)
  
    !! No valid observations
    IF(n==0) THEN
      cco1 = 0
      cca1 = 0
      cco2 = 0
      cca2 = 0
      RETURN
    ENDIF
    
    code = 3
    IF(code==1) THEN
    
      !! Extract individuals from each phase
      isS1 = ss==1_iks
      x1 = PACK(x, isS1)
      y1 = PACK(y, isS1)
      s1 = PACK(s, isS1)
      x2 = PACK(x, .NOT.isS1)
      y2 = PACK(y, .NOT.isS1)
      s2 = PACK(s, .NOT.isS1)
      
      !! Prepare for counting the two-locus genotypes
      isCo1 = s1==NumCo
      isCo2 = s2==NumCo
      
      !! Count the genotypes
      DO g1=0,2
        x1g1 = x1==INT(g1,iks)
        x2g1 = x2==INT(g1,iks)
        DO g2=0,2
          cco1(g1+1,g2+1) = COUNT(x1g1 .AND. isCo1 .AND. y1==INT(g2,iks))
          cca1(g1+1,g2+1) = COUNT(x1g1 .AND. .NOT.isCo1 .AND. y1==INT(g2,iks))
          cco2(g1+1,g2+1) = COUNT(x2g1 .AND. isCo2 .AND. y2==INT(g2,iks))
          cca2(g1+1,g2+1) = COUNT(x2g1 .AND. .NOT.isCo2 .AND. y2==INT(g2,iks))    
        ENDDO
      ENDDO
    
    ELSEIF(code==2) THEN
    
      !! Extract individuals from each phase
      isS1 = ss==1_iks
      x1 = PACK(x, isS1)
      y1 = PACK(y, isS1)
      s1 = PACK(s, isS1)
      x2 = PACK(x, .NOT.isS1)
      y2 = PACK(y, .NOT.isS1)
      s2 = PACK(s, .NOT.isS1)
      
      !! Initialize contingency table counts and sample size counters
      cco1 = 0
      cca1 = 0
    
      !! Separate data for pre-test and post-test
      DO j=1,n1
      
        !! Where to store the counts
        g1 = INT(x1(j)) + 1
        g2 = INT(y1(j)) + 1
        
        !! s contains only valid values 
        IF(s1(j)==NumCo) THEN
          cco1(g1,g2) = cco1(g1,g2) + 1
        ELSE
          cca1(g1,g2) = cca1(g1,g2) + 1
        ENDIF
    
      ENDDO
      
      cco2 = 0
      cca2 = 0
  
      !! Separate data for pre-test and post-test
      DO j=1,n2
      
        !! Where to store the counts
        g1 = INT(x2(j)) + 1
        g2 = INT(y2(j)) + 1
        
        !! s contains only valid values 
        IF(s2(j)==NumCo) THEN
          cco2(g1,g2) = cco2(g1,g2) + 1
        ELSE
          cca2(g1,g2) = cca2(g1,g2) + 1
        ENDIF
    
      ENDDO
      
    ELSEIF(code==3) THEN
    
      !n1 = COUNT(ss(1:nv)==1_iks)
      !n2 = nv - n1
      
      !! Initialize contingency table counts and sample size counters
      cco1 = 0
      cco2 = 0
      cca1 = 0
      cca2 = 0
      m1 = 0
      m2 = 0
    
      !! Separate data for pre-test and post-test
      DO j=1,n
      
        !! Where to store the counts
        g1 = INT(x(j)) + 1
        g2 = INT(y(j)) + 1
        
        !! s contains only valid values 
        SELECT CASE (s(j))
    
        !! CONTROL
        CASE(NumCo) 
    
          !! Determine whether to select or not for pretest
          IF(ss(j)==1_iks) THEN
            cco1(g1,g2) = cco1(g1,g2) + 1
            use_phase1 = .TRUE.
          ELSE
            cco2(g1,g2) = cco2(g1,g2) + 1
            use_phase1 = .FALSE.
          ENDIF
    
        !! CASE
        CASE(NumCa)
    
          !! Determine whether to select or not for pretest
          IF(ss(j)==1_iks) THEN
            cca1(g1,g2) = cca1(g1,g2) + 1
            use_phase1 = .TRUE.
          ELSE
            cca2(g1,g2) = cca2(g1,g2) + 1
            use_phase1 = .FALSE.
          ENDIF
          
        !! Neither control, nor case
        CASE DEFAULT
    
          CYCLE
          
        END SELECT
    
        !! If everything ok with this individual, save the data
        IF(use_phase1) THEN          
          m1 = m1 + 1 
          x1(m1) = x(j)
          y1(m1) = y(j)
          s1(m1) = s(j)
        ELSE
          m2 = m2 + 1
          x2(m2) = x(j)
          y2(m2) = y(j)
          s2(m2) = s(j)
        ENDIF
    
      ENDDO
      
    ELSE
      STOP "Wrong 'code' in GetCounts."
    ENDIF
    
    RETURN
    
  END SUBROUTINE GetCounts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END SUBROUTINE ScanData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetSubsample(subsamp, sts, n, d_co, d_ca, announce)
  IMPLICIT NONE
  INTEGER, INTENT(IN)           :: n
  INTEGER(iks), INTENT(OUT)      :: subsamp(n)
  INTEGER(iks), INTENT(IN)       :: sts(n) 
  REAL(dpp), INTENT(IN)         :: d_co, d_ca
  LOGICAL, INTENT(IN), OPTIONAL :: announce
  LOGICAL                       :: announce1
  INTEGER                       :: i, N_co, N_ca, N_co1, N_ca1, ss(n)
  
  !! Default values
  announce1 = .FALSE.
  IF(PRESENT(announce)) announce1 = announce
  
  !! Check for invalid d_co
  IF(isnan(d_co)) THEN
    CALL PrntE("Invalid d_co: NaN (GetSubsample)", Q=.TRUE.)
  ELSEIF(d_co<0 .OR. d_co>1) THEN
    CALL PrntE("Invalid d_co: "//TRIM(r2c(d_co))//" (GetSubsample)", Q=.TRUE.)
  ENDIF
  !! Check for invalid d_ca
  IF(isnan(d_ca)) THEN
    CALL PrntE("Invalid d_ca: NaN (GetSubsample)", Q=.TRUE.)
  ELSEIF(d_ca<0 .OR. d_ca>1) THEN
    CALL PrntE("Invalid d_ca: "//TRIM(r2c(d_ca))//" (GetSubsample)", Q=.TRUE.)
  ENDIF
  
  !! Check for allocated subsample and proper size
  IF(SIZE(ss)<=0) &
    CALL PrntE("Wrong vector size (GetSubsample)", Q=.TRUE.)
  IF(SIZE(ss)/=SIZE(sts)) &
    CALL PrntE("Vectors 'subsamp' and 'sts' must have the same size"//&
               " (GetSubsample)", Q=.TRUE.)
                    
  !! Announce the action
  IF(announce1) &
    CALL Prnt("Generating S1 subsample ... ", ADVANCE='NO')

  !! Count valid controls and cases
  N_co = 0
  N_ca = 0
  DO i=1,SIZE(sts)
    IF(sts(i) == NumCo) THEN
      N_co = N_co + 1
    ELSEIF(sts(i) == NumCa) THEN
      N_ca = N_ca + 1
    ENDIF
  ENDDO
  
  !! Subsample sizes
  N_co1 = round_int(d_co*N_co)
  N_ca1 = round_int(d_ca*N_ca)

  !! Determine subsamples of controls and cases
  ss = 0

  !! First controls, then cases
  IF(N_co1 > 0) ss(1:N_co1) = 1
  IF(N_ca1 > 0) ss(N_co+1:N_co+N_ca1) = 1

  !! Randomly permutate elements of the vector ss 
  IF(N_co1>0 .AND. N_co1<N_co) CALL Permutate(ss(1:N_co))
  IF(N_ca1>0 .AND. N_ca1<N_ca) CALL Permutate(ss(N_co+1:N_co+N_ca))

  !! Reorder ss so that it corresponds with individuals status
  N_co1 = 0
  N_ca1 = 0
  DO i=1,SIZE(sts)
    IF(sts(i) == NumCo) THEN
      N_co1 = N_co1 + 1
      subsamp(i) = INT(ss(N_co1), iks) 
    ELSEIF(sts(i) == NumCa) THEN
      N_ca1 = N_ca1 + 1
      subsamp(i) = INT(ss(N_co+N_ca1), iks) 
    ENDIF
  ENDDO
  
  !! Announce end
  IF(announce1) THEN
    CALL Prnt0("Finished.", log=.FALSE.)
    CALL Prnt("Finished generating S1 subsample.", screen=.FALSE.)
  ENDIF

  RETURN

END SUBROUTINE GetSubsample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetAutoLevel(plevel, ppower, lambda, slope, model, maf1, maf2, &
                        level2, level_min, level_max, K, b0, b1, b2, b3, &
                        ccco1, ccco2, ccca1, ccca2, disjoint, merge, ierr, &
                        MatOptimum, tol)
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)             :: plevel, ppower, lambda, slope
  INTEGER, INTENT(IN)                :: model
  REAL(dpp), INTENT(INOUT)           :: maf1, maf2
  REAL(dpp), INTENT(IN)              :: level2, level_min, level_max, b0, b1, &
                                        b2, b3, K
  LOGICAL, INTENT(IN)                :: disjoint, merge
  INTEGER, INTENT(INOUT)             :: ierr
  REAL(dpp), INTENT(IN), OPTIONAL    :: tol
  REAL(dpp), INTENT(INOUT), OPTIONAL :: MatOptimum(:,:,:)
  REAL(dpp), INTENT(IN)              :: ccco1(:,:), ccco2(:,:), &
                                        ccca1(:,:), ccca2(:,:)
  REAL(dpp)                          :: tol1, bleft, bright, p(3,1), q(1,3), &
                                        pi(3,3), vpi(dmV), p_co(3,1), &
                                        q_co(1,3), vp_co(3), vq_co(3),& 
                                        pq_co(3,3), pi_co(3,3), pi_ca(3,3), &
                                        vpi_co(dmV), vpi_ca(dmV), &
                                        beta0(3), beta1(4), KL0(dmV,dmF), &
                                        KL1(dmV,dmF), PsiKL0(dmV), PsiKL1(dmV),& 
                                        Fish(dmF,dmF), Fishco(dmF,dmF), &
                                        Fishca(dmF,dmF), N2, N_co, N_co1, &
                                        N_co2, N_ca, N_ca1, N_ca2, N_all, &
                                        a4(1,dmF), a4co(1,dmF), a4ca(1,dmF), &
                                        b(1,dmV), varS0(1,1), varS0co(1,1), &
                                        varS0ca(1,1), varT(1,1), V(dmV,dmV), &
                                        V_I(dmV,dmV), C(dmF,dmV), SD
  INTEGER                            :: ii, jj
  LOGICAL                            :: fresh_parameters, AS_failed
  
  !! Assign initial values
  AS_failed = .FALSE.
  fresh_parameters = .FALSE.
  lambda = NAneg
  slope  = NAneg
  plevel = NAneg
  ppower = NAneg
  
  !! If no mafs inputted (NAneg), then determine them from the pretest data
  IF(ANY((/maf1, maf2/) == NAneg)) CALL GetAlleleFreq(ccco1, maf1, maf2)

  !! Check whether lambda and slope were previously computed for the current 
  !! combination of mafs and if so, re-read it
  ii = 0
  jj = 0
  IF(PRESENT(MatOptimum)) THEN
    ii = MAX(1,INT(maf1*10**opt_ndig_mafs))  
    jj = MAX(1,INT(maf2*10**opt_ndig_mafs))
    !! Restore previously saved values
    plevel = MatOptimum(1,ii,jj)
    ppower = MatOptimum(2,ii,jj)
    lambda = MatOptimum(3,ii,jj)  
    slope  = MatOptimum(4,ii,jj)
    !! Increase the counter of encountered allele combinations
    MatOptimum(5,ii,jj) = MIN(MatOptimum(5,ii,jj) + one, one)
  ENDIF
  
  !! If lambda or slope still unknown at this point, compute them
  IF(ANY((/plevel, lambda, slope, ppower/) == NAneg)) THEN
    fresh_parameters = .TRUE.
  
    !! Set tolerances
    tol1 = epstol
    IF(PRESENT(tol)) tol1 = tol
  
    !! Get genotype probabilities for each locus based on mafs
    p(:,1) = (/ (one-maf1)**2, two*maf1*(one-maf1), maf1**2 /)
    q(1,:) = (/ (one-maf2)**2, two*maf2*(one-maf2), maf2**2 /)

    !! Get genotype probabilities for two-locus genotypes
    pi = MATMUL(p, q)
    
    !! Get the value of Psi at KL based on the values in beta1
    beta1 = (/b0, b1, b2, b3/)
    CALL GetKL(KL1, PsiKL1, model, beta1)
  
    !! Get two-locus genotypes deviated from equilibrium by presence of main
    !! effects or interactions (i.e. non-zero values in beta)
    pi_co = pi * ( one - TRANSPOSE(RESHAPE(PsiKL1, (/3,3/))) )
    pi_co = pi_co / SUM(pi_co)
    pi_ca = pi * TRANSPOSE(RESHAPE(PsiKL1, (/3,3/)))
    pi_ca = pi_ca / SUM(pi_ca)
  
    !! Get marginal probabilities in such population
    p_co = RESHAPE(SUM(pi_co, 2), (/3,1/))
    q_co = RESHAPE(SUM(pi_co, 1), (/1,3/))
    pq_co = MATMUL(p_co, q_co)
    
    !! Store probabilities as vectors
    vp_co = p_co(:,1)
    vq_co = q_co(1,:)
    vpi = RESHAPE(TRANSPOSE(pi), (/dmV/))
    vpi_co = RESHAPE(TRANSPOSE(pi_co), (/dmV/))
    vpi_ca = RESHAPE(TRANSPOSE(pi_ca), (/dmV/))
    
    !! Get pretest sample size
    N_co1 = SUM(ccco1)
    N_co2 = SUM(ccco2)
    N_ca1 = SUM(ccca1)
    N_ca2 = SUM(ccca2)
    N_co = N_co1 + N_co2
    N_ca = N_ca1 + N_ca2
    N_all = N_co1 + N_co2 + N_ca2
    
    !! ********************************************************************** !!
    !!    Compute the non-centrality parameter for the current population     !!
    !! ********************************************************************** !!
    lambda = N_co1 * SUM( (pi_co-pq_co)**2 * Invert(pq_co) )
  
    !! ********************************************************************** !!
    !!  Compute the slope parameter for disjoint or adjuste score statistic   !!
    !! ********************************************************************** !!
    beta0 = (/b0, b1, b2/)
    CALL GetKL(KL0, PsiKL0, model, beta0)
    
    !! Calculate projection a4 from all data or for cases and controls 
    !! separatelly
    print *,"merge=", merge
    IF(merge) THEN
      CALL GetFisher(Fish, PsiKL0, KL0, vpi, tol1)
      print *,"PsiKL0=", sum(abs(PsiKL0))
      print *,"KL0=", sum(abs(KL0))
      print *,"Fish=", sum(abs(Fish))
      CALL GetARM(a4, Fish, ierr, tol1)
      print *,"ierr=", ierr
      print *,"a4=", a4
      IF(ierr/=0) GOTO 10
      !! Get the slope of the disjoint score statistic
      varS0 = MATMUL(MATMUL(a4, Fish), TRANSPOSE(a4))
    ELSE
      CALL GetFisher(Fishco, PsiKL0, KL0, vpi_co, tol1)
      CALL GetFisher(Fishca, PsiKL0, KL0, vpi_ca, tol1)
      CALL GetARM(a4co, Fishco, ierr, tol1)
      IF(ierr/=0) GOTO 10
      CALL GetARM(a4ca, Fishca, ierr, tol1)
      IF(ierr/=0) GOTO 10
      varS0co = MATMUL(MATMUL(a4co, Fishco), TRANSPOSE(a4co))
      varS0ca = MATMUL(MATMUL(a4ca, Fishca), TRANSPOSE(a4ca))
      varS0 = (N_co * varS0co + N_ca * varS0ca) / N_all
      a4 = (N_co * a4co + N_ca * a4ca) / N_all
    ENDIF
    
    10 CONTINUE
    
    !! If there was an error, revert to default pretest level
    IF(ierr/=0) THEN
      plevel = def_level1
      CALL PrntE("Cannot compute optimal S1 level for the current test."//&
                 " Reverting to default level instead.")
      RETURN
    ENDIF

    !! Compute the slope for ADJUSTED score statistic
    IF(.NOT.disjoint) THEN
      IF(merge) THEN
        N2 = four * N_co * N_ca / N_all
      ELSE
        N2 = N_all
      ENDIF
      
      !! Compute variance of the pretest chisquare vector and its pseudoinverse
      CALL GetVarT(V, pi_co, vp_co, vq_co, vp_co, vq_co)
      CALL PseudoInv(V, I=V_I, ierr=ierr, tol=tol1)
      
      !! Check for errors occuring during the call to PseudoInv. If there 
      !! were any errors, do the disjoint test
      IF(ierr/=0) THEN
        AS_failed = .TRUE.
        CALL PrntE("Cannot compute pseudoinverse for the slope of AS."//&
                   " Will use the slope of DS instead.")
      ELSE
        !! Compute the covariance C, regression vector b and variance varT
        CALL GetCovST(C, beta0, N_all, N_co, pi, vp_co, vq_co, vp_co, vq_co, &
                      model)
        !C = C * SQRT(DBLE(N_co1)/N_co)
        print *,"a4=",a4
        b = MATMUL(MATMUL(a4, C), V_I)
        varT = MATMUL(MATMUL(b, V), TRANSPOSE(b))

        !! Now we have all to estimate the slope of AS
        print *,"N2=", N2
        print *,"varS0(1,1)=", varS0(1,1)
        print *,"varT(1,1)=", varT(1,1)
        SD = SQRT(DBLE(varS0(1,1) + varT(1,1)))
        IF(SD == zero) CALL PrntE("Zero variance! (GetAutoLevel)")
        slope = SQRT(DBLE(N2)) * varS0(1,1) / SD
      ENDIF
    ENDIF
    
    !! Compute the slope for DISJOINT score statistic
    IF(disjoint .OR. AS_failed) THEN
      N2 = N_co2 + N_ca2
      slope = SQRT(DBLE(N2 * varS0(1,1)))
    ENDIF

    !! Assign values to variables used by PowerFun
    PF_beta = b3
    PF_level = level2
    PF_K = K
    PF_lambda = lambda
    PF_slope = slope
    
    !! Determine the optimal pretest level given lambda and slope
    IF(K == one) THEN
      plevel = one
    ELSE
      bleft = MAX(level_min, MIN(half,one/K))
      bright = level_max
      plevel = FindMin(PowerFun, bleft, bright)
    ENDIF

    !! Compute the power for the pretest level
    ppower = -PowerFun(plevel)

  ENDIF
  
  !! Store the freshly calculated values unless AS_failed (in which case
  !! the slope for the current combination of maf1-maf2 will be computed again)
  IF(PRESENT(MatOptimum) .AND. fresh_parameters .AND. .NOT.AS_failed) THEN
    MatOptimum(1,ii,jj) = plevel
    MatOptimum(2,ii,jj) = ppower
    MatOptimum(3,ii,jj) = lambda
    MatOptimum(4,ii,jj) = slope
  ENDIF

  RETURN
  
END SUBROUTINE GetAutoLevel
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetChiTest(kind, ccT, T, Tp, ierr, df, Tn, ETn, V, VE, Z, model, &
                      ccE, no_exp, no_var, divide, var_indep, ccV, ccVE)
!! Calculate pretest statistic (T, Tp) based on kind, df, model (which type, 
!! whatdegrees of freedom are used if kind=4, which model if kind=1) using 
!! ccT and returns error code in ierr. Also returned are the pretest    
!! generating vector Tn (calculated from count) and the centerring vector ETn 
!! (calculated from ccE) unless no_exp is true. Additionally, if no_var is 
!! false, variances of Tn and ETn are returned in V and VE using ccV and
!! ccVE. The way Tn and ETn and their variances are calculated is determined   
!! by divide (division by sqrt(p_kq_l)) and var_indep (i.e. independence of  
!! genotypes is assumed if true). If kind=1 the scaling factors in the pretest 
!! with df=1 is returned in Z.
  IMPLICIT NONE
  INTEGER, INTENT(IN)              :: kind
  REAL(dpp), INTENT(IN)            :: ccT(dmS,dmS)
  REAL(dpp), INTENT(OUT), OPTIONAL :: T, Tp
  INTEGER, INTENT(OUT), OPTIONAL   :: ierr
  REAL(dpp), INTENT(OUT), OPTIONAL :: Tn(dmV), ETn(dmV), V(dmV,dmV), &
                                      VE(dmV,dmV), Z(dmV,dmV)
  INTEGER, INTENT(IN), OPTIONAL    :: df, model
  LOGICAL, INTENT(IN), OPTIONAL    :: no_exp, no_var, divide, var_indep
  REAL(dpp), INTENT(IN), OPTIONAL  :: ccE(dmS,dmS), ccV(dmS,dmS), ccVE(dmS,dmS)
  REAL(dpp) :: P(dmS,dmS), P_H0(dmS,dmS), mp(dmS), mq(dmS), Pe(dmS,dmS), & 
               Pe_H0(dmS,dmS), mpe(dmS), mqe(dmS), Pv(dmS,dmS), mpv(dmS), &
               mqv(dmS), Pve(dmS,dmS), mpve(dmS), mqve(dmS), denom(dmS,dmS), &
               mpd(dmS), mqd(dmS), T1, Tp1, s2, TnMat(dmV,1), ETnMat(dmV,1), &
               V1(dmV,dmV), Z1(dmV,dmV), z33(dmS,dmS), z9(dmV,1), NT
  INTEGER   :: model1
  LOGICAL   :: no_exp1, no_var1, divide1
  
  !print *,"P1"
  !! Default error indicator is "not done"
  IF(PRESENT(ierr)) ierr = err_not_done
  
  no_exp1 = .FALSE.
  no_var1 = .FALSE.
  divide1 = .FALSE.
  IF(PRESENT(no_exp)) no_exp1 = no_exp
  IF(PRESENT(no_var)) no_var1 = no_var
  IF(PRESENT(divide)) divide1 = divide
  
  !print *,"P2"
   !! Get the total sample size
  NT = REAL(SUM(ccT), dpp)
  
  !print *,"ccT=", ccT
  
  IF(NT<zero) CALL PrntE("Negative NT!", Q=.TRUE.)

  !! If sample size is 0, then return default values
  IF(NT==zero) THEN
    IF(PRESENT(T)) T = NAv
    IF(PRESENT(Tp)) Tp = NAp
    IF(PRESENT(ierr)) ierr = err_no_data
    RETURN
  ENDIF

  !print *,"P3"
  !! Calculate the test sample frequencies
  CALL GetGenoFreqs(ccT, P, mp, mq, P_H0)
  
  !! Calculate the centering sample frequencies
  IF(PRESENT(ccE)) THEN
    CALL GetGenoFreqs(ccE, Pe, mpe, mqe, Pe_H0)
  ELSE
    Pe = zero; mpe = zero; mqe = zero; Pe_H0 = zero
  ENDIF
  
  !! Check for missing degrees of freedom variable 
  IF(.NOT.PRESENT(df)) &
    CALL PrntE("Missing S1 degrees of freedom 'df'! (GetChiTest)", Q=.TRUE.)
      
  !print *,"P4"
  !! Calculate variance matrix V using 'ccV'
  IF(.NOT.no_var1) THEN
    
    IF(.NOT.PRESENT(ccV)) &
      CALL PrntE("'no_var' is false and 'ccV' is missing! (GetChiTest)", Q=.TRUE.)
    
    !! Calculate the test variance sample frequencies
    CALL GetGenoFreqs(ccV, Pv, mpv, mqv)
  
  !! Or simply use 'ccT' 
  ELSE
    Pv = P
    mpv = mp
    mqv = mq
  ENDIF
  
  !! Decide whether variance matrices are for the vector divided (or not) by
  !! marginal frequencies
  IF(divide1) THEN
    mpd = mp
    mqd = mq
  ELSE
    mpd = one
    mqd = one
  ENDIF
  
  !print *,"P5"
  !! Calculate variance matrix VE using 'ccVE'
  IF(.NOT.no_var1 .AND. .NOT.no_exp1) THEN
    
    IF(.NOT.PRESENT(ccVE)) &
      CALL PrntE("'no_var' and 'no_exp' are false and 'ccE' is missing!"//&
                 " (GetChiTest)", Q=.TRUE.)
    
    !! Calculate the centering variance sample frequencies
    CALL GetGenoFreqs(ccVE, Pve, mpve, mqve)
  
  !! Or simply use 'ccE' 
  ELSE
  
    Pve = Pe
    mpve = mpe
    mqve = mqe
  
  ENDIF 
  
  !print *,"P51"
  !print *,"kind=", kind
  !! Calculate the chisquare-4 test of independence statistic and p-value
  IF(kind == 4) THEN
  
    !print *,"P52"
    !! Calculate the generating and centering vectors
    IF(divide1) THEN
      denom = SQRT(DBLE(Invert(P_H0)))
    ELSE
      denom = one
    ENDIF
    
    TnMat = SQRT(NT) * RESHAPE(TRANSPOSE((P-P_H0) * denom), (/dmV,1/))
    ETnMat = SQRT(NT) * RESHAPE(TRANSPOSE((Pe-Pe_H0) * denom), (/dmV,1/))
    IF(PRESENT(Tn)) Tn = TnMat(:,1)
    IF(PRESENT(ETn)) ETn = ETnMat(:,1)
    
    !! Return the variance of the test generating and the centering vectors
    IF(PRESENT(V) .AND. .NOT.no_var1) &
      CALL GetVarT(V, Pv, mpv, mqv, mpd, mqd, var_indep)
    IF(PRESENT(VE) .AND. .NOT.no_var1) &
      CALL GetVarT(VE, Pve, mpve, mqve, mpd, mqd, var_indep)
      
    !! If that's all, jump to the end
    IF(.NOT.PRESENT(T) .AND. .NOT.PRESENT(Tp)) GOTO 10
    
    !print *,"P55"
    !! Calculate the standardized statistic
    denom = SQRT(DBLE(Invert(P_H0)))
    TnMat = SQRT(NT) * RESHAPE(TRANSPOSE((P-P_H0) * denom), (/dmV,1/))
    T1 = SUM(TnMat**2)
    !print *,"P56"
    
  !! Calculate the chisquare-1 test of independence statistic and p-value
  ELSEIF(kind==1) THEN

    !! Decide about the model
    model1 = 1
    IF(PRESENT(model)) model1 = model
    
      !! Get interaction coefficients
    CALL GetInteractionTable(z33, model1)
    z9 = RESHAPE(TRANSPOSE(z33), (/9,1/))

    !! Calculate the generating and centering vectors
    TnMat = SQRT(NT) * RESHAPE(TRANSPOSE(z33 * (P - P_H0)), (/dmV,1/))
    ETnMat = SQRT(NT) * RESHAPE(TRANSPOSE(z33 * (Pe - Pe_H0)), (/dmV,1/))
    
    !! Return the generating and centering vectors
    IF(PRESENT(Tn)) Tn = TnMat(:,1) 
    IF(PRESENT(ETn)) ETn = ETnMat(:,1)
    
    !! If that's all, jump to the end
    IF(.NOT.PRESENT(T) .AND. .NOT.PRESENT(Tp)) GOTO 10

    !! Calculate variance of the statistic
    CALL GetVarT(V1, Pv, mpv, mqv, one3, one3, var_indep)
    s2 = SUM(MATMUL(MATMUL(TRANSPOSE(z9), V1), z9))
    
    !! Return the variance of the test generating vector
    Z1 = diag(z9(:,1))
    IF(PRESENT(Z)) Z = Z1 
    IF(PRESENT(V)) V = MATMUL(MATMUL(Z1, V1), Z1)
    IF(PRESENT(VE) .AND. .NOT.no_var1) THEN
      CALL GetVarT(VE, Pve, mpve, mqve, one3, one3, var_indep)
      VE = MATMUL(MATMUL(Z1, VE), Z1)
    ENDIF
        
    !! Calculate the standardized statistic
    T1 = SUM(TnMat)**2 * Invert(s2)
    
  !! Unknown test selected
  ELSE
    CALL PrntE("Unknown value of 'kind'! (GetChiTest)", Q=.TRUE.)
  ENDIF

  !print *,"P6"
  if(.false.) then
    print *,""
    print *,"df=", df
    if(present(ETn)) then
      print *,"Tn=", Tn(1:3)
      print *,"Tn=", Tn(4:6)
      print *,"Tn=", Tn(7:9)
      print *,"ETn=", ETn(1:3)
      print *,"ETn=", ETn(4:6)
      print *,"ETn=", ETn(7:9)
      print *,"Tn-ETn=", Tn(1:3) - ETn(1:3)
      print *,"Tn-ETn=", Tn(4:6) - ETn(4:6)
      print *,"Tn-ETn=", Tn(7:9) - ETn(7:9)
      print *,""
      print *,"Pe", Pe(1,1:3)
      print *,"Pe", Pe(2,1:3)
      print *,"Pe", Pe(3,1:3)
      print *,"Pe_H0", Pe_H0(1,1:3)
      print *,"Pe_H0", Pe_H0(2,1:3)
      print *,"Pe_H0", Pe_H0(2,1:3)
    endif
  endif
  
  !! If the statistic is very small, make it exact zero
  IF(T1 < MINvalue) T1 = zero

  !! Calculate the p-value
  Tp1 = one
  IF(T1 > zero) Tp1 = one - pchisq(T1, df)
  
  !! Return the statistic and p-value
  IF(PRESENT(T)) T = T1
  IF(PRESENT(Tp)) Tp = Tp1

  10 CONTINUE

  IF(PRESENT(ierr)) ierr = 0

  !print *,"P7"
  RETURN
  
END SUBROUTINE GetChiTest
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE TestS1(ccco1, ccca1, ccco2, ccca2, WT1, T4, T4p, ierrT4, T4df, T1, &
                  T1p, ierrT1, doPT4, doPT1, doPO4, doPO1, doAS4, doAS1, doDS, &  
                  !! parameters are OPTIONAL from here
                  Tsp4, Tsp4p, ierrsp4, Tsp4df, Tsp1, Tsp1p, ierrsp1, &
                  Tpo4, Tpo4pval, ierrpo4, Tpo4df, Tpo1, Tpo1pval, ierrpo1, &
                  Tn4, V4, V4_I, V4_SqI, Tn1, V1, V1_I, V1_SqI, Z1, WT2, &
                  min_margin, zero_margs, cntrgrp, var_grouped, var_indep, &
                  divided, tol, model)!sts, x, y, n, 
!! Calculates pre-test quantities
!! Parameters are optional after doDS
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)              :: ccco1(dmS,dmS), ccca1(dmS,dmS), &
                                        ccco2(dmS,dmS), ccca2(dmS,dmS)
  INTEGER, INTENT(IN)                :: WT1, T4df
  REAL(dpp), INTENT(OUT)             :: T4, T1
  REAL(dpp), INTENT(INOUT)           :: T4p, T1p
  INTEGER, INTENT(OUT)               :: ierrT4, ierrT1
  LOGICAL, INTENT(IN), OPTIONAL      :: doPT4, doPT1, doPO4, doPO1, &
                                        doAS4, doAS1, doDS
  REAL(dpp), INTENT(INOUT), OPTIONAL :: Tsp4p, Tsp1p, Tpo4pval, Tpo1pval
  REAL(dpp), INTENT(OUT), OPTIONAL   :: Tn4(dmV), Tn1(dmV), Tsp4, Tsp1, Tpo4, Tpo1, &
                                        V4(dmV,dmV), V4_I(dmV,dmV), V4_SqI(dmV,dmV), &
                                        V1(dmV,dmV), V1_I(dmV,dmV), V1_SqI(dmV,dmV), &
                                        Z1(dmV,dmV)
  INTEGER, INTENT(OUT), OPTIONAL     :: ierrsp4, ierrsp1, ierrpo1, ierrpo4 
  !INTEGER(iks), INTENT(IN), OPTIONAL  :: sts(n), x(n), y(n)
  REAL(dpp), INTENT(IN), OPTIONAL    :: min_margin, tol
  LOGICAL, INTENT(IN), OPTIONAL      :: var_grouped, var_indep, divided
  LOGICAL, INTENT(OUT), OPTIONAL     :: zero_margs
  INTEGER, INTENT(IN), OPTIONAL      :: cntrgrp, Tsp4df, Tpo4df, WT2, model

  REAL(dpp) :: ccco(dmS,dmS), ccca(dmS,dmS), cc_all(dmS,dmS), ccT(dmS,dmS), &
               ccTvar(dmS,dmS), ccE(dmS,dmS), ccEvar(dmS,dmS), &
               ETn4(dmV), ETn1(dmV), TnMat(dmV,1), ETn4Mat(dmV,1), Tsp41, &
               Tsp4p1, Tsp11, Tsp1p1, Tpo41, Tpo4pval1, Tpo11, Tpo1pval1, &
               TnMat1(dmV,1), TnMat2(dmV,1), &
               P_co1(dmS,dmS), P_co2(dmS,dmS), &
               P_co(dmS,dmS), P_coH0(dmS,dmS), P_ca1(dmS,dmS), &
               P_ca2(dmS,dmS), P_ca(dmS,dmS), P_caH0(dmS,dmS), &
               P_co1H0(dmS,dmS), P_co2H0(dmS,dmS), &
               P_ca1H0(dmS,dmS), P_ca2H0(dmS,dmS), &
               mp_co(dmS), mq_co(dmS), mp_ca(dmS), mq_ca(dmS), &
               mp_co1(dmS), mq_co1(dmS), mp_co2(dmS), mq_co2(dmS), &
               mp_ca1(dmS), mq_ca1(dmS), mp_ca2(dmS), mq_ca2(dmS), &
               rho2, V_co(dmV,dmV), V_ca(dmV,dmV), &
               VTn4(dmV,dmV), VTn1(dmV,dmV), ZTn1(dmV,dmV), &
               VETn4(dmV,dmV), VETn1(dmV,dmV), &
               numer(dmS,dmS), num_p(dmS), num_q(dmS), &
               denom1(dmS,dmS), denom2(dmS,dmS), Tn41(dmV), &
               Tn11(dmV), N_co, N_co1, N_co2, N_ca, N_ca1, &
               N_ca2, N_all, NT, NE, min_mc, tol1
               !P_all(dmS,dmS), P_allH0(dmS,dmS), mp_all(dmS), mq_all(dmS), & 
               !P1(dmS,dmS), P1_H0(dmS,dmS), mp1(dmS), mq1(dmS), &
               !P2(dmS,dmS), P2_H0(dmS,dmS), mp2(dmS), mq2(dmS), &
               !denom(dmS,dmS), den_p(dmS), den_q(dmS), beta0(0:2), se_beta0(0:2),  
  INTEGER   :: cntrgrp1, WT21, Tsp4df1
  LOGICAL   :: both1, both2, zero_margs1, grp_var, doPT41, &
               doPT11, doPO41, doPO11, doAS41, doAS11, doDS1
  
  !print *," T1"
  
  !! Initialize error indicator
  IF(WT1>0) THEN
    ierrT4 = 0
    ierrT1 = 0
  ENDIF
  
  !! Check for "no data in S1" situation 
  IF(SUM(ccco1 + ccca1) == 0 .AND. WT1>0) THEN
    ierrT4 = err_no_cc 
    ierrT1 = err_no_cc 
    T4 = NAv
    T1 = NAv
    T4p = zero
    T1p = zero
    RETURN
  ENDIF
  
  !! "Initialize" optional variables
  doPT41 = .FALSE.
  doPT11 = .FALSE.
  doPO41 = .FALSE.
  doPO11 = .FALSE.
  doAS41 = .FALSE.
  doAS11 = .FALSE.
  doDS1 = .FALSE.
  min_mc = one
  grp_var = .TRUE.
  cntrgrp1 = 0
  tol1 = epstol
  WT21 = 0
  IF(PRESENT(doPT4))        doPT41 = doPT4 
  IF(PRESENT(doPT1))        doPT11 = doPT1 
  IF(PRESENT(doPO4))        doPO41 = doPO4 
  IF(PRESENT(doPO1))        doPO11 = doPO1 
  IF(PRESENT(doAS4))        doAS41 = doAS4 
  IF(PRESENT(doAS1))        doAS11 = doAS1 
  IF(PRESENT(doDS))         doDS1 = doDS 
  IF(PRESENT(min_margin))   min_mc = min_margin
  IF(PRESENT(var_grouped))  grp_var = var_grouped
  IF(PRESENT(cntrgrp))      cntrgrp1 = cntrgrp
  IF(PRESENT(tol))          tol1 = tol
  IF(PRESENT(WT2))          WT21 = WT2
  
  !! Set the returned variables to values indicating error; if all OK these
  !! values will be changed to return what they should return
  IF(PRESENT(Z1))     Z1 = zero
  IF(PRESENT(V4))     V4 = zero
  IF(PRESENT(V1))     V1 = zero
  IF(PRESENT(V4_I)) V4_I = zero
  IF(PRESENT(V1_I)) V1_I = zero

  !! Initialize local variables
  Tn41  = zero
  Tn11  = zero
  ETn4  = zero
  ETn1  = zero
  VTn4  = zero
  VTn1  = zero
  ZTn1  = zero
  VETn4 = zero
  VETn1 = zero
  numer = one
  num_p = one
  num_q = one
  !denom = one
  !den_p = one
  !den_q = one
  Tsp41 = NAv
  Tsp11 = NAv
  Tsp4p1 = NAp
  Tsp1p1 = NAp

  !print *," T2"
  !! Calculate total counts
  ccco  = ccco1 + ccco2
  ccca  = ccca1 + ccca2
  cc_all = ccco + ccca
  N_co1 = SUM(ccco1)
  N_co2 = SUM(ccco2)
  N_ca1 = SUM(ccca1)
  N_ca2 = SUM(ccca2)
  N_co  = N_co1 + N_co2
  N_ca  = N_ca1 + N_ca2
  N_all = N_co + N_ca
  
  !!! ************************************* !!!
  !!!        CALCULATE SCORE PRETEST        !!!  
  !!! ************************************* !!!
      
  !! SELECTED PRETEST : CONTROLS AND CASES, DO SCORE TEST
  IF(WT1 == T1sc) THEN
  
    CALL PrntE("S1 tests "//TRIM(i2cp(T1sc))//" has been disabled.", Q=.TRUE.)

    !IF(.NOT.(PRESENT(sts) .AND. PRESENT(x) .AND. PRESENT(y))) &
    !  CALL PrntE("All variables 'sts', 'x', 'y' must be present when"//&
    !             " WT1="//i2c(T1sc)//" (TestS1)!", Q=.TRUE.)
    !!! Get logistic regression parameters
    !CALL GetBeta(beta0, sts, x, y, ierrT4, se_beta0, model)
    !ierrT1 = ierrT4
    !
    !!! If error in GetBeta return with values indicating error
    !IF(ierrT4 /= 0) THEN
    !  T4 = NAv
    !  T1 = NAv
    !  IF(PRESENT(Tn4)) Tn4 = NAv
    !  IF(PRESENT(Tn1)) Tn1 = NAv
    !  RETURN
    !ENDIF
    !
    !!! Calculate score statistic
    !CALL GetScore4(T4, beta0, x, y, sts, ierrT4, model=model, tol=tol1, &
    !               var_bound=var_bound, min_ss=min_ss)
    !IF(ierrT4 /= 0) RETURN
    !
    !!! Obtain asymptotic p-value for score tests (if no error occured)
    !T4p = norm_pval(T4)
    !
    !!! Return the same in chisquare-1 test
    !T1 = T4
    !T1p = T4p
    
  ENDIF

  !print *," T3"
  !!! ************************************* !!!
  !!!       CALCULATE POOLED PRETESTS       !!!  
  !!! ************************************* !!!
  
  !! Select counts based on the pretest
  IF(WT1 == T1co) THEN
    ccT = ccco1
    ccTvar = ccT
    IF(grp_var) ccTvar = ccco
  ELSEIF(WT1 == T1ca) THEN
    ccT = ccca1
    ccTvar = ccT
    IF(grp_var) ccTvar = ccca
  ELSEIF(WT1 == T1po) THEN
    ccT = ccca1 + ccco1
    ccTvar = ccT
    IF(grp_var) ccTvar = ccca + ccco
  ENDIF

  !! Check for non-sensical input
  IF(WT1 == T1po .AND. cntrgrp1/=0) &  
    CALL PrntE("With WT1="//TRIM(i2cp(WT1))//" centering is not necessary."//&
               " This must be a programming error.")
                    
  SELECT CASE (cntrgrp1)

  !! Get centering sample counts
  CASE (0)
    
    ccE = zero
    ccEvar = zero
  
  CASE (1)
  
    IF(WT1 == T1co) THEN
      ccE = ccco2
      ccEvar = ccE
      IF(grp_var) ccEvar = ccco
    ELSEIF(WT1 == T1ca) THEN
      ccE = ccca2
      ccEvar = ccE
      IF(grp_var) ccEvar = ccca
    ENDIF
    
  CASE (2)
  
    IF(WT1 == T1co) THEN
      ccE = ccca1
      ccEvar = ccE
      IF(grp_var) ccEvar = ccca
    ELSEIF(WT1 == T1ca) THEN
      ccE = ccco1
      ccEvar = ccE
      IF(grp_var) ccEvar = ccco
    ENDIF

  CASE (3)

    IF(WT1 == T1co) THEN
      ccE = ccca + ccco2
      ccEvar = ccE
    ELSEIF(WT1 == T1ca) THEN
      ccE = ccco + ccca2
      ccEvar = ccE
    ENDIF
    IF(grp_var) ccEvar = ccca + ccco
    
  END SELECT
  
  !! Store sample sizes
  NT = SUM(ccT) 
  NE = SUM(ccE)
  
  !print *," T4"
  !! Calculate frequencies within the pooled sample
  !CALL GetGenoFreqs(cc_all, P_all, mp_all, mq_all, P_allH0)
  
  !!! Calculate frequencies within the centering sample
  !CALL GetGenoFreqs(ccT, P1, mp1, mq1, P1_H0)
  !CALL GetGenoFreqs(ccE, P2, mp2, mq2, P2_H0)

  !!! ************************************** !!!
  !!!     DO THE POOLED POPULATION TESTS     !!!
  !!! ************************************** !!!

  !! Calculate the pooled sample chisquare-4 test statistic and p-value
  IF(doPO41) THEN

    !! Get the statistic and p-value
    CALL GetChiTest(4, cc_all, T=Tpo41, Tp=Tpo4pval1, ierr=ierrpo4, &
                    df=Tpo4df, no_exp=.TRUE., no_var=.TRUE.)

    !! Return the calculated values
    IF(PRESENT(Tpo4)) Tpo4 = Tpo41
    IF(PRESENT(Tpo4pval)) Tpo4pval = Tpo4pval1
    
  ENDIF
  
  !! Calculate the pooled sample chisquare-1 test statistic and p-value (see Lewinger 2013)
  IF(doPO11) THEN

    !! Get the statistic and p-value
    CALL GetChiTest(1, cc_all, T=Tpo11, Tp=Tpo1pval1, ierr=ierrpo1, df=1, &
                    model=model, no_exp=.TRUE., no_var=.TRUE.)

    !! Return the calculated values
    IF(PRESENT(Tpo1)) Tpo1 = Tpo11
    IF(PRESENT(Tpo1pval)) Tpo1pval = Tpo1pval1
    
  ENDIF
  
  !!! ************************************** !!!
  !!!     DO THE SINGLE POPULATION TESTS     !!!
  !!! ************************************** !!!
  
  !print *," T5"
  !! Calculate the control sample chisquare-4 test statistic and p-value
  IF(ANY(WT1 == (/T1co,T1ca,T1po/)) .OR. doPT41) THEN
  
    !! Select degrees of freedom
    IF(ANY(WT1 == (/T1co,T1ca,T1po/))) THEN
      Tsp4df1 = T4df
    ELSEIF(PRESENT(Tsp4df)) THEN
      Tsp4df1 = Tsp4df
    ELSE
      Tsp4df1 = 4
    ENDIF
    
    !! Get the statistic and p-value
    IF(doAS41) THEN
      
      CALL GetChiTest(4, ccT, T=Tsp41, Tp=Tsp4p1, ierr=ierrsp4, df=Tsp4df1, &
                      Tn=Tn41, ETn=ETn4, V=VTn4, VE=VETn4, ccE=ccE, &
                      ccV=ccTvar, ccVE=ccEvar, divide=divided)

      !! Calculate inverse, pseudoinvere, square root pseudoinverse matrices
      IF(WT21>0) CALL PseudoInv(VTn4, I=V4_I, SqI=V4_SqI, tol=tol1, ierr=ierrT4)

    ELSE
      CALL GetChiTest(4, ccT, T=Tsp41, Tp=Tsp4p1, ierr=ierrsp4, df=Tsp4df1, &
                      no_exp=.TRUE., no_var=.TRUE.)
    ENDIF
    
    !! Return the calculated values
    IF(PRESENT(Tsp4))  Tsp4 = Tsp41
    IF(PRESENT(Tsp4p)) Tsp4p = Tsp4p1
    
    IF(ierrT4/=0) RETURN
    
  ENDIF
  
  !print *," T5a"
  !! Calculate the control sample chisquare-1 test statistic and p-value
  IF(ANY(WT1 == (/T1co,T1ca,T1po/)) .OR. doPT11) THEN

    !! Get the statistic and p-value
    IF(doAS11) THEN
      
      CALL GetChiTest(1, ccT, T=Tsp11, Tp=Tsp1p1, ierr=ierrsp1, df=1, Tn=Tn11, &
                      ETn=ETn1, V=VTn1, VE=VETn1, Z=ZTn1, model=model, &
                      ccE=ccE, ccV=ccTvar, ccVE=ccEvar)
                      
      !! Calculate inverse and square-invese matrices
      IF(WT21>0) CALL PseudoInv(VTn1, I=V1_I, SqI=V1_SqI, tol=tol1, ierr=ierrT1)
      !print *,"with AS"
      
    ELSE
    
      !print *,"ccT=", ccT(1,:); print *,"ccT=", ccT(2,:); print *,"ccT=", ccT(3,:)
      CALL GetChiTest(1, ccT, T=Tsp11, Tp=Tsp1p1, ierr=ierrsp1, df=1, &
                      model=model, no_exp=.TRUE., no_var=.TRUE.)
                      
      !print *,"without AS"
    ENDIF
    
    !! Return the calculated values
    IF(PRESENT(Tsp1))  Tsp1 = Tsp11
    IF(PRESENT(Tsp1p)) Tsp1p = Tsp1p1
    IF(PRESENT(Z1))    Z1 = ZTn1
    
    IF(ierrT1/=0) RETURN
  
  ENDIF

  !! Check if any other pretests wanted
  IF(WT1==0) THEN
    ierrT4 = err_not_done
    ierrT1 = err_not_done
    RETURN
  ENDIF    

  !! *************************************** !!
  !!       CHECK WHETHER COUNTS ARE OK       !!
  !! *************************************** !!
  !! Check for too low sample size in pretest

  both1 = ANY(cntrgrp1==(/1,3/))
  both2 = ANY(cntrgrp1==(/2,3/))

  !! Sample size checks (CURRENTLY DISABLED)
  IF(.FALSE.) THEN
    IF(WT21>0) THEN
      IF(WT1==T1cc .AND. MIN(N_co2,N_ca2)==zero .AND. & 
          (doDS1 .OR. cntrgrp1>0))                                  GOTO 102
      IF(WT1==T1co .AND. N_co2==zero .AND. doDS1)                   GOTO 102
      IF(WT1==T1co .AND. cntrgrp1>0 .AND. N_co2==zero .AND. both1)  GOTO 102
      IF(WT1==T1co .AND. cntrgrp1>0 .AND. N_ca2==zero .AND. both2)  GOTO 102
      IF(WT1==T1ca .AND. N_ca2==zero .AND. doDS1)                   GOTO 102
      IF(WT1==T1ca .AND. cntrgrp1>0 .AND. N_ca2==zero .AND. both1)  GOTO 102
      IF(WT1==T1ca .AND. cntrgrp1>0 .AND. N_co2==zero .AND. both2)  GOTO 102
    ENDIF
  ENDIF

  !! If score pretest selected skip the rest
  IF(WT1 == T1sc) RETURN

  !! *************************************** !!
  !!      CHECK WHETHER H0 COUNTS ARE OK     !!
  !! *************************************** !!

  !! Get necessary marginal counts and marginal and cell frequencies
  CALL GetGenoFreqs(ccco, P_co, mp_co, mq_co, P_coH0)
  CALL GetGenoFreqs(ccca, P_ca, mp_ca, mq_ca, P_caH0)
  CALL GetGenoFreqs(ccco1, P_co1, mp_co1, mq_co1, P_co1H0)
  CALL GetGenoFreqs(ccco2, P_co2, mp_co2, mq_co2, P_co2H0)
  CALL GetGenoFreqs(ccca1, P_ca1, mp_ca1, mq_ca1, P_ca1H0)
  CALL GetGenoFreqs(ccca2, P_ca2, mp_ca2, mq_ca2, P_ca2H0)

  zero_margs1 = .FALSE.
  IF(ANY(P_coH0*N_co < min_mc)) GOTO 201
  IF(ANY(P_caH0*N_ca < min_mc)) GOTO 201
  IF(WT1>0) THEN
  
    IF(WT1/=T1ca) THEN
      IF(ANY(P_co1H0*N_co1 == zero))    zero_margs1 = .TRUE.
      IF(ANY(P_co1H0*N_co1 < min_mc))   GOTO 201
    ENDIF
    IF(WT1/=T1co) THEN
      IF(ANY(P_ca1H0*N_ca1 == zero))    zero_margs1 = .TRUE.
      IF(ANY(P_ca1H0*N_ca1 < min_mc))   GOTO 201
    ENDIF
    IF(WT1==T1cc .AND. (cntrgrp1>0 .OR. doDS1)) THEN
      IF(ANY(P_co2H0*N_co2 == zero))    zero_margs1 = .TRUE.
      IF(ANY(P_ca2H0*N_ca2 == zero))    zero_margs1 = .TRUE.
      IF(ANY(P_co2H0*N_co2 < min_mc))   GOTO 201
      IF(ANY(P_ca2H0*N_ca2 < min_mc))   GOTO 201
    ENDIF
    IF(WT1==T1co .AND. cntrgrp1>0) THEN
      IF(both1.AND.ANY(P_co2H0*N_co2 < min_mc)) GOTO 201
      IF(both2.AND.ANY(P_ca2H0*N_ca2 < min_mc)) GOTO 201
    ENDIF
    IF(WT1==T1ca .AND. cntrgrp1>0) THEN
      IF(both1.AND.ANY(P_ca2H0*N_ca2 < min_mc)) GOTO 201
      IF(both2.AND.ANY(P_co2H0*N_co2 < min_mc)) GOTO 201
    ENDIF
    
  ENDIF
  IF(PRESENT(zero_margs)) zero_margs = zero_margs1

  !! ************************************ !!
  !!      GET THE VARIANCE OF T1 - T2     !!
  !! ************************************ !!
  
  !! Get the variance of Tn (CONTROLS AND CASES DIFFERENCE)
  IF(WT1 == T1cc) THEN
    NT = N_co1
    denom1 = SQRT((Invert(P_co1H0)))
    denom2 = SQRT((Invert(P_ca1H0)))
    IF(grp_var) THEN
      CALL GetVarT(V_co, P_co, mp_co, mq_co, one3, one3, var_indep)
      CALL GetVarT(V_ca, P_ca, mp_ca, mq_ca, one3, one3, var_indep)
    ELSE
      CALL GetVarT(V_co, P_co1, mp_co1, mq_co1, one3, one3, var_indep)
      CALL GetVarT(V_ca, P_ca2, mp_ca2, mq_ca2, one3, one3, var_indep)
    ENDIF
    rho2 = N_co / N_ca
    VTn4 = V_co + V_ca * rho2

    !!! remove this once the code is simplified, since VTn4 is already computed above !!!
    !CALL PseudoInv(VTn4, I=V4_I, SqI=V4_SqI, tol=tol1, ierr=ierrT4)
    !IF(ierrT4/=0) RETURN
  ENDIF

  !! *************************************** !!
  !!        GET THE CHISQUARE VECTOR         !!
  !! *************************************** !!

  T4p = one
  T1p = one
  IF(WT1 == T1cc) THEN
    !TnMat = SQRT((NT)) * RESHAPE(TRANSPOSE((P_co1-P_co1H0-P_ca1+P_ca1H0)*denom), (/dmV,1/))
    TnMat1 = SQRT((N_co1)) * RESHAPE(TRANSPOSE((P_co1-P_co1H0)*denom1), (/dmV,1/))
    TnMat2 = SQRT((N_ca1)) * RESHAPE(TRANSPOSE((P_ca1-P_ca1H0)*denom2), (/dmV,1/))
    TnMat = TnMat1 - TnMat2
    Tn41 = TnMat(:,1)
    T4 = SUM(Tn41**2)
    IF(T4>zero) T4p = one - pchisq(T4, T4df)
  ELSEIF(ANY(WT1 == (/T1co,T1ca,T1po/))) THEN
    T4 = Tsp41
    T4p = Tsp4p1
    T1 = Tsp11
    T1p = Tsp1p1
  ENDIF
  
  !! Return chisquare statistics and no errors if null .FALSE.
  IF(PRESENT(Tn4)) Tn4 = Tn41
  IF(PRESENT(Tn1)) Tn1 = Tn11
  IF(PRESENT(V4))  V4 = VTn4
  IF(PRESENT(V1))  V1 = VTn1
  
  !! ******************************************* !!
  !!   CALCULATE STUFF NECESSARY FOR ADJUSTING   !!
  !! ******************************************* !!
  
  !! Calculate centering variables (if centering wanted)
  IF((doAS41 .OR. doAS11) .AND. cntrgrp1 > 0) THEN

    !! Center the pretest vector for regression in the main test
    SELECT CASE (WT1)
    
    !! Using CONTROLS and CASES
    CASE (T1cc)
    
      !! Centering will be done with second portion of CONTROLS and CASES
      !ETn4Mat = SQRT((NT)) * RESHAPE(TRANSPOSE(P_co2-P_co2H0 - P_ca2-P_ca2H0), &
      !                                   (/dmV,1/))
      TnMat1 = SQRT((N_co1)) * RESHAPE(TRANSPOSE((P_co2-P_co2H0)*denom1), (/dmV,1/))
      TnMat2 = SQRT((N_ca1)) * RESHAPE(TRANSPOSE((P_ca2-P_ca2H0)*denom2), (/dmV,1/))
      ETn4Mat = TnMat1 - TnMat2
      
      rho2 = N_co1 / N_co2
      IF(grp_var) THEN
        VETn4 = VTn4 * rho2
      ELSE
        CALL GetVarT(V_co, P_co2, mp_co2, mq_co2, one3, one3, var_indep)
        CALL GetVarT(V_ca, P_ca2, mp_ca2, mq_ca2, one3, one3, var_indep)
        VETn4 = V_co + V_ca * rho2
      ENDIF

    !! Using ONLY CONTROLS
    CASE (T1co, T1ca)
      
      !print *,"NT=", NT
      !print *,"NE=", NE
      
      IF(NE==zero) CALL PrntE("Cannot divide by zero (NE).", Q=.TRUE.)
      
      !print *,"VETn4=", VETn4 
      !print *,"VETn1=", VETn1 
  
      !! Restandardize the variance matrices
      rho2 = NT / NE
      VETn4 = VETn4 * rho2
      VETn1 = VETn1 * rho2
       
    END SELECT

    !! Do the "centering" and get the corresponding variance matrices
    IF(PRESENT(Tn4)) Tn4 = Tn4 - ETn4
    IF(PRESENT(Tn1)) Tn1 = Tn1 - ETn1
    IF(PRESENT(V4))  V4 = VTn4 + VETn4
    IF(PRESENT(V1))  V1 = VTn1 + VETn1
    
  ENDIF

  !! If adjusting through regression not done, make sure the variances are 0
  IF(.NOT.doAS41) THEN
    IF(PRESENT(V4))     V4 = zero
    IF(PRESENT(V4_I))   V4_I = zero
    IF(PRESENT(V4_SqI)) V4_SqI = zero
  ENDIF
  
  IF(.NOT.doAS11) THEN
    IF(PRESENT(V1))     V1 = zero
    IF(PRESENT(V1_I))   V1_I = zero
    IF(PRESENT(V1_SqI)) V1_SqI = zero
  ENDIF
  
  !!! ***************************************** !!!
  !!!     END OF CALCULATING ADJUSTED STUFF     !!!  
  !!! ***************************************** !!!
  
  RETURN

  !! Error in number of cases or controls (data for STAGE 2) 
  102 CONTINUE
  ierrT4 = err_low_ss_S2
  ierrT1 = err_low_ss_S2
  RETURN
  
  !! Error in minimal MARGINAL counts 
  201 CONTINUE
  ierrT4 = err_low_margin
  ierrT1 = err_low_margin
  RETURN
      
END SUBROUTINE TestS1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE TestS2(sts, x, y, nv, n, model, CS, CSp, ierrCS, use_beta1, &
              !! parameters are OPTIONAL from here
              beta0, beta1, se_beta1, doAS4, doAS1, Debug, AS4, AS4p, &
              ierrAS4, AS1, AS1p, ierrAS1, min_ss, var_bound, var_cor_AS, & 
              tol, ccco1, ccco2, ccca1, ccca2, Tn4, V4, V4_I, V4_SqI, Tn1, &
              V1, V1_I, V1_SqI, Z1, varAS_sum, varAS_sum_n, WT1, cntrgrp, &
              cov_grouped, cov_divided, poststand, stand, nonpar_probs)
!! Calculates post-test statistics and p-values. Minimally, it calculates the
!! score statistic. Optionally, it also calculates AS statistic and p-value.
  IMPLICIT NONE
  !! Compulsory arguments
  INTEGER(iks), INTENT(IN)           :: sts(n), x(n), y(n)
  INTEGER, INTENT(IN)                :: nv, n, model
  REAL(dpp), INTENT(INOUT)           :: CS, CSp
  INTEGER, INTENT(OUT)               :: ierrCS
  LOGICAL, INTENT(IN)                :: use_beta1
  !! Optional arguments
  REAL(dpp), INTENT(OUT), OPTIONAL   :: beta0(0:2), beta1(0:3), se_beta1(0:3), &
                                        Tn4(dmV), Tn1(dmV)
  REAL(dpp), INTENT(IN), OPTIONAL    :: min_ss, tol, var_bound, &
                                        var_cor_AS, V4(dmV,dmV), V4_I(dmV,dmV), &
                                        V4_SqI(dmV,dmV), V1(dmV,dmV), V1_I(dmV,dmV), &
                                        V1_SqI(dmV,dmV), Z1(dmV,dmV), &
                                        ccco1(dmS,dmS), ccco2(dmS,dmS), &
                                        ccca1(dmS,dmS), ccca2(dmS,dmS)
  REAL(dpp), INTENT(INOUT), OPTIONAL :: Debug(sDebug), AS4, AS4p, AS1, AS1p, &
                                        varAS_sum, varAS_sum_n
  INTEGER, INTENT(OUT), OPTIONAL     :: ierrAS4, ierrAS1
  INTEGER, INTENT(IN), OPTIONAL      :: WT1, cntrgrp
  LOGICAL, INTENT(IN), OPTIONAL      :: doAS4, doAS1, poststand, &
                                        cov_grouped, cov_divided, nonpar_probs
  LOGICAL, INTENT(INOUT), OPTIONAL   :: stand

  REAL(dpp)                          :: beta01(0:2), se_beta0(0:2), &
                                        beta11(0:3), se_beta11(0:3), delta, &
                                        rr_ca(dmV), C4(dmF,dmV), C1(dmF,dmV), &
                                        C4var(dmF,dmV), C1var(dmF,dmV), tol1, &
                                        Debug1(sDebug-6)
  LOGICAL                            :: doAS41, doAS11
  CHARACTER(mmtl)                    :: text
  
  !! Check for any data
  IF(nv==0) THEN
    ierrCS = err_no_cc
    RETURN
  ENDIF
  
  !print *,"  a"
  
  tol1 = epstol
  IF(PRESENT(tol)) tol1 = tol
  
  !! Decide whether AS is produced
  doAS41 = .FALSE.
  doAS11 = .FALSE.
  IF(PRESENT(AS4)) doAS41 = .TRUE.
  IF(PRESENT(AS1)) doAS11 = .TRUE.
  IF(PRESENT(doAS4)) doAS41 = doAS4
  IF(PRESENT(doAS1)) doAS11 = doAS1
  IF(PRESENT(WT1)) THEN
    IF(WT1==0) THEN
      doAS41 = .FALSE.
      doAS11 = .FALSE.
    ENDIF
  ENDIF
  
  !! Get the null hypothesis estimate for beta
  !print *,""; print *,"n=", n; print *,"y=", size(y)
  CALL GetBeta(beta01, sts, x, y, n, ierrCS, se_beta0)

  !print *,"  c"
  IF(PRESENT(beta0))   beta0 = beta01
  IF(PRESENT(ierrAS4)) ierrAS4 = ierrCS
  IF(PRESENT(ierrAS1)) ierrAS1 = ierrCS
  IF(PRESENT(Debug))   Debug(1:6) = (/ beta01, se_beta0 /)
  IF(ierrCS/=0) RETURN
    
  !! Get also the general estimate for beta
  IF(use_beta1) THEN
    CALL GetBeta(beta11, sts, x, y, n, ierrCS, se_beta11, model)
  ELSE
    beta11 = zero
    se_beta11 = zero
  ENDIF

  IF(PRESENT(beta1)) beta1 = beta11
  IF(PRESENT(se_beta1)) se_beta1 = se_beta11
  
  !print *,"  d"
  !! If only CS wanted then get it and return
  IF(.NOT.doAS41 .AND. .NOT.doAS11) THEN
    
    !! Get the CS statistic
    CALL GetScore4(CS, beta01, x, y, sts, ierrCS, model=model, &
                   beta1=beta11, use_beta1=use_beta1, tol=tol1, &
                   var_bound=var_bound, min_ss=min_ss, Debug=Debug1)

  !print *,"  e"
    IF(PRESENT(Debug)) Debug(7:sDebug) = Debug1

    !! Get the CS p-value
    IF(ierrCS == 0) CSp = norm_pval(CS)
    IF(PRESENT(AS4)) AS4 = NAv
    IF(PRESENT(AS1)) AS1 = NAv
    IF(PRESENT(AS4p)) AS4p = NAp
    IF(PRESENT(AS1p)) AS1p = NAp
    IF(PRESENT(ierrAS4)) ierrAS4 = 0
    IF(PRESENT(ierrAS1)) ierrAS1 = 0
    
    RETURN
  
  ENDIF
  
  !! ************************************************************************ !!
  !!                          DO THE ADJUSTED SCORE                           !!
  !! ************************************************************************ !!
  
  text = ""
  IF(.NOT.PRESENT(AS4))          text = TRIM(text)//" AS4"
  IF(.NOT.PRESENT(AS4p))         text = TRIM(text)//" AS4p"
  IF(.NOT.PRESENT(ierrAS4))      text = TRIM(text)//" ierrAS4"
  IF(.NOT.PRESENT(AS1))          text = TRIM(text)//" AS1"
  IF(.NOT.PRESENT(AS1p))         text = TRIM(text)//" AS1p"
  IF(.NOT.PRESENT(ierrAS1))      text = TRIM(text)//" ierrAS1"
  IF(.NOT.PRESENT(Debug))        text = TRIM(text)//" Debug"
  IF(.NOT.PRESENT(min_ss))       text = TRIM(text)//" min_ss"
  IF(.NOT.PRESENT(var_bound))    text = TRIM(text)//" var_bound"
  IF(.NOT.PRESENT(Tn4))          text = TRIM(text)//" Tn4"
  IF(.NOT.PRESENT(Tn1))          text = TRIM(text)//" Tn1"
  IF(.NOT.PRESENT(V4))           text = TRIM(text)//" V4"
  IF(.NOT.PRESENT(V4_I))         text = TRIM(text)//" V4_I"
  IF(.NOT.PRESENT(V4_SqI))       text = TRIM(text)//" V4_SqI"
  IF(.NOT.PRESENT(V1))           text = TRIM(text)//" V1"
  IF(.NOT.PRESENT(V1_I))         text = TRIM(text)//" V1_I"
  IF(.NOT.PRESENT(V1_SqI))       text = TRIM(text)//" V1_SqI"
  IF(.NOT.PRESENT(Z1))           text = TRIM(text)//" Z1"
  IF(.NOT.PRESENT(WT1))          text = TRIM(text)//" WT1"
  IF(.NOT.PRESENT(cntrgrp))      text = TRIM(text)//" cntrgrp"
  IF(.NOT.PRESENT(varAS_sum))    text = TRIM(text)//" varAS_sum"
  IF(.NOT.PRESENT(varAS_sum_n))  text = TRIM(text)//" varAS_sum_n"
  IF(.NOT.PRESENT(ccco1))        text = TRIM(text)//" ccco1"
  IF(.NOT.PRESENT(ccco2))        text = TRIM(text)//" ccco2"
  IF(.NOT.PRESENT(ccca1))        text = TRIM(text)//" ccca1"
  IF(.NOT.PRESENT(ccca2))        text = TRIM(text)//" ccca2"
  IF(.NOT.PRESENT(poststand))    text = TRIM(text)//" poststand"
  IF(.NOT.PRESENT(stand))        text = TRIM(text)//" stand"
  IF(.NOT.PRESENT(cov_grouped))  text = TRIM(text)//" cov_grouped"
  IF(.NOT.PRESENT(cov_divided))  text = TRIM(text)//" cov_divided"
  IF(.NOT.PRESENT(nonpar_probs)) text = TRIM(text)//" nonpar_probs"

  IF(LEN_TRIM(text)>0) &
    CALL PrntE("Missing variables:"//TRIM(text)//" (TestS2).", Q=.TRUE.)

  !print *,"  e1"
  !! Get the covariance matrix of score and pretest chisquare-4 vector
  CALL GetCovariance(C4, C4var, delta, rr_ca, model, WT1, cntrgrp, &
                     cov_grouped, cov_divided, nonpar_probs, ccco1, &
                     ccco2, ccca1, ccca2, beta01, V4_SqI)

  !print *,"  e2"
  !! Get the covariance matrix of score and pretest chisquare-1 vector
  CALL GetCovariance(C1, C1var, delta, rr_ca, model, WT1, cntrgrp, &
                     cov_grouped, .FALSE., nonpar_probs, ccco1, &
                     ccco2, ccca1, ccca2, beta01, V1_SqI)
  
  !! Rescale the covariance matrices to account for scaling by z(u,v) in Tn1
  C1 = MATMUL(C1, Z1)        
  C1var = MATMUL(C1var, Z1)
  
  !print *,"  f"
  !! Get both AS and CS statistics
  CALL GetScore4(CS, beta01, x, y, sts, ierrCS=ierrCS, model=model, &
                 ierrAS4=ierrAS4, ierrAS1=ierrAS1, beta1=beta1, &
                 use_beta1=use_beta1, tol=tol1, var_bound=var_bound, &
                 var_cor_AS=var_cor_AS, min_ss=min_ss, &
                 Debug=Debug(7:sDebug), doAS=doAS41.OR.doAS11, &
                 AS4=AS4, Tn4=Tn4, C4=C4, C4var=C4var, V4_I=V4_I, V4=V4, &
                 AS1=AS1, Tn1=Tn1, C1=C1, C1var=C1var, V1_I=V1_I, V1=V1, &
                 delta=delta, cov_grouped=cov_grouped, WT1=WT1, &
                 varAS4_sum=varAS_sum, varAS4_sum_n=varAS_sum_n,&
                 stand=stand, rr_ca=rr_ca, nonparam=nonpar_probs)
  !print *,"  f2"
                 
  !! Check for errors
  IF(ierrAS4==err_low_var_CS .OR. ierrAS4==err_neg_var_CS) THEN
    ierrCS = ierrAS4
  ELSEIF(ierrAS4==err_low_var_AS .OR. ierrAS4==err_neg_var_AS) THEN
    IF(.NOT.poststand) AS4 = NAv
  ENDIF

  !! Get the p-values if no error
  IF(ierrCS==0 .AND. CS/=NAv) CSp = norm_pval(CS)
  IF(ierrAS4==0 .AND. AS4/=NAv .AND. stand) AS4p = norm_pval(AS4)
  IF(ierrAS1==0 .AND. AS1/=NAv .AND. stand) AS1p = norm_pval(AS1)

  !print *,"  g"
  RETURN
  
END SUBROUTINE TestS2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetBeta(beta, s, x, y, n, ierr, se_beta, model)
!! Calculates a parameter estimate for the parameter beta in logistic 
!! regression in 'beta'. Depending on the length of beta, which must be either
!! 3 or 4 it returns the null hypothesis estimate (interaction estimate 0) or 
!! the general estimate including the interaction effect.
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)           :: beta(:) 
  INTEGER(iks), INTENT(IN)          :: x(n), y(n), s(n)
  INTEGER, INTENT(IN)              :: n
  INTEGER, INTENT(OUT)             :: ierr
  INTEGER, INTENT(IN), OPTIONAL    :: model
  REAL(dpp), INTENT(OUT), OPTIONAL :: se_beta(:)     
  REAL(dpp)                        :: chisq, deviance, se_beta1(SIZE(beta)), &
                                      Z(SIZE(x),SIZE(beta)-1)
  INTEGER                          :: n1(SIZE(x)), ndf
    
  !! Check for proper lengths
  IF(n/=SIZE(y)) CALL PrntE("Different vector sizes (GetBeta).", Q=.TRUE.)
  IF(n<=0) CALL PrntE("Zero vector size (GetBeta).", Q=.TRUE.)

  !! Assign values for the regression matrix and n1
  Z(:,1) = REAL(x, dpp) 
  Z(:,2) = REAL(y, dpp) 
  
  !! If beta has size four, then calculate the general estimate (Z must have
  !! enough columns for it, so its size is 'size(beta) - 1') 
  IF(SIZE(beta)>3) THEN
    IF(.NOT.PRESENT(model)) &
      CALL PrntE("Missing variable 'model' (GetBeta)!", Q=.TRUE.)
    Z(:,3) = Interaction(Z(:,1), Z(:,2), model)
  ENDIF 
  n1 = 1

  !! Perform logistic regression estimation
  CALL logistic(n, Z, SIZE(Z,2), INT(s), n1, chisq, deviance, ndf, beta, &
                se_beta1, ierr)
                
  !! Announce logistic error
  !IF(ierr /= 0) ierr = 9000 + ierr

  !! Return standard errors
  IF(PRESENT(se_beta)) se_beta = se_beta1

  RETURN
  
END SUBROUTINE GetBeta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!                  PRIVATE FUNCTIONS AND SUBROUTINES                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetCovariance(C, Cvar, delta, rr_ca, model, WT1, cntrgrp, grouped, &
                         divided, nonpar_probs, ccco1, ccco2, ccca1, &
                         ccca2, beta, V_SqI)
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT) :: C(dmF,dmV), Cvar(dmF,dmV), delta, rr_ca(dmV)
  INTEGER, INTENT(IN)    :: model, WT1, cntrgrp
  LOGICAL, INTENT(IN)    :: grouped, divided, nonpar_probs
  REAL(dpp), INTENT(IN)  :: ccco1(dmS,dmS), ccco2(dmS,dmS), ccca1(dmS,dmS), &
                            ccca2(dmS,dmS), beta(:), V_SqI(dmV,dmV) 

  REAL(dpp)              :: P_co(dmS,dmS), P_co1(dmS,dmS), P_co2(dmS,dmS), &
                            P_ca(dmS,dmS), P_ca1(dmS,dmS), P_ca2(dmS,dmS), &
                            P_all(dmS,dmS), P2(dmS,dmS), &
                            mp_co(dmS), mq_co(dmS), mp_ca(dmS), mq_ca(dmS), &
                            mp_all(dmS), mq_all(dmS), & !mp2(dmS), mq2(dmS), &
                            mp_co1(dmS), mq_co1(dmS), mp_co2(dmS), mq_co2(dmS), &
                            mp_ca1(dmS), mq_ca1(dmS), mp_ca2(dmS), mq_ca2(dmS), &
                            C1(dmF,dmV), C2(dmF,dmV), C3(dmF,dmV), &
                            den_p(dmS), den_q(dmS), rho, &
                            ccco(dmS,dmS), ccca(dmS,dmS), &
                            cc_all(dmS,dmS), cc2(dmS,dmS), N1, N2, N, &
                            N_co, N_ca, N_co1, N_ca1, N_co2, N_ca2, N_all !, tau

  !! Number of observations
  N_co1 = SUM(ccco1)
  N_co2 = SUM(ccco2)
  N_ca1 = SUM(ccca1)
  N_ca2 = SUM(ccca2)
  N_co  = N_co1 + N_co2
  N_ca  = N_ca1 + N_ca2
  N_all = N_co + N_ca

  !! Get contingency table frequencies
  ccco = ccco1 + ccco2
  ccca = ccca1 + ccca2
  cc_all = ccco + ccca

  !! Get estimators of probabilities of genotype k,l in controls
  rr_ca = RESHAPE(TRANSPOSE(ccca * Invert(cc_all)), (/dmV/))

  !! Initialize the ratio or controls and cases
  rho  = zero
      
  !! Get marginal cell frequencies
  CALL GetGenoFreqs(ccco, P_co, mp_co, mq_co)
  CALL GetGenoFreqs(ccca, P_ca, mp_ca, mq_ca)
  CALL GetGenoFreqs(cc_all, P_all, mp_all, mq_all)
  CALL GetGenoFreqs(ccco1, P_co1, mp_co1, mq_co1)
  CALL GetGenoFreqs(ccco2, P_co2, mp_co2, mq_co2)
  CALL GetGenoFreqs(ccca1, P_ca1, mp_ca1, mq_ca1)
  CALL GetGenoFreqs(ccca2, P_ca2, mp_ca2, mq_ca2)

  !! If pretest vector not divided, denominators are equal to 1
  den_p = one
  den_q = one

  !! SELECTED PRETEST
  SELECT CASE (WT1)

  !! CONTROLS AND CASES WITH SPLITTING IN PRETEST          
  CASE (T1cc)
  
    delta = N_co1 / N_co
    !delta2 = N_ca1 / N_ca

      !! Compute covariance matrix of score and pretest vector
    IF(grouped) THEN
      CALL GetCovST(C1, beta, N_all, N_co, P_all, mp_co, mq_co, one3, &
                    one3, model, .TRUE., nonpar_probs, rr_ca)
      CALL GetCovST(C2, beta, N_all, N_ca, P_all, mp_ca, mq_ca, one3, &
                    one3, model, .FALSE., nonpar_probs, rr_ca)
    ELSE
      CALL GetCovST(C1, beta, N_all, N_co, P_all, mp_co1, mq_co1, one3, &
                    one3, model, .TRUE., nonpar_probs, rr_ca)
      CALL GetCovST(C2, beta, N_all, N_ca, P_all, mp_ca1, mq_ca1, one3, &
                    one3, model, .FALSE., nonpar_probs, rr_ca)
    ENDIF

    C = MATMUL(C1 - C2, V_SqI) * SQRT((delta))
    Cvar = zero

    !! If pretest NOT centered, then it is correlated with the score
    IF(cntrgrp==0) Cvar = C
  
  !! CONTROLS ONLY IN PRETEST
  CASE(T1co)

    !! Determine what the division is
    IF(divided) THEN
      den_p = mp_co1
      den_q = mq_co1
    ENDIF
    
    N1 = N_co1
    N = N_co

    !! Compute covariance matrix of score and NONCENTERED pretest vector
    IF(grouped) THEN
      CALL GetCovST(C1, beta, N_all, N, P_all, mp_co, mq_co, den_p, &
                    den_q, model, .TRUE., nonpar_probs, rr_ca)
    ELSE
      CALL GetCovST(C1, beta, N_all, N, P_all, mp_co1, mq_co1, den_p, &
                    den_q, model, .TRUE., nonpar_probs, rr_ca)
    ENDIF

    !! Apply factor to take into account that only part of the controls was used
    delta = N1 / N
    C = C1
    
    !! If pretest NOT centered, then it is correlated with the score
    !! Otherwise, get the covariance of score and centered pretest vector
    IF(cntrgrp==0) THEN
      
      Cvar = C1 * SQRT((delta))                             

    !! Compute covariance matrix of score and CENTERED pretest vector
    ELSEIF(cntrgrp==1) THEN
        
      !! Covariances computes from same sample => difference automatically 0
      IF(grouped) THEN
        Cvar = zero
      ELSE
        CALL GetCovST(C2, beta, N_all, N, P_all, mp_co2, mq_co2, den_p, den_q, &
                      model, .TRUE., nonpar_probs, rr_ca)
        Cvar = (C1 - C2) * SQRT((delta))
      ENDIF
      
    ELSEIF(cntrgrp==2) THEN
        
      IF(grouped) THEN
        !CALL GetCovST(C2, beta, N_all, N_ca, P_all, mp_ca, mq_ca, den_p, &
        !              den_q, model, .FALSE., nonpar_probs, rr_ca)
        CALL GetCovST(C2, beta, N_all, N_all-N, P_all, mp_ca, mq_ca, den_p, &
                      den_q, model, .FALSE., nonpar_probs, rr_ca)
      ELSE
        !CALL GetCovST(C2, beta, N_all, N_ca, P_all, mp_ca1, mq_ca1, den_p, &
        !              den_q, model, .FALSE., nonpar_probs, rr_ca)
        CALL GetCovST(C2, beta, N_all, N_all-N, P_all, mp_ca1, mq_ca1, den_p, &
                      den_q, model, .FALSE., nonpar_probs, rr_ca)
      ENDIF
      
      C2 = C2 * SQRT(N/(N_all-N))
      Cvar = (C1 - C2) * SQRT((delta))

    ELSEIF(cntrgrp==3) THEN
    
      cc2 = ccco2 + ccca1 + ccca2
      N2 = SUM(cc2)
      CALL GetGenoFreqs(cc2, P2) 

      !! Compute covariance matrix of score and CENTERED pretest vector
      CALL GetCovST(C2, beta, N_all, N_co, P_all, mp_co2, mq_co2, den_p, &
                    den_q, model, .TRUE., nonpar_probs, rr_ca)
      CALL GetCovST(C3, beta, N_all, N_ca, P_all, mp_ca, mq_ca, den_p, &
                    den_q, model, .FALSE., nonpar_probs, rr_ca)
      C2 = C2 + C3
      Cvar = C1 - C2 * SQRT((N1 / N2))

    ENDIF

  !! CASES ONLY IN PRETEST
  CASE(T1ca)

    !! Determine what the division is
    IF(divided) THEN
      den_p = mp_ca1
      den_q = mq_ca1
    ENDIF

    N1 = N_ca1
    N = N_ca

    !! Compute covariance matrix of score and NONCENTERED pretest vector
    IF(grouped) THEN
      CALL GetCovST(C1, beta, N_all, N, P_all, mp_ca, mq_ca, den_p, &
                    den_q, model, .FALSE., nonpar_probs, rr_ca)
    ELSE
      CALL GetCovST(C1, beta, N_all, N, P_all, mp_ca1, mq_ca1, den_p, &
                    den_q, model, .FALSE., nonpar_probs, rr_ca)
    ENDIF

    !! Apply factor to take into account that only part of cases was used
    delta = N1 / N
    C = C1
    
    !! If pretest NOT centered, then it is correlated with the score
    !! Otherwise, get the covariance of score and centered pretest vector
    IF(cntrgrp==0) THEN
    
      Cvar = C1 * SQRT((delta))

    !! Compute covariance matrix of score and CENTERED pretest vector
    ELSEIF(cntrgrp==1) THEN
        
      !! Covariances computes from same sample => difference automatically 0
      IF(grouped) THEN
        Cvar = zero
      ELSE
        CALL GetCovST(C2, beta, N_all, N, P_all, mp_ca2, mq_ca2, den_p, &
                      den_q, model, .FALSE., nonpar_probs, rr_ca)
        Cvar = (C1 - C2) * SQRT((delta))
      ENDIF

    ELSEIF(cntrgrp==2) THEN
        
      IF(grouped) THEN
        !CALL GetCovST(C2, beta, N_all, N_co, P_all, mp_co, mq_co, den_p, &
        !              den_q, model, .TRUE., nonpar_probs, rr_ca)
        CALL GetCovST(C2, beta, N_all, N_all-N, P_all, mp_co, mq_co, den_p, &
                      den_q, model, .TRUE., nonpar_probs, rr_ca)
      ELSE
        !CALL GetCovST(C2, beta, N_all, N_co, P_all, mp_co1, mq_co1, den_p, &
        !              den_q, model, .TRUE., nonpar_probs, rr_ca)
        CALL GetCovST(C2, beta, N_all, N_all-N, P_all, mp_co1, mq_co1, den_p, &
                      den_q, model, .TRUE., nonpar_probs, rr_ca)
      ENDIF
      
      C2 = C2 * SQRT(N/(N_all-N))
      Cvar = (C1 - C2) * SQRT((delta))

    ELSEIF(cntrgrp==3) THEN
    
      cc2 = ccca2 + ccco1 + ccco2
      N2 = SUM(cc2)
      CALL GetGenoFreqs(cc2, P2) 

      !! Compute covariance matrix of score and CENTERED pretest vector
      CALL GetCovST(C2, beta, N_all, N_ca, P_all, mp_ca2, mq_ca2, &
                    den_p, den_q, model, .FALSE., nonpar_probs, rr_ca)
      CALL GetCovST(C3, beta, N_all, N_co, P_all, mp_co, mq_co, den_p, &
                    den_q, model, .TRUE., nonpar_probs, rr_ca)
      C2 = C2 + C3
      Cvar = C1 - C2 * SQRT((N1 / N2))

    ENDIF

  END SELECT
  
  RETURN
  
END SUBROUTINE GetCovariance
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetScore4(S, beta0, x, y, sts, ierrCS, &
                     !! After this point the parameters are OPTIONAL from here
                     model, ierrAS4, ierrAS1, beta1, use_beta1, tol, &
                     var_bound, var_cor_AS, min_ss, Debug, doAS, AS4, Tn4, &
                     C4, C4var, V4_I, V4, AS1, Tn1, C1, C1var, V1_I, V1, &
                     delta, cov_grouped, WT1, varAS4_sum, varAS4_sum_n, &
                     varAS1_sum, varAS1_sum_n, stand, rr_ca, nonparam)
!! Computes the value of the CS (S) statistic for interaction and 
!! optionally also AS statistic for interaction.
!! After ierr the variables are OPTIONAL and can be left out. Except for model
!! which if left out gives always an error, but it is useful to have it as
!! optional because of the way it is called elsewhere, where model might be 
!! optional input there (and if it is missing there, this function will throw
!! the error for it).
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)              :: beta0(:)
  REAL(dpp), INTENT(OUT)             :: S
  INTEGER(iks), INTENT(IN)            :: sts(:), x(:), y(:)
  INTEGER, INTENT(INOUT)             :: ierrCS
  INTEGER, INTENT(INOUT), OPTIONAL   :: ierrAS4, ierrAS1
  REAL(dpp), INTENT(INOUT), OPTIONAL :: AS4, Tn4(:), AS1, Tn1(:), Debug(:)
  REAL(dpp), INTENT(IN), OPTIONAL    :: beta1(:), tol, var_bound, min_ss, &
                                        var_cor_AS, C4(:,:), C4var(:,:), &
                                        V4(:,:), V4_I(:,:), C1(:,:), &
                                        C1var(:,:), V1(:,:), V1_I(:,:), &
                                        delta, rr_ca(:)
  LOGICAL, INTENT(IN), OPTIONAL      :: use_beta1, doAS, nonparam, &
                                        cov_grouped
  LOGICAL, INTENT(INOUT), OPTIONAL   :: stand
  INTEGER, INTENT(IN), OPTIONAL      :: model, WT1
  REAL(dpp), INTENT(INOUT), OPTIONAL :: varAS4_sum, varAS4_sum_n, &
                                        varAS1_sum, varAS1_sum_n  
  REAL(dpp)                          :: S0, XY(SIZE(x),dmF), PsiXY(SIZE(x)), &
                                        PsiXYnp(SIZE(x)), varSmat(1,1), varS, &
                                        varS0, varAS4, varAS1, varT4(1,1), &
                                        varT1(1,1), b4(dmV), b4mat(1,dmV), &
                                        b1(dmV), a4(1,dmF), aC4var(1,dmV), &
                                        aC4(1,dmV), aC1(1,dmV), aC1var(1,dmV), & 
                                        b1mat(1,dmV), Fish(dmF,dmF), tol1, &
                                        var_bnd, var_cor, R4, R1, bb(dmV)
  CHARACTER(mmtl)                    :: text
  LOGICAL                            :: use_beta11, doAS1, var_nonparam, &
                                        cov_grouped1
  INTEGER                            :: ierrAS, k, N, sD
  
  !! Make sure model selection is not missing
  IF(.NOT.PRESENT(model)) &
    CALL PrntE("Model variable 'model' cannot be missing (GetScore4)!", Q=.TRUE.)
  
  !print *,"        1"
  !! Take care of some optional variables
  ierrAS     = 0
  varAS4     = -one
  varAS1     = -one
  use_beta11 = .FALSE.
  doAS1  = .FALSE.
  tol1       = epstol
  var_bnd    = epstol
  var_cor    = one
  cov_grouped1   = .FALSE.
  IF(PRESENT(use_beta1))   use_beta11   = use_beta1
  IF(PRESENT(doAS))    doAS1    = doAS
  IF(PRESENT(tol))         tol1         = tol  
  IF(PRESENT(var_bound))   var_bnd      = var_bound
  IF(PRESENT(var_cor_AS))  var_cor      = var_cor_AS
  IF(PRESENT(cov_grouped)) cov_grouped1 = cov_grouped
  
  !! Check for missing beta1
  IF(use_beta11 .AND. .NOT.PRESENT(beta1)) &
    CALL PrntE("Cannot use beta1 when beta1 is missing! (GetScore4)", Q=.TRUE.)  
  
  !print *,"        2"
  !! Store sample size in short name
  N = SIZE(x)
  
  !! Check for erroneous input
  IF(SIZE(y)/=N .OR. SIZE(sts)/=N) &
    CALL PrntE("Dimensions do not match (GetScore4).", Q=.TRUE.)
  
  !! Check for too small sample size
  IF(PRESENT(min_ss)) THEN
    IF(N < min_ss) THEN
      ierrCS = err_low_ss_S2
      RETURN
    ENDIF
  ENDIF
  
  !print *,"        3"
  !! ********************************************************** !!
  !!           START OF CLASSICAL SCORE CALCULATIONS            !!
  !! ********************************************************** !!

  !! Assign default values for the regression matrix XY and get PsiXY
  CALL GetKL(XY, PsiXY, model, beta0, REAL(x, dpp), REAL(y, dpp))

  !print *,"        4"
  !! Compute the CLASSICAL SCORE STATISTICS (4th component only)
  S0 = SUM((REAL(sts, dpp) - PsiXY) * XY(:,4)) / SQRT(REAL(N, dpp))

  !! If beta1 is to be used for calculating variance and covariance matrices,
  !! recalculate PsiXY using beta1
  IF(use_beta11) &
    CALL GetKL(XY, PsiXY, model, beta1, REAL(x, dpp), REAL(y, dpp))

  !print *,"        5"
  !! Compute the information matrix for given PsiXY and XY (based either on  
  !! beta0 or beta1 depending on whether use_beta11 is true or false)
  CALL GetFisher(Fish, PsiXY, XY)
  !print *,"        51"

  !! Get asymptotic representation matrix A using given Fish
  CALL GetARM(a4, Fish, ierrAS, tol1)
  !print *,"        52"
  IF(PRESENT(ierrAS4)) ierrAS4 = ierrAS
  IF(PRESENT(ierrAS1)) ierrAS1 = ierrAS
  !print *,"        53"
  IF(ierrAS/=0) RETURN
  
  !! Get the variance of the CS statistic
  varSmat = MATMUL(MATMUL(a4, Fish), TRANSPOSE(a4))
  !print *,"        54"
  varS = varSmat(1,1)
  
  !print *,"        6"
  !! Report error if variance NOT OK
  IF(varS < var_bnd) THEN
    ierrCS = err_low_var_CS
    IF(varS < zero) ierrCS = err_neg_var_CS
    GOTO 10
  ENDIF

  !! Standardize the CS statistic  
  S = S0 / SQRT((varS))
  ierrCS = 0
  
  !print *,"        7"
  !! If ONLY CLASSICAL SCORE statistics wanted skip the rest
  IF(.NOT.doAS1) THEN
    IF(PRESENT(AS4)) AS4 = S
    IF(PRESENT(AS1)) AS1 = S
    IF(PRESENT(ierrAS4)) ierrAS4 = 0  
    IF(PRESENT(ierrAS1)) ierrAS1 = 0  
    GOTO 10
  ENDIF
  
  !! ********************************************************** !!
  !!           START OF ADJUSTED SCORE CALCULATIONS             !!
  !! ********************************************************** !!

  !! Check for missing optional variables
  text = ""
  IF(.NOT.PRESENT(AS4))          text = TRIM(text)//" AS4"
  IF(.NOT.PRESENT(AS1))          text = TRIM(text)//" AS1"
  IF(.NOT.PRESENT(ierrAS4))      text = TRIM(text)//" ierrAS4"
  IF(.NOT.PRESENT(ierrAS1))      text = TRIM(text)//" ierrAS1"
  IF(.NOT.PRESENT(Tn4))          text = TRIM(text)//" Tn4"
  IF(.NOT.PRESENT(Tn1))          text = TRIM(text)//" Tn1"
  IF(.NOT.PRESENT(C4))           text = TRIM(text)//" C4"
  IF(.NOT.PRESENT(C4var))        text = TRIM(text)//" C4var"
  IF(.NOT.PRESENT(V4))           text = TRIM(text)//" V4"
  IF(.NOT.PRESENT(V4_I))         text = TRIM(text)//" V4_I"
  IF(.NOT.PRESENT(C1))           text = TRIM(text)//" C1"
  IF(.NOT.PRESENT(C1var))        text = TRIM(text)//" C1var"
  IF(.NOT.PRESENT(V1))           text = TRIM(text)//" V1"
  IF(.NOT.PRESENT(V1_I))         text = TRIM(text)//" V1_I"
  IF(.NOT.PRESENT(delta))        text = TRIM(text)//" delta"
  IF(.NOT.PRESENT(WT1))          text = TRIM(text)//" WT1"
  IF(.NOT.PRESENT(varAS4_sum))   text = TRIM(text)//" varAS4_sum"
  IF(.NOT.PRESENT(varAS4_sum_n)) text = TRIM(text)//" varAS4_sum_n"
  IF(.NOT.PRESENT(Debug))        text = TRIM(text)//" Debug"

  IF(LEN_TRIM(text)>0) &
    CALL PrntE("Missing variables"//TRIM(text)//" (GetScore4).", Q=.TRUE.)

  !print *,"        8"
  !! Indicator of nonparametric variance computation
  var_nonparam = .FALSE.
  IF(PRESENT(nonparam) .AND. PRESENT(rr_ca)) var_nonparam = nonparam 

  !! If var_nonparam is true, use the non-parametric estimator for PsiXY and 
  !! compute the score variance from it
  IF(var_nonparam) THEN
    DO k=1,N
      IF(3*x(k)+y(k)+1>SIZE(rr_ca)) CALL PrntE("In GetScore4!!!", Q=.TRUE.)
      PsiXYnp(k) = rr_ca(3*x(k)+y(k)+1)
    ENDDO
    !! Get nonparametric estimate of variance
    CALL GetFisher(Fish, PsiXYnp, XY)
    CALL GetARM(a4, Fish, ierrAS, tol1)
    IF(PRESENT(ierrAS4)) ierrAS4 = ierrAS
    IF(PRESENT(ierrAS1)) ierrAS1 = ierrAS
    IF(ierrAS/=0) RETURN
    varSmat = MATMUL(MATMUL(a4, Fish), TRANSPOSE(a4))
    varS0 = varSmat(1,1)
  ELSE
    varS0 = varS
  ENDIF
  
  !print *,"        9"
  !! Get covariance matrix AC
  aC4 = MATMUL(a4, C4)
  aC1 = MATMUL(a4, C1)
  aC4var = MATMUL(a4, C4var)
  aC1var = MATMUL(a4, C1var)

  !! Get regression vectors b4 and b1
  b4mat = SQRT(delta) * MATMUL(aC4, V4_I)
  b1mat = SQRT(delta) * MATMUL(aC1, V1_I)
  b4 = b4mat(1,:)
  b1 = b1mat(1,:)
  
  !print *,"        10"
  !! Compute AS
  R4 = SUM(b4*Tn4)
  R1 = SUM(b1*Tn1)
  !IF(R4*S0<=zero .OR. (R4*S0>zero .AND. R4>S0)) R4 = zero
  !IF(R1*S0<=zero .OR. (R1*S0>zero .AND. R1>S0)) R1 = zero
  
  if(.false.) then
    print *,""
    print *,""
    write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "aC4=", &
    aC4(1,1), aC4(1,2), aC4(1,3), aC4(1,4), aC4(1,5), aC4(1,6), aC4(1,7), aC4(1,8), aC4(1,9)
    write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "b4 =", &
    b4(1), b4(2), b4(3), b4(4), b4(5), b4(6), b4(7), b4(8), b4(9)
    write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "Tn4=", &
    Tn4(1), Tn4(2), Tn4(3), Tn4(4), Tn4(5), Tn4(6), Tn4(7), Tn4(8), Tn4(9)
    Tn1 = b4*Tn4
    write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "b4*Tn4=", &
    Tn1(1), Tn1(2), Tn1(3), Tn1(4), Tn1(5), Tn1(6), Tn1(7), Tn1(8), Tn1(9)
    print *,"R4=", R4
    print *,""
    
    do k=1,9
      write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "V4 =", &
      V4(k,1), V4(k,2), V4(k,3), V4(k,4), V4(k,5), V4(k,6), V4(k,7), V4(k,8), V4(k,9)
    enddo
  
    print *,""
  
    do k=1,9
      write(*,'(a,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)') "V4i=", &
      V4_I(k,1), V4_I(k,2), V4_I(k,3), V4_I(k,4), V4_I(k,5), V4_I(k,6), &
      V4_I(k,7), V4_I(k,8), V4_I(k,9)
    enddo
    
    print *,""
  endif
  
  !print *,"        11"
  !! Check whether projection onto Tn4 actually makes S0 larger (in abs. value)
  IF(.FALSE. .AND. ABS(S0 - R4) > ABS(S0)) THEN

    AS4 = S0
    varAS4 = varS0
    
  ELSE
    
    AS4 = S0 - R4
  
    !! Get the variance of AS (with possibly slightly inflated 
    !! variance of S0 when var_cor is larger than one)
    IF(.FALSE. .AND. cov_grouped1) THEN
      varT4 = MATMUL(MATMUL(b4mat, V4), TRANSPOSE(b4mat)) 
    ELSE 
      varT4 = MATMUL(MATMUL(b4mat, V4) - two*aC4var, TRANSPOSE(b4mat)) 
    ENDIF
    
    varAS4 = var_cor * varS0 + varT4(1,1)
  
    if(.false.) then
      print *,"S0=", S0
      print *,"R4=", R4
      print *,"S0-R4=", S0 - R4
      print *,"varT4=", varT4 
      print *,"varAS4=",varAS4
      stop
    endif 

  ENDIF
  
  !print *,"        12"
  !! Check whether projection onto Tn1 actually makes S0 larger (in abs. value)
  IF(.FALSE. .AND. ABS(S0 - R1) > ABS(S0)) THEN

    AS1 = S0
    varAS1 = varS0

  ELSE
    
    AS1 = S0 - R1
    
    !! Get the variance of AS (with possibly slightly inflated 
    !! variance of S0 when var_cor is larger than one)
    IF(.FALSE. .AND. cov_grouped1) THEN
      varT1 = MATMUL(MATMUL(b1mat, V1), TRANSPOSE(b1mat))
    ELSE 
      varT1 = MATMUL(MATMUL(b1mat, V1) - two*aC1var, TRANSPOSE(b1mat))
    ENDIF
    
    varAS1 = var_cor * varS0 + varT1(1,1)
  
  ENDIF
  
  if(.false.) then
    print *,""
    print *,"b4mat=", b4mat(1,1:9)
    print *,"Tn4=", Tn4(1:9)
    print *,""
    print *,"V4=", V4(9,9)
    print *,"V4_I=", V4_I(9,9)
    print *,""
    bb = b4*Tn4
    print *,"b4*Tn4=", bb(1:9)
    print *,""
    print *,"varS0=", varS0
    print *,"varT4=", varT4
    print *,""
    print *,"S0=", S0
    print *,"SUM(b4*Tn4)=", SUM(b4*Tn4)
    print *,"AS4=", AS4
    print *,"varAS4=", varAS4
    print *,""
    print *,"AS4/sqrt(varAS4)=", AS4 / SQRT(varAS4)
    print *,"S0/sqrt(varS0)=", S0 / SQRT(varS0)
    print *,""
  
    varT4 = MATMUL( MATMUL(b4mat, V4) - two*aC4var, TRANSPOSE(b4mat) ) 
    print *,"varT4=", varT4
    print *,""
    print *,"aC4*b4=", MATMUL(aC4var, TRANSPOSE(b4mat))
    print *,""
  endif

  IF(.FALSE.) THEN
    print *,""
    print *,"delta=", delta
    print *,"sum(aC1var)", sum(abs(aC1var))
    print *,"varS0=", varS0
    print *,"varT4=", varT4
    print *,"varAS4=", varAS4
  ENDIF
  
  !! If this point reached and stand missing, throw an error 
  IF(.NOT.PRESENT(stand)) & 
    CALL PrntE("Variable 'stand' missing (GetScore4).", Q=.TRUE.)

  !! Report error if variance of AS is NOT OK
  IF(varAS4 < var_bnd) THEN
    ierrAS = err_low_var_AS
    IF(varAS4 < zero) ierrAS = err_neg_var_AS
    IF(PRESENT(ierrAS4)) ierrAS4 = ierrAS
    varAS4 = NAv2
    stand = .FALSE.
  ENDIF
  IF(varAS1 < var_bnd) THEN
    ierrAS = err_low_var_AS
    IF(varAS1 < zero) ierrAS = err_neg_var_AS
    IF(PRESENT(ierrAS1)) ierrAS1 = ierrAS
    varAS1 = NAv2
    stand = .FALSE.
  ENDIF

  !print *,"        13"
  !! Variance problem => do not standardize 
  IF(MIN(varAS4,varAS1) < var_bnd) GOTO 10

  !! Standardize AS if variance OK
  AS4 = AS4 / SQRT((varAS4))
  AS1 = AS1 / SQRT((varAS1))
  !! Update the mean of variances of AS
  varAS4_sum_n = varAS4_sum_n + 1
  varAS4_sum = varAS4_sum + varAS4
  IF(PRESENT(varAS1_sum_n) .AND. PRESENT(varAS1_sum)) THEN
    varAS1_sum_n = varAS1_sum_n + 1
    varAS1_sum = varAS1_sum + varAS1
  ENDIF
  stand = .TRUE.
  
  !! ********************************************************** !!
  !!            END OF ADJUSTED SCORE CALCULATIONS              !!
  !! ********************************************************** !!

  !print *,"        14"
  10 CONTINUE
  
  !! Save the weighted sum of pretest vector and variances
  !! NOTE: If the number of return values changes (from 7 at this point) do not
  !! forget to change sDebug in epi_parameters.F90
  IF(PRESENT(Debug)) THEN
    sD = SIZE(Debug)
    IF(sD<1) RETURN
    Debug(1) = S0
    IF(sD<2) RETURN
    Debug(2) = varS
    IF(PRESENT(Tn4)) THEN
      IF(sD<3) RETURN
      Debug(3) = varT4(1,1)
      IF(sD<4) RETURN
      Debug(4) = varAS4
      IF(sD<5) RETURN
      Debug(5) = SUM(b4*Tn4)
      IF(sD<6) RETURN
      Debug(6) = S0 - SUM(b4*Tn4)
    ENDIF
    IF(PRESENT(Tn1)) THEN
      IF(sD<7) RETURN
      Debug(7) = varT1(1,1)
      IF(sD<8) RETURN
      Debug(8) = varAS1
      IF(sD<9) RETURN
      Debug(9) = SUM(b1*Tn1)
      IF(sD<10) RETURN
      Debug(10) = S0 - SUM(b1*Tn1)
    ENDIF
  ENDIF
  
  RETURN

END SUBROUTINE GetScore4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetGenoFreqs(cc, P, mp, mq, P_H0)
!! Routine to compute the marginal frequencies of P
!! Input variable cc is assummed to be a contingency table, while the output
!! variable P contains the cell frequencies and mp contains the row marginal 
!! frequencies and mq the column marginal frequencies
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: cc(dmS,dmS)
  REAL(dpp), INTENT(OUT)           :: P(dmS,dmS)
  REAL(dpp), INTENT(OUT), OPTIONAL :: mp(dmS), mq(dmS), P_H0(dmS,dmS)
  REAL(dpp)                        :: N, mp1(dmS), mq1(dmS) 
  
  P = zero
  N = SUM(cc)
  IF(N>zero) P = cc / N
  mp1 = SUM(P, 2)
  mq1 = SUM(P, 1)

  IF(PRESENT(mp)) mp = mp1
  IF(PRESENT(mq)) mq = mq1
  IF(PRESENT(P_H0)) P_H0 = MATMUL(RESHAPE(mp1, (/3,1/)), RESHAPE(mq1, (/1,3/)))

  RETURN
  
END SUBROUTINE GetGenoFreqs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetAlleleFreq(cc, maf1, maf2)
!! Computes minor allele frequency based on the contingency table cc.
!! 
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)  :: cc(:,:)
  REAL(dpp), INTENT(OUT) :: maf1, maf2
  REAL(dpp)              :: n
  
  !! Total number of alleles
  n = two * SUM(cc)
  
  !! Minor allele frequencies
  maf1 = SUM( SUM(cc, 2) * (/0,1,2/) ) / n
  maf2 = SUM( SUM(cc, 1) * (/0,1,2/) ) / n
  
  RETURN
  
END SUBROUTINE GetAlleleFreq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetCovST(Cov, beta, N_all, N, PQ_all, p, q, p2, q2, model, &
                    controls, nonparam, rr_ca)
!! Computes the covariance matrix of classical score and pretest generating
!! chisquare vector as a difference of two covariance matrices C and D. See
!! the thesis for details. The result is returned in Cov
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)          :: Cov(:,:)
  REAL(dpp), INTENT(IN)           :: beta(:), N_all, N
  REAL(dpp), INTENT(IN)           :: PQ_all(:,:), p(:), q(:), p2(:), q2(:)
  INTEGER, INTENT(IN)             :: model
  LOGICAL, INTENT(IN), OPTIONAL   :: controls, nonparam
  REAL(dpp), INTENT(IN), OPTIONAL :: rr_ca(:)
  LOGICAL                         :: controls1, nonparam1
  REAL(dpp)                       :: PQ(dmV), pp2(dmV), qq2(dmV), KL(dmV,dmF), &
                                     Psi_KL(dmV), C(SIZE(Cov,1),SIZE(Cov,2)), &
                                     D(SIZE(Cov,1),SIZE(Cov,2)), tau_sqrt, &
                                     d1(dmF,3), d2(dmF,3), PsiDif(dmV), &
                                     sp2(SIZE(p2)), sq2(SIZE(q2)), spp2(dmV), &
                                     sqq2(dmV), isp2(SIZE(p2)), &
                                     isq2(SIZE(q2)), rr_ca1(SIZE(Psi_KL))
  INTEGER                         :: i,j
    
  !! Assign values for the regression matrix KL to compute Cov and compute
  !! the values of Psi at KL
  CALL GetKL(KL, Psi_KL, model, beta)

  !! Ratio of subsample and sample
  tau_sqrt = REAL(SQRT((N / N_all)), dpp)
  
  !! Initialize vectors of probabilities
  pp2 = (/ (p2(i), p2(i), p2(i), i=1,3) /) 
  qq2 = (/ (q2, i=1,3) /)
  sp2 = REAL(SQRT((p2)), dpp) 
  sq2 = REAL(SQRT((q2)), dpp)
  spp2 = REAL(SQRT((pp2)), dpp) 
  sqq2 = REAL(SQRT((qq2)), dpp)
  isp2 = Invert(sp2) 
  isq2 = Invert(sq2) 

  !! Get nonparametric probablity of being a case (if nonparam is present and true)
  rr_ca1 = one
  nonparam1 = .FALSE.
  controls1 = .TRUE.
  IF(PRESENT(rr_ca))    rr_ca1 = rr_ca
  IF(PRESENT(nonparam)) nonparam1 = nonparam
  IF(PRESENT(controls)) controls1 = controls

  IF(nonparam1) Psi_KL = rr_ca1
  
  !! Store PQ_all as vector (row-wise)
  PQ = RESHAPE(TRANSPOSE(PQ_all), (/dmV/))
  
  PsiDif = Psi_KL * (one - Psi_KL)
  IF(controls1) PsiDif = - PsiDif
  
  !! Compute C as a sum of Psi' and z*z' over n
  DO j=1,SIZE(KL,2)
    C(j,:) = PQ * KL(:,j) * PsiDif * Invert(spp2*sqq2) / tau_sqrt 
  ENDDO

  !! Compute the elements which make up D
  d1(:,1) = sp2(1) * C(:,1) + sp2(2) * C(:,4) + sp2(3) * C(:,7)
  d1(:,2) = sp2(1) * C(:,2) + sp2(2) * C(:,5) + sp2(3) * C(:,8)
  d1(:,3) = sp2(1) * C(:,3) + sp2(2) * C(:,6) + sp2(3) * C(:,9)
                                                      
  d2(:,1) = sq2(1) * C(:,1) + sq2(2) * C(:,2) + sq2(3) * C(:,3)
  d2(:,2) = sq2(1) * C(:,4) + sq2(2) * C(:,5) + sq2(3) * C(:,6)
  d2(:,3) = sq2(1) * C(:,7) + sq2(2) * C(:,8) + sq2(3) * C(:,9)
                                                      
  !! Compute D from these elements                    
  D(:,1) = p(1) * isp2(1) * d1(:,1) + q(1) * isq2(1) * d2(:,1)
  D(:,2) = p(1) * isp2(1) * d1(:,2) + q(2) * isq2(2) * d2(:,1)
  D(:,3) = p(1) * isp2(1) * d1(:,3) + q(3) * isq2(3) * d2(:,1)
                                     
  D(:,4) = p(2) * isp2(2) * d1(:,1) + q(1) * isq2(1) * d2(:,2)
  D(:,5) = p(2) * isp2(2) * d1(:,2) + q(2) * isq2(2) * d2(:,2)
  D(:,6) = p(2) * isp2(2) * d1(:,3) + q(3) * isq2(3) * d2(:,2)
                                     
  D(:,7) = p(3) * isp2(3) * d1(:,1) + q(1) * isq2(1) * d2(:,3)
  D(:,8) = p(3) * isp2(3) * d1(:,2) + q(2) * isq2(2) * d2(:,3)
  D(:,9) = p(3) * isp2(3) * d1(:,3) + q(3) * isq2(3) * d2(:,3)
                                     
  !! Return the covariance matrix
  Cov = C - D
  
  RETURN    
  
END SUBROUTINE GetCovST

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetVarT(V, pq_in, p, q, p2, q2, LE)
!! Returns the variance of a chisquare vector with pq being the cell 
!! probabilities and p, q the marginals. p2, q2 contain the standardizing
!! probabilities.
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)        :: V(0:,0:)
  REAL(dpp), INTENT(IN)         :: pq_in(0:,0:), p(0:), q(0:), p2(0:), q2(0:)
  LOGICAL, INTENT(IN), OPTIONAL :: LE
  REAL(dpp)                     :: pq(0:2,0:2)
  REAL(dpd)                     :: p22(0:SIZE(p2)-1), q22(0:SIZE(q2)-1)
  INTEGER                       :: k, l, r, s, a, b

  !! Take observed cell probabilities LE assumption not wanted by user
  !pq = TRANSPOSE(RESHAPE(pq_in, (/3,3/)))
  pq = pq_in
  !! Take products of marginal probabilities if LE assumption wanted by user
  IF(PRESENT(LE)) THEN
    IF(LE) pq = MATMUL(RESHAPE(p, (/3,1/)), RESHAPE(q, (/1,3/)))   
  ENDIF
  
  !! Save variables as plain double precision
  p22 = DBLE(p2)
  q22 = DBLE(q2)

  !! Compute elements of the variance matrix
  DO k=0,2; DO l=0,2; DO r=0,2; DO s=0,2
  
    a = 3*k+l
    b = 3*r+s
    
    IF(k==r .AND. l==s) &
      V(a,b) = (pq(k,l)*(1-pq(k,l)) - 2*p(k)*pq(k,l)*(1-q(l)) &
                        - 2*q(l)*pq(k,l)*(1-p(k)) &
                        + 2*p(k)*q(l)*(pq(k,l)-p(k)*q(l)) & 
                        + p(k)*p(k)*q(l)*(1-q(l)) &
                        + q(l)*q(l)*p(k)*(1-p(k))) 
    IF(k==r .AND. l/=s) &
      V(a,b) = (-pq(k,l)*pq(k,s) + p(k)*pq(k,s)*q(l) &
                        + p(k)*pq(k,l)*q(s) - q(l)*pq(k,s)*(1-p(k)) &
                        - q(s)*pq(k,l)*(1-p(k)) &
                        + p(k)*q(l)*(pq(k,s)-p(k)*q(s)) &
                        + p(k)*q(s)*(pq(k,l)-p(k)*q(l)) &
                        - p(k)*p(k)*q(l)*q(s) + q(l)*q(s)*p(k)*(1-p(k))) 
    IF(k/=r .AND. l==s) &
      V(a,b) = (-pq(k,l)*pq(r,l) - p(k)*pq(r,l)*(1-q(l)) &
                        - p(r)*pq(k,l)*(1-q(l)) + p(k)*pq(r,l)*q(l) &
                        + p(r)*pq(k,l)*q(l) + p(r)*q(l)*(pq(k,l)-p(k)*q(l))&
                        + p(k)*q(l)*(pq(r,l)-p(r)*q(l)) &
                        + p(k)*p(r)*q(l)*(1-q(l)) - q(l)*q(l)*p(k)*p(r)) 
    IF(k/=r .AND. l/=s) &
      V(a,b) = (-pq(k,l)*pq(r,s) + 2*p(k)*pq(r,s)*q(l) &
                        + 2*p(r)*pq(k,l)*q(s) &
                        + p(r)*q(l)*(pq(k,s)-p(k)*q(s)) &
                        + p(k)*q(s)*(pq(r,l)-p(r)*q(l)) &
                        - 2*p(k)*p(r)*q(l)*q(s) )
                        
    V(a,b) = V(a,b) * Invert( REAL(SQRT(p22(k)*p22(r)*q22(l)*q22(s)), dpp) ) 

  ENDDO; ENDDO; ENDDO; ENDDO
  
  RETURN

END SUBROUTINE GetVarT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetARM(A, Fish, ierr, tol)
!! Computes the asymptotic representation matrix (ARM) A
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: Fish(:,:)
  REAL(dpp), INTENT(OUT)           :: A(:,:)
  INTEGER, INTENT(INOUT), OPTIONAL :: ierr
  REAL(dpp), INTENT(IN), OPTIONAL  :: tol
  REAL(dpp)                        :: Amat(dmF,dmF), Pi(dmF,dmF), P(dmF,dmF), &
                                      FishSq(dmF,dmF), FishSqI(dmF,dmF), &
                                      FishI(dmF,dmF), J0(dmF,dmF), v(dmF), &
                                      JFJ(dmF,dmF), JFJ_I(dmF,dmF), tol1
  !INTEGER :: i
  
  tol1 = epstol
  IF(PRESENT(tol)) tol1 = tol

  !! Create base matrix for H_{n,0} (the null hypothesis parameter space)
  v = one; v(dmF) = zero
  J0 = Diag(v)

  !! Get the product JFJ where F is the Fisher matrix
  JFJ = MATMUL(MATMUL(J0, Fish), J0)

  !! Compute pseudoinverse to JFJ
  CALL PseudoInv(JFJ, I=JFJ_I, tol=tol1, ierr=ierr)
  
  !! Check for error
  IF(PRESENT(ierr)) THEN
    IF(ierr/=0) RETURN
  ENDIF

  !! Get square root of Fish and square root of Fish^{-1}
  CALL PseudoInv(Fish, I=FishI, Sq=FishSq, SqI=FishSqI, tol=tol1, ierr=ierr)

  !! Check for error
  IF(PRESENT(ierr)) THEN
    IF(ierr/=0) RETURN
  ENDIF

  !! Compute Pi as a product of FishSq, J0 and the pseudoinverse
  P = MATMUL(FishSq, J0)
  Pi = MATMUL(MATMUL(P, JFJ_I), TRANSPOSE(P))

  !! Get matrix A that gives the of score statistics with estimated parameters
  Amat = MATMUL(MATMUL(FishSq, Id(4) - Pi), FishSqI)
  
  !! Return only the 4th row of Amat
  A(1,:) = Amat(4,:)

  RETURN

END SUBROUTINE GetARM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetFisher(Fish, Psi, XY, pi, tol)
!! Routine to compute the Fisher information matrix
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)          :: Fish(:,:)
  REAL(dpp), INTENT(IN)           :: Psi(:), XY(:,:)
  REAL(dpp), INTENT(IN), OPTIONAL :: pi(:), tol
  REAL(dpp)                       :: XY_k(SIZE(XY,2),1), tol1, pi1(SIZE(Psi)), &
                                     Psi2(SIZE(Psi))
  INTEGER                         :: k, N
  
  !! Store sample size locally
  N = SIZE(Psi)
  
  !! Check for proper size of pi and if OK assign it to 
  IF(PRESENT(pi)) THEN
    IF(SIZE(pi)/=N) CALL PrntE("Wrong size of pi (GetFisher)!")
    pi1 = pi
  ELSE
    pi1 = one / N
  ENDIF
  
  !! Set tolerance
  tol1 = EPSILON(zero)
  IF(PRESENT(tol)) tol1 = tol
  
  !! Compute the Fisher matrix as a sum of Psi' and z*z' over n
  Fish = zero
  Psi2 = Psi*(one-Psi)
  DO k=1,N
    XY_k(:,1) = XY(k,:)   !! k-th row of XY 
    Fish = Fish + Psi2(k) * MATMUL(XY_k, TRANSPOSE(XY_k)) * pi1(k)
    !Fish = Fish + Psi(k)*(one-Psi(k))*MATMUL(XY_k,TRANSPOSE(XY_k))*pi1(k)
    !Fish(1,1) = Fish(1,1) + Psi(k)*(1-Psi(k))*pi1(k)
  ENDDO

  !! Check for symmetry of Fish
  IF(ABS(SUM(Fish-TRANSPOSE(Fish))) > tol1) &
    CALL PrntE("Fish matrix is not symmetric (GetFisher).", Q=.TRUE.)

  RETURN
  
END SUBROUTINE GetFisher

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetKL(KL, PsiKL, model, beta, x, y)
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT)           :: KL(:,:), PsiKL(:)
  INTEGER, INTENT(IN)              :: model
  REAL(dpp), INTENT(IN)            :: beta(:)
  REAL(dpp), INTENT(IN), OPTIONAL  :: x(:), y(:)
  INTEGER                          :: i, j

  !! Check for proper sample sizes
  IF(SIZE(beta)>dmF) CALL PrntE("beta too long (GetKL).", Q=.TRUE.)
  IF(SIZE(KL,2)<dmF) CALL PrntE("KL too small (GetKL).", Q=.TRUE.)
  
  !! Assign values to KL all the way up to 4 columns even if beta is shorter, it
  !! is used elsewhere
  KL(:,1) = one 
  IF(PRESENT(x) .AND. PRESENT(y)) THEN
    KL(:,2) = REAL(x, dpp) 
    KL(:,3) = REAL(y, dpp) 
  ELSE
    IF(SIZE(KL,1)/=9) &
      CALL PrntE("Wrong dimensions of KL (GetKL).", Q=.TRUE.)
    KL(:,2) = REAL((/ (i,i,i, i=0,2) /), dpp) 
    KL(:,3) = REAL((/ ((i, i=0,2), j=0,2) /), dpp)  
  ENDIF
  KL(:,4) = Interaction(KL(:,2), KL(:,3), model)
  
  !! Compute the values of Psi at KL if PsiKL and beta present
  DO j=1,SIZE(KL,1)
    PsiKL(j) = Psi( SUM( beta * KL(j,1:SIZE(beta))) )
  ENDDO
  
  RETURN
END SUBROUTINE GetKL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Interaction(x, y, model) RESULT(xy)
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: x(:), y(:)
  INTEGER, INTENT(IN)   :: model
  REAL(dpp)             :: xy(SIZE(x)), coef(3,3)
  INTEGER               :: i, ix(SIZE(x)), iy(SIZE(y)) 
  
  !! Get the interaction table
  CALL GetInteractionTable(coef, model)

  !! Round x and y to integers
  ix = INT(x + half) + 1 
  iy = INT(y + half) + 1 
  
  !! Extract interaction coefficients
  DO i=1,SIZE(x)
    xy(i) = coef(ix(i),iy(i))
  ENDDO

  RETURN

END FUNCTION Interaction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION Psi(x)
!! Returns the value of the logistic function 1/(1+exp(-x))
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: x
  REAL(dpp)             :: Psi

  Psi = one / ( one + EXP(-x) )
  
  RETURN

END FUNCTION Psi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION PowerFun(plevel)
  IMPLICIT NONE
  INTEGER, PARAMETER    :: df1 = 4, df2 = 4
  REAL(dpd), INTENT(IN) :: plevel
  REAL(dpd)             :: PowerFun
  REAL(dpp)             :: plevel1, alpha2, QCH, QN, pow1, pow2a, pow2b, &
                           corlevel
  INTEGER               :: sts

  alpha2 = PF_level
  sts = 0
  plevel1 = REAL(plevel, dpp)
  
  !! Value -NAneg will be returned if error occurs
  PowerFun = NApos
  PF_level_unreliable = .TRUE.
  
  !! If plevel1 not positive, power is 0 (nothing can be rejected)
  IF(plevel1 < zero) CALL PrntE("Negative level in PowerFun!", Q=.TRUE.)
  IF(plevel1 == zero) THEN
    PowerFun = zero_d
    PF_level_unreliable = .FALSE.
    RETURN
  ENDIF
  
  !! Power function of the pretest (chisquare)
  !! If pretest level is one, no pretest is actually done and every pair goes
  !! trough the first phase, thus pow1 = 1. Otherwise compute it from the
  !! non-central chisquare distribution. 
  IF(plevel1==one) THEN
    pow1 = one
  ELSE
    QCH = qchisq(one - plevel1, df1, sts, silent=.TRUE.)
    IF(sts/=0) RETURN
    pow1 = one - pchisq(QCH, df2, PF_lambda, sts)
    IF(sts/=0) RETURN
  ENDIF

  !! Power function of the score test (normal)
  corlevel = alpha2 / MAX(PF_K*plevel1, one)
  IF(PF_two_sided) corlevel = corlevel / two

  QN = qnorm(corlevel)
  pow2a = pnorm(PF_beta*PF_slope + QN)

  pow2b = zero
  IF(PF_two_sided) pow2b = pnorm(-PF_beta*PF_slope + QN)

  !! Return the product of the two powers and change sign of the output, which 
  !! is then minimized by FindMin to find the maximum power
  PowerFun = - DBLE(pow1 * (pow2a + pow2b))
  PF_level_unreliable = .FALSE.
  
  RETURN

END FUNCTION PowerFun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_TESTPROC
