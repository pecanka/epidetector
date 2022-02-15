MODULE EPI_DATAIO

  USE EPI_MATH                !! MATHEMATICAL TOOLS
  USE EPI_INIT                !! DEFINITION AND INITIALIZATION OF VARIABLES
  USE EPI_SIMUL               !! DATA SIMULATION
  USE EPI_TESTPROC            !! TESTING ROUTINES (Needed for GetBeta and GetSubsample)
  USE EPI_REPORT              !! PRINTING OF REPORTS

  IMPLICIT NONE

  !! Default Accessibility
  PUBLIC

  PRIVATE :: LH, out_cnams, out_ncols, get_new_data, get_new_map

  CHARACTER(mstl), PARAMETER  :: analyzing = "analyzing files ..."
  INTEGER(1), PARAMETER       :: bed_byte1 = 108, bed_byte2 = 27, NAbed = 85
  
  CHARACTER(mstl), SAVE       :: out_cnams(1000)
  INTEGER, SAVE               :: out_ncols = -1, &
                                 rec_size_tmp = -1, &
                                 rec_size_out = -1, &
                                 out_cycle = 50000
  LOGICAL                     :: get_new_data, get_new_map
  
  !! Synchronizes minimalistic and regular temporary records
  !INTERFACE SyncTmp
  !  MODULE PROCEDURE :: SyncTmpFromTmpmin, SyncTmpToTmpmin
  !END INTERFACE 
  
  CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
!! ************************************************************************ !!  
!!                        PARAMETER DETERMINATION                           !!
!! ************************************************************************ !!  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SetDeltas(WT1, dAS, dDS, cgrp, dAS_co, dAS_ca, dDS_co, dDS_ca)
  IMPLICIT NONE  
  INTEGER, INTENT(IN)      :: WT1, cgrp
  REAL(dpp), INTENT(IN)    :: dAS, dDS
  REAL(dpp), INTENT(INOUT) :: dAS_co, dAS_ca, dDS_co, dDS_ca
  
  IF(WT1==0) THEN
  
    dAS_co = zero
    dAS_ca = zero

    dDS_co = zero
    dDS_ca = zero
    
  ELSEIF(WT1==T1co) THEN
    
    dAS_co = dAS
    dAS_ca = zero
    IF(cgrp==2) dAS_ca = dAS

    dDS_co = dDS
    dDS_ca = zero
  
  ELSEIF(WT1==T1ca) THEN
    
    dAS_ca = dAS
    dAS_co = zero
    IF(cgrp==2) dAS_co = dAS

    dDS_ca = dDS
    dDS_co = zero
    
  ELSEIF(WT1==T1cc) THEN
    
    dAS_co = dAS
    dAS_ca = dAS

    dDS_co = dDS
    dDS_ca = dDS
  
  ELSEIF(WT1==T1po) THEN
  
    dAS_co = dAS
    dAS_ca = dAS

    dDS_co = dDS
    dDS_ca = dDS
  
  ENDIF
  
  RETURN
  
END SUBROUTINE SetDeltas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetAutoDelta(delta, level1, simple_delta, simple_delta_frac, &
    min_delta, max_delta, ndelta, model, WT1, T1_df, OR, OR1, OR2, &
    prev, ntests, nco, nca, maf_fixed, maf1, maf2, min_level1, max_level1, &
    nlevel1, level2, K, correct_all, correction, min_count, min_ca, min_co) 
  IMPLICIT NONE
  REAL(dpp), PARAMETER            :: minimum_delta = 0.01
  REAL(dpp), INTENT(OUT)          :: delta, level1
  REAL(dpp), INTENT(IN)           :: min_delta, max_delta
  INTEGER, INTENT(IN)             :: ndelta, model, WT1, T1_df, nco, &
                                     nca, ntests, nlevel1
  REAL(dpp), INTENT(IN)           :: simple_delta_frac, OR, OR1, OR2, prev, &
                                     maf1, maf2, min_level1, max_level1, &
                                     level2, K
  LOGICAL, INTENT(IN)             :: simple_delta, maf_fixed 
  LOGICAL, INTENT(IN), OPTIONAL   :: correct_all 
  REAL(dpp), INTENT(IN), OPTIONAL :: correction, min_count, min_ca, min_co  
  CHARACTER, ALLOCATABLE          :: FAM(:,:)
  INTEGER(iks), ALLOCATABLE       :: X(:,:), sts(:), sex(:)
  INTEGER(iks)                    :: sts0(nco+nca), x0(nco+nca), y0(nco+nca), & 
                                     sts1(nco+nca), x1(nco+nca), y1(nco+nca), &
                                     sts2(nco+nca), x2(nco+nca), y2(nco+nca)
  INTEGER                         :: i, n, nd, id, ia, ierr0, ierrPT4, ierrPT1, &
                                     ierrDS, ierrCS, nv, nv1, nv2, wm(2), NR1, &
                                     NT1, NR2(ndelta,nlevel1), NT2(ndelta), nrep
  INTEGER(ikb)                    :: pct
  REAL(dpp)                       :: deltas(ndelta), delta1, a1, a2, adif, &
                                     alphas(nlevel1), tmp, tmp1(1), &
                                     c_co1(dmS,dmS), c_co2(dmS,dmS), &
                                     c_ca1(dmS,dmS), c_ca2(dmS,dmS), &
                                     PT(ntests), PTp(ntests), DS(ntests), &
                                     DSp(ntests), CS(ntests), CSp(ntests), &
                                     PT1, PT1p, RR1, RR2(ndelta,nlevel1)
  CHARACTER(mltl)                 :: text, stext
  LOGICAL                         :: comeback
  
  !! If delta should only be a simple excess, get it and return 
  IF(simple_delta) THEN
  
    CALL Prnt("Delta will be determined as a simple excess ratio"//&
                   " (fractioned by "//TRIM(r2c(simple_delta_frac))//").")
    !! Calculate delta differently for different pretests
    SELECT CASE(WT1)
      CASE(T1co)                      
        delta = REAL(nco - nca, dpp) * simple_delta_frac / nco
      CASE(T1ca)
        delta = REAL(nca - nco, dpp) * simple_delta_frac / nca
      CASE DEFAULT
        delta = def_delta
    END SELECT
    
    !! Check for too small delta
    delta = MAX(delta, minimum_delta) 
      
    text = "Automatically determined delta: "//r2c(delta)
    CALL Prnt(text)
    
    RETURN
  ENDIF
  
  !! Annouce optimization of delta
  CALL Prnt("Pretest sample size ratio delta will be optimized ...", skip1=1)
  
  !! Store number of delta candidates locally
  nd = ndelta
  !! If all deltas are the same, do just one of them
  IF(max_delta==min_delta) nd = 1
  !! Calculate candidates for delta based on user input
  IF(nd==1) THEN
    deltas = (min_delta + max_delta) / two
  ELSE
    deltas = min_delta + (/(i,i=0,nd-1)/) * (max_delta-min_delta)/(nd-1)
  ENDIF
  
  !! Calculate candidates for pretest level
  a1 = LOG10(min_level1)
  a2 = LOG10(max_level1)
  IF(a1>a2) THEN
    tmp = a1; a1 = a2; a2 = tmp
  ENDIF
  adif = (a2 - a1) / nlevel1
  alphas = ten**(a2 - (/(i,i=0,nlevel1-1)/)*adif)
  
  !! Announce candidate deltas
  text = "Optimization will be done assuming a total of "//TRIM(r2c(K))//" tests"
  CALL Prnt(text)
  IF(nd==1) THEN
    text = TRIM(i2c(nd))//" considered value of delta is "//&
           TRIM(r2c(deltas(1)))//" ..."
  ELSE
    text = TRIM(i2c(nd))//" considered values of delta range between "//&
           TRIM(r2c(deltas(1)))//" and "//TRIM(r2c(deltas(nd)))//" ..."
  ENDIF
  CALL Prnt(text)

  !! Announce candidate levels
  IF(nd==1) THEN
    text = TRIM(i2c(nlevel1))//" considered PHASE1 levels range between "//&
           TRIM(r2c(min_level1))//" and "//TRIM(r2c(max_level1))//"."
  ELSE
    text = TRIM(i2c(nlevel1))//" considered PHASE1 level is "//&
           TRIM(r2c(alphas(1)))//"."
  ENDIF
  CALL Prnt(text)

  !! Announce simulation settings
  text = "Optimization simulation settings: "//TRIM(i2c(ntests))//" tests"//&
         " with model "//TRIM(GetModelName(model))//", prevalence "//&
         TRIM(r2c(prev))//", OR1="//TRIM(r2c(OR1))//", OR2="//TRIM(r2c(OR2))//&
         " and OR="//TRIM(r2c(OR))//" with maf1="//TRIM(r2c(maf1))//&
         " and maf2="//TRIM(r2c(maf1))
  IF(maf_fixed) text = TRIM(text)//" (fixed)"
  IF(.NOT.maf_fixed) text = TRIM(text)//" (lower bound)"
  CALL Prnt(text)

  !! Nulify rejection and progress counters
  NR1 = 0
  NR2 = 0
  NT1 = 0
  NT2 = 0

  !! Point of repeat return
  nrep = 0
  comeback = .FALSE.
  17 CONTINUE
  
  !! Nulify statistic and p-value vectors
  PT = zero
  CS = zero
  DS = zero
  PTp = NAp
  CSp = NAp
  DSp = NAp
  
  !! Increase counter of repeats
  nrep = nrep + 1

  !! Simulate data
  CALL SimulData(X, sts, sex, FAM, .TRUE., model, OR, OR1, OR2, &
                 prev, LD, 2, ntests, nca, nco, maf_fixed, maf1, maf2)
                 
  n = nco+nca

  !! Loop over delta candidates
  DO id=1,nd

    delta1 = deltas(id)
    pct = -1
    stext = "Running optimization for delta equal to "//TRIM(r2c(delta1))//" ..."
    IF(.NOT.countdown) CALL Prnt(stext, advance='NO')
    
    !! Loop over simulated loci
    DO i=1,ntests
    
      !! Report progress
      CALL PrintCD(text, INT(i,ikb), INT(ntests,ikb), pct, finish=.TRUE., &
                   stext=stext)

      !! Recall original genotypes and status data
      CALL Bed2Ped(X(((i-1)*n+3)/4+1:(i*n+3)/4,1), x0)
      CALL Bed2Ped(X(((i-1)*n+3)/4+1:(i*n+3)/4,2), y0)
      sts0 = sts((i-1)*n+1:i*n)
      
      !! Create pre-test/post-test division, counts etc.
      CALL ScanData(n, sts0, x0, y0, nv, ierr0, delta1, delta1, WT1, sts1, x1, &
                    y1, nv1, sts2, x2, y2, nv2, c_co1, c_co2, c_ca1, c_ca2, &
                    correct_all, correction, min_count, min_ca, min_co)

      !!!!!!!!!!!!!!!!!!!!!!!!
      !!  SINGLE-STEP TEST  !!
      !!!!!!!!!!!!!!!!!!!!!!!!
      IF(id==1) THEN
        
        !! Calculate classical score
        CALL TestS2(sts0(1:nv), x0(1:nv), y0(1:nv), nv, nv, model, CS(i), &
                    CSp(i), ierrCS, var_use_beta1)

        !! Check for an error during pretest
        IF(ierrCS/=0) THEN
          CSp(i) = NAp
          CYCLE
        ENDIF

        !! If no errors, increase test counter
        NT1 = NT1 + 1

        !! Check for rejection by one-step test
        IF(CSp(i)<=level2/K) NR1 = NR1+1
        
      ENDIF
      
      !!!!!!!!!!!!!!!!!!!!!!!!
      !!    TWO-STEP TEST   !!
      !!!!!!!!!!!!!!!!!!!!!!!!
      !! For delta equal to 0 no need to redo the testing
      IF(delta1==zero) THEN
        PT(i) = zero
        PTp(i) = zero
        DS(i) = CS(i)
        DSp(i) = CSp(i)
      ELSE
                      
        !! Call phase 1 tests
        CALL TestS1(c_co1, c_ca1, c_co2, c_ca2, WT1, PT(i), PTp(i), &
                    ierrT4=ierrPT4, T4df=T1_df, T1=PT1, T1p=PT1p, &
                    ierrT1=ierrPT1, doDS=.TRUE., var_grouped=.FALSE.)
        
        !! Check for an error during pretest
        IF(ierrPT4/=0) THEN
          PTp(i) = NAp
          CYCLE
        ENDIF
        
        !! Calculate disjoint score
        CALL TestS2(sts2(1:nv2), x2(1:nv2), y2(1:nv2), nv2, nv2, model, DS(i), &
                    DSp(i), ierrDS, var_use_beta1)

        !! Check for an error during posttest
        IF(ierrDS/=0) THEN
          PTp(i) = NAp
          DSp(i) = NAp
          CYCLE
        ENDIF
        
        !! If no errors, increase test counter
        NT2(id) = NT2(id) + 1

      ENDIF

    ENDDO

    IF(.NOT.countdown) CALL Prnt0("")

    !! Optimize two-step over levels in alphas
    IF(delta1==zero) THEN
      NR2(id,:) = COUNT(DSp<=level2/K)
    ELSE
      DO ia=1,nlevel1
        NR2(id,ia) = COUNT(PTp<=alphas(ia) .AND. DSp<=level2/MAX(one,alphas(ia)*K))
      ENDDO
    ENDIF
  
  ENDDO
  
  CALL Prnt0(" ** RESULTS **", skip1=1)

  !! Find the maximizing combination of delta and level1
  wm = MAXLOC(NR2)
  delta = deltas(wm(1))
  level1 = alphas(wm(2))

  !! Announce results for one-step tests
  RR1 = zero
  IF(NT1>0) RR1 = REAL(NR1, dpp) / NT1
  text = "Rejection ratio by CS (rejections:tests): "//TRIM(r2c(RR1))//" ("//&
         TRIM(i2cp(NR1))//":"//TRIM(i2cp(NT1))//")"
  CALL Prnt(text)

  !! Announce results for two-step tests
  IF(comeback) THEN
    IF(NT2(1)>0) RR2(1,:) = REAL(NR2(1,:), dpp) / NT2(1)
    text = "Maximum rejection ratio by PT&DS for delta="//TRIM(r2c(delta))//&
           " (rejections:tests): "//TRIM(r2c(RR2(1,wm(2))))//" ("//&
           TRIM(i2cp(NR2(1,wm(2))))//":"//TRIM(i2cp(NT2(1)))//")"
    CALL Prnt(text)
    IF(ANY(NR2>0)) THEN
      tmp1 = alphas(MAXLOC(RR2(1,:)))
      text = "Corresponding maximizing level: "//r2c(tmp1(1))
      CALL Prnt(text)
    ENDIF
  ELSE
    text = "Maximum rejection ratios by PT&DS for increasing delta"//&
           " (rejections:tests):"
    RR2 = zero
    DO id=1,nd
      IF(NT2(id)>0) RR2(id,:) = REAL(NR2(id,:), dpp) / NT2(id)
      text = TRIM(text)//" "//TRIM(r2c(MAXVAL(RR2(id,:))))//" ("//&
             TRIM(i2cp(MAXVAL(NR2(id,:))))//":"//TRIM(i2cp(NT2(id)))//")"
    ENDDO 
    CALL Prnt(text)
    IF(ANY(NR2>0)) THEN
      text = "Corresponding maximizing levels:"
      DO id=1,nd
        tmp1 = alphas(MAXLOC(RR2(id,:))) 
        text = TRIM(text)//" "//TRIM(r2c(tmp1(1)))
      ENDDO 
      CALL Prnt(text)
    ENDIF
  ENDIF
  
  !! Check for absolutely no rejections
  IF(ALL(NR2==0)) THEN
    IF(nrep < def_auto_d_maxnreps) THEN
      text = "There was no rejection by DS or CS. I'll try a few more tests ..." 
      CALL PrntW(text)
      GOTO 17
    ELSE
      IF(.NOT.comeback) THEN
        IF(WT1==T1co) delta = MAX(zero, REAL(nco - nca, dpp) / nco)
        IF(WT1==T1ca) delta = MAX(zero, REAL(nca - nco, dpp) / nca)
        IF(WT1==T1cc) delta = half
        IF(WT1==T1po) delta = one
      ENDIF
      text = "There was still no rejection by DS or CS even after "//&
             TRIM(i2c(nrep*ntests))//" tests. Giving up on finding optimal"//&
             " PHASE1 sample size ratio and using value "//TRIM(r2c(delta))//"."//&
             " It is recommended to try different input parameters."
      CALL PrntW(text)
      level1 = one
    ENDIF
  !! Check for no difference in rejection counts
  ELSEIF(ALL(NR2==MINVAL(NR2))) THEN
    text = "All considered PHASE1 sample size ratios and levels yield the"//&
           " same number of rejections "//TRIM(i2c(NR2(1,1)))//". It is"//&
           " probably a good idea to try different input parameters."
    CALL PrntW(text, skip1=1, skip2=1)
  ENDIF

  !! Announce "optimal" delta
  text = "Optimized PHASE1 sample size ratio: "//TRIM(r2c(delta))
  IF(delta==MINVAL(deltas)) text = TRIM(text)//" (smallest considered value)"
  IF(delta==MAXVAL(deltas)) text = TRIM(text)//" (biggest considered value)"
  CALL Prnt(text)
  
  !! Check for delta close to zero
  IF(delta==zero .AND. .NOT.comeback) THEN
    IF(WT1==T1co) deltas = MAX(zero, REAL(nco - nca, dpp) / nco)
    IF(WT1==T1ca) deltas = MAX(zero, REAL(nca - nco, dpp) / nca)
    CALL PrntW("Replacing PHASE1 ratio 0.0 with the 'balancing' value "//TRIM(r2c(delta))//"!")
    NR1 = 0
    NR2 = 0
    comeback = .TRUE.
    nd = 1
    GOTO 17
  ENDIF
  
  !! Announce "optimal" pretest level for that delta, unless delta is zero, in
  !! which case make level1 equal to one
  IF(delta==zero) THEN
    level1 = one
  ELSE
    text = "Optimized PHASE1 level: "//TRIM(r2c(level1))
    IF(level1==MINVAL(alphas)) text = TRIM(text)//" (smallest considered value)"
    IF(level1==MAXVAL(alphas)) text = TRIM(text)//" (biggest considered value)"
    CALL Prnt(text)
  ENDIF
  CALL Prnt0("")
  
  RETURN
  
END SUBROUTINE GetAutoDelta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!! ************************************************************************ !!  
!!                                DATA INPUT                                !!
!! ************************************************************************ !!  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE RememberData()
!! This subroutine saves the input data from the previous run to be reused later 
  IMPLICIT NONE

    !! Save values
    prev_ped_nloci      = ped_nloci
    prev_ped_nsamp   = ped_nsamp
    prev_ped_ss = ped_ss
    prev_ped_nco  = ped_nco
    prev_ped_nca     = ped_nca
    prev_ped_nmale     = ped_nmale
    prev_ped_nfema   = ped_nfema
    
    !! Do not remember file info if data are to be simulated
    IF(simulate_data) THEN
      IF(ALLOCATED(prev_ped_file))    DEALLOCATE(prev_ped_file)
      IF(ALLOCATED(prev_bed_file))    DEALLOCATE(prev_bed_file)
      !IF(ALLOCATED(prev_ped_ncols))   DEALLOCATE(prev_ped_ncols)
      IF(ALLOCATED(prev_map_file))    DEALLOCATE(prev_map_file)
      IF(ALLOCATED(prev_sub_file)) DEALLOCATE(prev_sub_file)
    ENDIF
    
    !! Remember file info if data are to be read from PED file
    IF(ALLOCATED(ped_file)) THEN    
      CALL ResizeVar(prev_ped_file, SIZE(ped_file))
      prev_ped_file = ped_file
      !IF(ALLOCATED(ped_ncols)) THEN
      !  CALL ResizeVar(prev_ped_ncols, SIZE(ped_ncols))
      !  prev_ped_ncols = ped_ncols
      !ENDIF
    ENDIF

    !! Remember file info if data are to be read from BED file
    IF(ALLOCATED(bed_file)) THEN    
      CALL ResizeVar(prev_bed_file, SIZE(bed_file))
      prev_bed_file = bed_file
      !IF(ALLOCATED(ped_ncols)) THEN
      !  CALL ResizeVar(prev_ped_ncols, SIZE(ped_ncols))
      !  prev_ped_ncols = ped_ncols
      !ENDIF
    ENDIF
    !! Remember file info if data are to be read from MAP file
    IF(ALLOCATED(map_file)) THEN    
      CALL ResizeVar(prev_map_file, SIZE(map_file))
      prev_map_file = map_file
    ENDIF
    !! Remember file info if data are to be read from SUBMAP file
    IF(ALLOCATED(sub_file)) THEN    
      CALL ResizeVar(prev_sub_file, SIZE(sub_file))
      prev_sub_file = sub_file
    ENDIF
    
    RETURN
    
END SUBROUTINE RememberData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE CheckForSameInputFiles()
!! Checks if the current input (binary) ped files have the same names as
!! in previous run and if they are the old data will be reused
  IMPLICIT NONE

  !! If simulated data will not be saved, don't bother saving family info 
  IF(.NOT.save_input_data .AND. simulate_data) fam_ncols = 0

  !! Decide whether to re-read input ped files
  get_new_data = .TRUE.
  IF(ALLOCATED(prev_ped_file) .AND. ALLOCATED(ped_file)) THEN
    IF(SIZE(prev_ped_file) == SIZE(ped_file) .AND. &
    ALL(prev_ped_file == ped_file) .AND. ALL(ped_file/="")) &
      get_new_data = .FALSE.
  ENDIF
  IF(ALLOCATED(prev_bed_file) .AND. ALLOCATED(bed_file)) THEN
    IF(SIZE(prev_bed_file) == SIZE(bed_file) .AND. &
    ALL(prev_bed_file == bed_file) .AND. ALL(bed_file/="")) &
      get_new_data = .FALSE.
  ENDIF

  !! Assume mapping info will be re-read
  get_new_map = .TRUE.

  !! Decide whether to re-read input mapping/submapping files and redo exclusion
  IF(ALLOCATED(prev_map_file) .AND. ALLOCATED(map_file) .AND. &
  ALLOCATED(prev_sub_file) .AND. ALLOCATED(sub_file)) THEN
    IF(SIZE(prev_map_file) == SIZE(map_file) .AND. &
    SIZE(prev_sub_file) == SIZE(sub_file) .AND. &
    ALL(prev_map_file == map_file) .AND. ALL(map_file/="") .AND. &
    ALL(prev_sub_file == sub_file) .AND. ALL(sub_file/="")) &
      get_new_map = .FALSE.
  ENDIF
  
  !! If only tmp to out conversion is to happen, don't get new data
  IF(tmp_to_out_only) THEN
    get_new_data = .FALSE.
    get_new_map = .FALSE.
  ENDIF

  !! If data were simulated last run and should be reused, do so
  IF(simulate_data .AND. reuse_data .AND. ALLOCATED(X)) get_new_data = .FALSE.
  
  !! Nulify the counter of reusing data
  IF(get_new_data) THEN
    ped_nrepeat = 0
  !! Otherwise, If data should be reused, restore the variables
  ELSE  
    ped_nsamp   = prev_ped_nsamp
    ped_ss = prev_ped_ss
    ped_nco  = prev_ped_nco
    ped_nca     = prev_ped_nca
    ped_nmale     = prev_ped_nmale
    ped_nfema   = prev_ped_nfema
    ped_nloci      = prev_ped_nloci
    !IF(ALLOCATED(prev_ped_ncols)) THEN
    !  CALL ResizeVar(ped_ncols, SIZE(prev_ped_ncols))
    !  ped_ncols = prev_ped_ncols
    !ENDIF
    ped_nrepeat = ped_nrepeat + 1
  ENDIF
  
  RETURN
  
END SUBROUTINE CheckForSameInputFiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrepareDataVariables(exit_code)
!! Prepares data variables for input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_code

  !! If data should be simulated, do not try to read it from a file, 
  !! and set variables to proper values
  IF(simulate_data) THEN

    IF(get_new_data) THEN
      IF(ped_ss<=0) THEN
        CALL PrntE("Cannot simulate input data without user specified"//&
                     " sample size.")
        CALL Prnt("Hint: Use --samplesize or --ncontrols, --ncases.", &
                  Q=.TRUE., premature=.FALSE.)
      ENDIF
      ped_nloci = 2 + nneutralloci + 2*nLDpairs
    ENDIF

    !! If user gave just samplesize, split it half and half for controls/cases
    IF(ped_nca<0)    ped_nca = ped_ss/2
    IF(ped_nco<0) ped_nco = ped_ss - ped_ss/2
    test_same_chr_S1 = .TRUE.
    test_same_chr = .TRUE.

    !! Announce automatically determined allele frequency (in EPI_INITIALIZE)
    IF(auto_maf1 .OR. auto_maf2 .OR. auto_maf3) &
      CALL Prnt("Automatically determined minor allele frequency for"//&
                     " simulation is "//TRIM(r2c(auto_maf))//".")

  ELSEIF(get_new_data) THEN

    !! Announce problem with missing input files
    IF(ped_nfiles == 0 .AND. bed_nfiles == 0) THEN
      CALL PrntE("No (binary) ped file specified.", skip1=1)
      CALL Prnt("Use --ped or --bed to provide input ped file name.")
      premat_halt = .TRUE.
      exit_code = 3000
      RETURN
    ENDIF
    
    !! Check existence of all input files
    IF(ped_nfiles > 0) &
      CALL CheckFileExistence(fnlist=ped_file, msg="PED file")
    IF(bed_nfiles > 0) &
      CALL CheckFileExistence(fnlist=bed_file, msg="BED file")
    IF(map_nfiles > 0) &
      CALL CheckFileExistence(fnlist=map_file, msg="MAP file")

  ENDIF
  
  RETURN

END SUBROUTINE PrepareDataVariables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetData()
  IMPLICIT NONE
  REAL(dpp) :: level1, dAS_co, dAS_ca, dDS_co, dDS_ca
  LOGICAL   :: fam_ok

  !! Deallocate data arrays
  IF(get_new_data) THEN
    IF(ALLOCATED(X))   DEALLOCATE(X)
    IF(ALLOCATED(sts)) DEALLOCATE(sts)
    IF(ALLOCATED(sex)) DEALLOCATE(sex)
  ENDIF
  IF(get_new_map) THEN
    IF(ALLOCATED(MAPA))      DEALLOCATE(MAPA)
    IF(ALLOCATED(SUBMAPA))   DEALLOCATE(SUBMAPA)
    IF(ALLOCATED(MAPARS))    DEALLOCATE(MAPARS)
    IF(ALLOCATED(SUBMAPARS)) DEALLOCATE(SUBMAPARS)
  ENDIF
  
  !! Get the dimensions of ped files (number of lines, number of cols,...)
  !IF(get_new_data.AND..NOT.simulate_data.AND.ped_ss<=max_samplesize)THEN
  IF(get_new_data .AND. .NOT.simulate_data) THEN

    !! Get mapping info
    IF(map_nfiles>0) THEN
      IF(get_new_map) THEN
      
        !! Read the mapping file
        CALL ReadMAPfile(map_file, MAPA, MAPARS, map_ncols, map_nlines, &
                         map_nskip, map_cs)

        !! Read also the subset mapping info if file given
        IF(sub_nfiles>0) &
          CALL ReadMAPfile(sub_file, SUBMAPA, SUBMAPARS, sub_ncols, &
                           sub_nlines, sub_nskip, sub_cs, submap=.TRUE.)

      ELSE
        CALL Prnt("MAP file(s) will be kept from previous analysis.")
        IF(sub_nfiles>0) &
          CALL Prnt("SUBMAP file(s) will be kept from previous analysis.")
      ENDIF
    ENDIF

    !! Input is a binary ped file (mapping file has been read at this point,
    !! so we know what the sample size and the number of loci are)
    IF(bed_nfiles>0) THEN
    
      !! Read the family file information
      ped_nloci = SIZE(MAPA,1)
      CALL ReadFAMfile(fam_file, sts, sex, FAM, fam_ncols, fam_nskip, fam_cs, &
                       ped_nco, ped_nca, ped_nmale, ped_nfema, ped_nsamp)
      
      !! Calculate the total sample size (including possible NA records)
      IF(fam_single_sample) THEN
        ped_ss = SIZE(sts)
      ELSE
        ped_ss = SIZE(sts) / ped_nsamp
      ENDIF
      
      !! Check for problematic FAM file
      IF(fam_single_sample) THEN
        fam_ok = ped_ss == SIZE(sts)
      ELSE
        fam_ok = ped_ss * ped_nsamp == SIZE(sts)
      ENDIF 
      IF(.NOT.fam_ok) &
        CALL PrntE("FAM file seems inconsistent with "//TRIM(i2c(ped_nsamp))//&
                     " samples.", Q=.TRUE.)

      !! Check for missing cases and controls
      IF(ped_nca+ped_nco==0) &
        CALL PrntW("No valid status found in file "//TRIM(fam_file)//"!")
        
    ENDIF
    
  ENDIF
  
  !! Check for too big sample size
  !IF(ped_ss > max_samplesize) &
  !  CALL PrntE("Maximum sample size exceeded! Input sample size cannot"//&
  !               " be larger than "//i2cp(max_samplesize)//".",Q=.TRUE.)

  !********************************************************************!
  !!             GET NEW DATA FROM FILES OR SIMULATION                !!
  !********************************************************************!
  IF(get_new_data) THEN 

    !! Simulate input data and return binary formatted X 
    IF(simulate_data) THEN
    
      CALL Report(simul_input=.TRUE.)
      CALL SimulData(X, sts, sex, FAM, .TRUE., sim_model, OR, OR1, OR2, &
                     prevalence, LD, ped_nloci, ped_nsamp, ped_nca, ped_nco, &
                     nmale=ped_nmale, nfema=ped_nfema)
    
    ELSE
    
      !! Read in the plain text format ped data, return binary formatted X
      IF(ped_nfiles>0) THEN
        CALL ReadPEDfile(ped_file(1), uped, sts, sex, FAM, X, MAPA, ped_nloci, &
              ped_nskip, ped_cs, ped_nco, ped_nca, ped_nmale, ped_nfema)
        ped_ss = ped_nca + ped_nco                             
        IF(ped_ss==0) &
          CALL PrntW("No controls or cases found in file "//TRIM(fam_file)//"!")
      ENDIF
    
      !! Read binary ped data, return binary formatted X 
      IF(bed_nfiles>0) &
        CALL ReadBEDfile(bed_file, X, ped_ss, ped_nsamp, ped_nloci)
    
    ENDIF
    
  ENDIF

  !! If no mapping files given, define MAPA as if data from different files
  !! came from different chromosomes
  IF(map_nfiles == 0) THEN
  
    CALL SimulateMap(MAPA, ped_nloci, map_ncols, map_chcol, map_rscol, &
                     map_dscol, map_bpcol, map_a1col, map_a2col, &
                     binary=save_binary, ch_start=sim_chr_start, &
                     rs_start=sim_rs_start, sglchr=sim_single_chr)
    !CALL ResizeVar(MAPA, ped_nloci, map_ncols, "")
    !CALL SimulateMap(MAPA, ped_ncols, map_chcol, map_rscol, map_dscol, &
    !                 map_bpcol, map_a1col, map_a2col, bim=save_binary)
  
  ENDIF
  
  !! Set the pretest sample portion to an automatic value if selected
  IF(auto_delta) THEN

    !! Optimum delta for DS will be determined through simulation 
    IF(auto_d_optimize) THEN
      !! If oracle, then use the same setting as for simulation
      IF(auto_d_oracle) THEN
        auto_d_OR = OR
        auto_d_OR1 = OR1
        auto_d_OR2 = OR2
        auto_d_prev = prevalence
        auto_d_maf1 = allelefreqA
        auto_d_maf2 = allelefreqB
      ENDIF

      !! Determine optimal delta
      CALL GetAutoDelta(dDS, level1, simple_delta, auto_d_frac, auto_d_min, &
                        auto_d_max, auto_d_n, ana_model, WT1, T1_df, &
                        auto_d_OR, auto_d_OR1, auto_d_OR2, auto_d_prev, &
                        auto_d_ntests, ped_nco, ped_nca, auto_d_maf_fixed, &
                        auto_d_maf1, auto_d_maf2, auto_d_minlevel, &
                        auto_d_maxlevel, auto_d_nlevel, level_S2, auto_d_COR)
              
      !! If pretest level is to be automatically determined, use the value that
      !! is associated with the optimal delta 
      IF(auto_level_report) auto_level_optim = level1
                    
    ELSEIF(WT1==T1co .AND. ped_nco-ped_nca > min_ss) THEN
      dDS = round(REAL(ped_nco - ped_nca, dpp) / ped_nco, delta_precision)
    ELSEIF(WT1==T1ca .AND. ped_nca-ped_nco>min_ss) THEN
      dDS = round(REAL(ped_nca - ped_nco, dpp) / ped_nca, delta_precision) 
    ELSE
      auto_delta = .FALSE.
      auto_d_incapable = .TRUE.
      dDS = def_delta 
    ENDIF

    !! Use the same value for AS
    dAS = dDS
              
  ENDIF
  
  !! Allocate pretest subsamples pss_AS and pss_DS and set all values to -1
  CALL ResizeVar(pss_AS, ped_ss, allvalue=1_iks)
  CALL ResizeVar(pss_DS, ped_ss, allvalue=1_iks)
  
  !! Get the pretest selection status
  IF(fix_subsamp) THEN
   
    IF(pss_nfiles>0) THEN

      !! Read in the pretest selection status from
      CALL ReadPSSFile(pss_file, pss_AS, pss_nskip, dAS, pss=.TRUE.)
      dDS = dAS
      pss_DS = pss_AS
      
    ELSEIF(.NOT.no_testing) THEN

      !! Set the portions according to the selected pretest
      CALL SetDeltas(WT1, dAS, dDS, cntrgrp, dAS_co, dAS_ca, dDS_co, dDS_ca)
                           
      !! Generate random pretest subsample for AS
      CALL GetSubsample(pss_AS, sts(1:ped_ss), ped_ss, dAS_co, dAS_ca, announce=.TRUE.)
      
      !! Generate random pretest subsample for DS
      IF(dAS == dDS) THEN
        pss_DS = pss_AS
      ELSE
        CALL GetSubsample(pss_DS, sts(1:ped_ss), ped_ss, dDS_co, dDS_ca, announce=.TRUE.)
      ENDIF
      
    ENDIF

  ENDIF
  
  !! Read in case-control status from a separate file
  IF(sts_nfiles>0) CALL ReadPSSFile(sts_file, sts, sts_nskip)
    
  !! Check for conflict between pss and PHASE1 test
  IF(WT1==T1co .AND. COUNT(pss_AS==1_iks .AND. sts==NumCo)==0) &
    CALL PrntE("The selected PHASE1 tests require controls but none were selected.", Q=.TRUE.)
                             
  IF(WT1==T1ca .AND. COUNT(pss_AS==1_iks .AND. sts==NumCa)==0) &
    CALL PrntE("The selected PHASE1 tests require cases but none were selected.", Q=.TRUE.)

  RETURN
   
END SUBROUTINE GetData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ProcessData(exit_code)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_code
  INTEGER              :: nincl, n_incl_loci, i, ss
  LOGICAL, ALLOCATABLE :: is_co(:)

  !! Allocate InclLoc and set inclusion so that initially all loci included
  CALL ResizeVar(InclLoc, SIZE(X,2), .TRUE.)

  !! Exclude loci with excluded chromosome number, negative base-pair number
  !! and those that were not specified in subset mapping file
  IF(.NOT.simulate_data) &
    CALL CheckLociExcl(InclLoc, MAPA, SUBMAPA, MAPARS, SUBMAPARS, &
                       n_chr_excl, n_rs_excl, n_negbp_excl, n_sub_excl)
          
  !! Compute MAFs and exclude those loci that do not meet the MAF bound 
  IF(compute_maf) THEN
  
    CALL Prnt("Getting minor allele frequencies in the whole sample ...")
    CALL ComputeMAF(MAF, X, sex, MAPA, InclLoc, n_maf_excl, min_maf)
    
    !! Also get the MAFs in the two subsamples (cases/controls)
    IF(out_maf_coca) THEN
      
      CALL Prnt("Getting minor allele frequencies in the controls ...")
      IF(fam_single_sample) THEN
        ss = SIZE(sts)
      ELSE
        ss = SIZE(sts) / ped_nsamp
      ENDIF
      CALL ResizeVar(is_co, ss, .FALSE.)
      WHERE(sts(1:SIZE(is_co))==0) is_co = .TRUE.
      CALL ComputeMAF(MAFco, X, sex, MAPA, InclLoc, n_maf_excl, &
                      min_maf, is_co, warn_maf=.FALSE.)

      CALL Prnt("Getting minor allele frequencies in the cases ...")
      CALL ComputeMAF(MAFca, X, sex, MAPA, InclLoc, n_maf_excl, &
                      min_maf, .NOT.is_co, warn_maf=.FALSE.)
                      
    ENDIF
    
  ENDIF

  !! Skip first pair if no_causal_pair is true 
  IF(no_causal_pair .AND. ped_nloci>=2) InclLoc(1:2) = .FALSE.

  !! Change inclusion status so that the total number of loci does not exceed
  !! user given value of ntest_limit  
  IF(ntest_limit > 0) THEN
  
    CALL Prnt("Limiting the number of tests to (around) "//i2cp(ntest_limit)//&
              " ... ", ADVANCE='NO')
    
    n_incl_loci = CountTrue(InclLoc)
    nincl = 0 
    DO i=1,SIZE(InclLoc)
      IF(.NOT.InclLoc(i)) CYCLE
      nincl = nincl + 1
      IF(nincl*(nincl-1)/2 >= ntest_limit .AND. nincl < n_incl_loci) THEN
        n_maxntests_excl = CountTrue(InclLoc(i+1:))
        InclLoc(i+1:) = .FALSE.
        ntests_limit = .TRUE.
        EXIT
      ENDIF
    ENDDO
    
    IF(sub_make_pairs) &
      WHERE(Included_pairs>i) Included_pairs = -1  
  
    CALL Prnt0("Finished.", log=.FALSE.)
    CALL Prnt("Limiting finished.", screen=.FALSE.)

  ENDIF

  !! Count the included loci and announce the number
  n_incl_loci = CountTrue(InclLoc)
  CALL Prnt("After all checks there are "//i2cp(n_incl_loci)//" markers left.")

  !! Get the total number of included loci
  IF(n_incl_loci == 0) THEN
    CALL Report(all_excl=.TRUE.)
    no_testing = .TRUE.
    save_input_data = .FALSE.
    do_out_maf = .FALSE.
    exit_code = 3000
    RETURN
  ELSEIF(n_incl_loci == 1) THEN
    CALL Report(all_but_one_excl=.TRUE.)
    no_testing = .TRUE.
    RETURN
  ELSE
    !! Make sure the length of a cycle for writing output is not too small
    temp_cycle = MAX(NUMBER_OF_THREADS * n_incl_loci, temp_cycle)
  ENDIF

  !! Based on inclusion status remove from X, MAPA, MAF the loci that are excluded
  CALL FilterOutData(InclLoc, X, MAPA, MAF, MAFco, MAFca)
    
  RETURN
   
END SUBROUTINE ProcessData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SaveData()
  IMPLICIT NONE

  !! Save allele frequency data
  IF(do_out_maf) &
    CALL SaveMAF(out_maf_file, MAF, MAPA, InclLoc)
  
  !! Save genetical data
  IF(save_input_data) THEN
    
    IF(ped_nsamp>1) &
      CALL Report(partial_fam=.NOT.save_full_fam)
    
    CALL SaveInputData(save_ped_file, X, FAM, ped_nsamp*ped_ss, &
                       save_binary, save_fam_file, save_map_file, MAPA, &
                       map_cs(1))
    DEALLOCATE(FAM)
  
  ENDIF

  !! Save pretest selection status
  IF(do_out_pss) CALL SavePSS(out_pss_file, pss_AS)

  RETURN
  
END SUBROUTINE SaveData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SavePSS(fn, subsample)
!! Routine that writes the PHASE1 selection status into file
  IMPLICIT NONE
  INTEGER, PARAMETER          :: upss = 32 
  CHARACTER(*), INTENT(INOUT) :: fn
  INTEGER(iks), INTENT(IN)    :: subsample(:)
  CHARACTER(mmtl)             :: h(20,1)
  CHARACTER                   :: PSS(SIZE(subsample),1), c
  INTEGER                     :: i, ih
  
  c = comment
  !! Announce sorting
  CALL Prnt("Saving PHASE1 selection status into file ... ", ADVANCE='NO')
  
  CALL PrintProgramHeader(unit=upss, file=fn, commented=.TRUE.)
  
  ih = 1
  h(ih,1) = c//" This is a PHASE1 selection status (PSS) file generated"//&
               " at "//starttime_full
  ih = ih + 1
  h(ih,1) = c//" It contains information about which individuals were"//&
               " selected for the PHASE1 tests."
  ih = ih + 1
  h(ih,1) = c//" Possible values: 0 (not selected), 1 (selected)"
  ih = ih + 1
  h(ih,1) = c//" You can edit this file to change the PHASE1"//&
               " selection status of an individual and"
  ih = ih + 1
  h(ih,1) = c//" use it during analysis by adding --pss "//TRIM(fn)
  ih = ih + 1
  h(ih,1) = c//" If you edit this file keep in mind that any PSS"//&
               " file must have only 1 column!"
  ih = ih + 1
  h(ih,1) = c
                
  !! Convert to character
  PSS = "0"
  DO i=1,SIZE(subsample)
    IF(subsample(i)==1_iks) PSS(i,1) = "1"
  ENDDO

  !! Write the file header 
  CALL WriteFile(fn, upss, h(1:ih,:), attach=.TRUE.)

  !! Write the PHASE1 selection status into file 
  CALL WriteFile(fn, upss, PSS, attach=.TRUE.)

  !! Announce the end of writing
  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("PHASE1 selection saved into file ["//TRIM(fn)//"]", screen=.FALSE.)
  
  RETURN

END SUBROUTINE SavePSS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadMAPfile(fn, MAPA, NRS, ncols, nl, skip_nl, colsep, submap)
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: fu = umap
  
  CHARACTER(*), INTENT(OUT), ALLOCATABLE :: MAPA(:,:)
  INTEGER(ikb), INTENT(OUT), ALLOCATABLE :: NRS(:)

  !! Input variables
  CHARACTER(*), INTENT(INOUT)       :: fn(:)
  INTEGER, INTENT(IN)               :: skip_nl
  INTEGER, INTENT(INOUT)            :: ncols
  INTEGER, INTENT(OUT)              :: nl
  CHARACTER, INTENT(INOUT)          :: colsep(:)
  LOGICAL, INTENT(IN), OPTIONAL     :: submap
  !! Output variables
  CHARACTER(mml)                    :: chr, rs
  LOGICAL                           :: submap1
  INTEGER                           :: ios, total_file_nl, ifile, i, &
                                       j, file_nl, n_bad_chr, n_bad_rs, &
                                       chcol, rscol, ncols_expect
  CHARACTER(mmtl)                   :: text, text_data, colnames
  CHARACTER(LEN(MAPA)), ALLOCATABLE :: TMP(:,:)
  
  !! Determine whether current mapping file is a submap fil
  submap1 = .FALSE.
  IF(PRESENT(submap)) submap1 = submap
    
  !! Store the kind of data text i
  colnames = ""
  IF(submap1) THEN 
    text_data = "SUBMAP"
    IF(sub_make_pairs) THEN
      IF(sub_no_chr) THEN
        ncols = 2
        colnames = "RS1 RS2"
      ELSE
        ncols = 4
        colnames = "CHR1 RS1 CHR2 RS2"
      ENDIF
    ENDIF
  ELSE
    text_data = "MAP"
  ENDIF
  
  ncols_expect = ncols

  !! Announce what is being done
  CALL Prnt("Reading "//TRIM(text_data)//" file ("//TRIM(i2c(ncols))//&
            " columns expected) ... ", ADVANCE="NO")
  
  !! Initialize variables
  nl = 0
  ios = 0
  total_file_nl = 0
  
  !! Read mapping files
  DO ifile=1,SIZE(fn)

    file_nl = 0
    !! Read the input file
    CALL ReadFile(fn(ifile), fu, MAPA, nl, file_nl, ncols, colsep, skip_nl, &
                  cmt=comment, ios=ios, ignore_ncols=submap1, progress=.TRUE.)
                  
    !! Drop the redundant (empty) rows and columns
    CALL ResizeVar(MAPA, GetEffectiveNRow(MAPA), GetEffectiveNCol(MAPA))
                  
    !! Check for wrong number of columns
    IF(ios == 1005) THEN
      IF(submap1) THEN
        CALL Prnt("Hint: Use --submapsep to specify the correct subset"//&
                       " mapping file column separator and/or --submapncols"//&
                       " to specify the actual number of columns in the"//&
                       " subset mapping file.", Q=.TRUE.)
      ELSE
        CALL Prnt("Hint: Use --mapsep to specify the correct mapping"//&
                       " file column separator and/or --mapncols to specify"//&
                       " the actual number of columns in the mapping file.", &
                       Q=.TRUE.)
      ENDIF
    ENDIF
      
    !! If reading submap file determine the column with RS numbers unless
    !! sub_read_all is true in which case all cells will be considered
    IF(submap1) THEN
      IF(sub_read_all) THEN
        sub_chcol = 0
        sub_rscol = 1
      ELSEIF(.NOT.sub_make_pairs) THEN
        CALL GetRSCol(MAPA, sub_rscol, sub_chcol)
      ENDIF
    ENDIF

    !! Corruption error if no records were read from at least one of the files 
    IF(file_nl==0) &
      CALL PrntE("File ["//TRIM(fn(ifile))//"] apears to be empty.",&
                   ask_to_continue=.TRUE.)
                        
    !! Increase the total counter of lines
    total_file_nl = total_file_nl + file_nl

  ENDDO
  
  !! Check for improper number of columns
  IF(ncols /= ncols_expect) &
    CALL PrntE("The actual number of columns in the file is different"//&
                 " from the expected number.", Q=.TRUE.)
  
  !! Check for alternative names for chromosomes
  !! X    X chromosome                    -> 23
  !! Y    Y chromosome                    -> 24
  !! XY   Pseudo-autosomal region of X    -> 25
  !! MT   Mitochondrial                   -> 26
  IF(submap1) THEN
    chcol = sub_chcol
    rscol = sub_rscol
  ELSE
    chcol = map_chcol
    rscol = map_rscol
  ENDIF

  !! If all elements in submap are supposed to be RS numbers, put all of the
  !! read values into the first column 
  IF(submap1 .AND. sub_make_pairs) THEN
    CALL ResizeVar(TMP, SIZE(MAPA), ncols/2)
    TMP = TRANSPOSE(RESHAPE(TRANSPOSE(MAPA), (/ncols/2, SIZE(MAPA)/(ncols/2)/)))
    CALL ResizeVar(MAPA, SIZE(MAPA), ncols/2)
    MAPA = TMP
    DEALLOCATE(TMP)
  ENDIF
  
  !! Announce the end of map/submap file reading
  IF(.NOT.countdown) CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Finished reading "//TRIM(text_data)//" data!", screen=.FALSE.)
  
  !! Perform chromosome number checks and correction
  n_bad_chr = 0
  IF(chcol > 0 .AND. chcol < ncols) THEN
    
    text = "Validating chromosome information ..."
    CALL Prnt(TRIM(text)//" ", ADVANCE="NO")
    
    DO i=1,nl
      chr = MAPA(i,chcol)
      IF(chr=="X")  chr = "23"
      IF(chr=="Y")  chr = "24"
      IF(chr=="XY") chr = "25"
      IF(chr=="MT") chr = "26"
      
      !! Check for non-numerical representation
      IF(.NOT.IsInteger(chr)) THEN
        n_bad_chr = n_bad_chr + 1
        chr = "0"
      ELSEIF(c2i(chr)<0 .OR. c2i(chr)>26) THEN
        n_bad_chr = n_bad_chr + 1
        chr = "0"
      ENDIF

      MAPA(i,chcol) = chr
    ENDDO
    
    CALL Prnt0("Finished.", log=.FALSE.)
    CALL Prnt("Chromosome validation finished.", screen=.FALSE.)
  
  ELSE
    
    CALL PrntW("Chromosome validity check cannot be performed for "//&
               TRIM(text_data)//" data.")
  
  ENDIF 

  !! Perform RS number validity checks
  n_bad_rs = 0
  IF(rs_valid_chck) THEN

    IF(rscol > 0 .AND. rscol < ncols) THEN

      text = "Validating marker RS numbers ..."
      CALL Prnt(TRIM(text)//" ", ADVANCE="NO")
      CALL ResizeVar(NRS, MAX(1,nl), 0_ikb)
      j = 0

      DO i=1,nl
        rs = MAPA(i,rscol)

        !! Check for validity
        IF(upcasef(rs(1:2))/="RS" .OR. .NOT.IsInteger(rs(3:))) THEN
          CALL PrntW("Invalid RS number found: "//TRIM(rs))
          n_bad_rs = n_bad_rs + 1
          IF(exclude_invalid_rs) CYCLE
        ENDIF

        !! Remove non-numerical characters
        CALL KeepNumbers(rs)

        !! Store inside NRS
        j = j+1
        NRS(j) = c2i(rs)
      ENDDO

      CALL ResizeVar(NRS, MAX(1,j))
      CALL Prnt0("Finished.", log=.FALSE.)
      CALL Prnt("RS number validation finished.", screen=.FALSE.)

    ELSE

      CALL PrntW("RS number validity check could not be performed for"//&
                 " the "//TRIM(text_data)//" data.")

    ENDIF

  ENDIF

  !! Quit if no data was read
  IF(total_file_nl==0) THEN
    CALL PrntE("No "//TRIM(text_data)//" data was found!", Q=.TRUE.)
  ELSE
    CALL Prnt("Total of "//i2cp(total_file_nl)//" input lines found.")
  ENDIF

  !! Announce bad chromosome names
  IF(n_bad_chr > 0) &
    CALL PrntW("There were "//i2cp(n_bad_chr)//" (out of "//i2cp(nl)//&
               ") invalid chromosome names! These value were replaced"//&
               " with 0 (unplaced).")
  
  !! Announce bad chromosome rs numbers
  IF(n_bad_rs > 0) &
    CALL PrntW("There were "//i2cp(n_bad_rs)//" (out of "//i2cp(nl)//") invalid RS numbers!")
  
  RETURN

END SUBROUTINE ReadMAPfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadFAMfile(fn, sts, sex, FAM, ncols, nskip, sep, nco, nca, nmale, &
                       nfema, nsamp)

  IMPLICIT NONE
  INTEGER(iks), INTENT(OUT), ALLOCATABLE :: sts(:), sex(:)
  CHARACTER, INTENT(OUT), ALLOCATABLE    :: FAM(:,:)
  !! Input variables
  INTEGER, INTENT(IN)           :: nskip, nsamp
  INTEGER, INTENT(INOUT)        :: ncols, nco, nca, nmale, nfema
  CHARACTER, INTENT(INOUT)      :: sep(:)
  CHARACTER(*), INTENT(INOUT)   :: fn
  !! Output variables
  CHARACTER                     :: sts_read, sex_read
  CHARACTER(mhtl)               :: line, prev_line
  CHARACTER(mfml), ALLOCATABLE  :: words(:)
  INTEGER                       :: k, i, j, nlines, ss, nc, &
                                   prev_nc, ios1
  LOGICAL                       :: check_expected_ncol
  
  check_expected_ncol = .FALSE.
  
  !! If status is to be assumed, put first cases and then controls
  IF(assume_status) THEN
  
    CALL Prnt("No pedigree information file read. Affection status"//&
              " assumed as 'first cases, then controls'.")
                   
    IF(nco<=0 .OR. nca<=0) &
      CALL PrntE("Status cannot be assumed, numbers of controls or"//&
                      " cases missing or nonpositive.", Q=.TRUE.)        
    ss = nca + nco
    nlines = nsamp*ss
    CALL ResizeVar(sts, nlines, newvalue = NumNA)
    
    !! Fill status with cases and controls
    DO k=1,nsamp
      i = (k-1)*ss
      sts(i+1:i+nca) = NumCa
      sts(i+nca+1:i+ss) = NumCo
    ENDDO
    
    !! Assign male sex
    CALL ResizeVar(sex, nlines, newvalue = NumMa)
    nmale = ss
    nfema = 0

    RETURN
  ENDIF

  !! Otherwise, read the family data in
  CALL Prnt("Reading FAM file ... ", advance='NO')
  
  !! Read raw pedigree data without processing of column separators
  nlines = 0
  CALL ReadTable(fn, FAM)
  
  !! Resize status and sex vectors
  nlines = SIZE(FAM,1)
  CALL ResizeVar(sts, nlines, newvalue = NumNA)
  CALL ResizeVar(sex, nlines, newvalue = NumNA)
  
  !! Nulify counters
  nca = 0
  nco = 0
  nmale = 0
  nfema = 0
  
  !! Process the separators
  prev_nc = -1
  prev_line = ""
  DO i=1+nskip,nlines
    
    !! Store the current individuals's pedigree as a single line
    line = ""
    DO j=1,SIZE(FAM,2) 
      line(j:j) = FAM(i,j)
    ENDDO
    
    !! Separate the line into records (based on sep)
    CALL SplitString(TRIM(line), sep, words, nc, ios1)
    
    !! Check for unexpected number of columns
    IF(check_expected_ncol .AND. nc /= ncols) &
      CALL PrntE("Unexpected number of columns in FAM file (line "//i2cp(i)//").")

    !! Check for inconsistencies in the record count
    IF(i>1 .AND. nc /= prev_nc) THEN
      DO j=1,SIZE(FAM,2) 
        WRITE(*,'(2A)', ADVANCE='NO') FAM(i-1,j),":"
      ENDDO
      DO j=1,SIZE(FAM,2) 
        WRITE(*,'(2A)', ADVANCE='NO') FAM(i,j),":"
      ENDDO
      write (*,"(A)") ""
      CALL PrntE("Inconsistent FAM file (line "//i2cp(i)//").")
    ENDIF
    
    prev_nc = nc
    prev_line = line
    
    !! Process the status and sex info
    sts_read = TRIM(words(fam_stscol))
    sex_read = CharNAsex
    IF(fam_sexcol>0) sex_read = TRIM(words(fam_sexcol))

    !! Assign the values in sts_read to sts
    !! If sts is other than case or control, read next line
    IF(sts_read == CharCa) THEN
      sts(i) = NumCa
    ELSEIF(sts_read == CharCo) THEN
      sts(i) = NumCo
    ELSE
      CYCLE
    ENDIF

    !! Assign the values in sex_read to sex
    IF(assume_male .OR. sex_read == CharMa) sex(i) = NumMa
    IF(assume_fema .OR. sex_read == CharFe) sex(i) = NumFe

  ENDDO
  
  !! Count cases, control, males, females
  nca = COUNT(sts==NumCa)
  nco = COUNT(sts==NumCo)
  nmale = COUNT(sex==NumMa)
  nfema = COUNT(sex==NumFe)
  
  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Finished reading FAM file.", screen=.FALSE.)
  
  !IF(show_summary) THEN
    CALL Prnt("Sample size "//i2c(SIZE(sts))//" (controls: "//i2c(nco)//" | cases: "//&
              i2c(nca)//" || males: "//i2c(nmale)//" | females: "//i2c(nfema)//")")
  !ENDIF
  
  RETURN

END SUBROUTINE ReadFAMfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
SUBROUTINE ReadPEDfile(fn, fu, sts, sex, FAM, X, MAPA, nloci, &
              nskip, colsep, nco, nca, nmale, nfema)

  IMPLICIT NONE
  INTEGER(iks), INTENT(OUT), ALLOCATABLE    :: sts(:), sex(:), X(:,:)
  CHARACTER, INTENT(OUT), ALLOCATABLE      :: FAM(:,:)
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: MAPA(:,:)
  !! Input variables
  INTEGER, INTENT(IN)           :: fu
  INTEGER, INTENT(IN)           :: nskip
  INTEGER, INTENT(INOUT)        :: nloci, nco, nca, nmale, &
                                   nfema
  CHARACTER, INTENT(IN)         :: colsep(:)
  CHARACTER(*), INTENT(INOUT)   :: fn
  !! Output variables
  CHARACTER                     :: sts_read, sex_read, alleles(3), &
                                   sep(SIZE(colsep)) 
  INTEGER(iks), ALLOCATABLE      :: Y(:,:)
  INTEGER                       :: k, i, j, l, nlines, file_nlines, &
                                   nallele(2), nallele_miss, ncols, &
                                   min_map_ncols
  CHARACTER(2), ALLOCATABLE     :: GENES(:,:), genes1(:), G(:,:)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!             READ-IN THE PEDIGREE AND GENETIC INFORMATION              !!  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Otherwise, read the family data in
  CALL Prnt("Reading only pedigree data from the PED file ... ", advance='NO')

  nlines = 0
  file_nlines = 0
  sep = colsep
  
  CALL PrntE("Jakub, modify the code for reading FAM first (ReadPEDfile).")
  
  !! Read the input file (pedigree information first)
  CALL ReadFile(fn, fu, FAM, nlines, file_nlines, fam_ncols, sep, nskip, &
                read_ncols=.TRUE., cmt=comment, progress=.TRUE.)

  !! Announce the end of reading pedigree data
  CALL Prnt("Finished reading pedigree data!", screen=.FALSE.)
  
  !! Read the genetic data
  CALL Prnt("Reading genotype data from the PED file ... ", advance='NO')

  ncols = 0
  nlines = 0
  file_nlines = 0
  !! Learn the number of columns
  CALL ReadFile(fn, fu, G, nlines, file_nlines, ncols, sep, max_nl=1, &
                learn_ncols=.TRUE., cmt=comment)
                
  nlines = 0
  file_nlines = 0
  sep = colsep
  ncols = ncols - fam_ncols
  
  !! Read the genetic information from input file
  CALL ReadFile(fn, fu, GENES, nlines, file_nlines, ncols, sep, nskip, &
                skip_ncols=fam_ncols, cmt=comment, progress=.TRUE.)
            
  !! Check for missing values
  DO k=1,nlines
    IF(ANY(GENES(k,:)=="")) & 
      CALL PrntE("Missing data on line "//i2cp(k)//" of the input file.", &
                 Q=.TRUE.)
  ENDDO
            
  !! Announce the end of reading ped data
  CALL Prnt("Finished reading genotype data!", screen=.FALSE.)
  
  !! Check for improper number of columns found
  IF(IsOdd(ncols)) &
    CALL PrntE("An odd number of genotype columns found in the input file!", Q=.TRUE.)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!      PROCESS THE PEDIGREE INFORMATION INTO STATUS, SEX, ETC.          !!  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Otherwise, read the family data in
  CALL Prnt("Processing input data ... ", advance='NO')

  nloci = ncols / 2
  
  !! Resize sts and sex vectors
  CALL ResizeVar(genes1, nlines, allvalue = "")
  
  !! Nulify counters
  nca = 0
  nco = 0
  nmale = 0
  nfema = 0
  nallele_miss = ncols

  !! If MAPA has fewer columns than needed don't search for alleles
  min_map_ncols = MAX(map_a1col, map_a2col) 
  IF(.NOT.ALLOCATED(MAPA)) &
    CALL ResizeVar(MAPA, nloci, min_map_ncols, "")
  IF(SIZE(MAPA,1)<nloci .OR. SIZE(MAPA,2)<min_map_ncols) &
    CALL ResizeVar(MAPA, nloci, min_map_ncols, "")
  map_ncols = SIZE(MAPA,2)

  !! Find allele names
  DO k=1,ncols,2
    alleles = ""
    nallele = 0
    j = 1    
    DO l=0,1
      genes1 = GENES(:,k+l)
      DO i=1,nlines
        !! Count the alleles
        IF(genes1(i)==alleles(1)) nallele(1) = nallele(1) + 1
        IF(genes1(i)==alleles(2)) nallele(2) = nallele(2) + 1
        !! Check for a new allele
        IF(ALL(genes1(i)/=(/CharNAgen,alleles(1:MAX(1,j-1))/))) THEN
          alleles(j) = genes1(i)(1:1)
          
          IF(j<3) THEN
            nallele_miss = nallele_miss - 1
            j = j + 1
          ELSE
            IF(ALL(alleles(3)/=alleles(1:2))) &
              CALL PrntE("More than 2 alleles found at locus "//&
                           TRIM(MAPA((k+1)/2,map_rscol))//"!", Q=.TRUE.)
                      !" For individual ["//&
                      !TRIM(FAM(i,1))//" "//TRIM(FAM(i,2))//"] allele ["//&
                      !TRIM(alleles(3))//"] found when alleles ["//&
                      !TRIM(alleles(1))//","//TRIM(alleles(2))//&
                      !"] previously encountered for this locus.", Q=.TRUE.)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    
    !! Save the found alleles
    IF(nallele(1) > nallele(2)) THEN
      MAPA((k+1)/2,map_a1col) = alleles(2)
      MAPA((k+1)/2,map_a2col) = alleles(1)
    ELSE
      MAPA((k+1)/2,map_a1col) = alleles(1)
      MAPA((k+1)/2,map_a2col) = alleles(2)
    ENDIF

  ENDDO
  
  !! Allocate and resize
  CALL ResizeVar(sts, nlines, allvalue = NumNA)
  CALL ResizeVar(sex, nlines, allvalue = NumNA)
  CALL ResizeVar(genes1, ncols, allvalue = "")
  CALL ResizeVar(X, (nlines+3)/4, nloci)
  CALL ResizeVar(Y, nlines, nloci, allvalue = 0_iks)

  !! Assign sts and sex info and represent genes numerically
  DO i=1,nlines

    sts_read = TRIM(FAM(i,fam_stscol))
    sex_read = CharNAsex
    IF(fam_sexcol>0) sex_read = TRIM(FAM(i,fam_sexcol))

    !! Assign the values in sts_read to sts
    !! If sts is other than case or control, read next line
    IF(sts_read == CharCa) THEN
      sts(i) = NumCa
      nca = nca + 1
    ELSEIF(sts_read == CharCo) THEN
      sts(i) = NumCo
      nco = nco + 1
    ELSE
      CYCLE
    ENDIF

    !! Assign the values in sex_read to sex
    IF(sex_read == CharMa) THEN
      sex(i) = NumMa
      nmale = nmale + 1
    ELSEIF(sex_read == CharFe) THEN
      sex(i) = NumFe
      nfema = nfema + 1
    ELSEIF(assume_male) THEN
      sex(i) = NumMa 
      nmale = nmale + 1
    ELSEIF(assume_fema) THEN
      sex(i) = NumFe
      nfema = nfema + 1
    ENDIF

    !! Store locally the current individual's genetic information
    genes1 = GENES(i,:)

    !! Change to numerical representation of genes
    DO k=1,ncols,2
      j = (k+1)/2
      alleles(1:2) = (/ MAPA(j,map_a1col), MAPA(j,map_a2col) /)
      !! Compare odd/even columns with first line and assign either 1 (same) or 0 (different) or 9 (missing)
      IF(ALL(genes1(k)/=alleles(1:2)) .OR. ALL(genes1(k+1)/=alleles(1:2))) THEN  
        Y(i,j) = NumNA
      ELSE
        IF(genes1(k)==alleles(1))   Y(i,j) = Y(i,j) + 1_iks
        IF(genes1(k+1)==alleles(1)) Y(i,j) = Y(i,j) + 1_iks
      ENDIF
    ENDDO
    
  ENDDO
  
  !! Check whether minor allele is the one that is counted
  CALL CheckAlleleCount(Y, count_minor=.TRUE.)

  !! Convert the read ped info into binary format and return in X
  CALL Ped2Bed(Y, X, bed_major, ped_minor)

  !! Announce the end of processing
  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Finished processing input data!", screen=.FALSE.)
  
  IF(nallele_miss > 0) &
    CALL PrntE("Fewer than expected allele names were found in the PED file.", Q=.TRUE.)

  RETURN

END SUBROUTINE ReadPEDfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
SUBROUTINE ReadPSSFile(fn, status, nskip, d, pss)
!! Read the PHASE1 selection status file. This file must contain the same
!! number of lines as the ped_file (i.e. entire sample size) and the first
!! column (it shouldn't have more columns, but it's not a problem) contains
!! the PHASE1 selection status. That is 1 in n-th row if the n-th individual
!! in the sample is supposed to be selected for pretest and 0 otherwise.  
  IMPLICIT NONE
  INTEGER, PARAMETER               :: fu = 32
  CHARACTER(*), INTENT(INOUT)      :: fn
  INTEGER(iks), INTENT(INOUT)      :: status(:)
  INTEGER, INTENT(IN)              :: nskip
  REAL(dpp), INTENT(OUT), OPTIONAL :: d
  LOGICAL, INTENT(IN), OPTIONAL    :: pss
  CHARACTER(mfmll), ALLOCATABLE    :: S(:,:)
  INTEGER                          :: i, j, k, ncols, nlines, maxnitems, &
                                      nitems, file_nlines
  LOGICAL                          :: pss1
  
  !! Select what type of data is read
  pss1 = .FALSE.
  IF(PRESENT(pss)) pss1 = pss
  
  !! Announce what is being read
  IF(pss1) THEN
    CALL Prnt("Reading PHASE1 selection status data ... ", advance='NO')
  ELSE
    CALL Prnt("Reading status data ... ", advance='NO')
  ENDIF
  
  !! Initialize variables
  ncols = 1
  nlines = 0
  file_nlines = 0
  
  !! Read the input file
  CALL ReadFile(fn, fu, S, nlines, file_nlines, ncols, nskip=nskip, &
                cmt=comment, fixed_ncols=.TRUE., progress=.TRUE.)

  !! What is the maximum number of items read from the file
  maxnitems = SIZE(status)
  !! How many items were actually read from the file
  nitems = nlines*ncols
  
  !! Corruption error if no or too many items were read 
  IF(nitems==0) &
    CALL PrntE("File ["//TRIM(fn)//"] apears to be empty.", Q=.TRUE.)
  IF(nitems > maxnitems) &
    CALL PrntW("File ["//TRIM(fn)//"] contains too many lines"//&
               " ("//i2cp(nitems)//" lines found when "//i2cp(maxnitems)//&
               " lines expected). Only the first "//i2cp(maxnitems)//" lines"//&
               " will be used!", Q=.TRUE.)

  !! Process the read values
  k = 0
  DO i=1,nlines
    DO j=1,ncols
      k = k+1
      !! Process the next value
      IF(S(i,j)=="0") THEN
        IF(pss1) THEN
          status(k) = 0_iks
        ELSE
          status(k) = NumNA
        ENDIF
      ELSEIF(S(i,j)=="1") THEN
        IF(pss1) THEN
          status(k) = 1_iks
        ELSE
          status(k) = NumCo
        ENDIF
      ELSEIF(pss1) THEN
        CALL PrntE("Unknown value encountered! Value "//TRIM(S(i,j))//&
                   " found when only values 1 (selected) or 0 (not"//&
                   " selected) allowed.", Q=.TRUE.)
      ELSEIF(S(i,j)/="2") THEN
        CALL PrntE("Unknown value encountered! Value "//TRIM(S(i,j))//&
                   " found when only values 0 (NA), 1 (control) or 2 (case)"//&
                   " allowed.", Q=.TRUE.)
      ELSE
        status(k) = NumCa
      ENDIF
    ENDDO 
  ENDDO
  
  !! Set the PHASE1 sample size ratio
  IF(PRESENT(d)) d = Mean(status)

  !! Announce the end of file processing
  IF(pss1) THEN
    CALL Prnt("Finished reading PHASE1 selection status!", screen=.FALSE.)
  ELSE  
    CALL Prnt("Finished reading PSS file!", screen=.FALSE.)
  ENDIF
  
  RETURN

END SUBROUTINE ReadPSSFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! SUBROUTINE AnalyzeInputPEDFiles(fn, fu, ncols, maxncols, nlines, &
!               colsep, ss, nsamp, nloci, famncols, map_nlines)
!   IMPLICIT NONE
!   CHARACTER(*), INTENT(INOUT)       :: fn(:)
!   INTEGER, INTENT(IN)               :: fu, maxncols, famncols, map_nlines
!   INTEGER, INTENT(INOUT)            :: ss
!   INTEGER, INTENT(OUT)              :: ncols(:), nsamp, nloci, nlines 
!   CHARACTER(*), INTENT(INOUT)       :: colsep
!   CHARACTER(30)                     :: read1(MAX(1,famncols-1))
!   CHARACTER                         :: read2, read3(maxncols), sep(4)
!   CHARACTER(4*maxncols)             :: line
!   LOGICAL                           :: file_too_big
!   INTEGER                           :: k, nfound, next, previous, lc,&
!                                        sepcount(4), which_max, count, &
!                                        linelen, ios, ifile, col_counter, &
!                                        lines_read
!   INTEGER(ikb)                      :: filesize
! 
!   !! Check for array sizes
!   IF(SIZE(ncols) /= SIZE(fn)) &
!     CALL PrntE("Wrong array size. variables do not match.", Q=.TRUE.)
!   
!   !! Initialize variable with input files dimensions
!   lc = 0
! 
!   !! Print announcement of started reading
!   IF(SIZE(fn)==1) THEN
!     CALL Prnt("Reading ped data ... ", advance='NO')
!     !IF(countdown) WRITE(usto,'(A)', ADVANCE='NO') TRIM(analyzing)
!     IF(countdown) CALL Prnt0(analyzing, ADVANCE='NO', log=.FALSE.)
!   ELSE
!     CALL Prnt("Analyzing ped data ... ", advance='NO')
!   ENDIF
!       
!   !! Check for too big input files
!   file_too_big = .FALSE.
!   DO ifile=1,SIZE(fn)
!     CALL CheckFileExistence(fn=fn(ifile), delim=delim)
!     filesize = GetFileSize(fn(ifile)) 
!     IF(filesize > max_fs) THEN
!       CALL PrntE("PED file ["//TRIM(ped_file(ifile))//"] is too large!"//&
!               "Size of this file is "//i2c(filesize)//" bytes while the"//&
!               " maximum allowed file size if "//i2c(max_fs)//" bytes.")
!       CALL Prnt("Hint: Split this input file into several files and use"//&
!               " --pedlist to input it.")
!       file_too_big = .TRUE.
!     ENDIF
!   ENDDO
!   !! Terminate only now so that all too large files get announced
!   IF(file_too_big) CALL Terminate()
! 
!   DO ifile=1,SIZE(fn)
!   
!     !! Initialize values
!     nfound = 0
!     col_counter = 0
!     line = ""
!     read2 = CharNAgen
!     ios = 0
!     lines_read = 0
!  
!     !! Open the input file (only if it exists)
!     CALL OpenFile(UNIT=fu, FILE=fn(ifile), ACTION='R', &
!                   STATUS='O', FORM='F', POSITION="R", ACCESS='S')
! 
!     DO WHILE (read2/=CharCa .AND. read2/=CharCo .AND. ios==0)
! 
!       !! Read a line
!       !READ(fu,line_format,ADVANCE='NO',SIZE=linelen,IOSTAT=ios) line
!       READ(fu,'(A)', ADVANCE='NO', SIZE=linelen, IOSTAT=ios) line
!       
!       !! Get rid of trailing empty spaces
!       linelen = LEN_TRIM(line)
! 
!       !! Increase counter for error message purposes
!       lines_read = lines_read + 1
! 
!       !! If ios indicates error check whether it is just end-of-line
!       IF(ios==-2 .AND. linelen>0) THEN
!         ios = 0
!         IF(line(1:1)/=comment) lc = lc + 1
!       ENDIF
! 
!       IF(ios/=0) EXIT
! 
!       !! Get the number of separators if unknown
!       !! If separator was not specified by user, it is assumed to be 
!       !! either tab, space, comma or semicolon is the separator
!       !! horizontal tab .. 9, space .. 32, comma .. 44, semicolon .. 59 
!       !! Checks how many of each of these there are and chooses the most 
!       !! frequent one    
!       IF(colsep == NA_colsep) THEN
!       
!         sep = (/ tab, space, comma, semicolon /)
!         sepcount = 0
!         DO k=1,linelen
!           IF(line(k:k) == sep(1)) sepcount(1) = sepcount(1) + 1
!           IF(line(k:k) == sep(2)) sepcount(2) = sepcount(2) + 1
!           IF(line(k:k) == sep(3)) sepcount(3) = sepcount(3) + 1
!           IF(line(k:k) == sep(4)) sepcount(4) = sepcount(4) + 1
!         ENDDO
!         which_max = 1
!         DO k=2,size(sepcount)
!            IF(sepcount(k) > sepcount(which_max)) which_max = k        
!         ENDDO
!         colsep = sep(which_max)
!         !! The number of gene data cols is the number of separators minus 
!         !! fam_ncols plus 1
!         count = sepcount(which_max) - famncols + 1
!         !! Return the number of actual loci (1 locus = 2 gene data columns) 
!         col_counter = count / 2
!      
!       !! Get the number of separator if it is known
!       ELSE
!        
!         sepcount = 0
!         DO k=1,linelen
!           IF(line(k:k) == colsep) sepcount(1) = sepcount(1) + 1
!         ENDDO
!         !! The number of gene data cols is the number of separators minus 
!         !! fam_ncols plus 1
!         count = sepcount(1) - famncols + 1
!         !! Return the number of actual loci (1 locus = 2 gene data columns) 
!         col_counter = count / 2
!        
!       ENDIF
!       
!       !! Assign the read values into variables read1, read2, read3
!       read1=""; read2=""; read3=""
!       nfound = 1
!       previous = 0
!       next=1
!        
!       !! Split what was read according to colsep
!       DO WHILE(next<=linelen)
!         IF(line(next:next) /= colsep) THEN
!           k = nfound
!           IF(nfound < famncols .AND. next-previous<=LEN(read1(k)))&
!               read1(k) = TRIM(read1(k)) // line(next:next)
!           IF(nfound == famncols .AND. next-previous<=LEN(read2))  &
!              read2 = TRIM(read2) // line(next:next)
!           k = nfound - famncols
!           IF(nfound > famncols .AND. next-previous<=LEN(read3(k)) &
!             .AND. k<=size(read3)) read3(k)=TRIM(read3(k))//line(next:next)
!           IF(nfound > famncols .AND. k>size(read3)) &
!             CALL Prnt("Error in the input file on line "//&
!                            i2cp(lines_read)//".", Q=.TRUE.)
!             
!         ELSE
!           previous = next
!           nfound = nfound + 1
!         ENDIF
!         next = next+1
!       ENDDO
!     
!       !! Check for proper number of genes found
!       IF(count/=INT(count/2)*2 .OR. nfound/=famncols+2*col_counter .OR. &
!       count<2) THEN
!        CALL PrntE("Unexpected format of input ped data.")
!        IF(count<2) &
!         CALL Prnt("Number of gene data columns appears to be too low."//&
!                         " Only "//i2c(count)//" columns found.")
!        IF(count/=INT(count/2)*2) &
!         CALL Prnt("An odd number of gene data columns found!"//&
!                         " (Found "//i2c(count)//" columns.)")
!        IF(nfound/=famncols+2*col_counter) & 
!          CALL Prnt("Wrong number of columns in gene data.")
!        CALL Prnt("Suggestion: Check the input ped file.",Q=.TRUE.)
!       ENDIF
!       
!     ENDDO
!     
!     !! Count the number of remaining lines in the file
!     DO WHILE(ios==0)
!       READ(fu,*,IOSTAT=ios)
!       IF(ios==0) lc = lc + 1
!     ENDDO
!     
!     CLOSE(fu)
! 
!     IF(ifile==1) nlines = lc
!     ncols(ifile) = col_counter
! 
!   ENDDO ! for ifile
!   
!   IF(SIZE(fn)>1) THEN
!     CALL Prnt0("Finished.", skip2=1, log=.FALSE.)
!     CALL Prnt("Finished.", screen=.FALSE.)
!   ENDIF
! 
!   !! CHECK FOR ERRORS
!   !! Check whether anything was even read
!   IF(lc <= 0) &
!     CALL PrntE("No ped data found in the input file!", q=.true.)
! 
!   !! Total number of lines in all files must be a product of the number of files
!   !! and the number of lines in the first file (all files must have same number of lines)
!   IF(SIZE(fn)*nlines /= lc) &
!     CALL PrntE("All ped files must have the same number of lines!", &
!                     Q=.TRUE.)
! 
!   !! Get sample information
!   IF(ss <= 0) THEN
!     nsamp = 1
!     ss = nlines
!   ELSE
!     IF(MODULO(nlines, ss) /= 0) &
!       CALL PrntE("Number of lines in the input file is not a multiple"//&
!                       " of the given sample size.", Q=.TRUE.)
! 
!     nsamp = nlines/ss
!   ENDIF
!   
!   !! Maximum number of tests (gives as all possible tests in all samples)
!   nloci = SUM(ncols)
!   
!   !! Check whether the number of loci corresponds to the number in map file(s)
!   IF(map_nlines > 0 .AND. nloci /= map_nlines) &
!     CALL PrntE("Number of loci found in the ped file(s) does "//&
!                     "not match the number found in the mapping file(s).", &
!                     Q=.TRUE.)
! 
!   RETURN
! 
! END SUBROUTINE AnalyzeInputPEDFiles 
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! SUBROUTINE ReadPEDfile1(fn, unit, X, sts, sex, FAM, MAPA, ncases, &
!               ncontrols, nmale, nfema, ncols, nlines, nskip, famncols, &
!               map_ncols, read_alleles)
!   IMPLICIT NONE
!   INTEGER(iks), INTENT(OUT), ALLOCATABLE    :: X(:,:), sts(:), sex(:)
!   CHARACTER(*), INTENT(OUT), ALLOCATABLE   :: FAM(:,:)
!   CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: MAPA(:,:)
!   CHARACTER(*), INTENT(INOUT) :: fn(:)
!   INTEGER, INTENT(IN)         :: unit, ncols(:), nlines, famncols, &
!                                  nskip
!   INTEGER, INTENT(OUT)        :: ncontrols, ncases, nmale, nfema
!   INTEGER, INTENT(INOUT)      :: map_ncols
!   LOGICAL, INTENT(IN)         :: read_alleles
!   INTEGER(iks)                 :: Y(nlines, SUM(ncols))
!   CHARACTER(LEN(FAM(1,1)))    :: fam_line(famncols)
!   CHARACTER(mstl)             :: text
!   CHARACTER                   :: sts_read, sex_read, firstchar, &
!                                  genes_read(2*MAXVAL(ncols)), &
!                                  genes1(2*MAXVAL(ncols))
!   INTEGER                     :: lines_read, lc, ncol, pncol, &
!                                  i, j, k, ifile, fileline, ios, pct, &
!                                  total_lc, nloci, nallele_miss, &
!                                  min_map_ncols
!                                  
!   !! Save the total number of columns
!   nloci = SUM(ncols)
!   
!   !! Allocate arrays
!   CALL ResizeVar(sts, nlines)
!   CALL ResizeVar(sex, nlines)
!   CALL ResizeVar(FAM, nlines, famncols)
!   
!   !! Assign initial values to Y and counters
!   Y = 0
!   total_lc = 0
!   ncases = 0
!   ncontrols = 0
!   nmale = 0
!   nfema = 0
!   nallele_miss = 2*nloci
!   
!   !! If MAPA has fewer columns than needed don't search for alleles
!   min_map_ncols = MAX(map_a1col, map_a2col) 
!   IF(read_alleles) THEN
!     IF(.NOT.ALLOCATED(MAPA)) THEN
!       CALL ResizeVar(MAPA, nloci, min_map_ncols, "")
!     ELSEIF(SIZE(MAPA,1)<nloci .OR. SIZE(MAPA,2)<min_map_ncols) THEN
!       CALL ResizeVar(MAPA, nloci, min_map_ncols, "")
!     ENDIF
!     map_ncols = SIZE(MAPA,2)
!   ELSE
!     nallele_miss = 0
!   ENDIF
! 
!   !! Read all ped files
!   DO ifile=1,SIZE(fn)
!    
!     !! Announce reading and current filename
!     IF(SIZE(fn)>1) THEN
!       IF(ifile==1) CALL Prnt("Reading ped data ... ", skip2=1)
!       CALL Prnt("Reading file ["//TRIM(fn(ifile))//"] ... ",&
!                      advance='N')
!       CALL Prnt("", screen=.FALSE.)
!     ENDIF
! 
!     !! Number of columns in current file
!     ncol = ncols(ifile)
! 
!     !! Number of columns in previous files
!     pncol = 0
!     IF(ifile > 1) pncol = SUM(ncols(1:ifile-1))
!   
!     !! Open the input file (only if it exists)
!     CALL CheckFileExistence(fn=fn(ifile), delim=delim)
!     CALL OpenFile(UNIT=unit, FILE=fn(ifile), ACTION='R', &
!                   STATUS='O', FORM='F', POSITION="R", ACCESS='S')
!     
!     fileline = 0
!     lines_read = 0
!     ios = 0
!     
!     !! A cycle to skip specified number of lines 
!     !! Does not count skipping of lines that are commented!!!
!     k=0
!     DO WHILE (k < nskip)
!       READ(unit, '(A)', IOSTAT=ios) firstchar
!       IF(ios/=0) EXIT
!       fileline = fileline + 1
!       IF(firstchar/=comment) k = k + 1
!     ENDDO
! 
!     !! Set the variables for countdown
!     text = ""
!     IF(SIZE(fn)==1 .AND. countdown) text = analyzing
!     pct = -1
!   
!     !! Read the ped file into genes and sts_read      
!     DO WHILE(ios == 0 .AND. lines_read<=nlines)
! 
!       !! Print progress
!       CALL PrintCD(text, INT(lines_read,ikb), INT(nlines,ikb), &
!                           pct, finish=.TRUE.)
!     
!       !! Read the next line of the input file
!       READ(unit, *, IOSTAT=ios) fam_line, genes_read(:2*ncol)
!       firstchar = fam_line(1)
! 
!       !! If commented or empty line skip it
!       IF(ios==0) fileline = fileline + 1
!       IF(ios/=0 .OR. firstchar==comment) CYCLE
!                                          
!       !! If ios ok then increase counter lines_read
!       total_lc = total_lc + 1
!       lines_read = lines_read + 1
! 
!       !! Remember sex and sts
!       sts_read = fam_line(famncols)
!       sex_read = CharNAsex
!       IF(famncols>1) sex_read = fam_line(famncols-1)  
!       
!       !! Save lines_read into i to save space
!       i = lines_read
! 
!       !! Read sts from the first file only
!       IF(ifile==1) THEN 
!         
!         !! Remember the entire family part of the current line
!         FAM(i,:) = fam_line   
!         
!         !! Assign the values in sts_read to sts
!         !! If sts is other than case or control, read next line
!         IF(sts_read == CharCa) THEN
!           sts(i) = NumCa
!           ncases = ncases + 1
!         ELSEIF(sts_read == CharCo) THEN
!           sts(i) = NumCo
!           ncontrols = ncontrols + 1
!         ELSE
!           CYCLE
!         ENDIF
! 
!         !! Assign the values in sex_read to sex
!         IF(sex_read == CharMa) THEN
!           sex(i) = NumMa
!           nmale = nmale + 1
!         ELSEIF(sex_read == CharFe) THEN
!           sex(i) = NumFe
!           nfema = nfema + 1
!         ELSEIF(sex_read == CharNAsex) THEN
!           IF(assume_male) THEN
!             sex(i) = NumMa 
!             nmale = nmale + 1
!           ELSE
!             sex(i) = NumFe 
!             nfema = nfema + 1
!           ENDIF
!         ENDIF
!       ENDIF
! 
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       !!! THIS CODE CURRENTLY ASSUMES THAT THE FIRST LINE OF EACH INPUT !!!
!       !!! FILE DOES NOT CONTAIN ANY MISSING VALUES IN GENE COLUMNS      !!!
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       IF(lines_read==1) THEN
! 
!         !! Remember the first line
!         genes1 = genes_read
!         !! Check whether there are any NA or other wrong values on 1st line
!         IF(ANY(genes1=="" .OR. ANY(genes1==CharNAgen))) THEN 
!           CALL PrntE("Missing or NA data on the first input line (file"//&
!                   " line "//i2c(fileline)//"). First input line cannot"//&
!                   " contain any missing or NA values. Ignoring this line.")
!           lines_read = 0
!           CYCLE
!         ENDIF
! 
!         !! Transform read genetic data into numerical representation
!         Y(1,pncol+1:pncol+ncol) = 1
!         DO k=1,2*ncol,2
!           IF(genes_read(k) == genes_read(k+1)) &
!             Y(1,pncol+(k+1)/2) = Y(1,pncol+(k+1)/2) + 1_iks 
!         ENDDO
!         lc=1
! 
!       !! Read line only if sts is not NA, otherwise read next line
!       ELSEIF(sts(i) == NumCa .OR. sts(i) == NumCo) THEN
! 
!         !! Increase line counter
!         lc = lc+1
! 
!         !! Check whether genes has even number of columns, if not then stop
!         IF(ANY(genes_read=="")) & 
!           CALL PrntE("Missing data on line "//i2cp(k)//&
!                          " of the input file "//i2cp(ifile)//".", Q=.TRUE.)
! 
!         !! Change to numerical representation of genes
!         !! NOTE: at this point size(genes_read) is even
!         DO k=1,2*ncol,2
!           j = pncol+(k+1)/2
!           !! Compare odd/even columns with first line and assign either 1 (same) or 0 (different) or 9 (missing)
!           IF(genes_read(k)==CharNAgen .OR. genes_read(k+1)==CharNAgen) &  
!             Y(lc,j) = NumNA
!           IF(genes_read(k)==genes1(k)) &
!             Y(lc,j) = 1_iks
!           IF(genes_read(k+1)==genes1(k) .AND. Y(lc,j)/=NumNA) &
!             Y(lc,j) = Y(lc,j) + 1_iks
!         ENDDO
! 
!       ENDIF
!       
!       !! If allele names are not to be read or all have been found skip the rest
!       !! current line
!       IF(.NOT.read_alleles .OR. nallele_miss == 0) CYCLE
! 
!       !! If we still do not know all the allele names, look for them
!       DO k=1,2*ncol,2
!         j = pncol+(k+1)/2
!         !! Find the first allele
!         IF(LEN_TRIM(MAPA(k,map_a1col))==0) THEN
!           IF(genes_read(k)/=CharNAgen) THEN
!             MAPA(j,map_a1col) = genes_read(k)
!             nallele_miss = nallele_miss - 1
!           ELSEIF(genes_read(k+1)/=CharNAgen) THEN
!             MAPA(j,map_a1col) = genes_read(k+1)
!             nallele_miss = nallele_miss - 1
!           ENDIF
!         ENDIF 
!         !! Then find the other allele
!         IF(LEN_TRIM(MAPA(k,map_a2col))==0) THEN
!           IF(ALL(genes_read(k)/=(/CharNAgen,MAPA(j,map_a1col)/))) THEN
!             MAPA(j,map_a2col) = genes_read(k)
!             nallele_miss = nallele_miss - 1
!           ELSEIF(ALL(genes_read(k+1)/=(/CharNAgen,MAPA(j,map_a1col)/))) THEN
!             MAPA(j,map_a2col) = genes_read(k+1)
!             nallele_miss = nallele_miss - 1
!           ENDIF
!         ENDIF 
!       ENDDO
! 
!     ENDDO
!     
!     CLOSE(unit)
! 
!   ENDDO
!   
!   !! Check whether minor allele is the one that is counted
!   CALL CheckAlleleCount(Y, count_minor=.TRUE.)
! 
!   !! Convert the read ped info into binary format and return in X
!   CALL ResizeVar(X, (nlines+3)/4, nloci)
!   CALL Ped2Bed(Y, X, bed_major, ped_minor)
! 
!   !! Announce the end of input file reading
!   IF(SIZE(fn)>1) THEN
!     CALL Prnt("All input files read!", skip1=1)
!   ELSE
!     CALL Prnt0("Finished.", screen=.NOT.countdown,log=.FALSE.)
!     CALL Prnt("Finished reading ped data!", screen=.FALSE.)
!   ENDIF
!   
!   IF(total_lc==0) &
!     CALL PrntE("No data on input! Pedding information missing!", Q=.TRUE.)
! 
!   IF(read_alleles .AND. nallele_miss > 0) &
!     CALL PrntE("Less than expected allele names were found in the"//&
!                     " ped file!", Q=.TRUE.)
! 
!   RETURN
!     
! END SUBROUTINE ReadPEDfile1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadBEDfile(fn, X, ss, nsamp, nloci)
  IMPLICIT NONE
  INTEGER(iks), INTENT(OUT), ALLOCATABLE :: X(:,:)
  INTEGER, INTENT(IN)                    :: ss, nsamp, nloci
  CHARACTER(*), INTENT(INOUT)            :: fn(:)
  INTEGER(iks), ALLOCATABLE              :: Y(:,:)
  INTEGER                                :: nrec, ifile, iostat, nrowX, &
                                            ncolX, nrowY, ncolY
  INTEGER(ikb)                           :: filesize
  INTEGER(1)                             :: b1, b2, b3
  CHARACTER(mltl)                        :: text
  
  !! Loop over possibly more than 1 input files
  DO ifile=1,SIZE(fn)
   
    CALL Prnt("Reading BED data ... ", ADVANCE='NO')
    
    !! Open the input file (only if it exists)
    CALL CheckFileExistence(fn=fn(ifile), delim=delim)
    
    !!! Get size of the bedding file
    IF(check_bed_filesize) THEN
      filesize = GetFileSize(fn(ifile))
      IF(filesize==0) GOTO 13
    ENDIF
    
    !! If things ok until now, then open the file for reading first 3 bytes
    CALL OpenFile(UNIT=ubed, FILE=fn(ifile), ACTION='R', STATUS='O', &
                  FORM='U', POSITION="R", ACCESS='D', RECL=3)
    
    !! Read first 3 bytes from bed file
    READ(ubed, REC=1, IOSTAT=iostat) b1, b2, b3
    
    !! Close the file
    CLOSE(ubed)

    !! Check for proper values of these bytes
    IF(b1/=bed_byte1 .OR. b2/=bed_byte2 .OR. (b3/=0 .AND. b3/=1)) &
      CALL PrntE("Specified file is not a PLINK v1.00 BED file!", Q=.TRUE.)
    
    !! Determine whether bed file is SNP-major or individual-major
    nrowX = ((ss + 3) / 4) * nsamp
    ncolX = nloci
    CALL ResizeVar(X, nrowX, ncolX)
    X = NAbed

    IF(b3 == 1) THEN
      bed_data_SNPmajor = .TRUE.
      nrec = nrowX * ncolX
    ELSE
      bed_data_SNPmajor = .FALSE.
      nrowY = (nloci + 3) / 4 
      ncolY = nsamp * ss  
      CALL ResizeVar(Y, nrowY, ncolY)
      Y = NAbed
      nrec = nrowY * ncolY
    ENDIF

    !!! Check whether filesize fits the mapping file data
    IF(check_bed_filesize .AND. filesize/=nrec+3) GOTO 11
    
    !! Reopen the file with different RECL
    CALL OpenFile(UNIT=ubed, FILE=fn(ifile), POSITION='R', ACTION='R', &
                  STATUS='O', FORM='U', ACCESS='D', RECL=nrec+3)
                  
    !! Read the entire file at once
    IF(bed_data_SNPmajor) THEN
      !! Read into output variable X, no conversion of formats necessary
      READ(ubed, REC=1, IOSTAT=iostat) b1, b2, b3, X
      IF(iostat /= 0) GOTO 12
    ELSE
      !! Read into local variable Y, conversion to SNP-major format necessary
      READ(ubed, REC=1, IOSTAT=iostat) b1, b2, b3, Y
      IF(iostat /= 0) GOTO 12
      !! Do the conversion from indiv-major to snp-major and remove Y
      CALL ConvertBEDmajor(Y,X)
      DEALLOCATE(Y)      
    ENDIF
    
    CLOSE(ubed)

  ENDDO

  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Finished reading BED file!", screen=.FALSE.)

  !! Check for correct binary ped file format (only SNP-major supported)
  IF(bed_data_SNPmajor) THEN
    CALL Prnt("Input BED file has the Plink v1.00 SNP-major format.")
  ELSE
    CALL PrntE("Input BED file has the individual-major format,"//&
                 " which is currently not supported.")
  ENDIF
  
  RETURN

  !! Error messages
  11 CONTINUE
  
  text = "BED and FAM files are inconsistent! (BED file size is "//&
         i2c(filesize)//" bytes when size of "//i2c(nrec+3)//&
         " bytes expected based on the FAM file (given a total sample size "//&
         i2c(nsamp*ss)//" over "//i2c(nsamp)//" samples)."
#ifdef __INTEL_COMPILER
  text = TRIM(text)//" (If this error occurred due to a miscalculated"//&
         " size of the BED file make sure you used '-assume byterecl' flag"//&
         " during compilation by Intel Fortran compiler.)"
#endif

  CALL PrntE(text, Q=.TRUE.)
  
  12 CONTINUE
  
  text = "Unknown error occurred when reading the BED file (IOSTAT "//i2cp(iostat)//")."
  CALL PrntE(text, Q=.TRUE.)
    
  13 CONTINUE
  
  CALL PrntE("The file ["//TRIM(fn(ifile))//"] appears to be empty.", Q=.TRUE.)

END SUBROUTINE ReadBEDfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SavePrev(prev_file, prev)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)           :: prev_file
  REAL(dpp), INTENT(IN)                 :: prev

  !! Announce what's happening
  CALL Prnt("Saving actual observed prevalence into file ... ", ADVANCE='NO')

  !! Open the file for writing
  CALL OpenFile(UNIT=umaf, FILE=prev_file, ACTION='W', STATUS='U', &
                FORM='F', ACCESS='S', chck_exist=.TRUE., &
                overwrite=overwrite_files)
                
  !! Write the header: chromosome, snp, allele1, allele2, maf, non-missing allele count
  WRITE(umaf,'(A)') "Actual observed prevalence of cases: "//TRIM(r2c(prev))

  CLOSE(umaf)
  
  !! Announce the end of writing
  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Prevalence info saved into ["//TRIM(prev_file)//"]")

  RETURN
    
END SUBROUTINE SavePrev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SaveMAF(maf_file, MAF, MAPA, InclLoc)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)           :: maf_file
  REAL(dpp), INTENT(IN)                 :: MAF(:,:)
  CHARACTER(*), INTENT(IN)              :: MAPA(:,:)
  LOGICAL, INTENT(IN), OPTIONAL         :: InclLoc(:)
  LOGICAL                               :: InclLoc1(SIZE(MAF))
  INTEGER                               :: k
  CHARACTER(10)                         :: fm

  !! Check for sizes
  IF(SIZE(MAPA,1) /= SIZE(MAF,1)) &
    CALL PrntE("MAPA and MAF must have same number of rows (SaveMAF)!", Q=.TRUE.)
  IF(SIZE(MAPA,2)<MAX(map_chcol, map_rscol, map_a1col, map_a2col)) &
    CALL PrntE("Missing some columns in the mapping array (SaveMAF)!", Q=.TRUE.)

  !! If filename missing throw an error
  IF(LEN_TRIM(maf_file)==0) &
    CALL PrntE("Missing name for minor allele frequency file. When"//&
               " --saveinput is present and the data is not simulated,"//&
               " the root name of the output ped file must be"//&
               " specified (use --savefile).", Q=.TRUE.)

  !! ******************** !!
  !!     SAVE MAF DATA    !!
  !! ******************** !!
  
  !! Announce what's happening
  CALL Prnt("Saving allele frequency data into file ... ", ADVANCE='NO')

  InclLoc1 = .TRUE.
  IF(PRESENT(InclLoc)) InclLoc1 = InclLoc

  !! Open the file for writing
  CALL OpenFile(UNIT=umaf, FILE=maf_file, ACTION='W', STATUS='U', &
                FORM='F', ACCESS='S', chck_exist=.TRUE., &
                overwrite=overwrite_files)
                
  fm = "F0."//i2cp(ndig_mafs)

  !! Write the header: chromosome, snp, allele1, allele2, maf, non-missing 
  !! allele count
  WRITE(umaf,'(A)') "Chr"//out_cs//"SNP"//out_cs//"A1"//&
                    out_cs//"A2"//out_cs//&
                    "MAF"//out_cs//"NCHROBS"//out_cs//&
                    "MAFco"//out_cs//"NCHROBSco"//out_cs//&
                    "MAFca"//out_cs//"NCHROBSca"
  !! Write the mapping info
  DO k=1,SIZE(MAF,1)
    
    IF(.NOT.InclLoc1(k)) CYCLE
    
    !! Write location info
    WRITE(umaf,'(A)', ADVANCE='NO') TRIM(MAPA(k,map_chcol))//out_cs
    WRITE(umaf,'(A)', ADVANCE='NO') TRIM(MAPA(k,map_rscol))//out_cs
    WRITE(umaf,'(A)', ADVANCE='NO') TRIM(MAPA(k,map_a1col))//out_cs
    WRITE(umaf,'(A)', ADVANCE='NO') TRIM(MAPA(k,map_a2col))//out_cs
    
    !! Write MAF (both values and the non-missing allele counts in whole data
    WRITE(umaf,'('//TRIM(fm)//',A)', ADVANCE='NO') MAF(k,1), out_cs
    WRITE(umaf,'(A)', ADVANCE='NO') i2cp(INT(MAF(k,2)))//out_cs
    
    !! Write MAF (both values and the non-missing allele counts in controls
    WRITE(umaf,'('//TRIM(fm)//',A)', ADVANCE='NO') MAFco(k,1), out_cs
    WRITE(umaf,'(A)', ADVANCE='NO') i2cp(INT(MAFco(k,2)))//out_cs
    
    !! Write MAF (both values and the non-missing allele counts in cases
    WRITE(umaf,'('//TRIM(fm)//',A)', ADVANCE='NO') MAFca(k,1), out_cs
    WRITE(umaf,'(A)') i2cp(INT(MAFca(k,2)))
    
  ENDDO
  
  CLOSE(umaf)
  
  !! Announce the end of writing
  CALL Prnt0("Finished.", log=.FALSE.)
  CALL Prnt("Allele frequencies saved into ["//TRIM(maf_file)//"]")

  RETURN
  
END SUBROUTINE SaveMAF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SaveInputData(ped_file, X, FAM, ss, save_binary, fam_file, &
                         map_file, MAPA, map_sep)
  IMPLICIT NONE
  CHARACTER, PARAMETER                  :: LF = CHAR(10), CR = CHAR(13)
  CHARACTER(*), INTENT(INOUT)           :: ped_file
  INTEGER(iks), INTENT(IN)              :: X(:,:)
  CHARACTER, INTENT(IN)                 :: FAM(:,:)
  INTEGER, INTENT(IN)                   :: ss
  LOGICAL, INTENT(IN)                   :: save_binary
  CHARACTER(*), INTENT(INOUT), OPTIONAL :: fam_file, map_file
  CHARACTER(*), INTENT(IN), OPTIONAL    :: MAPA(:,:), map_sep
  CHARACTER(mfl)                        :: fam_file1, map_file1
  CHARACTER(mifll)                      :: line
  INTEGER                               :: nrec, k, l, i, j, ifam, nloci, &
                                           save_ncols, nlines !, ios
  INTEGER(iks), ALLOCATABLE             :: Y(:,:)
  INTEGER(1)                            :: b3
  CHARACTER(3)                          :: text
  CHARACTER                             :: Alleles(SIZE(X,2),2), A1, A2, sep
  
  !! If filename missing throw an error
  IF(LEN_TRIM(ped_file)==0) &
    CALL PrntE("Missing name for ped data output file. When"//&
                    " --saveinput is present and the data is not simulated,"//&
                    " the root name of the output ped file must be"//&
                    " specified (use --savefile).", Q=.TRUE.)

  !! Get other filenames
  fam_file1 = "" 
  map_file1 = ""
  IF(PRESENT(fam_file)) fam_file1 = fam_file
  IF(PRESENT(map_file)) map_file1 = map_file
  
  nloci = SIZE(X,2)
  nlines = SIZE(FAM,1)
  
  sep = fam_cs(1)
  
  !! ************************************** !!
  !!   SAVE PEDDING DATA IN BINARY FORMAT   !!
  !! ************************************** !!
  IF(save_binary) THEN

    !! Announce what's happening
    IF(out_bed_data_SNPmajor) THEN
      CALL Prnt("Saving ped data into file in a binary SNP-major"//&
                     " format ... ", ADVANCE='NO')
    ELSE
      CALL Prnt("Saving ped data into file in a binary individual"//&
                     "-major format ... ", ADVANCE='NO')
    ENDIF    


    !! Save BED data into a file
    IF(out_bed_data_SNPmajor) THEN
      b3 = 1
      !! Open binary file for writing with the correct number of records
      !! Direct access needed to avoid those DAMN HEADERS!
      nrec = SIZE(X,1)*SIZE(X,2) 
      CALL OpenFile(UNIT=ubed, FILE=ped_file, ACTION='W', STATUS='R', &
                    FORM='U', ACCESS='D', RECL=nrec+3, chck_exist=.TRUE., &
                    overwrite=overwrite_files)
                    
      !! Write the entire array into file at once
      WRITE(ubed, REC=1) bed_byte1, bed_byte2, b3, X
    ELSE
      b3 = 0
      CALL ResizeVar(Y, (nloci+3)/4, ss)
      CALL ConvertBEDmajor(X,Y)

      !! Open binary file for writing with the correct number of records
      !! Direct access needed to avoid those DAMN HEADERS!
      nrec = SIZE(Y,1)*SIZE(Y,2) 
      CALL OpenFile(UNIT=ubed, FILE=ped_file, ACTION='W', STATUS='R', &
                    FORM='U', ACCESS='D', RECL=nrec+3, chck_exist=.TRUE., &
                    overwrite=overwrite_files)
                    
      !! Write the entire array into file at once
      WRITE(ubed, REC=1) bed_byte1, bed_byte2, b3, Y
    ENDIF

    CLOSE(ubed)

    !! Check for errors during reading
    !IF(ios /= 0) GOTO 12

    !! Announce the end of writing
    CALL Prnt0("Finished.", log=.FALSE.)
    
    !! ************************************************* !!
    !!  WRITE THE FAMILY INFORMATION INTO SEPARATE FILE  !!
    !! ************************************************* !!

    !! Check for proper family file name 
    IF(LEN_TRIM(fam_file1)==0) &
      CALL PrntE("Output pedigree file name missing!", Q=.TRUE.)

    !! Announce what's happening
    CALL Prnt("Saving pedigree data into file ... ", ADVANCE='NO')

    !! Open the file for writing
    CALL OpenFile(UNIT=ufam, FILE=fam_file1, ACTION='W', STATUS='U', FORM='F', &
                  ACCESS='S', chck_exist=.TRUE., overwrite=overwrite_files)
    
    !! Write the pedigree (family) info to a file
    DO k=1,nlines
      
      !! Concatenate the raw pedigree data into a single line
      line = ""
      DO l=1,SIZE(FAM,2)
        IF(FAM(k,l)==CR .OR. FAM(k,l)==LF) EXIT
        line(l:l) = FAM(k,l)
      ENDDO
      
      !! Write the line
      WRITE(ufam,'(A)') TRIM(line)
      
    ENDDO
    
    CLOSE(ufam)

    !! Announce the end of writing
    CALL Prnt0("Finished.", log=.FALSE.)

  !! ************************************* !!
  !!   SAVE PEDDING DATA IN PLAIN FORMAT   !!
  !! ************************************* !!
  ELSE

    !! Announce what's happening
    CALL Prnt("Saving ped data into file in a plain format ... ", ADVANCE='NO')

    !! Open the file for writing
    CALL OpenFile(UNIT=uped, FILE=ped_file, ACTION='W', STATUS='U', FORM='F', &
                  ACCESS='S', chck_exist=.TRUE., overwrite=overwrite_files)
                  
    !! Write the plain format ped file
    i = 0
    ifam = 0
    CALL ResizeVar(Y, 4, nloci)

    !! Set either generic allele names or get them from MAPA if large enough
    Alleles(:,1) = "G" 
    Alleles(:,2) = "A"
    IF(PRESENT(MAPA)) THEN
      IF(SIZE(MAPA,2) >= MAX(map_a1col,map_a2col)) THEN
        Alleles(:,1) = MAPA(:,map_a1col) 
        Alleles(:,2) = MAPA(:,map_a2col)
      ENDIF
    ENDIF

    !! Loop over all rows of X (four individuals per cell)
    DO k=1,SIZE(X,1)

      !! Convert for each row all cells into ped data 
      DO l=1,nloci
        CALL Bed2Ped(X(k:k,l), Y(:,l), bed_major, ped_minor)
      ENDDO
      
      !! Print all loci records for each of the current (at most) 4 individuals
      DO j=1,4
      
        i = i+1
        ifam = ifam+1
      
        !! If all individuals were written, exit the inner loop
        !! NOTE: This can happen only on the last row of X and if the sample 
        !!       size is not a multiple of 4 
        IF(i > ss) EXIT
        
        !! First check if ifam would exceed bounds of FAM, which could happen
        !! if ped_nsamp>1 and save_full_fam is false
        IF(ifam > nlines) ifam = 1
        
        !! Concatenate the raw pedigree data into a single line
        line = ""
        DO l=1,SIZE(FAM,2)
          IF(FAM(k,l)==CR .OR. FAM(k,l)==LF) EXIT
          line(l:l) = FAM(k,l)
        ENDDO
        
        !! Write the line
        WRITE(uped,'(A)', ADVANCE='NO') TRIM(line)
      
        !!! Then write the family info
        !DO l=1,fam_ncols
        !  IF(l>1) WRITE(uped,'(A)', ADVANCE='NO') sep
        !  WRITE(uped,'(A)', ADVANCE='NO') TRIM(FAM(ifam,l))
        !ENDDO

        !! Then the genes of each individual (minor allele being "G")
        DO l=1,nloci
        
          IF(fam_ncols>0 .OR. l>1) WRITE(uped,'(A)', ADVANCE='NO') sep

          !! Get the allele names
          IF(ped_minor) THEN
            A1 = Alleles(l,1) 
            A2 = Alleles(l,2)
          ELSE
            A1 = Alleles(l,2) 
            A2 = Alleles(l,1)
          ENDIF 

          !! Write the two alleles
          SELECT CASE (Y(j,l))
          CASE (0_iks)
            text = A2//sep//A2
          CASE (1_iks)
            text = A1//sep//A2
          CASE (2_iks)
            text = A1//sep//A1
          CASE DEFAULT
            text = def_CharNAgen//sep//def_CharNAgen
          END SELECT

          WRITE(uped,'(A)', ADVANCE='NO') TRIM(text)
           
        ENDDO

        !! And finally write the end-of-line symbol
        WRITE(uped,'(A)') ""
        
      ENDDO

    ENDDO
    
    CLOSE(uped)

    !! Check for errors during reading
    !IF(ios /= 0) GOTO 12

    !! Announce the end of writing
    CALL Prnt0("Finished.", log=.FALSE.)
  ENDIF

  !! ******************** !!
  !!   SAVE MAPPING DATA  !!
  !! ******************** !!
  
  !! Write also the mapping information into a file
  IF(PRESENT(MAPA)) THEN
    !! Check for proper mapping file name 
    IF(LEN_TRIM(map_file1)==0) &
      CALL PrntE("Output mapping file name missing!", Q=.TRUE.)
    
    !! Check for proper mapping file name 
    IF(.NOT.PRESENT(map_sep)) &
      CALL PrntE("Missing mapping file separator!", Q=.TRUE.)

    !! Announce what's happening
    CALL Prnt("Saving mapping data into file ... ", ADVANCE='NO')

    !! Open the file for writing
    CALL OpenFile(UNIT=umap, FILE=map_file1, ACTION='W', STATUS='U', &
                  FORM='F', ACCESS='S', chck_exist=.TRUE., &
                  overwrite=overwrite_files)
                  
    !! Decide how many columns will be saved
    IF(save_binary) THEN
      save_ncols = MIN(def_bim_ncols, SIZE(MAPA,2))
    ELSE
      save_ncols = MIN(def_map_ncols, SIZE(MAPA,2))
    ENDIF
    
    !! Write the mapping info
    DO k=1,SIZE(MAPA,1)
      line = ""
      DO l=1,save_ncols
        IF(l>1) line = TRIM(line)//map_sep
        line = TRIM(line)//MAPA(k,l) 
      ENDDO
      WRITE(umap,'(A)') TRIM(line)
    ENDDO
    
    CLOSE(umap)
    
    !! Announce the end of writing
    CALL Prnt0("Finished.", log=.FALSE.)
    
  ENDIF

  !! *************************************** !!
  !!   ANNOUNCE DETAIL INFO ON OUTPUT FILES  !!
  !! *************************************** !!

  CALL Prnt("Data were saved into the files")
  IF(save_binary) THEN
    CALL Prnt0("  ["//TRIM(ped_file)//"]", lead=2)
    CALL Prnt0("  ["//TRIM(fam_file1)//"]", lead=2)
  ELSE
    CALL Prnt0("  ["//TRIM(ped_file)//"]", lead=2)
  ENDIF
  
  IF(PRESENT(MAPA)) THEN
    IF(save_binary) THEN
      CALL Prnt0("  ["//TRIM(map_file1)//"]", lead=2)
      IF(SIZE(MAPA,2) < def_bim_ncols) &
        CALL PrntW("Simulated BIM file does not have the standard number of"//&
                   " columns ("//TRIM(i2cp(def_bim_ncols))//")!")
    ELSE
      CALL Prnt0("  ["//TRIM(map_file1)//"]", lead=2)
      IF(SIZE(MAPA,2) < def_map_ncols) &
        CALL PrntE("Simulated MAP file does not have the standard number of"//&
                    " columns ("//TRIM(i2cp(def_map_ncols))//")!", Q=.TRUE.)
    ENDIF
  ELSE
    CALL PrntW("No MAP file produced")
  ENDIF

  RETURN

END SUBROUTINE SaveInputData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ConvertBEDmajor(IN,OUT)
!! This routine converts binary ped data (BED) from individual-major format
!! to SNP-major format and vice-versa depending on what the format of IN is
!! No distinction of these two formats in the code is needed because the 
!! conversion between the two formats is symmetric.
  IMPLICIT NONE
  INTEGER(iks), PARAMETER    :: NAbed = 85
  INTEGER(iks), INTENT(IN)  :: IN(:,:)
  INTEGER(iks), INTENT(OUT) :: OUT(:,:)
  INTEGER(iks)              :: Z(SIZE(OUT,2),0:3)
  INTEGER                  :: nrowOUT, ncolOUT, nrowIN, ncolIN, i, j, k, up

  !! Arrays sizes
  nrowOUT = SIZE(OUT,1)
  ncolOUT = SIZE(OUT,2)
  nrowIN  = SIZE(IN,1)
  ncolIN  = SIZE(IN,2)
  
  !! Check for incompatible dimensions of input and output
  IF(nrowIN/=(ncolOUT+3)/4 .OR. nrowOUT/=(ncolIN+3)/4) &
    CALL PrntE("Incompatible dimension (ConvertBEDmajor)", Q=.TRUE.)

  !! Do the conversion
  OUT = NAbed 
  i = 1
  DO j=1,ncolIN-1,4
    up = MIN(3,ncolIN-j)
    !! First decompress the next (up to) 4 individuals 
    DO k=0,up
      CALL Bed2Ped(IN(:,j+k),Z(:,k), bed_major, ped_minor)
    ENDDO
    !! Recode the loci of the 4 indivuals
    CALL Ped2Bed(TRANSPOSE(Z(:,0:up)),OUT(i:i,:), bed_major, ped_minor)
    i=i+1
  ENDDO

  RETURN
  
END SUBROUTINE ConvertBEDmajor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SimulData(X, sts, sex, FAM, bed_data, model, OR, OR1, OR2, prev, &
                     LD, nloci, nsamp, nca, nco, maf_fixed, maf1, maf2, maf3, &
                     nmale, nfema, silent)
  IMPLICIT NONE
  CHARACTER, INTENT(OUT), ALLOCATABLE    :: FAM(:,:)
  INTEGER(iks), INTENT(OUT), ALLOCATABLE  :: X(:,:), sts(:), sex(:)
  REAL(dpp), INTENT(IN)                  :: OR, OR1, OR2, prev, LD
  INTEGER, INTENT(IN)                    :: nloci, nsamp, nca, nco, model
  LOGICAL, INTENT(IN)                    :: bed_data
  INTEGER, INTENT(OUT), OPTIONAL         :: nmale, nfema
  REAL(dpp), INTENT(IN), OPTIONAL        :: maf1, maf2, maf3
  LOGICAL, INTENT(IN), OPTIONAL          :: maf_fixed, silent
  INTEGER(iks)                            :: X1((nca+nco+3)/4, nloci)
  INTEGER(iks), ALLOCATABLE               :: Y(:,:)
  CHARACTER(mmtl)                        :: text, line
  INTEGER                                :: i, j, k, l, ss, nsamp1, rc, nc, &
                                            fnlines, nindiv_total
  INTEGER(ikb)                           :: pct
  LOGICAL                                :: fam_ncols_decreased, silent1, mfix
  REAL(dpp)                              :: mafA, mafB, mafC, coef(3,3), &
                                            prev_actual
  
  !! Select announcement mode
  silent1 = .FALSE.
  IF(PRESENT(silent)) silent1 = silent
  
  !! Initialize actual observed prevalence variable and total simulated 
  !! individuals counter
  prev_actual = -one
  nindiv_total = -1
  
  !! Select whether maf will be fixed: use global variable fixed_MAF by default
  !! and override it if maf_fixed present
  mfix = fixed_MAF
  IF(PRESENT(maf_fixed)) mfix = maf_fixed 

  !! Select minor allele frequencies
  mafA = allelefreqA
  mafB = allelefreqB
  mafC = allelefreqC
  IF(PRESENT(maf1)) mafA = maf1
  IF(PRESENT(maf2)) mafB = maf2
  IF(PRESENT(maf3)) mafC = maf3
  
  !! Announce mafs used in simulation
  text = "MAFs used in simulation: "//TRIM(r2c(mafA))//" "//TRIM(r2c(mafB))
  IF(nneutralloci>0) text = TRIM(text)//" "//TRIM(r2c(mafC)) 
  CALL Prnt(text)                   
    
  !! Announcement
  CALL Prnt("Simulating input data ... ", ADVANCE='NO', screen=.NOT.silent1)

  !! Get number of samples
  ss = nca + nco
  rc = (ss+3)/4
  
  fam_single_sample = .FALSE.
  
  CALL ResizeVar(X, nsamp * rc, nloci)
  CALL ResizeVar(sts, nsamp * ss, NumNA)
  CALL ResizeVar(sex, nsamp * ss, NumNA)

  !! If full fam is to be saved, then use all nsamp. Otherwise use only the
  !! first sample in sts
  IF(save_full_fam) THEN
    nsamp1 = nsamp
  ELSE
    nsamp1 = 1
  ENDIF 
  fnlines = ss*nsamp1

  !! Check for invalid input values
  IF(OR<=zero) &
    CALL PrntE("Invalid OR ("//r2c(OR)//").", Q=.TRUE.)
  IF(OR1<=zero) &
    CALL PrntE("Invalid OR1 ("//r2c(OR1)//").", Q=.TRUE.)
  IF(OR2<=zero) &
    CALL PrntE("Invalid OR2 ("//r2c(OR2)//").", Q=.TRUE.)
  !IF(prev<=zero) &
  !  CALL PrntE("Invalid prevalence ("//r2c(prev)//").", Q=.TRUE.)
  
  !! Set the variables for countdown
  text = ""
  pct = -1
  
  !! Check for too many family file columns for memory reasons
  !! If too big, try to decrease it to minimum
  fam_ncols_decreased = .FALSE.
  IF(save_input_data) THEN
    
    IF(INT(fnlines,ikb)*mfml*fam_ncols > mem_warn_limit) THEN
      fam_ncols = 1
      fam_stscol = 1
      fam_sexcol = 0
      fam_ncols_decreased = .TRUE. 
    ENDIF
    
    CALL ResizeVar(FAM, fnlines, 30, "")
    
  ELSE
    fam_ncols = 0
  ENDIF
  
  !! Generate case/control sts
  DO k=1,nsamp1
    sts(1+(k-1)*ss:k*ss) = (/(NumCa,i=1,nca), (NumCo,j=1,nco)/)
  ENDDO
  
  !! Generate sex sts
  sex = NumNA
  IF(assume_male) THEN
    sex = NumMa
    IF(PRESENT(nmale)) nmale = nco + nca
    IF(PRESENT(nfema)) nfema = 0
  ELSEIF(assume_fema) THEN
    sex = NumFe
    IF(PRESENT(nfema)) nfema = nco + nca
    IF(PRESENT(nmale)) nmale = 0
  ENDIF
  
  !! Generate the rest of family info
  IF(fam_ncols>0) THEN
  
    DO i=1,fnlines
    
      line = ""
      DO j=1,fam_ncols
        
        !! Start with NAs for all
        text = "0"
        
        !! Write Family and Individual ID columns, sex, sts.
        IF(j==fam_fidcol) THEN
          text = i2cp(sim_fid_start + i - 1)
  
        ELSEIF(j==fam_iidcol) THEN
          text = i2cp(sim_iid_start + i - 1)
  
        ELSEIF(j==fam_stscol) THEN
          IF(sts(i)==NumCo) THEN
            text = CharCo
          ELSE
            text = CharCa
          ENDIF
  
        ELSEIF(j==fam_sexcol) THEN
          IF(sex(i)==NumMa) THEN
            text = CharMa
          ELSE
            text = CharFe
          ENDIF
        ENDIF
        
        if(j>1) line = TRIM(line)//fam_cs(1)
        line = TRIM(line)//TRIM(text)      
          
      ENDDO
  
      !! Store the line into FAM
      nc = LEN_TRIM(line)
      IF(nc>SIZE(FAM,2)) CALL ResizeVar(FAM, fnlines, nc, "")
      DO l=1,nc
        FAM(i,l) = line(l:l)
      ENDDO
       
    ENDDO
  
  ENDIF

  !! Assign initial values to X
  X = 0_iks
  text = ""
  pct = -1

  !! Get the interaction table required just below
  CALL GetInteractionTable(coef, model)
    
!!*******************************************************!!
!!  TWO VERSIONS OF CODE: PEDDING OR BEDDING SIMULATION  !!
!!*******************************************************!!
IF(1==2) THEN  
  
  !! If multiple samples, then generate in 
  IF(nsamp>1) THEN

    !! Allocate array for genotypes in ped format
    CALL ResizeVar(Y, nsamp*ss, nloci)

    !! Generate all of the samples at once
    CALL GenerateCaCo(nsamp*nca, nsamp*nco, mafA, mafB, mafC, &
                      mfix, OR, OR1, OR2, LD, coef, prev, nneutralloci, &
                      nLDpairs, HWE=simulate_HWE, binary=.FALSE., &
                      X=Y, prev_actual=prev_actual, nindiv_total=nindiv_total)
                
    !! Convert the ped samples into binary representation sample by sample
    CALL Prnt0("")
    CALL Prnt("Converting to binary format ... ", ADVANCE='NO', screen=.NOT.silent1)
    DO k = 1,nsamp
      CALL Ped2Bed(Y(1+(k-1)*ss:k*ss,:), X(1+(k-1)*rc:k*rc,:), .TRUE., .TRUE.)
    ENDDO

  ELSE
  
    !! Generate the single sample already in bed format
    CALL GenerateCaCo(nca, nco, mafA, mafB, mafC, fixed_MAF, OR, OR1, & 
                      OR2, LD, coef, prev, nneutralloci, nLDpairs, &
                      HWE=simulate_HWE, binary=bed_data, X=X, &
                      prev_actual=prev_actual, nindiv_total=nindiv_total)                        

  ENDIF

  CALL Prnt0("")
  CALL Prnt("Data simulation finished.")

ELSE

  !!***********************!!
  !! START PARALLEL REGION !!
  !!***********************!!
  !!$OMP PARALLEL DEFAULT(NONE) &
  !!$OMP FIRSTPRIVATE(nsamp, nca, nco, prev) &
  !!$OMP FIRSTPRIVATE(OR, OR1, OR2, LD, model, fixed_MAF, min_MAF) &
  !!$OMP FIRSTPRIVATE(nneutralloci, nLDpairs, countdown, text) &
  !!$OMP FIRSTPRIVATE(pct) &
  !!$OMP PRIVATE(k, genes, genes1, X1) &
  !!$OMP SHARED(X)
  
   
    !! Simulate samples      
    DO k = 1,nsamp
  
  !!    IF(OMP_GET_THREAD_NUM()==0) &
      !IF(nsamp>1) &
        IF(.NOT.silent1) &
          CALL PrintCD(text, k-1, nsamp, pct, finish=.TRUE., adv=.FALSE.)
    
      !! Generate the next sample
      CALL GenerateCaCo(nca, nco, mafA, mafB, mafC, fixed_MAF, OR, &
                        OR1, OR2, LD, coef, prev, nneutralloci, &
                        nLDpairs, HWE=simulate_HWE, binary=bed_data, X=X1, &
                        prev_actual=prev_actual, nindiv_total=nindiv_total)
  
  !!$OMP CRITICAL
      X(1+(k-1)*rc:k*rc,:) = X1
  !!$OMP END CRITICAL
  
  !!    IF(OMP_GET_THREAD_NUM()==0) &
      !IF(nsamp>1) &
    ENDDO
  
  !!$OMP END DO
  !!$OMP END PARALLEL 
    
  IF(.NOT.silent1) THEN
    CALL PrintCD(text, nsamp, nsamp, pct, finish=.TRUE.)
    CALL Prnt0("Finished.", screen=.NOT.countdown, log=.FALSE.)
  ENDIF
  CALL Prnt("Data simulation finished.", screen=.FALSE.)

ENDIF

  CALL Prnt("Total number of simulated individuals: "//i2c(nindiv_total))
  CALL Prnt("Actual observed prevalence of cases: "//r2c(prev_actual))
  
  !! Output the actual observed prevalence of cases into file
  IF(do_out_actual_prev) CALL SavePrev(out_prev_file, prev_actual)

!!*******************************************************!!
!!  TWO VERSIONS OF CODE: PEDDING OR BEDDING SIMULATION  !!
!!*******************************************************!!

  !! Warn about the dereased number of pedigree columns
  IF(fam_ncols_decreased) &
    CALL PrntW("The number of pedigree columns decreased to 1 to save memory.")

  RETURN
      
END SUBROUTINE SimulData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE AdjustChrNames(MAPA)
!! This routine does a check for non-numerical names of chromosomes and if it
!! finds any, changes them no numerical ones
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT), TARGET :: MAPA(:,:)
  CHARACTER(LEN(MAPA)), POINTER       :: chr
  INTEGER                             :: i  
  
  DO i=1,SIZE(MAPA,1)
    chr => MAPA(i,map_chcol)
    IF(chr == "X") THEN
      chr = i2cp(special_chr(1))
    ELSEIF(chr == "Y") THEN
      chr = i2cp(special_chr(2))
    ELSEIF(chr == "XY") THEN
      chr = i2cp(special_chr(3))
    ELSEIF(chr == "MT") THEN
      chr = i2cp(special_chr(4))
    ENDIF 
  ENDDO
  
  RETURN
  
END SUBROUTINE AdjustChrNames
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ComputeMAF(MAF, X, sex, MAPA, InclLoc, n_maf_excl, &
                      min_maf, maf_keep, warn_maf)
  IMPLICIT NONE
  REAL(dpp), INTENT(OUT), ALLOCATABLE :: MAF(:,:)
  INTEGER(iks), INTENT(IN)            :: X(:,:)
  INTEGER(iks), INTENT(IN), OPTIONAL  :: sex(:)
  CHARACTER(*), INTENT(IN), OPTIONAL  :: MAPA(:,:)
  LOGICAL, INTENT(INOUT), OPTIONAL, ALLOCATABLE :: InclLoc(:)
  INTEGER, INTENT(OUT), OPTIONAL      :: n_maf_excl
  REAL(dpp), INTENT(IN), OPTIONAL     :: min_maf
  LOGICAL, INTENT(IN), OPTIONAL       :: maf_keep(:), warn_maf
  REAL(dpp)                           :: maf1, avgmaf1
  INTEGER(iks), ALLOCATABLE            :: x_all(:), sex1(:)
  CHARACTER(mml), ALLOCATABLE         :: MAPAchr(:)
  LOGICAL, ALLOCATABLE                :: InclLoc1(:)
  INTEGER, ALLOCATABLE                :: chr(:)
  INTEGER                             :: i, j, n_maf_half, n_maf_chck, &
                                         NCHROBS, n_maf_excl1, nloci
  INTEGER(ikb)                        :: pct
  INTEGER(iks), ALLOCATABLE           :: X1(:)
  INTEGER(ikb)                        :: maxim, im
  CHARACTER(mmtl)                     :: text, stext
  LOGICAL                             :: warn_maf1

  nloci = SIZE(X,2)
  
  warn_maf1 = .TRUE.
  IF(PRESENT(warn_maf)) warn_maf1 = warn_maf
  
  !! Initialize MAF, where rows are the number of loci, columns are the number
  !! of samples times 2, where for each sample we store maf (1:ped_nsamp) and  
  !! the non-missing allele count (ped_nsamp+1:2*ped_nsamp)
  CALL ResizeVar(MAF, nloci, 2*ped_nsamp, NAneg)

  !! Initialize InclLoc1 to length equal to the number of loci
  CALL ResizeVar(InclLoc1, nloci, .TRUE.)

  !! Allocate InclLoc and set inclusion so that initially all loci included
  IF(PRESENT(InclLoc)) THEN
    IF(.NOT.ALLOCATED(InclLoc)) &
      CALL ResizeVar(InclLoc, nloci, .TRUE.)
    IF(SIZE(InclLoc)/=nloci) &
      CALL PrntE("InclLoc has wrong size (ComputeMAF)!", Q=.TRUE.)
    InclLoc1 = InclLoc
  ENDIF
  
  !! *********************************************************************** !!
  !!              MINOR ALLELE FREQUENCY COMPUTATION AND CHECK               !!
  !! *********************************************************************** !!

  !! Store chromosome information
  CALL ResizeVar(chr, nloci, 0)
  IF(PRESENT(MAPA)) THEN
    CALL ResizeVar(MAPAchr, nloci)
    MAPAchr = MAPA(:,map_chcol) 
    chr = c2i(MAPAchr)
  ELSE
    CALL PrntW("Minor allele frequency computation cannot take into"//&
            " account ANY CHROMOSOME INFORMATION because chromosome"//&
            " information is missing.", skip1=1, skip2=1)    
  ENDIF

  !! Allocate x_all and sex1
  CALL ResizeVar(x_all, ped_ss)
  CALL ResizeVar(sex1, ped_ss, 0_iks)
  CALL ResizeVar(X1, (ped_ss+3)/4)


  !! Check for presence of sex data
  IF(.NOT.PRESENT(sex)) &
    CALL PrntW("MAFs cannot be computed properly for SEX CHROMOSOMES due to"//&
               " missing sex information.", skip1=1, skip2=1)

  im = 0
  maxim = nloci

  !! Announce current action
  stext = "Calculating allele frequencies ..."
  !CALL Prnt(stext, screen=.FALSE.)

  !! Reset counter of loci checked and the counter of loci with maf above 0.5
  n_maf_chck = 0
  n_maf_half = 0
  n_maf_excl1 = 0

  !! Perform the check for low MAF for all loci
  DO i=1,nloci

    !! Print count down progress
    CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext, cdstep=5_ikb)

    !! Increase counter of checked loci
    im = im + 1

    !! Don't bother checking maf exclusion if loci already excluded
    IF(.NOT.InclLoc1(i)) CYCLE
    
    !!! Get ped data to compute MAF
    !CALL Bed2Ped(X(:,i), x_all, bed_major, ped_minor, chr(i), sex1, maf1, NCHROBS)

    !! Initialize counter of average maf
    avgmaf1 = zero
    
    IF(PRESENT(sex)) THEN
      IF(ALL(SIZE(sex)/=ped_ss*(/1,ped_nsamp/))) &
        CALL Prnt0("Invalid size of 'sex' (ComputeMAF).", Q=.TRUE.)
    ENDIF
    
    !! Calculated MAFs for different samples
    DO j=1,ped_nsamp
    
      !! Extract sex information
      IF(PRESENT(sex)) THEN
        IF(SIZE(sex)==ped_ss) THEN
          sex1 = sex
        ELSE
          sex1 = sex((j-1)*ped_ss+1:j*ped_ss)
        ENDIF
      ENDIF
      
      !! Extract gene data for i-th locus and j-th sample
      X1 = X(((j-1)*ped_ss+3)/4+1:(j*ped_ss+3)/4,i)  
  
      !! Get ped data to compute MAF
      CALL Bed2Ped(X1, x_all, bed_major, ped_minor, chr(i), sex1, maf1, NCHROBS, maf_keep)

      !! Check for erroneous maf
      n_maf_chck = n_maf_chck + 1 
      IF(maf1 > half) n_maf_half = n_maf_half + 1
  
      !! Save computed MAF and the non-missing allele count (NCHROBS)
      MAF(i,j) = maf1
      MAF(i,ped_nsamp+j) = REAL(NCHROBS, dpp)
      
      !! Add the current maf to the average maf counter 
      avgmaf1 = avgmaf1 + maf1 / ped_nsamp
  
    ENDDO
    
    !! Check for too low MAF and exclude
    IF(PRESENT(min_maf)) THEN
      IF(avgmaf1 < min_maf) THEN
        n_maf_excl1 = n_maf_excl1 + 1 
        InclLoc1(i) = .FALSE.
      ENDIF
    ENDIF
    
  ENDDO
  
  !! Announce the end
  CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext)
  CALL Prnt("Finished calculating allele frequencies.", screen=.FALSE.)
  
  !! Return the number of excluded loci based on MAF
  IF(PRESENT(n_maf_excl)) n_maf_excl = n_maf_excl1
  IF(PRESENT(InclLoc)) InclLoc = InclLoc1 

  !! Check for MAFs larger than 0.5
  IF(n_maf_half > 0 .AND. warn_maf1) &
    CALL PrntW("There were "//i2c(n_maf_half)//" loci ("//&
            TRIM(r2cp(hundred*round(DBLE(n_maf_half)/n_maf_chck, 5)))//&
            "% of "//i2cp(n_maf_chck)//" loci checked) with"//&
            " MAF above 0.5!")
  
  !! Check for min_MAF larger than 0.5
  IF(min_maf==half) & 
    CALL PrntW("MAF lower bound is set to 0.5! Only loci with both"//&
            " allele frequency equal will pass the minimum MAF test!")
  
  !!! Announce warning if multiple samples present
  !IF(ped_nsamp>1) &
  !  CALL PrntW("For multiple samples MAFs are calculated and checked"//&
  !          " only for loci of the first sample and these values are"//&
  !          " assummed and reported for all loci of the remaining samples!")
  
  RETURN
  
END SUBROUTINE ComputeMAF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE CheckLociExcl(InclLoc, MAPA, SUBMAPA, MAPARS, SUBMAPARS, &
                         n_chr_excl, n_rs_excl, n_negbp_excl, n_sub_excl)
  IMPLICIT NONE
  LOGICAL, INTENT(OUT), ALLOCATABLE     :: InclLoc(:)
  CHARACTER(*), INTENT(INOUT), OPTIONAL, ALLOCATABLE :: MAPA(:,:), SUBMAPA(:,:)
  INTEGER(ikb), INTENT(IN), OPTIONAL, ALLOCATABLE    :: MAPARS(:), SUBMAPARS(:)
  INTEGER, INTENT(OUT), OPTIONAL        :: n_chr_excl, n_negbp_excl, &
                                           n_rs_excl, n_sub_excl

  LOGICAL                               :: Found, numerical_compare
  CHARACTER(mml), ALLOCATABLE           :: SUBrs(:)
  INTEGER(iks), ALLOCATABLE              :: is_loc(:)
  INTEGER, ALLOCATABLE                  :: order(:), keep(:) 
  CHARACTER(mml)                        :: MAPArs_i, MAPAbp_i
  INTEGER                               :: i, j, k, IncP1(2), order1(2), &
                                           whr(2), nfound, nloci, &
                                           nloci_sub, npairs, &
                                           n_sub_excl1, n_chr_excl1, &
                                           n_sub_incl1, n_rs_excl1, &
                                           n_negbp_excl1, n_incl_loci, &
                                           max_sub_incl, MAPAchr_i, &
                                           n_chr_excl_bd(-2:SIZE(excluded_chr)-1)
  INTEGER(ikb)                          :: pct, maxim, im
  CHARACTER(mmtl)                       :: text, stext
  
  nloci = SIZE(X,2)
  
  CALL Prnt("Processing input data ...")

  !! Allocate InclLoc and set inclusion so that initially all loci included
  CALL ResizeVar(InclLoc, nloci, .TRUE.)

  !! If MAPA not present of unallocated, skip this part and go directly to mafs
  IF(.NOT.PRESENT(MAPA)) GOTO 20
  IF(.NOT.ALLOCATED(MAPA)) GOTO 20

  !! *********************************************************************** !!
  !!                 CHROMOSOME AND BP EXCLUSION CHECK                       !!
  !! *********************************************************************** !!

  IF(.NOT.(PRESENT(n_chr_excl) .OR. PRESENT(n_negbp_excl) .OR. &
  PRESENT(n_rs_excl))) GOTO 10

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!    CHECK FOR EXCLUSION BASED ON CHROMOSOMES AND NEGATIVE BP   !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize
  im = 0
  nfound = 0
  n_chr_excl1 = 0 
  n_negbp_excl1 = 0
  n_rs_excl1 = 0
  maxim = nloci

  n_incl_loci = CountTrue(InclLoc)
  CALL Prnt("Currently "//i2cp(n_incl_loci)//" loci are included.")

  !! Exclude/include based on chromosome and bp numbers
  IF(n_incl_loci>0) THEN
  
    n_chr_excl_bd = 0

    !! Announce the next action
    stext = "Excluding loci based on chromosome, BP and RS numbers ..."
    !CALL Prnt(stext, screen=.FALSE.)
    
    !! Disable negative BP check if BP numbers are not available
    IF(map_bpcol > SIZE(MAPA,2)) keep_all_bp = .FALSE.
  
    !! Loop over loci
    DO i=1,nloci
  
      !! Print count down progress
      CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext)
  
      !! Increase counter of checked loci
      im = im + 1
  
      !! If current loci already excluded (for some reason) skip DO cycle
      IF(.NOT.InclLoc(i)) CYCLE
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!              CHROMOSOME EXCLUSION              !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(.NOT.keep_all_chr) THEN
        !! Check for non-numerical chromosome number
        IF(.NOT.IsInteger(MAPA(i,map_chcol))) THEN
          n_chr_excl_bd(-2) = n_chr_excl_bd(-2) + 1
          InclLoc(i) = .FALSE.
          n_chr_excl1 = n_chr_excl1 + 1
          CYCLE
        ENDIF
        MAPAchr_i = c2i(MAPA(i,map_chcol))
        
        !! Check whether the chromosome number is too low or too high
        IF(MAPAchr_i<0 .OR. MAPAchr_i>(SIZE(excluded_chr)-1)) THEN
          n_chr_excl_bd(-1) = n_chr_excl_bd(-1) + 1
          InclLoc(i) = .FALSE.
          n_chr_excl1 = n_chr_excl1 + 1
          CYCLE
        
        !! Check whether the chromosome is on the list of excluded chromosomes
        ELSEIF(ANY(excluded_chr)) THEN
          IF(excluded_chr(MAPAchr_i)) THEN
            n_chr_excl_bd(MAPAchr_i) = n_chr_excl_bd(MAPAchr_i) + 1
            InclLoc(i) = .FALSE.
            n_chr_excl1 = n_chr_excl1 + 1
            CYCLE
          ENDIF
        ENDIF
      ENDIF
  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!          NEGATIVE BASE PAIR EXCLUSION          !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(.NOT.keep_all_bp) THEN
        !! Check for negative base-pair number and if found, exclude the loci
        MAPAbp_i = MAPA(i,map_bpcol)
        IF(MAPAbp_i(1:1) == "-") THEN
          InclLoc(i) = .FALSE.
          n_negbp_excl1 = n_negbp_excl1 + 1
          CYCLE
        ENDIF
      ENDIF
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!           INVALID RS NUMBER EXCLUSION          !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Check for negative base-pair number and if found, exclude the loci
      IF(rs_valid_chck .AND. exclude_invalid_rs) THEN
        MAPArs_i = MAPA(i,map_rscol)
        IF(upcasef(MAPArs_i(1:2)) /= "RS" .OR. .NOT.IsInteger(MAPArs_i(3:))) THEN
          InclLoc(i) = .FALSE.
          n_rs_excl1 = n_rs_excl1 + 1
          CYCLE
        ENDIF
        
        IF(IsInteger(MAPArs_i(3:))) THEN
          IF(c2i(MAPArs_i(3:))<=0) THEN
            InclLoc(i) = .FALSE.
            n_rs_excl1 = n_rs_excl1 + 1
            CYCLE
          ENDIF
        ENDIF
      ENDIF
  
    ENDDO
  
    !! Announce the end
    CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext)
    CALL Prnt("Finished excluding loci based on chromosome, BP and RS numbers.")
    
    !! Announce chromosome exclusion breakdown
    IF(ANY(n_chr_excl_bd>0)) THEN
      CALL Prnt("Chromosome exclusion breakdown:")
      IF(n_chr_excl_bd(-2)>0) &
        CALL Prnt0("     Non-numerical chromosome number: "//&
                        i2cp(n_chr_excl_bd(-1)))
      IF(n_chr_excl_bd(-1)>0) &
        CALL Prnt0("     Numerical chromosome number out of range: "//&
                        i2cp(n_chr_excl_bd(-1)))
      DO i=0,UBOUND(excluded_chr,1)
        IF(n_chr_excl_bd(i)==0) CYCLE
        CALL Prnt0("     Chromosome number "//TRIM(i2cp(i))//": "//&
                        i2cp(n_chr_excl_bd(i)))
      ENDDO
    ENDIF
    
    !! Announce how many still included
    n_incl_loci = CountTrue(InclLoc)
    CALL Prnt("Currently "//i2cp(n_incl_loci)//" loci are included.")

  ENDIF

  !! Return the number of excluded loci if counters present
  IF(PRESENT(n_chr_excl))   n_chr_excl = n_chr_excl1 
  IF(PRESENT(n_negbp_excl)) n_negbp_excl = n_negbp_excl1 
  IF(PRESENT(n_rs_excl))    n_rs_excl = n_rs_excl1 

  10 CONTINUE
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!         CHECK FOR PRESENCE OF MAPA RECORDS IN SUBMAPA         !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(.NOT.PRESENT(SUBMAPA))     GOTO 20
  IF(.NOT.ALLOCATED(SUBMAPA))   GOTO 20

  !! SUBMAPA size and max number of checks
  maxim = nloci

  !! Ititialize counters
  im = 0
  n_sub_incl1 = 0
  n_sub_excl1 = 0

  !! If SUBMAPARS missing or unallocated or sub_read_all is true, skip to 
  !! string comparison implementation of the exclusion
  numerical_compare = .TRUE.
  IF(sub_read_all)              numerical_compare = .FALSE.
  IF(.NOT.PRESENT(SUBMAPARS))   numerical_compare = .FALSE.
  IF(.NOT.ALLOCATED(MAPARS))    numerical_compare = .FALSE.
  IF(.NOT.ALLOCATED(SUBMAPARS)) numerical_compare = .FALSE.
  
  !! Allocate SUBrs and fill it with RS numbers from from SUBMAPA
  IF(sub_read_all) THEN
    nloci_sub = SIZE(SUBMAPA)
    CALL ResizeVar(SUBrs, nloci_sub)
    SUBrs = RESHAPE(SUBMAPA, (/nloci_sub/))
  ELSEIF(.NOT.numerical_compare) THEN
    nloci_sub = SIZE(SUBMAPA, 1)
    CALL ResizeVar(SUBrs, nloci_sub)
    IF(sub_rscol<1 .OR. sub_rscol>SIZE(SUBMAPA,2)) &
      CALL PrntE("sub_rscol="//i2cp(sub_rscol)//" is out of range,"//&
                      " which means that something is not properly"//&
                      " implemented.", Q=.TRUE.)
    SUBrs = RESHAPE(SUBMAPA(:,sub_rscol), (/nloci_sub/))
  ELSE
    nloci_sub = SIZE(SUBMAPARS)
  ENDIF
  
  !! Exclude/include based on submap file
  IF(n_incl_loci>0) THEN
  
    max_sub_incl = nloci
    
    !! Make sure there is not a conflict between 'subset mapping file being 
    !! paired up' and whether 'markers listed in it are included or excluded'
    IF(sub_make_pairs) THEN
      IF(IsOdd(nloci_sub)) &
        CALL PrntE("With 'pairing' the number of items in the subset"//&
                        " mapping file must be even!", Q=.TRUE.)
      sub_include = .TRUE.
    ENDIF

    !! Announce the next action
    IF(sub_include) THEN
      stext = "Excluding loci absent from the submap file ..."
    ELSE
      stext = "Excluding loci present in the submap file ..."
    ENDIF
    !CALL Prnt(stext, screen=.FALSE.)
  
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!
    !! DO THE STRING OR NUMERICAL COMPARISON BASED EXCLUSION !!
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!
    
    DO i=1,nloci
  
      !! Print count down progress
      CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext)
  
      !! Increase counter of checked loci
      im = im + 1
  
      !! Assume exclusion until the current loci is found in SUBMAPA
      Found = .FALSE.
  
      !! Do the comparison
      IF(.NOT.numerical_compare) THEN
        
        !! Remember current loci mapping info
        MAPArs_i = MAPA(i,map_rscol)
        
        !! Compare current MAPA record with all SUBMAPA records until it's found
        DO j=1,nloci_sub
          IF(MAPArs_i == SUBrs(j)) THEN
            Found = .TRUE.                                                
            nfound = nfound + 1
            EXIT
          ENDIF
        ENDDO
      
      ELSE
      
        !! Compare current MAPA record with all SUBMAPA records until it's found
        DO j=1,nloci_sub
          IF(MAPARS(i) == SUBMAPARS(j)) THEN
            Found = .TRUE.
            nfound = nfound + 1
            EXIT
          ENDIF
        ENDDO
      
      ENDIF
  
      !! If current pair already excluded, skip it (NOTE: the comparison above is 
      !! performed to determine how many pairs from SUBMAPA are (not) found in 
      !! MAPA, i.e. to get the correct value for counter 'nfound').
      IF(.NOT.InclLoc(i)) CYCLE
  
      !! i-th loci should be excluded/included if found 
      !IF((sub_include .AND. .NOT.Found) .OR. (.NOT.sub_include .AND. Found)) THEN
      IF(sub_include .NEQV. Found) THEN
        InclLoc(i) = .FALSE.
        n_sub_excl1 = n_sub_excl1 + 1
      ELSE
        n_sub_incl1 = n_sub_incl1 + 1
      ENDIF
  
      !! Perform a check for the maximum number of included loci
      !! If the maximum number of included records reached, stop checking
      IF(n_sub_incl1 >= max_sub_incl) THEN
        n_sub_excl1 = n_sub_excl1 + CountTrue(InclLoc(i+1:))
        InclLoc(i+1:) = .FALSE.
      ENDIF
      
    ENDDO
    
    !! Announce the end
    CALL PrintCD(text, im, maxim, pct, finish=.TRUE., stext=stext)
    CALL Prnt("Finished excluding loci based on the submap file.", screen=.FALSE.)

    !! Announce how many loci were found
    CALL Prnt("Total of "//i2cp(nfound)//" loci ("//&
                   TRIM(r2cp(hundred*round(DBLE(nfound)/nloci_sub,5)))//&
                   " % of "//i2cp(nloci_sub)//") from the file ["//&
                   TRIM(sub_file(1))//"] were found in the mapping file.")
  
    !! Announce how many still included
    n_incl_loci = CountTrue(InclLoc)
    CALL Prnt("Currently "//i2cp(n_incl_loci)//" loci are included.")

  ENDIF ! IF(n_incl_loci>0)
  
  !! If tests are only to be done for the pairs of loci in subset mapping file,
  !! create the pairs and sort them according their positions in MAPA to 
  !! accommodate easier looping over the pairs when testing
  IF(n_incl_loci>0 .AND. sub_make_pairs) THEN
    
    stext = "Extracting valid pairs ..."
    !CALL Prnt(stext, screen=.FALSE.)
  
    !! Allocate arrays to store the included pairs and the position indicator
    npairs = nloci_sub / 2
    CALL ResizeVar(Included_pairs, npairs, 2, 0)
    CALL ResizeVar(order, npairs)
    CALL ResizeVar(keep, npairs, 0)
    CALL ResizeVar(is_loc, SIZE(MAPA, 1))
    
    !! Loop over loci in SUBMAP and find their positions in MAPA's rs column
    k = 0 
    DO i=1,npairs

      !! Print progress
      CALL PrintCD(text, i-1, npairs, pct, finish=.TRUE., stext=stext, adv=.FALSE.)

      !! Locate the marker names in SUBrs within MAPA
      whr = 0
      DO j=1,2
        is_loc = 0
        WHERE(SUBrs(2*(i-1)+j) == MAPA(:,map_rscol)) is_loc = 1
        IF(SUM(is_loc)==0) EXIT
        whr(j) = MAXLOC(is_loc, 1)
      ENDDO
      
      !! Only if both loci found keep them
      IF(ALL(whr/=0)) THEN
        k = k+1
        keep(k) = i
        Included_pairs(i,1) = COUNT(InclLoc(1:whr(1)))
        Included_pairs(i,2) = COUNT(InclLoc(1:whr(2)))
      ENDIF
      
    ENDDO
    
    !! Announce the end
    CALL PrintCD(text, npairs, npairs, pct, finish=.TRUE., stext=stext)

    !! If none of the pairs are valid, quit
    IF(k==0) &
      CALL PrntE("There are no valid pairs of loci in the data!", Q=.TRUE.)
    
    !! Extract only the rows of Included_pairs with valid positions and drop
    !! the extra rows
    CALL ResizeVar(keep, k)
    Included_pairs(1:k,:) = Included_pairs(keep,:)
    CALL ResizeVar(Included_pairs, k, 2)
    npairs = k
    
    CALL Prnt("Total of "//i2cp(k)//" valid pairs found.")
    
    !! Sort the loci based on position in MAPA
    CALL Prnt("Sorting the pairs ...")
    
    !! Sort each pair
    DO i=1,npairs
      order1 = (/1, 2/)
      IncP1 = Included_pairs(i,:)                             
      CALL qsortInt(IncP1, 2, order1)
      Included_pairs(i,:) = IncP1
      !CALL QuickSort(Included_pairs(i,:))
    ENDDO
    
    !! Perform multi-column sort of all pairs
    CALL MultiQSort(Included_pairs)
    
    CALL Prnt("Pairs sorted.")
    
  ENDIF

  !! Return the number of excluded loci if counters present
  IF(PRESENT(n_sub_excl)) n_sub_excl = n_sub_excl1 
 
  20 CONTINUE
  
  RETURN

END SUBROUTINE CheckLociExcl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetRSCol(MAPA, rscol, chcol, ios)
  IMPLICIT NONE
  INTEGER, PARAMETER             :: maxncheck = 10
  CHARACTER(*), INTENT(IN)       :: MAPA(:,:)
  INTEGER, INTENT(INOUT)         :: rscol, chcol
  INTEGER, INTENT(OUT), OPTIONAL :: ios
  INTEGER                        :: i, j, col, maxcount, ncheck, &
                                    count(SIZE(MAPA,2))
  CHARACTER(mml)                 :: cvalue
  CHARACTER(mltl)                :: text
  
  !! Default status is 0 (no change)
  IF(PRESENT(ios)) ios = 0
  
  !! Initial and input value
  col = rscol
  
  ncheck = MIN(maxncheck, SIZE(MAPA,1))
  
  !! Look for the column with most values of form 'rs*'
  count = 0
  DO j=1,SIZE(MAPA,2)
    DO i=1,ncheck
      cvalue = MAPA(i,j)
      IF(upcasef(cvalue(1:2)) == "RS") count(j) = count(j) + 1
    ENDDO
  ENDDO

  !! Check if any RS numbers were found
  IF(SUM(count) == 0) &
    CALL PrntE("No RS numbers found in the submap file.", skip1=1, Q=.TRUE.)
                    
  !! Check if there are more than 1 columns with RS numbers 
  IF(SUM(MIN(count,1)) > 1) THEN
    IF(sub_make_pairs) THEN
      chcol = -9
      rscol = -9
    ELSE
      text = "Multiple columns of the submap file seem to contain marker"//&
             " names (RS numbers). There are 2 options how to process this file. If you"//&
             " want all of the values in the file to be treated as marker names use"//&
             " --submap-all-names. Alternatively, if you want to test markers in pairs,"//&
             " use --submap-pairs. Then the file is expected to contain marker"//&
             " information for two markers on each line which make up the pairs to test."
      CALL PrntE(text, skip1=1, Q=.TRUE.)
    ENDIF

  !! Select the column with largest number of RS numbers
  ELSE
    maxcount = 0
    DO i=1,SIZE(count)
      IF(count(i) > maxcount) THEN
        maxcount = count(i)
        rscol = i
      ENDIF
    ENDDO
  ENDIF
  
  !! Make sure rscol and chcol do not colide
  IF(rscol == chcol) chcol = -1
  
  !! Warn user if the identified rs column is different from the value on input
  IF(PRESENT(ios) .AND. col /= rscol) ios = 1
  
  RETURN

END SUBROUTINE GetRSCol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE FilterOutData(InclLoc, X, MAPA, MAF, MAFco, MAFca)
!! Filters out the loci that are not have false inclusion status in InclLoc
  IMPLICIT NONE
  LOGICAL, INTENT(INOUT), ALLOCATABLE      :: InclLoc(:)
  INTEGER(iks), INTENT(INOUT), ALLOCATABLE  :: X(:,:)
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: MAPA(:,:)
  REAL(dpp), INTENT(INOUT), ALLOCATABLE    :: MAF(:,:), MAFco(:,:), MAFca(:,:)
  INTEGER(iks), ALLOCATABLE                 :: X1(:,:)
  CHARACTER(LEN(MAPA)), ALLOCATABLE        :: MAPA1(:,:)
  REAL(dpp), ALLOCATABLE                   :: MAF1(:,:), MAF2(:,:), MAF3(:,:)
  INTEGER                                  :: ninc, j
  INTEGER(ikb)                             :: i, nloci, pct
  CHARACTER(mmtl)                          :: text, stext
  
  !! Count the included loci
  nloci = SIZE(InclLoc)
  ninc = CountTrue(InclLoc)
  
  !! If all loci are included, just leave
  IF(ninc == nloci) RETURN
  
  !! Otherwise temp arrays for storing included loci
  CALL ResizeVar(X1, SIZE(X,1), ninc)
  CALL ResizeVar(MAPA1, ninc, SIZE(MAPA,2))
  IF(compute_maf) THEN
    IF(.NOT.ALLOCATED(MAF)) CALL PrntE("MAF not allocated.", Q=.TRUE.)
    CALL ResizeVar(MAF1, ninc, SIZE(MAF,2))
    IF(out_maf_coca) THEN
      IF(.NOT.ALLOCATED(MAFco)) CALL PrntE("MAFco not allocated.", Q=.TRUE.)
      IF(.NOT.ALLOCATED(MAFca)) CALL PrntE("MAFca not allocated.", Q=.TRUE.)
      CALL ResizeVar(MAF2, ninc, SIZE(MAF,2))
      CALL ResizeVar(MAF3, ninc, SIZE(MAF,2))
    ENDIF
  ENDIF

  !! Perform inclusion check
  stext = "Removing data of excluded loci ..."
  j = 0
  DO i=1,nloci

    !! Announce progress
    CALL PrintCD(text, i-1, nloci, pct, finish=.TRUE., stext=stext, &
                 adv=.FALSE.)
    
    !! If the current pair is excluded, skip to the next one
    IF(.NOT.InclLoc(i)) CYCLE
    
    !! Otherwise save the current loci
    j = j + 1
    X1(:,j) = X(:,i)
    MAPA1(j,:) = MAPA(i,:)
    IF(compute_maf) THEN
      MAF1(j,:) = MAF(i,:)
      IF(out_maf_coca) THEN
        MAF2(j,:) = MAFco(i,:)
        MAF3(j,:) = MAFca(i,:)
      ENDIF
    ENDIF
        
  ENDDO  
  
  !! Shrink X
  CALL ResizeVar(X, SIZE(X,1), ninc)
  X = X1

  !! Shrink MAPA
  CALL ResizeVar(MAPA, ninc, SIZE(MAPA,2))
  MAPA = MAPA1

  !! Shrink MAF
  IF(compute_maf) THEN
    CALL ResizeVar(MAF, ninc, SIZE(MAF,2))
    MAF = MAF1
    IF(out_maf_coca) THEN
      CALL ResizeVar(MAFco, ninc, SIZE(MAFco,2))
      CALL ResizeVar(MAFca, ninc, SIZE(MAFca,2))
      MAFco = MAF2
      MAFca = MAF3
    ENDIF
  ENDIF

  !! Shrink InclLoc
  CALL ResizeVar(InclLoc, ninc)
  InclLoc = .TRUE.

  !! Announce progress
  CALL PrintCD(text, nloci, nloci, pct, finish=.TRUE., stext=stext)

  RETURN
  
END SUBROUTINE FilterOutData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!! ************************************************************************ !!  
!!                              DATA OUTPUT                                 !!
!! ************************************************************************ !!  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE TmpRecordSize(recsize)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: recsize
  
  IF(out_minimal) THEN
    INQUIRE(IOLENGTH=recsize) NAtmpmin
  ELSE
    INQUIRE(IOLENGTH=recsize) NAtmp
  ENDIF
  
  RETURN

END SUBROUTINE TmpRecordSize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SyncTmpFromTmpmin(tmp, tmpmin)
  IMPLICIT NONE
  TYPE(TMPTYPE), INTENT(OUT)   :: tmp
  TYPE(TMPTYPEmin), INTENT(IN) :: tmpmin
  
  tmp          = NAtmp
  tmp%PTAS4p   = tmpmin%PTAS4p
  tmp%PTAS1p   = tmpmin%PTAS1p
  tmp%PTDS4p   = tmpmin%PTDS4p
  tmp%PTDS1p   = tmpmin%PTDS1p
  tmp%PT4p     = tmpmin%PT4p
  tmp%PT1p     = tmpmin%PT1p
  tmp%PTPO4p   = tmpmin%PTPO4p
  tmp%PTPO1p   = tmpmin%PTPO1p
  tmp%AS4p     = tmpmin%AS4p
  tmp%AS1p     = tmpmin%AS1p
  tmp%DSp      = tmpmin%DSp
  tmp%CSp      = tmpmin%CSp
  tmp%errPTAS4 = tmpmin%errPTAS4 
  tmp%errPTAS1 = tmpmin%errPTAS1 
  tmp%errPTDS4 = tmpmin%errPTDS4 
  tmp%errPTDS1 = tmpmin%errPTDS1 
  tmp%errPT4   = tmpmin%errPT4 
  tmp%errPT1   = tmpmin%errPT1 
  tmp%errPTPO4 = tmpmin%errPTPO4 
  tmp%errPTPO1 = tmpmin%errPTPO1  
  tmp%errAS4   = tmpmin%errAS4  
  tmp%errAS1   = tmpmin%errAS1 
  tmp%errDS    = tmpmin%errDS  
  tmp%errCS    = tmpmin%errCS 
  tmp%PLevel   = tmpmin%PLevel
  tmp%Pos      = tmpmin%Pos  

  RETURN
  
END SUBROUTINE SyncTmpFromTmpmin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SyncTmpToTmpmin(tmpmin, tmp)
  IMPLICIT NONE
  TYPE(TMPTYPEmin), INTENT(OUT) :: tmpmin
  TYPE(TMPTYPE), INTENT(IN)     :: tmp
  
  tmpmin = NAtmpmin
  tmpmin%PTAS4p = tmp%PTAS4p
  tmpmin%PTAS1p = tmp%PTAS1p
  tmpmin%PTDS4p = tmp%PTDS4p
  tmpmin%PTDS1p = tmp%PTDS1p
  tmpmin%PT4p = tmp%PT4p
  tmpmin%PT1p = tmp%PT1p
  tmpmin%PTPO4p = tmp%PTPO4p
  tmpmin%PTPO1p = tmp%PTPO1p
  tmpmin%AS4p = tmp%AS4p
  tmpmin%AS1p = tmp%AS1p
  tmpmin%DSp = tmp%DSp
  tmpmin%CSp = tmp%CSp
  tmpmin%errPTAS4 = tmp%errPTAS4 
  tmpmin%errPTAS1 = tmp%errPTAS1 
  tmpmin%errPTDS4 = tmp%errPTDS4 
  tmpmin%errPTDS1 = tmp%errPTDS1 
  tmpmin%errPT4 = tmp%errPT4 
  tmpmin%errPT1 = tmp%errPT1 
  tmpmin%errPTPO4 = tmp%errPTPO4 
  tmpmin%errPTPO1 = tmp%errPTPO1  
  tmpmin%errAS4 = tmp%errAS4  
  tmpmin%errAS1 = tmp%errAS1 
  tmpmin%errDS = tmp%errDS  
  tmpmin%errCS = tmp%errCS 
  tmpmin%PLevel = tmp%PLevel
  tmpmin%Pos = tmp%Pos  

  RETURN
  
END SUBROUTINE SyncTmpToTmpmin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetNumRecords(nrec, fn)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(OUT)   :: nrec
  CHARACTER(*), INTENT(INOUT) :: fn(:)
  INTEGER                     :: i
  INTEGER(ikb)                :: fsize
  
  !! Make sure the size of 1 tmp record is known
  IF(rec_size_tmp<=0) CALL TmpRecordSize(rec_size_tmp)
  
  !! Do the output writing
  nrec = 0
  DO i=1,SIZE(fn)
  
    !! Get the file size
    fsize = MAX(0_ikb, GetFileSize(fn(i)))
    
    !! Make sure it is not out of range
    IF(fsize > HUGE(i)) &
      CALL PrntE("File '"//TRIM(fn(i))//"' is too big! (GetNumRecords)", Q=.TRUE.)
    
    !! Make sure it is compatible with the current record size
    IF(.NOT.IsMultiple(fsize, rec_size_tmp)) &
      CALL PrntE("Wrong size of file ["//TRIM(fn(i))//"]! File"//&
                   " size of "//TRIM(i2c(fsize))//" bytes is not a"//&
                   " multiple of the record size "//&
                   TRIM(i2c(rec_size_tmp))//" bytes.", Q=.TRUE.)
      
    !! Calculate the number of records 
    nrec = nrec + fsize / rec_size_tmp 
  
  ENDDO
  
  RETURN

END SUBROUTINE GetNumRecords  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadTempFile(fn, tmp, maxnrec, announce, ios, &
                        closefile)
  IMPLICIT NONE
  INTEGER, PARAMETER                         :: chunk = 10000
  CHARACTER(*), INTENT(INOUT)                :: fn
  TYPE(TMPTYPE), INTENT(OUT), ALLOCATABLE    :: tmp(:)
  INTEGER, INTENT(IN), OPTIONAL              :: maxnrec
  LOGICAL, INTENT(IN), OPTIONAL              :: announce, closefile
  INTEGER, INTENT(OUT), OPTIONAL             :: ios
  TYPE(TMPTYPE)                              :: tmp1
  TYPE(TMPTYPEmin)                           :: tmpmin1
  LOGICAL                                    :: announce1, closefile1, &
                                                file_opened, read_at_once
  INTEGER                                    :: nl, tmp_size, ios1, &
                                                maxnrec1, i
  INTEGER(ikb)                               :: fsize
  TYPE(TMPTYPEmin), ALLOCATABLE              :: tmpmin(:)
  
  maxnrec1 = HUGE(1)
  closefile1   = .TRUE.
  announce1    = .FALSE.
  !read_at_once = .FALSE. !! If this was false, sequential needs to be changed to direct
  read_at_once = .TRUE.
  
  IF(PRESENT(maxnrec)) THEN
    maxnrec1 = maxnrec
    read_at_once = .FALSE.
  ENDIF
  IF(PRESENT(announce))  announce1 = announce
  IF(PRESENT(closefile)) closefile1 = closefile
  
  !! Announce sorting
  IF(announce1) CALL Prnt("Reading temporary output file(s) ...")

  !! Open output file and WRITE the report
  CALL CheckFileExistence(fn, delim=delim)
  INQUIRE(FILE=fn, OPENED=file_opened)

  !! Read the temporary file
  IF(read_at_once) THEN

    !! Make sure the size of 1 tmp record is known
    IF(rec_size_tmp<=0) CALL TmpRecordSize(rec_size_tmp)
    
    !! Prepare tmp
    fsize = GetFileSize(fn)
    IF(fsize > HUGE(tmp_size)) &
      CALL PrntE("File '"//TRIM(fn)//"' is too big! (ReadTempFile)", Q=.TRUE.)
    tmp_size = INT(fsize)
    nl = tmp_size / rec_size_tmp 
    CALL ResizeVar(tmp, nl)
    
    !! Open file for reading
    IF(.NOT.file_opened) &
      CALL OpenFile(UNIT=utmp, FILE=fn, ACTION='R', STATUS="O", &
                    FORM='U', RECL=tmp_size, ACCESS='D')
  
    !! Read the entire temporary file at once
    IF(out_minimal) THEN
      
      CALL ResizeVar(tmpmin, nl)
      
      READ(utmp, REC=1, IOSTAT=ios1) tmpmin

      DO i=1,SIZE(tmp)
        CALL SyncTmpFromTmpmin(tmp(i), tmpmin(i))
      ENDDO

    ELSE
    
      READ(utmp, REC=1, IOSTAT=ios1) tmp
    
    ENDIF
    
  ELSE

    !! Prepare tmp
    tmp_size = MAX(1, MIN(maxnrec1, chunk)) 
    CALL ResizeVar(tmp, tmp_size)
  
    IF(.NOT.file_opened) &
      CALL OpenFile(UNIT=utmp, FILE=fn, ACTION='R', STATUS="O", &
                    FORM='U', ACCESS='S')
  
    !! Read from temporary file(s) and write into output file(s)
    nl = 0
    DO WHILE (nl<=maxnrec1)
  
      !! Read the next line from temporary file
      IF(out_minimal) THEN

        READ(utmp, IOSTAT=ios1) tmpmin1

        CALL SyncTmpFromTmpmin(tmp1, tmpmin1)

      ELSE

        READ(utmp, IOSTAT=ios1) tmp1

      ENDIF
      
      !! Exit if an error occurred
      IF(ios1 /= 0) EXIT
      
      !! Increase record counter
      nl = nl + 1
      
      !! Make sure tmp is big enough
      IF(tmp_size < nl) THEN
        tmp_size = tmp_size + chunk
        CALL ResizeVar(tmp, tmp_size)
      ENDIF 
      
      !! Assign the read value
      tmp(nl) = tmp1
      
    ENDDO
    !! Remove the extra items
    CALL ResizeVar(tmp, nl)

  ENDIF
  
  !! Close the file
  IF(closefile1 .OR. ios1/=0) CLOSE(utmp)
  
  !! Check if all ok
  IF(ALL(ios1/=err_eof)) GOTO 100
  
  !! Erase the announcement of sorting
  IF(announce1) THEN
    !! If no results read, just quit
    IF(nl==0) THEN
      CALL Prnt("Temporary output file(s) read but it contains data.")
      ios1 = -55
    ELSE
      CALL Prnt("Temporary output file(s) read.")
    ENDIF
  ENDIF

  100 CONTINUE
  
  !! Return ios if present
  IF(PRESENT(ios)) THEN
    ios = ios1
    RETURN
  ENDIF
  
  !! If error occurred report it (if no ios present)
  IF(ios1/=0) CALL Prnt("Temporary output file(s) could not be read.")
  
  RETURN

END SUBROUTINE ReadTempFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE WriteTempFile(tmp_file, tmp, tmp_nl, tmp_cor, save_all, &
                         announce, ios)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: tmp_file(:)
  TYPE(TMPTYPE), INTENT(IN)                :: tmp(:)
  INTEGER(ikb), INTENT(INOUT), OPTIONAL    :: tmp_nl, tmp_cor
  LOGICAL, INTENT(IN), OPTIONAL            :: save_all, announce
  INTEGER, INTENT(OUT), OPTIONAL           :: ios
  LOGICAL                                  :: announce1, save_all1, out_it
  INTEGER                                  :: ios1, Err(10)
  INTEGER(ikb)                             :: i, temp_size, tmp_nl1, &
                                              tmp_cor1, cur_nl, init_size
  REAL(dpp)                                :: PTAS4p, PTAS1p, PTDS4p, PTDS1p, &
                                              PTPO4p, PTPO1p, PT4p, PT1p, &
                                              AS4p, AS1p, DSp, CSp, min_p1, &
                                              min_p2, min_pC
  TYPE(TMPTYPEmin), ALLOCATABLE            :: tmpmin(:)
                                              
  !! Initial values
  save_all1 = .FALSE.
  tmp_nl1 = 0
  tmp_cor1 = 0
  ios1 = 0
  announce1 = .FALSE.

  !! Change variables according to optional input
  IF(PRESENT(save_all))   save_all1 = save_all
  IF(PRESENT(tmp_nl)) tmp_nl1 = MAX(0_ikb, tmp_nl)
  IF(PRESENT(tmp_cor))    tmp_cor1 = MAX(0_ikb, tmp_cor)
  IF(PRESENT(announce))   announce1 = announce
  
  !! Announce sorting
  IF(announce1) CALL Prnt("Saving temporary output file(s) ... ")

  !! Save the number of temp files
  tmp_nfiles = SIZE(tmp_file)
  
  !! Increase mean sample size counter
  DO i=1,SIZE(tmp)
    IF(tmp(i)%Empty) CYCLE
    IF(tmp(i)%ErrPTAS4<=0 .AND. tmp(i)%ErrPTAS4/=err_not_done) & 
      MeanSS(0) = MeanSS(0) + tmp(i)%PTASNco + tmp(i)%PTASNca
    IF(tmp(i)%ErrAS4<=0 .AND. tmp(i)%ErrAS4/=err_not_done) &
      MeanSS(1) = MeanSS(1) + tmp(i)%ASNco + tmp(i)%ASNca
    IF(tmp(i)%ErrAS1<=0 .AND. tmp(i)%ErrAS1/=err_not_done) &
      MeanSS(2) = MeanSS(2) + tmp(i)%ASNco + tmp(i)%ASNca
    IF(tmp(i)%ErrDS<=0 .AND. tmp(i)%ErrDS/=err_not_done) &
      MeanSS(3) = MeanSS(3) + tmp(i)%DSNco + tmp(i)%DSNca
    IF(tmp(i)%ErrCS<=0 .AND. tmp(i)%ErrCS/=err_not_done) &
      MeanSS(4) = MeanSS(4) + tmp(i)%CSNco + tmp(i)%CSNca
  ENDDO
  
  !! Save the content of tmp into tmpmin (only the relevant fields)
  IF(out_minimal) THEN
    CALL ResizeVar(tmpmin, SIZE(tmp))
    DO i=1,SIZE(tmp)
      CALL SyncTmpToTmpmin(tmpmin(i), tmp(i))
    ENDDO
  ENDIF
  
  !! Make sure the size of 1 tmp record is known
  IF(rec_size_tmp<=0) CALL TmpRecordSize(rec_size_tmp)
  
  !! WRITE THE RESULTS INTO THE TEMPORARY OUTPUT FILE
  IF(save_all1) THEN

    !! Open temp file and store the results all at once
    CALL OpenFile(UNIT=utmp, FILE=tmp_file(tmp_nfiles), ACTION='W', &
                  STATUS="U", FORM='U', POSITION="A", ACCESS='D', &
                  RECL=rec_size_tmp*SIZE(tmp))

    !! Write the temporary data
    IF(out_minimal) THEN
      WRITE(utmp, REC=tmp_nl1+1, IOSTAT=ios1) tmpmin
    ELSE
      WRITE(utmp, REC=tmp_nl1+1, IOSTAT=ios1) tmp
    ENDIF
    
    !! Increase temp record counter
    tmp_nl1 = tmp_nl1 + SIZE(tmp)

  ELSE

    !! Store initial file size and number of records
    init_size = MAX(0_ikb, GetFileSize(tmp_file(tmp_nfiles)))
    cur_nl = init_size / rec_size_tmp
    
    !! Check if the file size is a multiple of the record size
    IF(.NOT.IsMultiple(init_size, rec_size_tmp)) &
      CALL PrntE("File ["//TRIM(tmp_file(tmp_nfiles))//"] is corrupted!", Q=.TRUE.)    

    !! Open temp file and store the results one by one
    CALL OpenFile(UNIT=utmp, FILE=tmp_file(tmp_nfiles), ACTION='W', &
                  STATUS="U", FORM='U', POSITION="A", ACCESS='D', &
                  RECL=rec_size_tmp)
                  
    !! Print if all p-values meet the bounds or report_all is .TRUE.
    !! Note that negative error indicator is not really considered an error
    DO i=1,SIZE(tmp)
      
      !! Check if the current record is empty and if so, just skip it
      IF(tmp(i)%Empty) CYCLE

      !! Now perform report bound checks
      
      Err(1) = tmp(i)%ErrPTAS4
      Err(2) = tmp(i)%ErrPTAS1
      Err(3) = tmp(i)%ErrPTDS4
      Err(4) = tmp(i)%ErrPTDS1
      Err(5) = tmp(i)%ErrPTPO4
      Err(6) = tmp(i)%ErrPTPO1
      Err(7) = tmp(i)%ErrAS4
      Err(8) = tmp(i)%ErrAS1
      Err(9) = tmp(i)%ErrDS
      Err(10) = tmp(i)%ErrCS
      
      PTAS4p = tmp(i)%PTAS4p
      PTAS1p = tmp(i)%PTAS1p
      PTDS4p = tmp(i)%PTDS4p
      PTDS1p = tmp(i)%PTDS1p
      PTPO4p = tmp(i)%PTPO4p
      PTPO1p = tmp(i)%PTPO1p
      PT4p = tmp(i)%PT4p
      PT1p = tmp(i)%PT1p
      min_p1 = MIN(PTAS4p,PTAS1p,PTDS4p,PTDS1p,PTPO4p,PTPO1p,PT4p,PT1p)
      
      AS4p = tmp(i)%AS4p
      AS1p = tmp(i)%AS1p
      DSp = tmp(i)%DSp
      CSp = tmp(i)%CSp
      min_p2 = MIN(AS4p,AS1p,DSp,CSp)
      
      !! Increase the approximate counter of tests in the second step
      IF(min_p1 < tmp(i)%PLevel(1)) tmp_cor1 = tmp_cor1 + 1
      min_pC = min_p2 * MAX(1_ikb, tmp_cor1)
      
      !! Decide whether to save the current record
      out_it = .FALSE.
      IF(report_all .AND. ANY(Err<=0))     out_it = .TRUE.
      IF(report_errs .AND. ANY(Err>0))     out_it = .TRUE.
      IF(min_p1 < olim1 .AND. olim2 > one) out_it = .TRUE.
      IF(min_p2 < olim2 .AND. olim3 > one) out_it = .TRUE.
      IF(min_pC < olim3)                   out_it = .TRUE.
                                                               
      !! Don't print if nothing was computed
      IF(MIN(min_p1,min_p2)==NAp .AND. ALL(Err==0)) out_it = .FALSE.
  
      !! If none of the above conditions are met, don't write this record
      IF(.NOT.out_it) CYCLE
    
      !! Increase counter of currently outputted tests
      tmp_nl1 = tmp_nl1 + 1
      cur_nl = cur_nl + 1 
  
      !! Calculate the temp file size after writing the current record
      temp_size = cur_nl*rec_size_tmp
      
      !! Check for filesize whether it exceeds maximum allowed size and if so,  
      !! run a routine that closes current file and reopens a new one with _N 
      !! appended to the file name (where N is the number of output files with  
      !! name matching the current file's pattern)
      IF(temp_size>=max_out_size) THEN

        CALL ExcessFileSize(tmp_file, utmp, temp_size, max_out_size, &
                            cur_nl, tmp_ext, check_size=.FALSE., form='U', &
                            access='D', rec_size=rec_size_tmp)

        tmp_nfiles = SIZE(tmp_file) 
        cur_nl = 1

      ENDIF
  
      !! Save the current record into the temp file
      !WRITE(utmp, REC=cur_nl, IOSTAT=ios1) tmp(i)
      IF(out_minimal) THEN
        WRITE(utmp, REC=cur_nl, IOSTAT=ios1) tmpmin(i)
      ELSE
        WRITE(utmp, REC=cur_nl, IOSTAT=ios1) tmp(i)
      ENDIF
      
      !! Check for an error
      IF(ios1/=0) EXIT
      
    ENDDO
  ENDIF
  
  !! Remove end-of-file status
  IF(ANY(ios1==err_eof)) ios1 = 0

  !! Close the temp file and delete it if nothing written into it (yet)
  IF(cur_nl<=0) THEN
    CLOSE(utmp, STATUS='DELETE')
  ELSE
    CLOSE(utmp)    
  ENDIF

  !! Return modified counters
  IF(PRESENT(tmp_nl)) tmp_nl = tmp_nl1
  IF(PRESENT(tmp_cor)) tmp_cor = tmp_cor1

  !! Return ios if present
  IF(PRESENT(ios)) THEN
    ios = ios1
    RETURN
  ENDIF
  
  !! If error occurred report it (if no ios present)
  IF(ios1/=0) THEN
    CALL Prnt("Temporary output file(s) might not be saved (IOSTAT "//&
                   i2cp(ios1)//").", skip1=1, skip2=1)
  ELSEIF(announce1) THEN  
    CALL Prnt("Temporary output file(s) saved.")
  ENDIF

  RETURN

END SUBROUTINE WriteTempFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SaveTempFile(tmp, tmp_nl, tmp_cor, last)

!! Writes the computed p-values into a specified files
!! Quality checks of the p-values are performed
  IMPLICIT NONE
  TYPE(TMPTYPE), INTENT(IN)   :: tmp(:)
  INTEGER(ikb), INTENT(INOUT) :: tmp_nl, tmp_cor
  LOGICAL, INTENT(IN)         :: last
  INTEGER(ikb)                :: k
  CHARACTER(mstl)             :: text
  INTEGER                     :: ierr, err_low, err_neg, j, chr1, chr2
  
  !! Count OK and problematic tests
  DO k=1,SIZE(tmp)

    !! Count error-free tests
    IF(tmp(k)%AS4p /= NAp) NumS2AS4OK = NumS2AS4OK + 1 
    IF(tmp(k)%AS1p /= NAp) NumS2AS1OK = NumS2AS1OK + 1 
    IF(tmp(k)%DSp /= NAp)  NumS2DSOK = NumS2DSOK + 1 
    IF(tmp(k)%CSp /= NAp)  NumS2CSOK = NumS2CSOK + 1 

    !! Check the problem with linear dependence in x and y
    IF(tmp(k)%errAS4==err_lin_dep) NumLinDep(1) = NumLinDep(1) + 1
    IF(tmp(k)%errAS1==err_lin_dep) NumLinDep(2) = NumLinDep(2) + 1
    IF(tmp(k)%errDS==err_lin_dep) NumLinDep(3) = NumLinDep(3) + 1
    IF(tmp(k)%errCS==err_lin_dep) NumLinDep(4) = NumLinDep(4) + 1
    !IF(ANY((/tmp(k)%errAS4,tmp(k)%errAS1,tmp(k)%errDS,tmp(k)%errCS/)==err_lin_dep)) &
    !  NumLinDep = NumLinDep + 1
      
    !! Count tests that include a sex chromosome
    IF(ALLOCATED(MAPA)) THEN
      chr1 = c2i(MAPA(tmp(k)%Pos(1),map_chcol))
      chr2 = c2i(MAPA(tmp(k)%Pos(2),map_chcol))
      IF(ANY(chr1 == sex_chr) .OR. ANY(chr2 == sex_chr)) &
        NumS2_sex = NumS2_sex + 1
    ENDIF 

    !! Check the problem with low or negative variances
    ierr = 0
    err_low = 0
    err_neg = 0
    DO j=1,4
    
      IF(j==1) ierr = tmp(k)%errAS4
      IF(j==2) ierr = tmp(k)%errAS1
      IF(j==3) ierr = tmp(k)%errDS
      IF(j==4) ierr = tmp(k)%errCS
    
      IF(j==1) err_low = err_low_var_AS
      IF(j==2) err_low = err_low_var_AS
      IF(j==3) err_low = err_low_var_DS
      IF(j==4) err_low = err_low_var_CS
      
      IF(j==1) err_neg = err_neg_var_AS
      IF(j==2) err_neg = err_neg_var_AS
      IF(j==3) err_neg = err_neg_var_DS
      IF(j==4) err_neg = err_neg_var_CS

      !! Check error code
      IF(ierr > 0) THEN
        IF(j==1) NumS2AS4NotOK = NumS2AS4NotOK + 1
        IF(j==2) NumS2AS1NotOK = NumS2AS1NotOK + 1
        IF(j==3) NumS2DSNotOK = NumS2DSNotOK + 1
        IF(j==4) NumS2CSNotOK = NumS2CSNotOK + 1 
        IF(ANY(ierr==(/err_low, err_neg/))) &
          NumVarErr(j) = NumVarErr(j) + 1
        IF(ierr==err_low) &
          NumLowVar(j) = NumLowVar(j) + 1 
        IF(ierr==err_neg) &
          NumNegVar(j) = NumNegVar(j) + 1 
      ENDIF
    ENDDO 

  ENDDO
  
  !! Count total number of tests
  NumS2AS4 = NumS2AS4OK + NumS2AS4NotOK
  NumS2AS1 = NumS2AS1OK + NumS2AS1NotOK
  NumS2DS = NumS2DSOK + NumS2DSNotOK
  NumS2CS = NumS2CSOK + NumS2CSNotOK

  !! Make sure the size of 1 tmp record is known
  IF(rec_size_tmp<=0) CALL TmpRecordSize(rec_size_tmp)

  !! Announce current action
  IF(last) THEN
    CALL Prnt("Saving results ... ", ADVANCE='NO')
  ELSEIF(countdown) THEN
    text = " (saving results ...)"
    WRITE(usto,'(A)', ADVANCE='NO') TRIM(text)
  ENDIF

  !! Save the results
  CALL WriteTempFile(tmp_file, tmp, tmp_nl, tmp_cor)
                        
  !! Announce end of current action
  IF(last) THEN
    CALL Prnt0("Finished.", log=.FALSE.)
    CALL Prnt("Results saved.", screen=.FALSE.)
  ELSEIF(countdown) THEN
    CALL EraseText(LEN_TRIM(text))
  ENDIF
  
  !! Sort the temporary file only if no more than 1 tmp files produced
  IF(out_sort .AND. last .AND. (SIZE(tmp_file)==1 .OR. out_sort_multiple_files)) THEN
    CALL SortFile(tmp_file, temp=.TRUE., secondary=out_sort_second)
    tmp_sorted = .TRUE.
    out_sort = .FALSE.
  ENDIF

  RETURN
  
END SUBROUTINE SaveTempFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE WriteOutputFile(tmp_file, tmp_nl, append)
!! Writes the computed p-values into a specified files
!! Quality checks of the p-values are performed
  IMPLICIT NONE
  INTEGER, PARAMETER                    :: chunk = 1000
  CHARACTER(*), INTENT(INOUT)           :: tmp_file(:)
  INTEGER(ikb), INTENT(INOUT), OPTIONAL :: tmp_nl
  LOGICAL, INTENT(IN), OPTIONAL         :: append
  TYPE(TMPTYPE)                         :: tmps(1:chunk), tmp
  TYPE(TMPTYPEmin)                      :: tmpmins(1:chunk)
  REAL(dpp)                             :: PTAS4p, PTAS1p, PTDS4p, PTDS1p, &
                                           PTPO4p, PTPO1p, PT4p, PT1p, AS4p, &
                                           AS1p, DSp, CSp, pCor(7), min_p1, &
                                           min_p2, min_pC, startt, stopt
  LOGICAL                               :: out_it, is_err, write_errfree !, signif                                           
  CHARACTER(mfl), ALLOCATABLE           :: fn(:)
  CHARACTER                             :: hc, POS
  CHARACTER(mmtl)                       :: text
  CHARACTER(mstl)                       :: text1
  CHARACTER(mttl)                       :: fv, fp
  CHARACTER(mtsl)                       :: startt_full, stopt_full, runtime
  INTEGER                               :: ifile, ios, i, &
                                           j, record_length, next_chunk
  INTEGER(ikb)                          :: pct, ij, tmp_irec, &
                                           tmp_nrecs, nl_now, &
                                           nl_all, nl_val, nl_err, &
                                           nerr, out_size, &
                                           tmp_nl1, tmp_size  
  LOGICAL                               :: append1 
    
  !! Initialize counters of lines and files
  nl_all = 0
  nl_now = 0
  nl_val = 0
  nl_err = 0
  tmp_nl1 = -1
  out_size = 0
  out_cnams = ""
  out_nfiles = 1
  write_errfree = .TRUE.

  !! Set the format to print variables in output file
  fv = "F0."//i2cp(ndig_stat)
  fp = "ES"//i2cp(6+ndig_pval)//"."//i2cp(ndig_pval)
  
  IF(PRESENT(tmp_nl)) tmp_nl1 = tmp_nl
  
  !! Make sure the size of 1 tmp record is known
  IF(rec_size_tmp<=0) CALL TmpRecordSize(rec_size_tmp)

  !! Determine whether file is overwritten or not
  POS = 'R'
  append1 = .FALSE.
  IF(PRESENT(append)) append1 = append
  IF(append1) POS = 'A'

  !! Make sure output file is defined
  IF(LEN_TRIM(out_file)==0) THEN
    out_file = TRIM(tmp_file(1))//out_ext
    keep_temp = .TRUE.
  ENDIF

  !! Save out_file locally so that it will remain unchanged for later use
  CALL ResizeVar(fn, out_nfiles, out_file)
  
  !! Set the out_cycle
  out_cycle = MAX( INT(HUGE(1) * (one - mfsf) / rec_size_out), out_cycle )
      !! Determines after how many output lines file size of output file will be 
      !! checked. It should not be too large to make sure the file size does not
      !! exceed HUGE(1)=2^31. it's determined by max_out_size and the maximum
      !! potential size of output record. Currently the value of out_cycle is 
      !! about 400,000, which comes from the rec_size_out of 520 (max output 
      !! record size is currently estimated to be 400 bytes, thus the value 
      !! of 520 seems enough) because 2,100,000,000 * (1-0.9) / 520 ~= 400,000, 
      !! where (2,100,000,000 ~= HUGE(1))
      
  CALL Prnt("Output file size check every "//TRIM(i2c(out_cycle))//" records.")
      
  !! Inform about interaction effect computation
  IF(out_epi_effect) &
    CALL Prnt("Interaction effects will be computed while output file is produced."// &
              " If there are many results to output this process can be very slow!")

  !! Open the latest output file for writing
  CALL OpenFile(UNIT=uout, FILE=fn(out_nfiles), ACTION='W', STATUS="U",&
              FORM='F', POSITION=POS, ACCESS='S', ask=.TRUE., &
              overwrite=overwrite_files, chck_exist=check_fexist)
              
  !! Get start time
  CALL GetCurrentTime(startt_full, startt)

  !! Get the total number of records in all temp files
  IF(tmp_nl1<0) CALL GetNumRecords(tmp_nl1, tmp_file)

  !! If temp files are empty, then warn and return
  IF(tmp_nl1<=0) THEN
    IF(PRESENT(tmp_nl)) tmp_nl = tmp_nl1
    CALL Prnt("There are no testing results to output.")
    CALL LH("No test results to output (see log file).", out_anyway=.TRUE.)
    CLOSE(uout)
    RETURN
  ELSE
    CALL Prnt("Records in temporary result files: "//TRIM(i2cp(tmp_nl1)))
  ENDIF
  
  333 CONTINUE
  
  !! Reset counters
  ij = 0
  nerr = 0
  pct = -1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!          WRITE THE NON-ERRONEOUS RESULTS INTO THE OUTPUT FILE          !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Otherwise announce the end of writing to temp file and start of writing
  !! to output file
  IF(write_errfree) THEN
    text = "error-free"
  ELSE
    text = "erroneous"
  ENDIF
  
  text = "Writing "//TRIM(text)//" results into output file(s) ..."
  CALL Prnt(TRIM(text)//" ", ADVANCE='NO')
  
  !! Do the output writing
  text = ""
  DO ifile=1,SIZE(tmp_file)
  
    !! Make sure the file exists
    CALL CheckFileExistence(tmp_file(ifile), delim=delim)

    !! Store initial file size and number of records
    tmp_size = MAX(0_ikb, GetFileSize(tmp_file(ifile)))
    tmp_nrecs = tmp_size / rec_size_tmp
    
    !! Check if the file size is a multiple of the record size
    IF(.NOT.IsMultiple(tmp_size, rec_size_tmp)) &
      CALL PrntE("File ["//TRIM(tmp_file(tmp_nfiles))//"] is corrupted!", Q=.TRUE.)    

    !! Open output file and WRITE the report
    CALL OpenFile(UNIT=utmp, FILE=tmp_file(ifile), ACTION='R', STATUS="O", FORM='U', &
                  POSITION="R", ACCESS='D', RECL=rec_size_tmp)

    tmp_irec = 0
    !! Read from temporary file(s) and write into output file(s)
    DO

      !! Print progress
      CALL PrintCD(text, ij, tmp_nl1, pct, finish=.TRUE., adv=.FALSE.)
      
      !! Quit if all has been read
      IF(tmp_irec==tmp_nrecs) EXIT
      
      !! Calculate how many records to read next
      next_chunk = INT(MIN(INT(chunk, ikb), tmp_nrecs - tmp_irec))
      
      !! Increase the progress counter of read lines
      ij = ij + next_chunk
      
      !! Read the next group of records
      IF(out_minimal) THEN
        DO j=1,next_chunk
          READ(utmp, REC=tmp_irec+j, IOSTAT=ios) tmpmins(j)
        ENDDO
      ELSE
        DO j=1,next_chunk
          READ(utmp, REC=tmp_irec+j, IOSTAT=ios) tmps(j)
        ENDDO
      ENDIF
       
      !! Move the counter of read lines to the end of the last group of records
      tmp_irec = tmp_irec + next_chunk
      
      !! If error or end of file occurred then exit
      IF(ios /= 0) EXIT
      
      !! Loop over the read records
      DO j=1,next_chunk
      
        IF(out_minimal) THEN
          CALL SyncTmpFromTmpmin(tmp, tmpmins(j))
        ELSE
          tmp = tmps(j)
        ENDIF

        !! If the current record is not an error, skip to the next one
        is_err = .FALSE.
        IF(doAS4 .AND. ANY((/tmp%ErrPTAS4, tmp%ErrAS4/)>0)) is_err = .TRUE.
        IF(doAS1 .AND. ANY((/tmp%ErrPTAS1, tmp%ErrAS1/)>0)) is_err = .TRUE.
        IF(doDS4 .AND. ANY((/tmp%ErrPTDS4, tmp%ErrDS/)>0))  is_err = .TRUE.
        IF(doDS1 .AND. ANY((/tmp%ErrPTDS1, tmp%ErrDS/)>0))  is_err = .TRUE.
        IF(doPO4 .AND. ANY((/tmp%ErrPTPO4, tmp%ErrCS/)>0))  is_err = .TRUE.
        IF(doPO1 .AND. ANY((/tmp%ErrPTPO1, tmp%ErrCS/)>0))  is_err = .TRUE.
        IF(doPT4 .AND. ANY((/tmp%ErrPT4,   tmp%ErrAS4/)>0)) is_err = .TRUE.
        IF(doPT1 .AND. ANY((/tmp%ErrPT1,   tmp%ErrAS1/)>0)) is_err = .TRUE.
        IF(doCS .AND. tmp%ErrCS>0)                          is_err = .TRUE.

        !! Decide whether it should be printed or not
        IF(write_errfree) THEN
  
          !! If writing error-free results and the current record has an error, skip it
          IF(is_err) THEN
            nerr = nerr + 1
            CYCLE
          ENDIF
          
          !! If AS was not standardized then standardize it with mean of varAS
          IF(.NOT.out_minimal .AND. var_poststand .AND. .NOT.tmp%ASstd &
          .AND. MAX(var_bound, varAS_sum)>zero) THEN
            IF(poststand_by_varbound .AND. var_bound>zero) THEN
              tmp%AS4 = tmp%AS4 / SQRT(var_bound)
            ELSEIF(varAS_sum>0) THEN
              tmp%AS4 = tmp%AS4 * SQRT(varAS_sum_n / varAS_sum)
            ELSE
              GOTO 111
            ENDIF
            NumS2AS4OK = NumS2AS4OK + 1
            NumS2AS4PostStand = NumS2AS4PostStand + 1
            tmp%AS4p = norm_pval(tmp%AS4)
            IF(tmp%AS4p < level_S2_nom) NumEpi(1) = NumEpi(1) + 1
            111 CONTINUE
          ENDIF
    
        ELSE
        
          !! If writing error results and the current record has no error, skip it
          IF(.NOT.is_err .OR. (is_err .AND. .NOT.report_errs)) CYCLE
        
        ENDIF
        
        !! Save p-values, error and significance information locally
        PTAS4p = tmp%PTAS4p
        PTAS1p = tmp%PTAS1p
        PTDS4p = tmp%PTDS4p
        PTDS1p = tmp%PTDS1p
        PTPO4p = tmp%PTPO4p
        PTPO1p = tmp%PTPO1p
        PT4p = tmp%PT4p
        PT1p = tmp%PT1p
        AS4p = tmp%AS4p
        AS1p = tmp%AS1p
        DSp = tmp%DSp
        CSp = tmp%CSp
  
        !! Get multiple testing corrected p-values
        pCor = NAp
        IF(AS4p /= NAp) pCor(1) = MAX(AS4p, MIN(one, MTC(1)*AS4p))
        IF(AS1p /= NAp) pCor(2) = MAX(AS1p, MIN(one, MTC(2)*AS1p))
        IF(DSp /= NAp)  pCor(3) = MAX(DSp,  MIN(one, MTC(3)*DSp))
        IF(DSp /= NAp)  pCor(4) = MAX(DSp,  MIN(one, MTC(4)*DSp))
        IF(CSp /= NAp)  pCor(5) = MAX(CSp,  MIN(one, MTC(5)*CSp))
        IF(CSp /= NAp)  pCor(6) = MAX(CSp,  MIN(one, MTC(6)*CSp))
        IF(CSp /= NAp)  pCor(7) = MAX(CSp,  MIN(one, MTC(7)*CSp))
        
        !! Count significant multiple testing corrected statistics
        DO i=1,7
          IF(pCor(i) < level_S2) NumEpiCorr(i) = NumEpiCorr(i) + 1
        ENDDO

        !! Print if all p-values meet the bounds or report_all is true
        !! Note that Errs(i)<0 is not really considered an error
        min_p1 = MIN(PTAS4p,PTAS1p,PTDS4p,PTDS1p,PTPO4p,PTPO1p)
        min_p2 = MIN(AS4p,AS1p,DSp,CSp)
        min_pC = MINVAL(pCor)
  
        !! Start by assuming the current record will not be printed      
        out_it = .FALSE.
        
        !! Change out_it if reporting conditions are met
        IF(report_all) out_it = .TRUE.
        IF(MIN(olim2,olim3)>one .AND. olim1<=one .AND. min_p1<olim1) out_it = .TRUE.
        IF(MIN(olim1,olim3)>one .AND. olim2<=one .AND. min_p2<olim2) out_it = .TRUE.
        IF(MIN(olim1,olim2)>one .AND. olim3<=one .AND. min_pC<olim3) out_it = .TRUE.
           
        !! Skip to next cycle if the line should not be printed
        IF(.NOT.out_it) CYCLE
        
        !! Increase counters of lines since last file size check
        nl_now = nl_now + 1
        nl_all = nl_all + 1
        IF(write_errfree) THEN
          nl_val = nl_val + 1 
        ELSE
          nl_err = nl_err + 1 
        ENDIF 
  
        !! Check for exceedingly large output file
        IF(out_size>=max_out_size) THEN

          IF(out_zipped) CALL ZipFile(fn, announce=.TRUE.)

          CALL ExcessFileSize(fn, uout, out_size, max_out_size, nl_now, out_ext, check_size=.FALSE., &
                              form='F', access='S', rec_size=rec_size_tmp)

          out_nfiles = SIZE(fn)
          out_size = 0

        ENDIF
                            
        !! Write header if this is the first result of the current file
        IF(out_header .AND. nl_now==1 .AND. (.NOT.append1 .OR. out_nfiles>1)) THEN
          hc = ""
          IF(out_header_comment) hc = comment
          CALL PrintHeader(hc, starttime_full, out_nfiles)
        ENDIF
  
        !! Write the record
        !!print *,'ij=',ij
        CALL WriteOutputRecord(nl_all, tmp, pCor, MTC, fv, fp, record_length)
        
        out_size = out_size + record_length 

      ENDDO ! DO j=1,next_chunk
 
    ENDDO
    
    CALL PrintCD(text, tmp_nl1, tmp_nl1, pct, finish=.TRUE.)
      
    !! Close the current temp file
    IF(write_errfree) THEN
    
      !! For error-free results, if no errors present or they should not be 
      !! reported, close and delete (unless keep_temp is true)
      IF(nerr == 0 .OR. .NOT.report_errs) THEN
        IF(keep_temp) THEN
          CLOSE(utmp)
        ELSE
          CLOSE(utmp, STATUS='DELETE')
        ENDIF
      !! Otherwise just close
      ELSE
        CLOSE(utmp)
      ENDIF
    
    ELSE
      
      !! For erroneous results, close and delete unless keep_temp is true
      IF(keep_temp) THEN
        CLOSE(utmp)
      ELSE
        CLOSE(utmp, STATUS='DELETE')
      ENDIF
    
    ENDIF
    
    IF(.NOT.countdown) CALL Prnt0("", log=.FALSE.)

  ENDDO
  

  !! Check errors are still to be printed and if so, go to the top
  IF(nerr>0 .AND. report_errs .AND. write_errfree) THEN
    write_errfree = .FALSE.
    GOTO 333
  ENDIF 
  
  !! Divide mean sample size counter by test counts
  MeanSS = MeanSS / (/MAX(1_ikb, NumS1OK), &
                      MAX(1_ikb, NumS2AS4OK), &
                      MAX(1_ikb, NumS2AS1OK), &
                      MAX(1_ikb, NumS2DSOK), &
                      MAX(1_ikb, NumS2CSOK)/)
  IF(NumS1OK == 0)    MeanSS(0) = 0
  IF(NumS2AS4OK == 0) MeanSS(1) = 0
  IF(NumS2AS1OK == 0) MeanSS(2) = 0
  IF(NumS2DSOK == 0)  MeanSS(3) = 0
  IF(NumS2CSOK == 0)  MeanSS(4) = 0
  
  !! Get the time
  CALL GetCurrentTime(stopt_full, stopt)
  runtime = GetTimeDifference(startt, stopt)
  
  !! Announce how many output files were produced
  IF(nl_all>0) THEN

    IF(out_zipped) CALL ZipFile(fn, announce=.TRUE.)

    out_available = .TRUE.
    CALL Prnt(i2c(out_nfiles)//" output file(s) produced with "//i2c(nl_all)//&
              " results ("//i2c(nl_val)//" error-free, "//i2c(nl_err)//&
              " with errors, writing time "//TRIM(runtime)//"):")
    DO ifile=1,SIZE(fn)
      CALL Prnt0("["//TRIM(fn(ifile))//"]", lead=3)
    ENDDO
    
  ELSE
  
    CALL LH("None of the test results satisfied the p-value reporting bounds!")
    CALL Prnt("No results were outputted.")
  
  ENDIF
  
  CLOSE(uout)
  
  !! Announce the number of significant MT corrected statistics
  text = ""
  DO j=1,7
    IF(j==1 .AND. .NOT.doAS4) CYCLE
    IF(j==2 .AND. .NOT.doAS1) CYCLE
    IF(j==3 .AND. .NOT.doDS4) CYCLE
    IF(j==4 .AND. .NOT.doDS1) CYCLE
    IF(j==5 .AND. .NOT.doPO4) CYCLE
    IF(j==6 .AND. .NOT.doPO1) CYCLE
    IF(j==7 .AND. .NOT.doCS) CYCLE
    IF(j==1) text1 = "AS4"
    IF(j==2) text1 = "AS1"
    IF(j==3) text1 = "DS4"
    IF(j==4) text1 = "DS1"
    IF(j==5) text1 = "PO4"
    IF(j==6) text1 = "PO1"
    IF(j==7) text1 = "CS"
    IF(LEN_TRIM(text)>0) text = TRIM(text)//","
    text = TRIM(text)//" "//TRIM(text1)//": "//i2c(NumEpiCorr(j))
  ENDDO  
  CALL Prnt("Results summary (significant corrected p-values):"//text)
    
  IF(out_sort) sort_file = out_file
  IF(PRESENT(tmp_nl)) tmp_nl = tmp_nl1
  
  RETURN
    
END SUBROUTINE WriteOutputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadOutputFile(fn, header, body, ncols, announce, warn_empty, ios)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)            :: fn
  CHARACTER(*), INTENT(OUT), ALLOCATABLE :: header(:,:)
  CHARACTER(*), INTENT(OUT), ALLOCATABLE :: body(:,:)
  INTEGER, INTENT(INOUT)                 :: ncols
  LOGICAL, INTENT(IN), OPTIONAL          :: announce, warn_empty
  INTEGER, INTENT(OUT), OPTIONAL         :: ios
  LOGICAL                                :: announce1
  INTEGER                                :: read_nlines, file_nlines, hncols, &
                                            ios1
  CHARACTER                              :: sep(1)
  
  announce1 = .FALSE.
  IF(PRESENT(announce)) announce1 = announce
  
  !! Announce sorting
  IF(announce1) CALL Prnt("Reading output file(s) ... ", ADVANCE='NO')

  IF(PRESENT(ios)) ios = 0

  !! Read the output file header
  hncols = 1
  read_nlines = 0
  file_nlines = 0
  ios1 = 0
  sep = out_cs
  CALL ReadFile(fn, uout, header, read_nlines, file_nlines, hncols, sep, &
                cmt=comment, skip_cmted=.FALSE., max_nl=100, &
                only_cmted=.TRUE., warn_empty=warn_empty, ios=ios1)
                
  !! Check if all ok
  IF(ios1/=0) GOTO 100
  
  !! If column headers not commented read one more line              
  IF(.NOT.out_header_comment) &
    CALL ReadFile(fn, uout, header, read_nlines, file_nlines, hncols, &
                  sep, attach=.TRUE., max_nl=1, warn_empty=warn_empty, ios=ios1)

  !! Check if all ok
  IF(ios1/=0) GOTO 100

  !! Read the output file result data
  read_nlines = 0
  file_nlines = 0
  CALL ReadFile(fn, uout, body, read_nlines, file_nlines, ncols, sep, &
                skip_cmted=.FALSE., nskip=SIZE(header,1), &
                warn_empty=warn_empty, progress=announce1, ios=ios1)
  
  !! Check if all ok
  IF(ios1/=0) GOTO 100

  !! Erase the announcement of sorting
  IF(announce1) THEN
    !! If no results read, just quit
    IF(read_nlines==0) THEN
      CALL Prnt("Output file(s) read but it contains no results.")
      ios1 = -55
    ELSE
      CALL Prnt("Output file(s) read.")
    ENDIF
  ENDIF

  100 CONTINUE
  
  !! Return ios if present
  IF(PRESENT(ios)) THEN
    ios = ios1
    RETURN
  ENDIF
  
  !! If error occurred report it (if no ios present)
  IF(ios1/=0) CALL Prnt("Output file(s) could not be read.")
  
  RETURN

END SUBROUTINE ReadOutputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SaveOutputFile(fn, header, body, announce, ios)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: fn
  CHARACTER(*), INTENT(IN)       :: header(:,:)
  CHARACTER(*), INTENT(IN)       :: body(:,:)
  LOGICAL, INTENT(IN), OPTIONAL  :: announce
  INTEGER, INTENT(OUT), OPTIONAL :: ios
  LOGICAL                        :: announce1
  INTEGER                        :: ios1
  
  !! Set announcement
  announce1 = .FALSE.
  IF(PRESENT(announce)) announce1 = announce
  
  !! Announce sorting
  IF(announce1) CALL Prnt("Saving output file(s) ... ")

  !! Assume everythink will go ok
  IF(PRESENT(ios)) ios = 0

  !! Write the header
  CALL WriteFile(fn, uout, header, sep=out_cs, ios=ios)

  !! Check if all ok
  IF(ios/=0) GOTO 100

  !! Now write the body
  CALL WriteFile(fn, uout, body, sep=out_cs, attach=.TRUE., ios=ios1)

  !! Check if all ok
  IF(ios/=0) GOTO 100

  !! Erase the announcement of sorting
  IF(announce1) CALL Prnt("Output file(s) saved.")

  RETURN

  100 CONTINUE
  
  !! Return ios if present
  IF(PRESENT(ios)) THEN
    ios = ios1
    RETURN
  ENDIF
  
  !! If error occurred report it (if no ios present)
  IF(ios1/=0) CALL Prnt("Output file(s) could not be read.")
  
  RETURN

END SUBROUTINE SaveOutputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ZipFile(fn, announce)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: fn(:)
  LOGICAL, INTENT(IN), OPTIONAL :: announce
  CHARACTER(mfl)                :: zipname
  CHARACTER(mcl)                :: cmd
  INTEGER                       :: i, ios
  LOGICAL                       :: file_exists, announce1
  
  announce1 = .TRUE.
  IF(PRESENT(announce)) announce1 = announce
  
  !! Loop over all files
  DO i=1,SIZE(fn)
  
    !! Announce what's happening
    IF(announce1) CALL Prnt("Zipping file ["//TRIM(fn(i))//"] ...")
    
    !! Make sure the file exists
    INQUIRE(FILE=fn(i), EXIST=file_exists)
    IF(.NOT.file_exists) THEN
      
      CALL PrntW("File ["//TRIM(fn(i))//"] cannot be zipped"//&
                        " because it does not exist.")
    
    !! If it does, make a system call to a zip command
    ELSE
      
      zipname = TRIM(fn(i))//".zip"
      cmd = "zip -r9Xm "//TRIM(zipname)//" "//TRIM(fn(i)) 
      CALL SYSTEM(cmd, ios)
      
      !! Check the result code of the system call for success or failure
      IF(ios==0) THEN
        fn(i) = zipname
        IF(announce1) CALL Prnt("File zipped.")
      ELSE
        CALL PrntW("System call '"//TRIM(cmd)//"'' failed with code "//i2cp(ios)//&
                   ". Please make sure a system call to 'zip' is possible.")
      ENDIF
      
    ENDIF
    
  ENDDO
    
  RETURN

END SUBROUTINE ZipFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE SortFile(fn, temp, secondary, announce, skip1)
  IMPLICIT NONE
  INTEGER, PARAMETER            :: announce_size = 10000000 !! 10 MB
  CHARACTER(*), INTENT(INOUT)   :: fn(:)
  LOGICAL, INTENT(IN)           :: temp
  LOGICAL, INTENT(IN), OPTIONAL :: secondary, announce, skip1
  CHARACTER(mfl), ALLOCATABLE   :: fn_temp(:)
  TYPE(TMPTYPE), ALLOCATABLE    :: tmp(:)
  CHARACTER(mltl), ALLOCATABLE  :: header(:,:)
  CHARACTER(moil), ALLOCATABLE  :: body(:,:)
  REAL(dpp), ALLOCATABLE        :: column(:)
  INTEGER, ALLOCATABLE          :: order(:)
  INTEGER                       :: nresults, nresults2, i, ios, ifile
  CHARACTER(mtsl)               :: startt_full, stopt_full, runtime
  REAL(dpp)                     :: startt, stopt
  LOGICAL                       :: secondary1, announce1
  
  !! If no sorting wanted, exit
  IF(out_sort_sel=="0" .OR. out_sort_col==0) RETURN

  !! Otherwise check for valid selection
  IF(temp) THEN
    IF(ALL(out_sort_sel/=sort_selections)) &
      CALL PrntE("Missing data specifier by which to sort the temporary"//&
                      " output file. Hint: Use --out-sort-by x where x can"//&
                      " take one of the following values: A, D, C, P, PA,"//&
                      " PB, PC, PD, 0.", Q=.TRUE.)    
  ELSEIF(out_sort_col<=0) THEN
    IF(out_sort_col<0) &
      CALL PrntE("Invalid sorting column selection. Use --out-sort-by x"//&
                      " to specify the sorting column.", Q=.TRUE.)
  ENDIF

  !! Decide whether to announce progress
  announce1 = .FALSE.
  IF(SIZE(fn)>1) THEN
    announce1 = .TRUE.
  ELSE
    IF(GetFileSize(fn(1)) > announce_size) announce1 = .TRUE.
  ENDIF
  IF(PRESENT(announce))  announce1 = announce

  secondary1 = .FALSE.
  IF(PRESENT(secondary)) secondary1 = secondary

  CALL GetCurrentTime(startt_full, startt)

  !! Skip 1 line
  IF(PRESENT(skip1)) THEN 
    IF(skip1) CALL Prnt0("")
  ENDIF
  
  !! Loop over files to sort
  DO ifile = 1,SIZE(fn)
    
    !! Announce sorting
    IF(temp) THEN
      !! Check for missing fn
      IF(LEN_TRIM(fn(ifile))==0) &
        CALL PrntE("Missing temporary results file to sort.", Q=.TRUE.)
        
      !! Announce sorting
      CALL Prnt("Sorting temporary output file(s) ... ")
  
      !! Read the output file header and body
      CALL ReadTempFile(fn(ifile), tmp, announce=announce1, ios=ios)
                          
      nresults = SIZE(tmp)
      
      IF(nresults==0) CYCLE
  
    ELSE
    
      !! Check for missing filename
      IF(LEN_TRIM(fn(ifile))==0) &
        CALL PrntE("Missing results file to sort.", Q=.TRUE.)
        
      !! Announce sorting
      CALL Prnt("Sorting output file(s) ... ")
  
      !! Read the output file header and body
      CALL ReadOutputFile(fn(ifile), header, body, out_ncols, announce1, ios=ios)
                        
      nresults = SIZE(body,1)
  
      !! Check for an empty file
      IF(nresults==1 .AND. ALL(body=="")) GOTO 50
       
      !! Check the given value of out_sort_col
      IF(SIZE(body,2) < out_sort_col) THEN
        ios = -99 
        GOTO 100
      ENDIF
  
    ENDIF
    
    !! If no results read, just quit
    IF(ios==-55) RETURN
  
    !! Check for other errors
    IF(ios/=0) GOTO 100
    
    !! ----------------------------------------------------------- !!
    !!            FIND THE (ASCENDING) ORDER OF P-VALUES           !!
    !! ----------------------------------------------------------- !!
    
    !! Announce sorting phase
    IF(announce1) CALL Prnt("Sorting the values ... ") 
    
    !! Store values in column
    CALL ResizeVar(column, nresults, NAp)
    
    IF(temp) THEN
      IF(out_sort_sel=="A4") column = tmp%AS4p
      IF(out_sort_sel=="A1") column = tmp%AS1p
      IF(out_sort_sel=="D")  column = tmp%DSp
      IF(out_sort_sel=="C")  column = tmp%CSp
      IF(out_sort_sel=="P4") column = tmp%PTAS4p
      IF(out_sort_sel=="P1") column = tmp%PTAS1p
      IF(out_sort_sel=="PA") column = tmp%PT4p
      IF(out_sort_sel=="PB") column = tmp%PT1p
      IF(out_sort_sel=="PC") column = tmp%PTPO4p
      IF(out_sort_sel=="PD") column = tmp%PTPO1p
    ELSE
      DO i=1,nresults
        IF(body(i,out_sort_col) == NA) CYCLE
        column(i) = c2r(body(i,out_sort_col))
      ENDDO
    ENDIF
    
    !! Sort the file only if any results available
    IF(ALL(column == NAp)) CYCLE
    
    !! Prepare the vector with indices
    CALL ResizeVar(order, nresults)
    order = (/(i, i=1,nresults)/)

    !! Call the QUICKSORT subroutine
    CALL qsort(column, nresults, order)
    !CALL QuickSort(column, order)
    
    IF(temp) THEN
      !! Reorder the rows of tmp
      tmp = tmp(order)
    ELSE
      !! Reorder the rows of body
      body = body(order,:)
      !! Change the column with line number
      body(:,1) = i2c((/(i, i=1,nresults)/), c=.FALSE., s=.FALSE.)
    ENDIF
    
    !! Do secondary sort
    IF(secondary1) THEN
      
      !! Count the non-NA pvalues
      nresults2 = CountOccurrence(column, NAp)
      IF(nresults2 > 0) THEN
      
        !! Reinitialize order and column
        order = (/(i, i=1,nresults)/)
        column = NAp
        IF(temp) THEN
          column = tmp%PTAS4p
        ELSE
          DO i=1,nresults
            IF(body(i,out_sort_col)==NA) CYCLE
            column(i) = c2r(body(i,out_sort_col))
          ENDDO
        ENDIF
        
        !! Perform the secondary sort
        CALL qsort(column(nresults-nresults2+1:), nresults2, &
                   order(nresults-nresults2+1:))
  
        !! Reorder values
        IF(temp) THEN
          !! Reorder the rows of tmp
          tmp = tmp(order)
        ELSE
          !! Reorder the rows of body
          body = body(order,:)
          !! Change the column with line number
          body(:,1) = i2c((/(i, i=1,nresults)/), c=.FALSE., s=.FALSE.)
        ENDIF
    
      ENDIF
      
    ENDIF
      
    !! Save (temporary) the output file
    CALL ResizeVar(fn_temp, 1)
    fn_temp = fn(ifile)
    IF(temp) THEN
      CALL RemoveSubstrEnd(fn_temp, ".sorted"//tmp_ext)
      CALL RemoveSubstrEnd(fn_temp, tmp_ext)
      fn_temp = TRIM(fn_temp(1))//".sorted"//tmp_ext
      CALL WriteTempFile(fn_temp, tmp, save_all=.TRUE., announce=announce1,&
                         ios=ios)
    ELSE
      CALL RemoveSubstrEnd(fn_temp, ".sorted"//out_ext)
      CALL RemoveSubstrEnd(fn_temp, out_ext)
      fn_temp = TRIM(fn_temp(1))//".sorted"//out_ext
      !! Save the output file
      CALL SaveOutputFile(fn_temp(1), header, body, announce=announce1, ios=ios)
    ENDIF
  
    !! Check if all ok
    IF(ios/=0) GOTO 100
  
    !! When file sorted delete the old file and change the name of the temp file
    !! to the newly created sorted file 
    IF(out_discard_unsorted) THEN
      CALL DeleteFile(fn(ifile))
      fn(ifile) = fn_temp(1)
      !!! After successful writing, rename the file
      !CALL RenameFile(fn_temp(1), fn(ifile), ios)
      !!! Check if all ok
      IF(ios/=0) GOTO 100
    ENDIF
    
  ENDDO

  50 CONTINUE

  CALL GetCurrentTime(stopt_full, stopt)
  runtime = GetTimeDifference(startt, stopt)

  !! Announce the end of sorting
  IF(temp) THEN
    CALL Prnt("Temporary output file(s) was sorted in "//&
                   TRIM(runtime)//".")
  ELSE
    CALL Prnt("Output file(s) was sorted after "//TRIM(runtime)//".")
  ENDIF

  RETURN

  100 CONTINUE
  
  IF(temp) THEN
    CALL Prnt("TMP file(s) could not be sorted (IOSTAT "//i2c(ios)//").")
  ELSE
    CALL Prnt("Output file(s) could not be sorted (IOSTAT "//i2c(ios)//").")
  ENDIF
  
  RETURN

END SUBROUTINE SortFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! SUBROUTINE SortOutputFile(fn, announce)
!   IMPLICIT NONE
!   INTEGER, PARAMETER            :: announce_size = 100000000 !! 100 MB
!   CHARACTER(*), INTENT(INOUT)   :: fn
!   LOGICAL, INTENT(IN), OPTIONAL :: announce
!   CHARACTER(mfl)                :: fn_temp
!   CHARACTER(mltl), ALLOCATABLE  :: header(:,:)
!   CHARACTER(mstl), ALLOCATABLE  :: body(:,:)
!   REAL(dpp), ALLOCATABLE        :: column(:)
!   INTEGER, ALLOCATABLE          :: order(:)
!   INTEGER                       :: nresults, i, ios
!   CHARACTER(mtsl)               :: startt_full, stopt_full, runtime_total
!   REAL(dpp)                     :: startt, stopt
!   LOGICAL                       :: announce1
!   
!   !! Decide whether to announce progress
!   announce1 = .FALSE.
!   IF(GetFileSize(fn)>announce_size) announce1 = .TRUE.
!   IF(PRESENT(announce)) announce1 = announce
! 
!   !! Check for missing filename
!   IF(LEN_TRIM(fn)==0) &
!     CALL PrntE("Missing results file to sort.", Q=.TRUE.)
!   
!   CALL GetCurrentTime(startt_full, startt)
! 
!   !! Announce sorting
!   CALL Prnt("Sorting output file(s) ... ")
! 
!   !! Read the output file header and body
!   CALL ReadOutputFile(fn, header, body, out_ncols, announce=announce1, &
!                       warn_empty=.FALSE., ios=ios)
!                       
!   !! If no results read, just quit
!   IF(ios==-55) RETURN
! 
!   !! Check for other errors
!   IF(ios/=0) GOTO 100
!   
!   nresults = SIZE(body,1)
!   
!   !! Check for an empty file
!   IF(nresults==1 .AND. ALL(body=="")) GOTO 50 
! 
!   !! Check the given value of out_sort_col
!   IF(SIZE(body,2) < out_sort_col .OR. out_sort_col<=0) THEN
!     ios = -99 
!     GOTO 100
!   ENDIF
! 
!   !! ----------------------------------------------------------- !!
!   !!            FIND THE (ASCENDING) ORDER OF RESULTS            !!
!   !! ----------------------------------------------------------- !!
!   
!   !! Announce sorting phase
!   IF(announce1) CALL Prnt("Sorting the values ... ") 
!   
!   !! Store values in column
!   CALL ResizeVar(column, nresults, NAp)
!   DO i=1,nresults
!     IF(body(i,out_sort_col)/=NA) column(i) = c2r(body(i,out_sort_col))
!   ENDDO
!   
!   !! Sort the file only if any results available
!   IF(ANY(column/=NAp)) THEN
!     CALL ResizeVar(order, nresults)
!     order = (/(i, i=1,nresults)/)
!     CALL qsort(column, nresults, order)
!   
!     !! Reorder the rows of body
!     body = body(order,:)
!   
!     !! Change the column with line number
!     body(:,1) = i2c((/(i, i=1,nresults)/), c=.FALSE., s=.FALSE.) 
!     
!     fn_temp = TRIM(fn)//".sorted"//out_ext
!     
!     !! Save the output file
!     CALL SaveOutputFile(fn_temp, header, body, announce=.FALSE., &
!                         ios=ios)
! 
!     !! Check if all ok
!     IF(ios/=0) GOTO 100
! 
!     !! After successful writing, rename the file
!     !CALL RenameFile(fn_temp, fn, ios)
!     
!     !!! Check if all ok
!     !IF(ios/=0) GOTO 100
!   ENDIF
! 
!   50 CONTINUE
! 
!   CALL GetCurrentTime(stopt_full, stopt)
!   runtime_total = GetTimeDifference(startt, stopt)
! 
!   !! Announce the end of sorting
!   CALL Prnt("Output file(s) was sorted after "//TRIM(runtime_total)//".")
! 
!   RETURN
! 
!   100 CONTINUE
!   
!   CALL Prnt("Output file(s) could not be sorted (IOSTAT "//i2c(ios)//&
!                  ").")
!   
!   RETURN
! 
! END SUBROUTINE SortOutputFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE LH(string, head, nc, nc2, advance, skip1, skip2, cmt, no_space, &
              no_colnames, out_anyway)
!! Write the output file header
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)           :: string
  CHARACTER(*), INTENT(IN), OPTIONAL :: head, advance, cmt
  INTEGER, INTENT(INOUT), OPTIONAL   :: nc
  INTEGER, INTENT(IN), OPTIONAL      :: nc2, skip1, skip2
  LOGICAL, INTENT(IN), OPTIONAL      :: no_space, no_colnames, out_anyway
  CHARACTER(mltl)                    :: line
  CHARACTER                          :: advance1, cmt1
  INTEGER                            :: i
  LOGICAL                            :: colnames, space, out_anyway1
  
  advance1 = 'Y'
  colnames = .TRUE.
  space = .TRUE.
  out_anyway1 = .FALSE.
  IF(PRESENT(advance)) advance1 = upcasef(advance(1:1))
  IF(PRESENT(no_colnames)) colnames = .NOT.no_colnames
  IF(PRESENT(no_space)) space = .NOT.no_space
  IF(PRESENT(out_anyway)) out_anyway1 = out_anyway
  
  cmt1 = comment
  IF(PRESENT(cmt)) cmt1 = cmt
  
  line = ""
  
  !! Attach head at the beginning of line
  IF(PRESENT(head) .AND. out_header_names) THEN
    IF(LEN_TRIM(head)>0) line = TRIM(head)//out_cs
    IF(PRESENT(nc) .AND. colnames) THEN
      IF(PRESENT(nc2)) THEN
        out_cnams(nc:nc2) = TRIM(head)
      ELSE
        out_cnams(nc) = TRIM(head)
      ENDIF
    ENDIF
  ENDIF
  
  !! Generate the line and increase line counter
  IF(PRESENT(nc)) THEN
    IF(.NOT.PRESENT(nc2)) THEN
      IF(space) THEN
        line = i2cp(nc)//". "//TRIM(line)
      ELSE
        line = i2cp(nc)//"."//TRIM(line)
      ENDIF
      nc = nc+1
    ELSE
      IF(space) THEN
        line = i2cp(nc)//".-"//i2cp(nc2)//"."//TRIM(line)
      ELSE
        line = i2cp(nc)//".-"//i2cp(nc2)//". "//TRIM(line)
      ENDIF
      nc = nc2+1
    ENDIF
  ENDIF

  !! If no legend wanted, don't print it
  IF(.NOT.out_legend .AND. .NOT.out_anyway1) RETURN
  
  !! Attach the string (ADJUSTL in case line was empty)
  IF(space) THEN
    line = ADJUSTL(TRIM(line)//" "//TRIM(string))
  ELSE
    line = ADJUSTL(TRIM(line)//TRIM(string))
  ENDIF
  
  !! Skip lines before writing
  IF(PRESENT(skip1)) THEN
    DO i=1,skip1; WRITE(uout, '(A)') TRIM(cmt1); ENDDO
  ENDIF 

  !! Write the line
  IF(LEN_TRIM(cmt1)>0) WRITE(uout, '(A)', ADVANCE='NO') cmt1//"  " 
  WRITE(uout, '(A)', ADVANCE='NO') TRIM(line) 
  IF(advance1=='Y') WRITE(uout, '(A)')

  !! Skip lines after writing
  IF(PRESENT(skip2)) THEN
    DO i=1,skip2; WRITE(uout, '(A)') TRIM(cmt1); ENDDO
  ENDIF 
  
  RETURN

END SUBROUTINE LH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrintHeader(hc, starttime, nfile)
  IMPLICIT NONE
  CHARACTER, INTENT(IN)    :: hc
  CHARACTER(*), INTENT(IN) :: starttime
  INTEGER, INTENT(IN)      :: nfile
  INTEGER                  :: k, kk, nc, df
  CHARACTER(mttl)          :: Tkind, Tname, Tpname, Tpcname, Tmtcname, &
                              Terrname, Tncaname, Tnconame, Tdf
  CHARACTER(mmtl)          :: Tfullname, Tfull
  
  !! Write the detailed header only into the first file
  IF(nfile == 1) THEN
  
    nc = 1
    CALL LH(TRIM(program_fullname)//" output file created at "//starttime)
    CALL LH("Filename : ["//TRIM(out_file)//"]")
    CALL LH("For details on input and settings see the accompanying log file")
    CALL LH("Log file : ["//TRIM(log_file)//"]")
    IF(ndig_pval<0) &
      CALL LH("P-values compressed to 3 digit format 'XYZ' means 'X.0E-YZ'.") 
    CALL LH("Column guide :")
    
    IF(out_testid) CALL LH("Line number", "Id", nc)

    !! Global error indicator
    IF(out_errcode .AND. out_errcode_glob) &
      CALL LH("Global error code", "GErr", nc)

    !! ******************************************************************** !!
    IF(do_all_pretests) THEN
      IF(out_statistic) CALL LH("Pre-test statistics self centered", "PTX", nc)
      CALL LH("Pre-test statistics p-value self centered", "PTpcoX", nc)
      IF(out_statistic) CALL LH("AS self centered", "AScoX", nc)
      CALL LH("AS p-value self centered", "ASpcoX", nc)
      IF(out_statistic) CALL LH("Disjoint score self centered", "DScoX", nc)
      CALL LH("Disjoint score p-value self centered", "DSpcoX", nc)

      IF(out_statistic) CALL LH("Pre-test statistics cases centered", "PTY", nc)
      CALL LH("Pre-test statistics p-value cases centered", "PTpcoY", nc)
      IF(out_statistic) CALL LH("AS cases centered", "AScoY", nc)
      CALL LH("AS p-value cases centered", "ASpcoY", nc)
      IF(out_statistic) CALL LH("Disjoint score cases centered", "DScoY", nc)
      CALL LH("Disjoint score p-value cases centered", "DSpcoY", nc)

      IF(out_statistic) CALL LH("Pre-test statistics noncentered", "PT0", nc)
      CALL LH("Pre-test statistics p-value noncentered", "PTpco0", nc)
      IF(out_statistic) CALL LH("AS noncentered", "ASco0", nc)
      CALL LH("AS p-value noncentered", "ASpco0", nc)
      IF(out_statistic) CALL LH("Disjoint score noncentered", "DSco0", nc)
      CALL LH("Disjoint score p-value noncentered", "DSpco0", nc)

      IF(out_statistic) CALL LH("Pre-test statistics for difference-based pretest noncentered", "PTcc", nc)
      CALL LH("Pre-test statistics p-value for difference-based pretest noncentered", "PTpcc", nc)
      IF(out_statistic) CALL LH("AS for difference-based pretest noncentered", "AScc", nc)
      CALL LH("AS p-value for difference-based pretest noncentered", "ASpcc", nc)
      IF(out_statistic) CALL LH("Disjoint score for difference-based pretest noncentered", "DScc", nc)
      CALL LH("Disjoint score p-value for difference-based pretest noncentered", "DSpcc", nc)
    ENDIF
    !! ******************************************************************** !!
    
    !! Choose name tag for pre-test and post-test
    Tkind = ""
    Tfull = ""
    Tdf = ""
    df = 0
    SELECT CASE (WT1)
      CASE (T1co)
        Tkind = "co"
        Tfull = "control-only"
      CASE (T1ca)
        Tkind = "ca"
        Tfull = "case-only"
      CASE (T1cc) 
        Tkind = "cc"
        Tfull = "control+case"
      CASE (T1po)
        Tkind = "po"
        Tfull = "control+case (pooled)"
      CASE (T1sc) 
        Tkind = "cs"
        Tfull = "classical score"
    END SELECT

    !! PRE-TEST info
    DO k=0,7
    
      IF(k==0 .AND. (.NOT.doAS4 .OR. .NOT.doPTAS4)) CYCLE
      IF(k==1 .AND. (.NOT.doAS1 .OR. .NOT.doPTAS1)) CYCLE
      IF(k==2 .AND. (.NOT.doDS4 .OR. .NOT.doPTDS4)) CYCLE
      IF(k==3 .AND. (.NOT.doDS1 .OR. .NOT.doPTDS1)) CYCLE
      IF(k==4 .AND. .NOT.doPO4) CYCLE
      IF(k==5 .AND. .NOT.doPO1) CYCLE
      IF(k==6 .AND. .NOT.doPT4) CYCLE
      IF(k==7 .AND. .NOT.doPT1) CYCLE

      df = 1
      IF(ANY(k==(/0,2,4,6/))) df = T1_df

      IF(k==0 .OR. k==1) THEN
        Tname = "PTAS"
        Tdf = "-"//i2cp(df)
        Tfullname = "AS-"//TRIM(Tfull)
        IF(Tkind/="cs") Tfullname = TRIM(Tfullname)//Tdf 
      ELSEIF(k==2 .OR. k==3) THEN
        Tname = "PTDS"
        Tdf = "-"//i2cp(df)
        Tfullname = "DS-"//TRIM(Tfull)
        IF(Tkind/="cs") Tfullname = TRIM(Tfullname)//Tdf 
      ELSEIF(k==4 .OR. k==5) THEN
        Tname = "PTPO"
        Tdf = "-"//i2cp(df)
        Tfullname = "pooled"//Tdf 
      ELSEIF(k==6 .OR. k==7) THEN
        Tname = "PTASco"
        Tdf = "-"//i2cp(df)
        Tfullname = Tfull
        IF(Tkind/="cs") Tfullname = TRIM(Tfullname)//Tdf 
      ENDIF
        
      IF(WT1>0) Tname = TRIM(Tname)//i2cp(df)
      Tpname = TRIM(Tname)//"p"
      Terrname = TRIM(Tname)//"err"
      Tncaname = TRIM(Tname)//"nca"
      Tnconame = TRIM(Tname)//"nco"

      IF(out_statistic) &
        CALL LH("Pre-test "//TRIM(Tfullname)//" statistics", Tname, nc)
      
      IF(k==0 .AND. out_sort_sel=="PA") out_sort_col = nc
      IF(k==2 .AND. out_sort_sel=="PB") out_sort_col = nc
      IF(k==3 .AND. out_sort_sel=="PC") out_sort_col = nc
      IF(k==4 .AND. out_sort_sel=="PD") out_sort_col = nc
      CALL LH(TRIM(Tname)//" p-values (raw)", Tpname, nc)
      
      IF(out_errcode .AND. .NOT.out_errcode_glob) &
        CALL LH(TRIM(Tname)//" error code", Terrname, nc)
    
      IF(out_ss) THEN
        CALL LH(TRIM(Tname)//" number of cases", Tncaname, nc)
        CALL LH(TRIM(Tname)//" number of controls", Tnconame, nc)
      ENDIF
      
    ENDDO

    !! POST-TEST info
    DO k=1,7
    
      IF(k==1 .AND. .NOT.doAS4) CYCLE
      IF(k==2 .AND. .NOT.doAS1) CYCLE
      IF(k==3 .AND. .NOT.doDS4) CYCLE
      IF(k==4 .AND. .NOT.doDS1) CYCLE
      IF(k==5 .AND. (.NOT.doPO4 .OR. .NOT.doPO4CS)) CYCLE
      IF(k==6 .AND. (.NOT.doPO1 .OR. .NOT.doPO1CS)) CYCLE
      IF(k==7 .AND. .NOT.doCS) CYCLE

      df = 1
      IF(ANY(k==(/1,3,5/))) df = T1_df

      IF(k<=2) THEN
        Tname = "AS"!//TRIM(Tkind)
        Tfullname = "Post-test AS-"//TRIM(Tfull)
        IF(WT1>0) THEN
          Tname = TRIM(Tname)//TRIM(i2cp(df))
          Tfullname = TRIM(Tfullname)//"-"//TRIM(i2cp(df))
        ENDIF
      ELSEIF(k>=3 .AND. k<=4) THEN
        Tname = "DS"!//TRIM(Tkind) 
        Tfullname = "Post-test DS-"//TRIM(Tfull)
        IF((k==3 .AND. doDS1) .OR. k==4 .AND. WT1>0) THEN
          Tname = TRIM(Tname)//TRIM(i2cp(df))
          Tfullname = TRIM(Tfullname)//"-"//TRIM(i2cp(df))
        ENDIF
      ELSEIF(k>=5 .AND. k<=6) THEN
        Tname = "CSpo"
        IF(WT1>0) Tname = TRIM(Tname)//i2cp(df)
        Tfullname = "CS with pooled-"//TRIM(i2cp(df))
      ELSEIF(k==7) THEN
        Tname = "CS"
        Tfullname = "Classical score"
      ENDIF
      
      Tpname = TRIM(Tname)//"p"
      Tpcname = TRIM(Tname)//"pC"
      Tmtcname = TRIM(Tname)//"mtc"
      Terrname = TRIM(Tname)//"err"
      Tncaname = TRIM(Tname)//"nca"
      Tnconame = TRIM(Tname)//"nco"

      IF(out_statistic) &
        CALL LH(TRIM(Tfullname)//" statistics", Tname, nc)

      IF(k==1 .AND. out_sort_sel=="A4") out_sort_col = nc
      IF(k==2 .AND. out_sort_sel=="A1") out_sort_col = nc
      IF(k==3 .AND. out_sort_sel=="D")  out_sort_col = nc
      IF(k==7 .AND. out_sort_sel=="C")  out_sort_col = nc
      CALL LH(TRIM(Tname)//" p-values (raw)", Tpname, nc)

      IF(out_pvalcor) &
        CALL LH(TRIM(Tname)//" p-values (Bonferroni corrected by "//&
                TRIM(r2c(MTC(k)))//")", Tpcname, nc)

      IF(out_mtc) &
        CALL LH(TRIM(Tname)//" Bonferroni correction ", Tmtcname, nc)

      IF(out_errcode .AND. .NOT.out_errcode_glob) &
        CALL LH(TRIM(Tname)//" error code", Terrname, nc)

      IF(out_ss) THEN
        CALL LH(TRIM(Tname)//" number of cases", Tncaname, nc)
        CALL LH(TRIM(Tname)//" number of controls", Tnconame, nc)
      ENDIF

    ENDDO
    
    !! Interaction effect
    IF(out_epi_effect) THEN
      CALL LH("Estimated interaction effect", "beta3", nc)
      CALL LH("Odds ratio of interaction effect (exp(beta3))", "OR", nc)
    ENDIF
    
    !! Error code and loci id info
    IF(out_loc_info) THEN
      CALL LH("Chromosome of the 1st locus in a test", "Chr1", nc)
      CALL LH("RS number of the 1st loci in a test", "RS1", nc)
      CALL LH("Chromosome of the 2nd locus in a test", "Chr2", nc)
      CALL LH("RS number of the 2nd loci in a test", "RS2", nc)
    ENDIF

    !! Minor allele frequencies (entire sample - used for automatic PHASE1 level)
    IF(out_maf) THEN
      CALL LH("Minor allele frequency of 1st locus (all) (based on the"//&
              " entire sample or all samples if the value given by"//&
              " --nsamples is larger than 1, same for all MAF below)", &
              "maf1", nc)
      CALL LH("Minor allele frequency of 2nd locus (all)", "maf2", nc)
      IF(out_maf_coca) THEN
        CALL LH("Minor allele frequency of 1st locus (controls only)", "mafco1", nc)
        CALL LH("Minor allele frequency of 2nd locus (controls only)", "mafco2", nc)
        CALL LH("Minor allele frequency of 1st locus (cases only)", "mafca1", nc)
        CALL LH("Minor allele frequency of 2nd locus (cases only)", "mafca2", nc)
      ENDIF
      IF(auto_level_report) THEN
        CALL LH("Allele frequency 1 used for PHASE1 level determination ", "af1o", nc)
        CALL LH("Allele frequency 2 used for PHASE1 level determination ", "af2o", nc)
      ENDIF
    ENDIF
    
    !! Used and automatic PHASE1 level
    IF(out_level_S1) THEN
      CALL LH("Used PHASE1 level", "lev1u", nc)
      IF(auto_level_report) THEN
        CALL LH("Automatic PHASE1 level", "lev1a", nc)
        CALL LH("Estimated non-centrality parameter (determines"//&
                " automatic PHASE1 level)", "ncp", nc)
        CALL LH("Estimated slope parameter (determines automatic"//&
                " PHASE1 level)", "slope", nc)
      ENDIF
    ENDIF
    
    !! Output allele names for th two loci
    IF(out_alleles) CALL LH("Allele names", "alleles", nc)

    !! Print user supplied value of the actual prevalence
    IF(res_prev>=zero) CALL LH("Actual prevalence of cases", "aprev", nc)

    !! Debugging info
    IF(out_debug_info) THEN

      IF(out_epi_effect) &
        CALL LH("Estimate of standard error of beta3", "sebeta3", nc)
      IF(sDebug<1) GOTO 110
      CALL LH("Estimate of prevalence", "beta0", nc)
      IF(sDebug<2) GOTO 110
      CALL LH("Estimate of first main effect", "beta1", nc)
      IF(sDebug<3) GOTO 110
      CALL LH("Estimate of second main effect", "beta2", nc)
      IF(sDebug<4) GOTO 110
      CALL LH("Estimate of standard error of beta0", "sebeta0", nc)
      IF(sDebug<5) GOTO 110
      CALL LH("Estimate of standard error of beta1", "sebeta1", nc)
      IF(sDebug<6) GOTO 110
      CALL LH("Estimate of standard error of beta2", "sebeta2", nc)
      
    ENDIF
    
    IF(out_debug_info .OR. out_variance) THEN
    
      IF(sDebug<7) GOTO 110
      CALL LH("Non-standardized classical score statistics", "S0", nc)
      IF(sDebug<8) GOTO 110
      CALL LH("Variance of the non-standardized classical score"//&
              " statistics", "varS0", nc)
      !! AS4
      IF(sDebug<9) GOTO 110
      IF(doAS4) CALL LH("Variance of the regression vector b4T4", "varb4T4", nc)
      IF(sDebug<10) GOTO 110
      IF(doAS4) CALL LH("Variance of AS4 (varAS4 = varS0 + varb4T4)", "varAS40", nc)
      IF(sDebug<11) GOTO 110
      IF(doAS4) THEN
        CALL LH("Regression term for AS4", "b4T4", nc)
        IF(sDebug<12) GOTO 110
        CALL LH("Non-standardized AS4 statistics", "AS40", nc)
      ENDIF

      !! AS1
      IF(sDebug<13) GOTO 110
      IF(doAS1) CALL LH("Variance of the regression vector b1T1", "varb1T1", nc)
      IF(sDebug<14) GOTO 110
      IF(doAS1) CALL LH("Variance of AS10 (varAS10 = varS0 + varb1T1)", "varAS10", nc)
      IF(sDebug<15) GOTO 110
      IF(doAS1) THEN
        CALL LH("Regression term for AS1", "b1T1", nc)
        IF(sDebug<16) GOTO 110
        CALL LH("Non-standardized AS1 statistics", "AS10", nc)
      ENDIF

    ENDIF
    
    110 CONTINUE
    
    !! Pretest vector coordinates
    IF(out_S1_vector) THEN
      CALL LH("Pre-test generating vector for AS4", "b4T4", nc, nc+dmV-1)  
      CALL LH("Pre-test generating vector for AS1", "b1T1", nc, nc+dmV-1)
    ENDIF  
    
    CALL LH("An error code '0' means no error; -1 mean the test"//&
            " was not performed.", skip1=1)
    CALL LH("Value '"//NA//"'' indicates that the corresponding"//&
            " statistic or p-value was not computed.", skip2=1)
    
    !! Negate the last increase of nc
    nc = nc - 1
    out_ncols = nc
      
  ELSE
    
    !! Write short initial information into the new file
    WRITE(uout,'(A)') &
      comment//" Continuation of the output file created at "//starttime
    WRITE(uout,'(A)') comment//" Output file part number "//i2cp(nfile)
    
  ENDIF
  
  !! Print column numbers
  IF(.FALSE.) THEN
    WRITE(uout, '(A)', ADVANCE='NO') TRIM(comment)
    DO k=1,out_ncols
      kk = k
      WRITE(uout, '(A)', ADVANCE='NO') i2cp(kk)
      IF(k < out_ncols) WRITE(uout, '(A)', ADVANCE='NO') out_cs
    ENDDO
    WRITE(uout, '(A)')
  ENDIF
  
  !! Print the column names (if there are any)
  IF(ANY(out_cnams/="")) THEN
    WRITE(uout, '(A)', ADVANCE='NO') TRIM(hc)
    DO k=1,out_ncols
      WRITE(uout, '(A)', ADVANCE='NO') TRIM(out_cnams(k))
      IF(k < out_ncols) WRITE(uout, '(A)', ADVANCE='NO') out_cs
    ENDDO
    WRITE(uout, '(A)')
  ENDIF
  
  RETURN

END SUBROUTINE PrintHeader

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
FUNCTION FmtPv(pval, ndig) RESULT(cpval)
!! Function that convert p-value into a scientific format and optionally 
!! compresses it into a 3 digit format if ndig==-1 in the sence that X.E-YZ
!! is coded as XYZ  
  IMPLICIT NONE
  REAL(dpp), INTENT(IN) :: pval
  INTEGER, INTENT(IN)   :: ndig
  CHARACTER(mstl)       :: cpval, fmt
    
  !! Convert the p-value into scientific format
  fmt = "(ES"//TRIM(i2cp(6+MAX(0,ndig)))//"."//TRIM(i2cp(MAX(0,ndig)))//")"
  WRITE(cpval, fmt) pval

  !! Make sure there are no trialing spaces at the beginning 
  cpval = ADJUSTL(cpval)
  
  !! Compress the p-value if ndig==0
  IF(ndig==0) cpval = cpval(1:1)//cpval(3:)
  IF(ndig==-1) cpval = cpval(1:1)//cpval(5:)
  
  RETURN
END FUNCTION FmtPv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE FmtPvR(CP, P, ndig)
!! Subroutine that does the same as the function FmtPv above
  IMPLICIT NONE
  CHARACTER(mstl), INTENT(OUT) :: CP 
  REAL(dpp), INTENT(IN)        :: P
  INTEGER, INTENT(IN)          :: ndig
  REAL(dpp)                    :: PP
  CHARACTER(mttl)              :: FMT
    
  !! Store p-value locally
  PP = P
  
  !! Make sure it's not too small to fit within the 2 digit scientific format
  IF(PP < 1E-98_dpp) PP = zero
  
  !! Convert the p-value into scientific format
  FMT = "(ES"//TRIM(i2cp(6+MAX(0,ndig)))//"."//TRIM(i2cp(MAX(0,ndig)))//")"
  WRITE(CP, FMT) PP
  
  !! Make sure there are no trialing spaces at the beginning 
  CP = ADJUSTL(CP)

  !! Compress the p-value if ndig==-1
  IF(ndig==-1) CP = CP(1:1)//CP(5:6)
  
  RETURN
END SUBROUTINE FmtPvR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE WriteOutputRecord(n, tmp, pCor, mtc, fv, fp, rec_len)
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN)  :: n
  TYPE(TMPTYPE), INTENT(IN) :: tmp
  REAL(dpp), INTENT(IN)     :: pCor(:), mtc(:)
  CHARACTER(*), INTENT(IN)  :: fp, fv
  INTEGER, INTENT(OUT)      :: rec_len
  CHARACTER(mltl)           :: LL
  CHARACTER(mstl)           :: L1, L2, fmv, fmv3, fmv9, fmp, fmm
  CHARACTER(mml)            :: ch1, ch2, rs1, rs2
  CHARACTER                 :: cs, a1, a2, b1, b2
  INTEGER                   :: ii, jj, k, j, ierr, nv, icycle, &
                               Err, Nca, Nco, GErrN, GErrP
  REAL(dpp)                 :: mafij(2), mafcoij(2), mafcaij(2), beta1(0:3), &
                               se_beta1(0:3), S, Sp, T, Tp
  INTEGER(iks)              :: sts1(SIZE(sts)), xi(SIZE(sts)), yj(SIZE(sts))
  LOGICAL                   :: all_pval_NA

  cs = out_cs
  
  !! Set the fortmat to print variables in output file
  fmv = '('//TRIM(fv)//')'
  fmv3 = '(3'//TRIM(fv)//')'
  fmv9 = '(9'//TRIM(fv)//')'
  fmp = '('//TRIM(fp)//')'
  fmm = "F0."//i2cp(ndig_mafs)
  
  ii = tmp%Pos(1)
  jj = tmp%Pos(2)
  
  !! Assign default (NA) mapping info
  ch1 = NA; rs1 = NA; ch2 = NA; rs2 = NA
  a1 = NA; a2 = NA; b1 = NA; b2 = NA
  GErrN = 0; GErrP = 0

  !! If non-NA mapping info known, assign that
  IF(ALLOCATED(MAPA) .AND. MIN(ii,jj)>0) THEN
    ch1 = MAPA(ii, map_chcol)
    rs1 = MAPA(ii, map_rscol)
    ch2 = MAPA(jj, map_chcol)
    rs2 = MAPA(jj, map_rscol)
    IF(SIZE(MAPA,2)>=MAX(map_a1col, map_a2col)) THEN
      a1 = MAPA(ii, map_a1col)(1:1)
      a2 = MAPA(ii, map_a2col)(1:1)
      b1 = MAPA(jj, map_a1col)(1:1)
      b2 = MAPA(jj, map_a2col)(1:1)
    ENDIF
  ENDIF

  !! If allele frequency calculated, use it
  mafij = NAv
  mafcoij = NAv
  mafcaij = NAv
  IF(ALLOCATED(MAF))   mafij = MAF((/ii,jj/),1)
  IF(ALLOCATED(MAFco)) mafcoij = MAFco((/ii,jj/),1)
  IF(ALLOCATED(MAFca)) mafcaij = MAFca((/ii,jj/),1)

  !! Decide whether to compute interaction effect
  all_pval_NA = ALL((/tmp%AS4p, tmp%AS1p, tmp%DSp, tmp%CSp/) == NAp)  
  
  !! Get estimate of the interaction effect
  IF(out_epi_effect .AND. .NOT.all_pval_NA) THEN
    IF(var_use_beta1) THEN
      beta1 = tmp%beta1
      se_beta1 = tmp%se_beta1
    ELSE
      !! Check if indices are ok
      IF(MAX(ii,jj)>SIZE(X,2)) &
        CALL PrntE("Index out of bound (in WriteOutputRecord)!", Q=.TRUE.)
  
      !! Get genetic data  
      CALL Bed2Ped(X(:,ii), xi, bed_major, ped_minor, c2i(ch1), sex)
      CALL Bed2Ped(X(:,jj), yj, bed_major, ped_minor, c2i(ch2), sex)
      sts1 = sts
      CALL ScanData(SIZE(xi), sts1, xi, yj, nv, ierr, d_co=one, d_ca=one)
      !! Compute the interaction effect
      CALL GetBeta(beta1, sts1(1:nv), xi(1:nv), yj(1:nv), nv, ierr, se_beta1, ana_model)
    ENDIF
  ENDIF

  LL = ""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Report all of the pretests that we done
  IF(do_all_pretests) THEN
    DO icycle=1,12
      SELECT CASE (icycle)
        CASE(1);  S = tmp%PT_co_X; Sp = tmp%PTp_co_X; 
        CASE(2);  S = tmp%AS_co_X; Sp = tmp%ASp_co_X; 
        CASE(3);  S = tmp%DS_co_X; Sp = tmp%DSp_co_X; 
        CASE(4);  S = tmp%PT_co_Y; Sp = tmp%PTp_co_Y; 
        CASE(5);  S = tmp%AS_co_Y; Sp = tmp%ASp_co_Y; 
        CASE(6);  S = tmp%DS_co_Y; Sp = tmp%DSp_co_Y; 
        CASE(7);  S = tmp%PT_co_0; Sp = tmp%PTp_co_0; 
        CASE(8);  S = tmp%AS_co_0; Sp = tmp%ASp_co_0; 
        CASE(9);  S = tmp%DS_co_0; Sp = tmp%DSp_co_0; 
        CASE(10); S = tmp%PT_cc;   Sp = tmp%PTp_cc; 
        CASE(11); S = tmp%AS_cc;   Sp = tmp%ASp_cc; 
        CASE(12); S = tmp%DS_cc;   Sp = tmp%DSp_cc;
      END SELECT 
      IF(Sp==NAp) THEN
        LL = NA//cs//TRIM(LL)
        IF(out_statistic) LL = TRIM(LL)//cs//NA
      ELSE
        WRITE(L1, fmp) Sp
        LL = TRIM(L1)//cs//TRIM(LL)
        IF(out_statistic) THEN
          WRITE(L1, fmv) S
          LL = TRIM(L1)//cs//TRIM(LL)
        ENDIF
      ENDIF
    ENDDO
  ENDIF 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Print pretest generating vector (also for debugging purposes)
  IF(out_S1_vector) THEN
    DO k=2,1,-1
      DO j=9,1,-1
        IF(k==1) WRITE(L1, fmv) tmp%Tn4(j)
        IF(k==2) WRITE(L1, fmv) tmp%Tn1(j)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDDO
    ENDDO
  ENDIF
                      
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print also debugging information
  IF(out_debug_info .OR. out_variance) THEN

    DO k=sDebug,7,-1
      IF(k==9 .AND. .NOT.doS1) CYCLE 
      IF(ANY(k==(/10,11,12,14,15,16/)) .AND. .NOT.doAS4 .AND. .NOT.doAS1) CYCLE 
      IF(ANY(tmp%Debug(k)==(/NAv,NAv2/))) THEN
        LL = NA//cs//TRIM(LL)
      ELSE
        WRITE(L1, fmv) tmp%Debug(k)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDDO
    
  ENDIF

  IF(out_debug_info) THEN

    DO k=MIN(6,sDebug),1,-1
      IF(ANY(tmp%Debug(k)==(/NAv,NAv2/))) THEN
        LL = NA//cs//TRIM(LL)
      ELSE
        WRITE(L1, fmv) tmp%Debug(k)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDDO

    IF(out_epi_effect) THEN
      IF(se_beta1(3)==NAneg) THEN
        LL = NA//cs//TRIM(LL)
      ELSE
        WRITE(L1, fmv) se_beta1(3)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDIF

  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print the actual prevalence
  IF(res_prev>=zero) THEN
    WRITE(L1, '('//fm_lvl//')') res_prev
    LL = TRIM(L1)//cs//TRIM(LL)
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print also allele names
  IF(out_alleles) THEN
    IF(a1==NA) THEN
      LL = NA//cs//TRIM(LL)
    ELSE
      LL = a1//a2//b1//b2//cs//TRIM(LL)
    ENDIF
  ENDIF  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print (automatic) pretest level used for the current pair
  IF(out_level_S1) THEN

    IF(auto_level_report) THEN
      IF(tmp%PLevel(2) == NAneg) THEN
        LL = NA//cs//NA//cs//NA//cs//TRIM(LL)
      ELSE

        IF(tmp%Slope == NAneg) THEN
          LL = NA//cs//TRIM(LL)
        ELSEIF(tmp%Slope == zero) THEN
          LL = "0.0"//cs//TRIM(LL)
        ELSE
          WRITE(L1, fmv) tmp%Slope
          LL = TRIM(L1)//cs//TRIM(LL)
        ENDIF

        IF(tmp%Lambda == NAneg) THEN
          LL = NA//cs//TRIM(LL)
        ELSEIF(tmp%Lambda == zero) THEN
          LL = "0.0"//cs//TRIM(LL)
        ELSE
          WRITE(L1, fmv) tmp%Lambda
          LL = TRIM(L1)//cs//TRIM(LL)
        ENDIF

        IF(tmp%PLevel(2) >= 0.99999_dpp) THEN
          LL = "1.0"//cs//TRIM(LL)
        ELSE
          WRITE(L1, '('//fm_lvl//')') tmp%PLevel(2)
          LL = TRIM(L1)//cs//TRIM(LL)
        ENDIF

      ENDIF
    ENDIF

    IF(tmp%PLevel(1)==one) THEN
      LL = "1.0"//cs//TRIM(LL)
    ELSE
      WRITE(L1, '('//fm_lvl//')') tmp%PLevel(1)
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDIF

  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !! Print minor allele frequencies
  IF(out_maf) THEN
  
    IF(out_maf_coca) THEN
      DO k=2,1,-1
        L1 = NA
        IF(mafcaij(k) /= NAneg) WRITE(L1,'('//TRIM(fmm)//')') mafcaij(k)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDDO
      DO k=2,1,-1
        L1 = NA
        IF(mafcoij(k) /= NAneg) WRITE(L1,'('//TRIM(fmm)//')') mafcoij(k)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDDO
    ENDIF
    
    DO k=2,1,-1
      L1 = NA
      IF(mafij(k) /= NAneg) WRITE(L1,'('//TRIM(fmm)//')') mafij(k)
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDDO
    
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print mapping info (only if map files were specified)
  IF(out_loc_info) THEN
    !!IF(map_nfiles>0) THEN
      LL = TRIM(ch1)//cs//TRIM(rs1)//cs//TRIM(ch2)//cs//TRIM(rs2)//cs//TRIM(LL)
    !!ELSE
    !!  LL = NA//cs//NA//cs//NA//cs//NA//cs//TRIM(LL)
    !!ENDIF
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print the interaction effect
  IF(out_epi_effect) THEN
    IF(se_beta1(3)==NAneg) THEN
      LL = NA//cs//NA//cs//TRIM(LL)
    ELSE
      WRITE(L1, fmv) EXP(beta1(3))
      WRITE(L2, fmv) beta1(3)
      LL = TRIM(L1)//cs//TRIM(L2)//cs//TRIM(LL)
    ENDIF
  ENDIF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Print POST-TEST results
  DO k=7,1,-1
  
    IF(k==1 .AND. .NOT.doAS4) CYCLE
    IF(k==2 .AND. .NOT.doAS1) CYCLE
    IF(k==3 .AND. .NOT.doDS4) CYCLE
    IF(k==4 .AND. .NOT.doDS1) CYCLE
    IF(k==5 .AND. (.NOT.doPO4 .OR. .NOT.doPO4CS .OR. out_minimal)) CYCLE
    IF(k==6 .AND. (.NOT.doPO1 .OR. .NOT.doPO1CS .OR. out_minimal)) CYCLE
    IF(k==7 .AND. .NOT.doCS) CYCLE
    
    !!print *,'k=',k

    IF(k==1) THEN
      T = tmp%AS4; Tp = tmp%AS4p; Err = tmp%ErrAS4; Nca = tmp%ASNca; Nco = tmp%ASNco
      !!print *,'================'
      !!print *,'T=',T
      !!print *,'Tp=',Tp
    ELSEIF(k==2) THEN
      T = tmp%AS1; Tp = tmp%AS1p; Err = tmp%ErrAS1; Nca = tmp%ASNca; Nco = tmp%ASNco
    ELSEIF(k==3) THEN
      T = tmp%DS;  Tp = tmp%DSp;  Err = tmp%ErrDS;  Nca = tmp%DSNca; Nco = tmp%DSNco
    ELSEIF(k==4) THEN
      T = tmp%DS;  Tp = tmp%DSp;  Err = tmp%ErrDS;  Nca = tmp%DSNca; Nco = tmp%DSNco
    ELSEIF(k>=5) THEN
      T = tmp%CS;  Tp = tmp%CSp;  Err = tmp%ErrCS;  Nca = tmp%CSNca; Nco = tmp%CSNco
    ENDIF
    
    !! Process the error codes for possible printing of the global error
    IF(Err<0) THEN
      IF(.NOT.(GErrN==-1 .AND. Err==-1)) GErrN = GErrN + Err
    ELSE
      GErrP = GErrP + Err
    ENDIF 

    !! Print sample size information
    IF(out_ss) THEN
      WRITE(L1, '(I0,A,I0)') Nca,cs,Nco
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDIF

    !! Print error codes
    IF(out_errcode .AND. .NOT.out_errcode_glob) THEN
      WRITE(L1, '(I0)') Err
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDIF

    !! Print effective multiple testing correction
    IF(out_mtc) THEN
      IF(Tp==NAp .OR. pCor(k)==NAp) THEN
        LL = NA//cs//TRIM(LL)
      ELSE
        WRITE(L1, '(I0)') INT(mtc(k), ikb)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDIF
  
    !! Print p-values
    IF(out_pvalcor) THEN
      IF(Tp==NAp .OR. pCor(k)==NAp) THEN
        LL = NA//cs//TRIM(LL)
      ELSEIF(pCor(k)==one) THEN
        LL = "1.0"//cs//TRIM(LL)
      ELSEIF(pCor(k)==zero) THEN
        LL = "0.0"//cs//TRIM(LL)
      ELSE
        WRITE(L1, fmp) pCor(k)
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDIF
  
    !! Print test statistics
    IF(Tp==NAp) THEN
      LL = NA//cs//TRIM(LL)
      IF(out_statistic) LL = NA//cs//TRIM(LL)
    ELSE
      WRITE(L1, fmp) Tp
      LL = TRIM(L1)//cs//TRIM(LL)
      IF(out_statistic) THEN
        WRITE(L1, fmv) T
        LL = TRIM(L1)//cs//TRIM(LL)
      ENDIF
    ENDIF
      
  ENDDO
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print PRETEST results
  DO k=7,0,-1
  
    IF(k==0 .AND. (.NOT.doAS4 .OR. .NOT.doPTAS4)) CYCLE
    IF(k==1 .AND. (.NOT.doAS1 .OR. .NOT.doPTAS1)) CYCLE
    IF(k==2 .AND. (.NOT.doDS4 .OR. .NOT.doPTDS4)) CYCLE
    IF(k==3 .AND. (.NOT.doDS1 .OR. .NOT.doPTDS1)) CYCLE
    IF(k==4 .AND. .NOT.doPO4)  CYCLE
    IF(k==5 .AND. .NOT.doPO1)  CYCLE
    IF(k==6 .AND. .NOT.doPT4) CYCLE
    IF(k==7 .AND. .NOT.doPT1) CYCLE
    
    !! Retrieve results
    IF(k==0) THEN
      T = tmp%PTAS4; Tp = tmp%PTAS4p; Err = tmp%ErrPTAS4; 
      Nca = tmp%PTASNca; Nco = tmp%PTASNco
    ELSEIF(k==1) THEN
      T = tmp%PTAS1; Tp = tmp%PTAS1p; Err = tmp%ErrPTAS1; 
      Nca = tmp%PTASNca; Nco = tmp%PTASNco
    ELSEIF(k==2) THEN
      T = tmp%PTDS4; Tp = tmp%PTDS4p; Err = tmp%ErrPTDS4; 
      Nca = tmp%PTDSNca; Nco = tmp%PTDSNco
    ELSEIF(k==3) THEN
      T = tmp%PTDS1; Tp = tmp%PTDS1p; Err = tmp%ErrPTDS1; 
      Nca = tmp%PTDSNca; Nco = tmp%PTDSNco
    ELSEIF(k==4) THEN
      T = tmp%PTPO4; Tp = tmp%PTPO4p; Err = tmp%ErrPTPO4; 
      Nca = tmp%CSNca; Nco = tmp%CSNco
    ELSEIF(k==5) THEN
      T = tmp%PTPO1; Tp = tmp%PTPO1p; Err = tmp%ErrPTPO1; 
      Nca = tmp%CSNca; Nco = tmp%CSNco
    ELSEIF(k==6) THEN
      T = tmp%PT4; Tp = tmp%PT4p; Err = tmp%ErrPT4; 
      Nca = tmp%PTASNca; Nco = tmp%PTASNco
    ELSEIF(k==7) THEN
      T = tmp%PT1; Tp = tmp%PT1p; Err = tmp%ErrPT1; 
      Nca = tmp%PTASNca; Nco = tmp%PTASNco
    ENDIF
    
    !! Process the error codes for possible printing of the global error
    IF(Err<0) THEN
      IF(.NOT.(GErrN==-1 .AND. Err==-1)) GErrN = GErrN + Err
    ELSE
      GErrP = GErrP + Err
    ENDIF 

    !! Print sample size information
    IF(out_ss) THEN
      IF(Err/=-1) THEN
        WRITE(L1, '(I0,A,I0)') Nca,cs,Nco
      ELSE
        WRITE(L1, '(3A)') NA,cs,NA 
      ENDIF
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDIF

    !! Print error codes
    IF(out_errcode .AND. .NOT.out_errcode_glob) THEN
      WRITE(L1, '(I0)') Err
      LL = TRIM(L1)//cs//TRIM(LL)
    ENDIF

    !! Print p-value and test statistics
    IF(Tp==NAp .OR. Tp<zero) THEN
      LL = NA//cs//TRIM(LL)
      IF(out_statistic) L1 = NA
    ELSEIF(Tp==one) THEN
      LL = "1.0"//cs//TRIM(LL)
      IF(out_statistic) L1 = "0.0"
    ELSEIF(Tp==zero) THEN
      LL = "0.0"//cs//TRIM(LL)
      IF(out_statistic) WRITE(L1, fmv) T
    ELSE
      WRITE(L1, fmp) Tp
      LL = TRIM(L1)//cs//TRIM(LL)
      IF(out_statistic) WRITE(L1, fmv) T
    ENDIF
    
    IF(out_statistic) LL = TRIM(L1)//cs//TRIM(LL)
    
    !!print *,'LL=', LL(1:30)
    !!print *,'NA=', NA

  ENDDO
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IF(out_errcode .AND. out_errcode_glob) THEN
    IF(GErrP>0) THEN
      WRITE(L1, '(I0)') GErrP
    ELSE
      WRITE(L1, '(I0)') GErrN
    ENDIF
    LL = TRIM(L1)//cs//TRIM(LL)
  ENDIF  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Print output if one of the above conditions is met
  IF(out_testid) THEN
    WRITE(L1,'(I0)') n
    LL = TRIM(L1)//cs//TRIM(LL)
  ENDIF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Drop the (possible) trailing column separator
  IF(LL(LEN_TRIM(LL):LEN_TRIM(LL))==cs) LL(LEN_TRIM(LL):LEN_TRIM(LL)) = ""  
  
  !! Write the line into the file
  WRITE(uout,'(A)') TRIM(ADJUSTL(LL))
  
  !! Return record length (+1 for EOL symbol)
  rec_len = LEN_TRIM(LL) + 1

END SUBROUTINE WriteOutputRecord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_DATAIO
