  MODULE EPI_TEST

  !USE INITIALIZE        !! INITIALIZATION OF VARIABLES
  USE EPI_UTILS         !! VARIOUS TOOLS
  USE EPI_TESTPROC      !! TESTING ROUTINES
  USE EPI_REPORT        !! REPORTING ROUTINES
  USE EPI_DATAIO        !! INPUT/OUTPUT ROUTINES
  
  IMPLICIT NONE

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE DoTests(X, sts, sex, n, nco, nca, nsamp, dAS, dDS, ssAS, ssDS, &
                   nrow, MAPA, MAF, InclLoc, model_S1, model_S2)

  IMPLICIT NONE

  !! Input / Output variables
  !INTEGER(iks), INTENT(INOUT), ALLOCATABLE :: X(:,:), sts(:), sex(:), ssAS(:), &
  !                                            ssDS(:)
  INTEGER(iks), INTENT(INOUT)           :: X(:,:), sts(:), sex(:), ssAS(:), &
                                           ssDS(:)
  INTEGER, INTENT(IN)                   :: nrow, n, nco, nca, nsamp, model_S1, &
                                           model_S2
  LOGICAL, INTENT(IN)                   :: InclLoc(:)
  CHARACTER(*), INTENT(INOUT)           :: MAPA(:,:)
  REAL(dpp), INTENT(IN)                 :: dAS, dDS
  REAL(dpp), INTENT(INOUT), ALLOCATABLE :: MAF(:,:)

  !! Local variables
  TYPE(TMPTYPE)   :: tmp(nrow)
  INTEGER(iks)    :: XX((n+3)/4), x_all(n), y_all(n), sts_all(n), sex_all(n), &
                     ssAS_all(n), ssDS_all(n), x_o(n), sts_o(n), sex_o(n), &
                     ssAS_o(n), ssDS_o(n), sts1AS(n), sts2AS(n), sts1DS(n), &
                     sts2DS(n), x1AS(n), y1AS(n), x2AS(n), y2AS(n), x1DS(n), &
                     y1DS(n), x2DS(n), y2DS(n)
  CHARACTER(mml)  :: MAPA1(SIZE(MAPA,1),2)
  CHARACTER(mtsl) :: time1f, time2f, runtime 
  LOGICAL         :: is_same_chr, ASstd, zermarAS, zermarDS, only_CS, signif, &
                     rejected, err_allzero_x, excl_male1, excl_fema1, doAS, &
                     doDS, doPO, now_sex_chr, skip_S1, out_it, any_neg_err, &
                     any_pos_err
  REAL(dpp)       :: time1, time2, Tn4(dmV), Tn1(dmV), PTAS4, PTDS4, PTAS1, &
                     PTDS1, AS4, AS1, DS, CS, PTAS4p, PTDS4p, PTAS1p, PTDS1p, &
                     AS4p, AS1p, DSp, CSp, PT4, PT1, PTPO4, PT4p, PT1p, &
                     PTPO4p, PTPO1, PTPO1p, V4(dmV,dmV), V4_I(dmV,dmV), &
                     V4_SqI(dmV,dmV), V1_I(dmV,dmV), V1_SqI(dmV,dmV), min_p1, &
                     V1(dmV,dmV), Z1(dmV,dmV), mafii, mafjj, mafii_opt, &
                     mafjj_opt, lambda, slope, beta0(0:2), beta1(0:3), min_p2, &
                     se_beta1(0:3), beta0DS(0:2), beta1DS(0:3), optpower, &
                     DebugAS(sDebug), levels(SIZE(MatOpt,2),SIZE(MatOpt,3)), &
                     weights(SIZE(MatOpt,2),SIZE(MatOpt,3))
  INTEGER(ikb)    :: ij, tmp_nl, tmp_cor, ntested, nrec, prcnt, maxntests, &
                     NumTests, k, ntests_cd
  INTEGER         :: ii, iii, jj, xrow, nthreads, nloci, ieAPL, ie0AS, ie0DS, &
                     ie1(4), ie1AS(2), ie1DS(2), ie1a4, ie1a1, ie1PO(2), &
                     ie2(4), ie2AS(2), ie2DS, ieCS, isamp, i_lbnd, i_ubnd, &
                     j_lbnd, j_ubnd, j_len, nv1AS, nv2AS, nv1DS, nv2DS, nv, i, &
                     imaf, jmaf, n_incl_loci, chrii, chrjj, rsii, rsjj, &
                     PTASNco, PTASNca, PTDSNco, PTDSNca, ASNco, ASNca, DSNco, &
                     DSNca, CSNco, CSNca, tt, cntrgrp1
  INTEGER(ikb)    :: posii, posjj  
  REAL(dpp)       :: cc_co1(dmS,dmS), cc_co2(dmS,dmS), cc_ca1(dmS,dmS), &
                     cc_ca2(dmS,dmS), ccAS_co1(dmS,dmS), ccAS_co2(dmS,dmS), &
                     ccAS_ca1(dmS,dmS), ccAS_ca2(dmS,dmS), ccDS_co1(dmS,dmS), &
                     ccDS_co2(dmS,dmS), ccDS_ca1(dmS,dmS), ccDS_ca2(dmS,dmS)
#if _OPENMP
  INTEGER         :: OMP_GET_NUM_THREADS !, OMP_GET_THREAD_NUM
#endif
  CHARACTER(mmtl) :: text, stext
  
  INTEGER         :: icycle, WT11
  REAL(dpp)       :: dASco, dASca, dDSco, dDSca, AS_co_X, ASp_co_X, DS_co_X, &
                     DSp_co_X, PT_co_X, PTp_co_X, AS_co_Y, ASp_co_Y, DS_co_Y, &
                     DSp_co_Y, PT_co_Y, PTp_co_Y, AS_co_0, ASp_co_0, DS_co_0, &
                     DSp_co_0, PT_co_0, PTp_co_0, AS_cc, ASp_cc, DS_cc, &
                     DSp_cc, PT_cc, PTp_cc
 

  !! Initialize counters
  CALL InitCounters()

  !! Print testing options report
  CALL Report(test_opts=.TRUE.)
  
  !! Make sure that there is no problem if MAF is not allocated
  IF(.NOT.ALLOCATED(MAF)) compute_maf = .FALSE.

  !! Get the total number of (included) loci
  nloci = SIZE(InclLoc)
  n_incl_loci = CountTrue(InclLoc)

  !! Count the tests
  IF(sub_make_pairs) THEN
    ntests = SIZE(Included_pairs,1)
    maxntests = ntests
    ntests_cd = nloci * (nloci-1)/2                         
  ELSE
    ntests = INT(n_incl_loci,ikb)*(INT(n_incl_loci,ikb)-1)/2 
    maxntests = nsamp * nloci * (nloci-1)/2
    ntests_cd = maxntests
  ENDIF
  
  !! Increase the number of tests to include all samples
  ntests = nsamp * ntests
  IF(ntests_assumed <= zero) THEN
    IF(MTC(1)>=one) THEN
      ntests_assumed = MTC(1)
    ELSE
      ntests_assumed = maxntests
    ENDIF
  ENDIF

  !! Set the number of threads for parallel computing
#ifdef _OPENMP
  nthreads = OMP_GET_NUM_THREADS()
#else
  nthreads = 1  
#endif

  !! Set the maf variables to NA values
  mafii = NAneg
  mafjj = NAneg
  mafii_opt = NAneg
  mafjj_opt = NAneg
   
  !! Set random generation seed
  !CALL SetSeed(put=.TRUE., seed=seed, nthreads=nthreads)
      
  !print *,"size(seed)=", size(seed)
  !print *,"seed=", seed
  !stop
  
  !! ***************************************************** !!
  !!            DETERMINE OPTIMAL PRETEST LEVEL            !!  
  !! ***************************************************** !!
  IF((auto_level_report .OR. auto_level_use) .AND. auto_level_single .AND. &
  auto_level_optim<zero) THEN

    !! Announce that S1 level is automatically determined (if selected)
    CALL Report(auto_level=auto_level_report .OR. auto_level_use)

    stext = "Determining optimal S1 level ..."
    ij = 0
    
    !! Loop over samples 
    DO isamp=1,nsamp
    
      !! Loop over loci
      DO ii = 1, nloci-1
        
        !! Print countdown
        CALL PrintCD(text, ij, ntests_cd, prcnt, stext=stext, finish=.TRUE., &
                     adv=.FALSE.)     
                            
        !! Set bounds
        j_lbnd = ii + 1
        j_len = nloci - ii
        
        !! Skip this loci if excluded
        IF(.NOT.sub_make_pairs) THEN 
          IF(.NOT.InclLoc(ii)) THEN
            ij = ij + j_len
            CYCLE
          ENDIF
        ENDIF
  
        !! Store the ii-th loci maf
        IF(compute_maf) mafii = MAF(ii,isamp)
  
        !! Loop over loci
        DO jj = ii+1, nloci
  
          !! Increase counter
          ij = ij + 1

          !! Skip this pair if excluded 
          IF(.NOT.sub_make_pairs) THEN
            IF(.NOT.InclLoc(jj)) CYCLE
          ENDIF
          
          !! Store the jj-th loci maf 
          IF(compute_maf) mafjj = MAF(jj,isamp)
          
          !! Fill cc_?? with such values that the sample sizes fit
          cc_co1 = zero
          cc_ca1 = zero
          cc_co2 = zero
          cc_ca2 = zero
          cc_co1(1,1) = dAS * nco
          cc_ca2(1,1) = nca
          cc_co2(1,1) = nco
          IF(auto_level_DS) cc_co2(1,1) = (one-dAS) * nco
          
          ieAPL = 0
          !! Determine optimal S1 level using mafs and other parameters
          CALL GetAutoLevel(auto_level_optim, optpower, lambda, slope, &
                            model_S2, mafii, mafjj, level_S2, auto_level_min, &
                            auto_level_max, ntests_assumed, auto_level_prev, &
                            LOG(auto_level_OR1), LOG(auto_level_OR2), &
                            LOG(auto_level_OR), cc_co1, cc_co2, cc_ca1, &
                            cc_ca2, auto_level_DS, auto_level_merge, ieAPL, &
                            MatOpt)
                  
        ENDDO
        
      ENDDO
      
    ENDDO
    
    !! Get the optimal level by weighing the individual optimal S1 levels
    !! by optimal power
    levels = LOG10(MAX(MIN(ABS(MatOpt(1,:,:)), one), minimum_test_level))
    weights = MAX(MatOpt(2,:,:), zero)

    !! Weight also by frequency of occurence of maf pairs
    IF(auto_level_weight_by_counts) &
      weights = weights * MAX(MatOpt(5,:,:), zero)
    
    !! If positive weights, compute the optimal level as weighted average
    IF(SUM(weights)>zero) &
      auto_level_optim = 10**(SUM(levels*weights)/SUM(weights)) &
                                * auto_level_factor
    
    !! Print countdown finish
    CALL PrintCD(text, ij, ntests_cd, prcnt, finish=.TRUE.)
  
    IF(auto_level_use .AND. auto_level_optim>0) level_S1 = auto_level_optim 

  ENDIF
  
  !! Announce the automatically determined S1 level
  IF(auto_level_report .OR. auto_level_use) &
    CALL Report(auto_level_value=.TRUE.)
    
  !! If auto level should not be used and no user input, set it to default 
  IF(.NOT.auto_level_use .AND. level_S1<zero) level_S1 = def_level1 

  !! ***************************************************** !!
  !!        END OF DETERMINE OPTIMAL PRETEST LEVEL         !!  
  !! ***************************************************** !!
  
  !! Store chr and rs numbers locally and numeralize rs information
  MAPA1 = MAPA(:,(/map_chcol, map_rscol/))
  CALL upcase(MAPA1(:,2))
  CALL DelAllSubstr(MAPA1(:,2), "RS")
  DO i=1,SIZE(MAPA1,1)
    IF(.NOT.IsInteger(MAPA1(i,2))) MAPA1(i,2) = "0"
  ENDDO

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!               PRINT WARNINGS AND ANNOUNCEMENTS               !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Get start time of testing stage
  CALL GetCurrentTime(time1f, time1)

  !! Announce beginning of testing stage
  CALL Report(test_start=.TRUE., starttime=time1f)

  !! Print a warning all loci were marked as excluded and exit
  IF(n_incl_loci == 0) THEN
    out_sort = .FALSE.
    CALL Report(all_excl=.TRUE.)
    RETURN
  ENDIF

  !! Announce that the number of tests is limited by user input
  CALL Report(ntests_limit=ntests_limit)
  
  !! Announce that delta is automatically determined (if selected)
  CALL Report(delta=.TRUE.)
  
  !! Check and warn that the settings make all pairs pass S1
  CALL Report(all_pretests_rejected=reject_all_pretests)
  
  !! Announce the length of a temp cycle
  CALL Report(cycle_length=.TRUE.)
  
  !! Check for too many threads or no parallelism
  CALL Report(check_too_many_threads=.TRUE.)
  
  !! Warn about the setting that only pairs within Included_pairs will be tested 
  CALL Report(test_only_neighbours=sub_make_pairs)

  !! Set deltas so that it can be evaluate whether DS pretest needs to be
  !! reported 
  CALL SetDeltas(WT1, dAS, dDS, cntrgrp, dASco, dASca, dDSco, dDSca)
  IF((dASco /= dDSco .OR. dASca /= dDSca) .AND. doS1) THEN
    IF(doDS4) doPTDS4 = .TRUE.
    IF(doDS1 .OR. doDS4) THEN
      doDS1 = .TRUE.
      doPTDS1 = .TRUE.
    ENDIF
  ENDIF
                           
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                  START THE TESTING PROCESS                   !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Set bounds for do cycles
  IF(sub_make_pairs) THEN
    i_lbnd = 1
    i_ubnd = SIZE(Included_pairs,1)
  ELSE
    i_lbnd = 1
    i_ubnd = nloci - 1
  ENDIF
  
  !! Reset the counters
  prcnt = -1
  tmp_nl  = -1
  ij      = 0
  nrec    = 0
  tmp_cor = 0
  ntested = 0
  n_samechr_excl = 0
  n_samechr_incl = 0

  !! Get the number of lines in X that contain current sample data
  xrow = (n + 3)/4

  !! Nulify variables
  text = ""
  tmp  = NAtmp

  !! Set indicators of which two-stage procedures should be done
  doAS = doAS4 .OR. doAS1
  doDS = doDS4 .OR. doDS1
  doPO = doPO4 .OR. doPO1
          
  !! Announce the start of testing stage
  CALL Prnt("Number of tests to perform: "//i2c(maxntests))
  IF(ntests < maxntests)& 
    CALL Prnt("User set limit on the number of tests: "//i2c(ntests))
  CALL Prnt("Testing in progress ...")

  !! ******************************* !!
  !! START TESTING LOOP OVER SAMPLES !!
  !! ******************************* !!
  
  DO isamp=1,nsamp  !! Announce the start of testing stage

    !! Reread case-control status, sex and S1 selection status info
    IF(SIZE(sts)==n) THEN
      sts_o = sts
      sex_o = sex
      ssAS_o = ssAS
      ssDS_o = ssDS
    ELSE
      sts_o = sts(1+(isamp-1)*n:isamp*n)
      sex_o = sex(1+(isamp-1)*n:isamp*n)
      ssAS_o = ssAS(1+(isamp-1)*n:isamp*n)
      ssDS_o = ssDS(1+(isamp-1)*n:isamp*n)
    ENDIF
    
    !! ************************************** !!
    !!   START FIRST TESTING LOOP OVER LOCI   !!
    !! ************************************** !!
    DO iii = i_lbnd, i_ubnd
    
      !! Print countdown
      CALL PrintCD(text, ij, maxntests, prcnt, testing=.TRUE., show_pct=.FALSE., &
                   adv=.FALSE.)
      
      !! Check if the maximum number of tests has been reached
      IF(ij > maxntests) EXIT
      
      !! Set the bounds for the inter loop and calculate the length of the loop
      !! differently for when only pairs are to be tested
      IF(sub_make_pairs) THEN
        ii = Included_pairs(iii,1)
        j_lbnd = Included_pairs(iii,2)
        j_ubnd = Included_pairs(iii,2)
        j_len = 1
      ELSE
        ii = iii
        j_lbnd = ii + 1
        j_ubnd = nloci
        j_len = j_ubnd - j_lbnd + 1
      ENDIF
      
      !! Skip all remaining tests if the cap on the number of tests was reached
      IF(ntested >= ntests) THEN
        ij = ij + j_len
        CALL Prnt
        CYCLE
      ENDIF

      !! Skip this loci if excluded
      IF(.NOT.sub_make_pairs) THEN 
        IF(.NOT.InclLoc(ii)) THEN
          ij = ij + j_len
          CYCLE
        ENDIF
      ENDIF

      !! Remember ii-th loci's mapping information
      IF(compute_maf) mafii = MAF(ii,isamp)
      chrii = c2i(MAPA1(ii,1))
      rsii = c2i(MAPA1(ii,2)) 

      !! Covert data from a binary format into 0-1-2 format (current loci only)
      XX = X(1+(isamp-1)*xrow:isamp*xrow,ii)
      CALL Bed2Ped(XX, x_o, bed_major, ped_minor, chrii, sex_o)
      
      !! Check for all-zero-values error in input data
      err_allzero_x = ALL(x_o==0)
      
!!***********************************!!
!!       START PARALLEL REGION       !!
!!***********************************!!

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(sex, x_o, sts_o, sex_o, n_incl_loci, ii, iii) &
!$OMP FIRSTPRIVATE(nrow, ssAS_o, ssDS_o, mafii, mafii_opt) &
!$OMP FIRSTPRIVATE(err_allzero_x, bed_major, ped_minor) &
!$OMP FIRSTPRIVATE(isamp, nsamp, i_lbnd, i_ubnd, j_lbnd, j_ubnd) &
!$OMP FIRSTPRIVATE(dAS, dDS, dASco, dASca, dDSco) &
!$OMP FIRSTPRIVATE(dDSca, var_use_beta1, var_indep, silent, chrii) &
!$OMP FIRSTPRIVATE(rsii, SUBMAPARS, replication_mode, rs_valid_chck) &
!$OMP FIRSTPRIVATE(level_S1, auto_level_optim, auto_level_min) &
!$OMP FIRSTPRIVATE(ntests_assumed, auto_level_max, auto_level_report) &
!$OMP FIRSTPRIVATE(auto_level_DS, auto_level_merge, auto_level_use) &
!$OMP FIRSTPRIVATE(auto_level_prev, auto_level_OR1, auto_level_OR2) &
!$OMP FIRSTPRIVATE(auto_level_OR, auto_level_maf, MatOpt) &
!$OMP FIRSTPRIVATE(level_S2, level_S2_nom, regres_divided) &
!$OMP FIRSTPRIVATE(doS1, doS2, test_same_chr_S1, test_same_chr) &
!$OMP FIRSTPRIVATE(no_causal_pair, only_S1, olim1, olim2) &
!$OMP FIRSTPRIVATE(xrow, InclLoc, no_output, sub_make_pairs) &
!$OMP FIRSTPRIVATE(WT1, T1_df, T1_df_co, T1_df_po, WT2, T2_df, cntrgrp) &
!$OMP FIRSTPRIVATE(do_all_pretests, countdown, map_chcol, map_rscol) &
!$OMP FIRSTPRIVATE(var_poststand, var_bound, var_cor_AS, min_marg_cnt) &
!$OMP FIRSTPRIVATE(correct_all_cells, cell_count_cor, min_cell_cnt) &
!$OMP FIRSTPRIVATE(fix_subsamp, nonpar_probs, timestamp, sex_chr_chck) &
!$OMP FIRSTPRIVATE(excl_male, excl_fema, report_all, report_errs, doAS, doAS4) &
!$OMP FIRSTPRIVATE(doAS1, doDS, doDS4, doDS1, doPO, doPO4, doPO1, doCS) &
!$OMP FIRSTPRIVATE(doCS_all, reject_all_pretests, always_pretest) &
!$OMP FIRSTPRIVATE(min_ss, min_co, min_ca, ntest_limit2) &
!$OMP FIRSTPRIVATE(var_grouped, cov_grouped, compute_maf) &
!$OMP PRIVATE(WT11, jj, mafjj, mafjj_opt, imaf, jmaf, lambda, slope) &
!$OMP PRIVATE(chrjj, rsjj, x1AS, y1AS, x2AS, y2AS, sts1AS, sts2AS) &
!$OMP PRIVATE(x1DS, y1DS, x2DS, y2DS, sts1DS, sts2DS, optpower, now_sex_chr) &
!$OMP PRIVATE(sts_all, sex_all, x_all, y_all, XX, posii, posjj, ssAS_all) &
!$OMP PRIVATE(ssDS_all, is_same_chr, cntrgrp1, beta0, beta1, se_beta1) &
!$OMP PRIVATE(beta0DS, beta1DS, k, DebugAS, V4, V4_I, V4_SqI, V1, V1_I) &
!$OMP PRIVATE(V1_SqI, Z1, PTASNco, PTASNca, PTDSNco, PTDSNca, ASNco) &
!$OMP PRIVATE(ASNca, DSNco, DSNca, CSNco, CSNca, Tn4, Tn1) &
!$OMP PRIVATE(PTAS4, PTDS4, PTAS1, PTDS1, PTAS4p, PTAS1p, PTDS4p, PTDS1p) &
!$OMP PRIVATE(PT4, PT1, PTPO4, PTPO1, CS, AS4, AS1, DS, PT4p, PT1p, PTPO4p) &
!$OMP PRIVATE(PTPO1p, CSp, AS4p, AS1p, DSp, ASstd, zermarAS, zermarDS, ie0AS) & 
!$OMP PRIVATE(rejected, signif, ie0DS, ie1, ie1AS, ie1DS, ie1a4, ie1a1, ie1PO) &
!$OMP PRIVATE(ie2, ie2AS, ie2DS, ieCS, ieAPL, only_CS, cc_co1, cc_co2) &
!$OMP PRIVATE(cc_ca1, cc_ca2, ccAS_co1, ccAS_co2, ccAS_ca1, ccAS_ca2) &
!$OMP PRIVATE(ccDS_co1, ccDS_co2, ccDS_ca1, ccDS_ca2, nv, nv1AS, nv2AS) &
!$OMP PRIVATE(nv1DS, nv2DS, min_p1, min_p2, out_it, any_neg_err, any_pos_err) &
!$OMP PRIVATE(excl_male1, excl_fema1, skip_S1, j_len) &
!$OMP PRIVATE(AS_co_X, ASp_co_X, DS_co_X, DSp_co_X, PT_co_X, PTp_co_X) &
!$OMP PRIVATE(AS_co_Y, ASp_co_Y, DS_co_Y, DSp_co_Y, PT_co_Y, PTp_co_Y) &
!$OMP PRIVATE(AS_co_0, ASp_co_0, DS_co_0, DSp_co_0, PT_co_0, PTp_co_0) &
!$OMP PRIVATE(AS_cc, ASp_cc, DS_cc, DSp_cc, PT_cc, PTp_cc) &
!$OMP SHARED(tmp, model_S1, model_S2, X, MAF, MAPA1, ij, nrec) &
!$OMP SHARED(NumS1, NumS1OK, NumS1NotOK, NumS2) &
!$OMP SHARED(NumLD, NumEpi, NumErrs, NumLowCell) &
!$OMP SHARED(NumLowCellNumCor, NumLowSS, NumLowSS1) &
!$OMP SHARED(NumLowSS2, NumLowMarg, varAS_sum, varAS_sum_n) &
!$OMP SHARED(ntested, n_samechr_incl, n_samechr_excl, nthreads) &
!$OMP SHARED(Included_pairs, n)

#ifdef _OPENMP
!$OMP MASTER
      IF(isamp==1 .AND. ii==i_lbnd) nthreads = OMP_GET_NUM_THREADS()
      !print *,"OMP_GET_NUM_THREADS()=", OMP_GET_NUM_THREADS()
!$OMP END MASTER
#endif
      
!$OMP DO SCHEDULE(DYNAMIC) REDUCTION(+:NumS1, NumS1OK, NumS1NotOK) &
!$OMP REDUCTION(+:NumS2, NumErrs, NumLD, NumEpi, NumLowMarg) &
!$OMP REDUCTION(+:NumLowCell, NumLowCellNumCor) &
!$OMP REDUCTION(+:NumLowSS, NumLowSS1, NumLowSS2) & 
!$OMP REDUCTION(+:varAS_sum, varAS_sum_n) &
!$OMP REDUCTION(+:n_samechr_excl, n_samechr_incl, ij, ntested)
! !$OMP REDUCTION(+:nrec) 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! START SECOND TESTING LOOP OVER LOCI !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO jj = j_lbnd, j_ubnd
      
        !! Store locally if pretest is to be done
        skip_S1 = .NOT.doS1
        
        !! Increase total counter of tests
        ij = ij + 1
        
        !! Skip this pair if excluded
        IF(.NOT.sub_make_pairs) THEN
          IF(.NOT.InclLoc(JJ)) CYCLE
        ENDIF

        !! Count the included loci
        ntested = ntested + 1

        !! Remember jj-th loci's mapping information
        mafjj = NAneg
        IF(compute_maf) mafjj = MAF(jj,isamp)
        chrjj = c2i(MAPA1(jj,1))
        rsjj = c2i(MAPA1(jj,2))
        
        !! If replication mode is on, skip this pair if not numbered according
        !! to 1-2, 3-4, ..., (2k-1)-2k, etc. 
        IF(replication_mode .AND. rs_valid_chck .AND. sub_make_pairs) THEN
          posii = MAXLOC(-ABS(SUBMAPARS-rsii),1)
          posjj = MAXLOC(-ABS(SUBMAPARS-rsjj),1)
          IF(posjj - posii /= 1 .OR. IsEven(posii)) CYCLE
        ENDIF
        
        !! Skip to the next cycle if both loci are on the same chromosome
        is_same_chr = chrii == chrjj
        IF(is_same_chr) THEN
          IF(test_same_chr_S1) THEN
            n_samechr_incl = n_samechr_incl + 1
          ELSE
            n_samechr_excl = n_samechr_excl + 1
            skip_S1 = .TRUE.  
            IF(.NOT.test_same_chr .OR. only_S1 .OR. (.NOT.doS2 .AND. .NOT.doCS_all)) THEN
              ntested = ntested - 1
              CYCLE
            ENDIF
          ENDIF
        ENDIF
        
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!
        !! Loop over 4 kinds of pre-tests
        DO icycle=4,1,-1
        
          !! Skip all of the extra pretests if "do all pretests" is false
          IF(icycle>1 .AND. (.NOT.do_all_pretests .OR. skip_S1)) CYCLE
          
        !write(70,'(A)') "B"

          !! Select which test is to be done now
          SELECT CASE (icycle)
          CASE(1)

            IF(skip_S1) THEN
              WT11 = 0
            ELSE
              WT11 = WT1
            ENDIF
            
            cntrgrp1 = cntrgrp
            
            CALL SetDeltas(WT1, dAS, dDS, cntrgrp1, dASco, dASca, dDSco, dDSca)

          CASE(2)
            WT11 = T1co
            dASco = one
            dASca = zero
            dDSco = one
            dDSca = zero
            cntrgrp1 = 2
          CASE(3)
            WT11 = T1co
            dASco = one
            dASca = zero
            dDSco = one
            dDSca = zero
            cntrgrp1 = 0
          CASE(4)
            WT11 = T1cc
            dASco = dAS
            dASca = dAS
            dDSco = dDS
            dDSca = dDS
            cntrgrp1 = 0
          END SELECT
        !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!

          !! Initialize test statistics, p-values and other variables
          CALL SetVal(NAv, PTAS4, PTAS1, PTDS4, PTDS1, PT4, PT1, PTPO4, PTPO1, &
                      CS, AS4, AS1, DS)
          CALL SetVal(NAp, PTAS4p, PTAS1p, PTDS4p, PTDS1p, PT4p, PT1p, PTPO4p, &
                      PTPO1p, CSp, AS4p, AS1p, DSp)
          CALL SetVal(NAv, Tn4, Tn1, beta0, beta1)
          CALL SetVal(NAv2, DebugAS)
          CALL SetVal(.FALSE., signif, rejected, only_CS, zermarAS, zermarDS)
          CALL SetVal(err_not_done, ie0AS, ie0DS, ie2DS, ieCS, ie1a4, ie1a1, &
                      ieAPL)
          CALL SetVal(err_not_done, ie1, ie2, ie1AS, ie1DS, ie2AS, ie1PO)
          CALL SetVal(NAneg, mafii_opt, mafjj_opt, lambda, slope)
          CALL SetVal(0, nv, nv1AS, nv1DS, nv2AS, nv2DS, PTASNco, PTASNca, &
                      PTDSNco, PTDSNca, ASNco, ASNca, DSNco, DSNca, CSNco, &
                      CSNca)
          
          CALL SetVal(NAv, PT_cc, PT_co_x, PT_co_y, PT_co_0, AS_cc, AS_co_x, &
                      AS_co_y, AS_co_0, DS_cc, DS_co_x, DS_co_y, DS_co_0)
          CALL SetVal(NAp, PTp_cc, PTp_co_x, PTp_co_y, PTp_co_0, ASp_cc, &
                      ASp_co_x, ASp_co_y, ASp_co_0, DSp_cc, DSp_co_x, &
                      DSp_co_y, DSp_co_0)
                      
          CALL SetVal(zero, ccAS_co1, ccAS_co2, ccAS_ca1, ccAS_ca2)
          CALL SetVal(zero, ccDS_co1, ccDS_co2, ccDS_ca1, ccDS_ca2)

          se_beta1 = -one
          ASstd = .TRUE.

          !print *,OMP_GET_THREAD_NUM(),"--------------"
          !! Reread status info, sex info, x_all, etc. since these could have 
          !! been modified by ScanData on the previous loop (ii x jj)
          sts_all = sts_o
          sex_all = sex_o
          x_all = x_o
          ssAS_all = ssAS_o
          ssDS_all = ssDS_o
          
        !write(70,'(A)') "C"
          !! Assign values to temporary gene data storage y_all
          XX = X(1+(isamp-1)*xrow:isamp*xrow, jj) 
          CALL Bed2Ped(XX, y_all, bed_major, ped_minor, chrjj, sex_all)
          
          !! Check for all-zero-values error in input data          
          IF(err_allzero_x .OR. ALL(y_all==0)) THEN
            CALL PrntE("Data for current pair contains only zero values.")
            CALL Prnt("NOTE: Skipping the current pair of markers.")
            WRITE(usto,*) "ii=", ii, "jj=", jj, "SUM(y_all)=", SUM(y_all)
            CALL SetVal(err_all_zeros, ie1AS, ie1DS, ie2AS, ie1PO)
            CALL SetVal(err_all_zeros, ie2DS, ieCS, ie1a4, ie1a1)
            GOTO 24
          ENDIF
          
          !! Set gender exclusion status
          excl_male1 = excl_male
          excl_fema1 = excl_fema
          now_sex_chr = .FALSE.
          IF(sex_chr_chck) THEN
            IF(ANY(chrii == sex_chr) .OR. ANY(chrjj == sex_chr)) THEN
              now_sex_chr = .TRUE.
              IF(ANY((/chrii,chrjj/)==sex_chr(1))) excl_male1 = .TRUE.
              IF(ANY((/chrii,chrjj/)==sex_chr(2))) excl_fema1 = .TRUE.
            ENDIF          
          ENDIF
          
          !! Clean up the data in terms of invalid values and NAs and get counts
          !! for AS
          IF(doAS .OR. doPO .OR. doCS) THEN
          
            CALL ScanData(n, sts_all, x_all, y_all, nv, ie0AS, dASco, dASca, &
                          WT11, sts1AS, x1AS, y1AS, nv1AS, sts2AS, x2AS, y2AS, &
                          nv2AS, ccAS_co1, ccAS_co2, ccAS_ca1, ccAS_ca2, &
                          correct_all_cells, cell_count_cor, min_cell_cnt, &
                          min_ca, min_co, ssAS_all, fix_subsamp, sex_all, & 
                          excl_male1, excl_fema1)
                          
            !! Check for any valid data
            IF(nv==0) THEN
              CALL SetVal(err_no_cc, ie1AS, ie1DS, ie1PO)
              CALL SetVal(err_no_cc, ie1a4, ie1a1)
              GOTO 24
            ENDIF
            
          ENDIF
          
          !print *,"D"
          !! Get counts for DS. If delta's are the same for AS and DS the counts
          !! for AS will be used by DS. Otherwise, rescan the data
          IF(doDS) THEN
          
            IF(doAS .AND. dASco==dDSco .AND. dASca==dDSca) THEN
               
               ie0DS = ie0AS
               sts1DS = sts1AS
               sts2DS = sts2AS
               x1DS = x1AS
               x2DS = x2AS
               y1DS = y1AS
               y2DS = y2AS
               nv1DS = nv1AS
               nv2DS = nv2AS
               ccDS_co1 = ccAS_co1
               ccDS_ca1 = ccAS_ca1
               ccDS_co2 = ccAS_co2
               ccDS_ca2 = ccAS_ca2
               
            ELSE
              
              CALL ScanData(n, sts_all, x_all, y_all, nv, ie0DS, dDSco, dDSca, &
                            WT11, sts1DS, x1DS, y1DS, nv1DS, sts2DS, x2DS, y2DS, &
                            nv2DS, ccDS_co1, ccDS_co2, ccDS_ca1, ccDS_ca2, &
                            correct_all_cells, cell_count_cor, min_cell_cnt, &
                            min_ca, min_co, ssDS_all, fix_subsamp, sex_all, &
                            excl_male1, excl_fema1)
                            
            ENDIF
            
          ENDIF
          
          !print *,"D2"
          !! Save samplesizes
          PTASNco = INT(SUM(ccAS_co1))               ! PTAS controls
          PTASNca = INT(SUM(ccAS_ca1))               ! PTAS cases
          PTDSNco = INT(SUM(ccDS_co1))               ! PTDS controls
          PTDSNca = INT(SUM(ccDS_ca1))               ! PTDS cases
          ASNco = INT(SUM(ccAS_co1+ccAS_co2))    ! AS4&AS1 controls
          ASNca = INT(SUM(ccAS_ca1+ccAS_ca2))    ! AS4&AS1 cases
          DSNco = INT(SUM(ccDS_co2))                 ! DS4&DS1 controls
          DSNca = INT(SUM(ccDS_ca2))                 ! DS4&DS1 cases
          CSNco = INT(SUM(ccAS_co1+ccAS_co2))    ! full sample controls (CS, PTPO4, PTPO1)
          CSNca = INT(SUM(ccAS_ca1+ccAS_ca2))    ! full sample cases (CS, PTPO4, PTPO1)
            
          !print *,"D2"
        !write(70,'(A)') "E"
          !! If counts corrected register it and remove error indicator
          IF(icycle==1 .AND. ANY((/ie0AS,ie0DS/) == err_low_cell_cor)) THEN
            NumLowCellNumCor = NumLowCellNumCor + 1
            ie0AS = 0
            ie0DS = 0
          ENDIF
  
          !! If counts too low (or other error occured) then do no perform testing
          IF(icycle==1) THEN
            DO tt=1,2
              IF(tt==1 .AND. .NOT.doAS .AND. .NOT.doPO) CYCLE
              IF(tt==2 .AND. .NOT.doDS) CYCLE
              IF((tt==1 .AND. ie0AS/=0) .OR. (tt==2 .AND. ie0DS/=0)) THEN
                IF(tt==1) ie1AS = ie0AS
                IF(tt==2) ie1DS = ie0DS
                IF(ANY(ie0AS==errs_low_ss) .OR. ANY(ie0DS==errs_low_ss)) THEN
                  NumLowSS = NumLowSS + 1
                ELSE
                  IF(ANY((/ie0AS,ie0DS/)==err_low_cell)) &
                    NumLowCell = NumLowCell + 1
                  IF(.NOT.skip_S1) NumS1NotOK = NumS1NotOK + 1
                ENDIF
                IF(.NOT.(doAS .OR. doDS .OR. doPO .OR. doCS)) GOTO 23
              ENDIF
            ENDDO
          ENDIF
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!    DETERMINE OPTIMAL PRETEST LEVEL     !!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
          !! Get the optimal S1 level if it has not been determined already
          IF(auto_level_report .AND. auto_level_optim < zero) THEN
            
            !! If pretest is based on LD in controls, perform power maximization
            IF(WT11 == T1co) THEN
              IF(.NOT.auto_level_maf) THEN
                mafii_opt = mafii 
                mafjj_opt = mafjj 
              ENDIF
              
              !! Determine optimal pretest level using mafs and other parameters
              IF(auto_level_DS .AND. doDS) THEN
                cc_co1 = ccDS_co1
                cc_ca1 = ccDS_ca1
                cc_co2 = ccDS_co2
                cc_ca2 = ccDS_ca2
              ELSEIF(doAS) THEN
                cc_co1 = ccAS_co1
                cc_ca1 = ccAS_ca1
                cc_co2 = ccAS_co2
                cc_ca2 = ccAS_ca2
              ELSE
                CALL PrntE("Automatic level supported only for AS/DS.")
              ENDIF
              CALL GetAutoLevel(auto_level_optim, optpower, lambda, slope, &
                                model_S2, mafii_opt, mafjj_opt, level_S2, &
                                auto_level_min, auto_level_max, ntests_assumed, &
                                auto_level_prev, LOG(auto_level_OR1), &
                                LOG(auto_level_OR2), LOG(auto_level_OR), &
                                cc_co1, cc_co2, cc_ca1, cc_ca2, auto_level_DS, &
                                auto_level_merge, ieAPL, MatOpt)
            !! Otherwise return what "may be a reasonable value"
            ELSE
              auto_level_optim = SQRT(level_S2 / ntests_assumed)
            ENDIF
            
          ELSEIF(compute_maf) THEN
          
            !! Save the mafs, lambda and slope to be printed in output
            mafii_opt = mafii 
            mafjj_opt = mafjj
            imaf = MAX(1,INT(mafii*10**opt_ndig_mafs))  
            jmaf = MAX(1,INT(mafjj*10**opt_ndig_mafs))
            lambda = MatOpt(3,imaf,jmaf)
            slope = MatOpt(4,imaf,jmaf)
              
          ENDIF
  
          !! Set the actual pretest level to the optimal value
          IF(auto_level_use) level_S1 = auto_level_optim 
  
          !! If pretest level implies no pretest, just skip it
          IF(level_S1 >= one .AND. .NOT.always_pretest) skip_S1 = .TRUE.
  
          !! Check for negative S1 level (we are being VERY cautious)
          IF(level_S1 < zero) level_S1 = def_level1
  
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!    CALL PRETEST SUBROUTINES     !!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
        !write(70,'(A)') "E2"
          !print *,OMP_GET_THREAD_NUM(),"a"
          !! Do S1 only if level_S1 is smaller than 1
          IF(doS1) THEN
          
            !! Perform selected S1 and/or S2 test(s)
            !print *,"A"
            IF(doAS .OR. doPO) &
              CALL TestS1(ccAS_co1, ccAS_ca1, ccAS_co2, ccAS_ca2, WT1=WT11, &
                          T4=PTAS4, T4p=PTAS4p, ierrT4=ie1AS(1), T4df=T1_df, &
                          T1=PTAS1, T1p=PTAS1p, ierrT1=ie1AS(2), &
                          doPT4=.TRUE., doPT1=.TRUE., doPO4=doPO4, &
                          doPO1=doPO1, doAS4=doAS4, doAS1=doAS1, &
                          doDS=doDS, Tsp4=PT4, Tsp4p=PT4p, &
                          ierrsp4=ie1a4, Tsp4df=T1_df_co, Tsp1=PT1, &
                          Tsp1p=PT1p, ierrsp1=ie1a1, Tpo4=PTPO4, Tpo4pval=PTPO4p, &
                          ierrpo4=ie1PO(1), Tpo4df=T1_df_po, Tpo1=PTPO1, &
                          Tpo1pval=PTPO1p, ierrpo1=ie1PO(2), Tn4=Tn4, V4=V4, &
                          V4_I=V4_I, V4_SqI=V4_SqI, Tn1=Tn1, V1=V1, V1_I=V1_I, &
                          V1_SqI=V1_SqI, Z1=Z1, WT2=WT2, min_margin=min_marg_cnt, &
                          zero_margs=zermarAS, cntrgrp=cntrgrp1, &
                          var_grouped=var_grouped, var_indep=var_indep, &
                          divided=regres_divided, model=model_S1)
            !print *,"B"

            !print *,OMP_GET_THREAD_NUM(),"b"
            IF(doDS) THEN
            
              IF(.FALSE. .AND. XOR(doAS4, doDS4) .AND. XOR(doAS1, doDS1) .AND. dASco==dDSco &
               .AND. dASca==dDSca) THEN
              
                ie1DS(1:2) = ie1AS(1:2)
                PTDS4 = PTAS4
                PTDS4p = PTAS4p
                PTDS1 = PTAS1
                PTDS1p = PTAS1p
                
              ELSE
              
  if(.false.) then
              
  print *,""              
  print *,"cco=", ccDS_co1(1,:)+ccDS_co2(1,:)
  print *,"cco=", ccDS_co1(2,:)+ccDS_co2(2,:)
  print *,"cco=", ccDS_co1(3,:)+ccDS_co2(3,:)
  print *,""              
  print *,"ccDS_co1=", ccDS_co1(1,:)
  print *,"ccDS_co1=", ccDS_co1(2,:)
  print *,"ccDS_co1=", ccDS_co1(3,:)
  print *,""              
  print *,"ccDS_co2=", ccDS_co2(1,:)
  print *,"ccDS_co2=", ccDS_co2(2,:)
  print *,"ccDS_co2=", ccDS_co2(3,:)
  print *,""              
  print *,"cca=", ccDS_ca1(1,:)+ccDS_ca2(1,:)
  print *,"cca=", ccDS_ca1(2,:)+ccDS_ca2(2,:)
  print *,"cca=", ccDS_ca1(3,:)+ccDS_ca2(3,:)
  print *,""              
  print *,"ccDS_ca1=", ccDS_ca1(1,:)
  print *,"ccDS_ca1=", ccDS_ca1(2,:)
  print *,"ccDS_ca1=", ccDS_ca1(3,:)
  print *,""              
  print *,"ccDS_ca2=", ccDS_ca2(1,:)
  print *,"ccDS_ca2=", ccDS_ca2(2,:)
  print *,"ccDS_ca2=", ccDS_ca2(3,:)
  print *,""
  
  endif              

            !print *,OMP_GET_THREAD_NUM(),"b1"
              !print *,"Aa"
                CALL TestS1(ccDS_co1, ccDS_ca1, ccDS_co2, ccDS_ca2, WT11, &
                            T4=PTDS4, T4p=PTDS4p, ierrT4=ie1DS(1), &
                            T4df=T1_df, T1=PTDS1, T1p=PTDS1p, &
                            ierrT1=ie1DS(2), doDS=.TRUE., WT2=WT2, &
                            min_margin=min_marg_cnt, zero_margs=zermarDS, &
                            cntrgrp=cntrgrp1, var_grouped=var_grouped, &
                            var_indep=var_indep, model=model_S1)
              !print *,"Bb"

                !print *,"B PTDS1=", PTDS1
                !print *,"B PTDS1p=", PTDS1p
                                
            !print *,OMP_GET_THREAD_NUM(),'b2'
              ENDIF
              
            ENDIF
            
            !! Erase results that are not to be reported
            IF(.NOT.doAS4) THEN 
              PTAS4 = NAv
              PTAS4p = NAp
            ENDIF
            IF(.NOT.doAS1) THEN
              PTAS1 = NAv
              PTAS1p = NAp
            ENDIF
            IF(.NOT.doDS4) THEN
              PTDS4 = NAv
              PTDS4p = NAp
            ENDIF
            IF(.NOT.doDS1) THEN
              PTDS1 = NAv
              PTDS1p = NAp
            ENDIF
            
          ENDIF
          
          !! Process the results of pre-testing
          IF(.NOT.skip_S1) THEN
          
            !! Increase stage 1 test counter
            IF(icycle==1) NumS1 = NumS1 + 1
    
            !! Error occurred
            ie1 = (/ie1AS, ie1DS/)
            IF(ANY(ie1 > 0)) THEN
            
              !! Increase the erroneous S1 test counter
              IF(icycle==1) NumS1NotOK = NumS1NotOK + 1
              
              !! Check for low sample size errors
              IF(ANY(ie1AS(1)==errs_low_ss) .OR. ANY(ie1DS(1)==errs_low_ss)) THEN
                IF(icycle==1) THEN 
                  NumLowSS = NumLowSS + 1
                  !! Sample size for stage 1 too small (prioritized error)
                  IF(ANY(ie1 == err_low_ss_S1)) THEN
                    NumLowSS1 = NumLowSS1 + 1
                  !! Sample size for stage 2 too small
                  ELSEIF(ANY(ie1 == err_low_ss_S2)) THEN
                    NumLowSS2 = NumLowSS2 + 1
                  ENDIF
                ENDIF
                GOTO 24
              ENDIF
  
              !! Count low marginal count errors
              IF(icycle==1 .AND. ANY(ie1 == err_low_margin)) &
                NumLowMarg = NumLowMarg + 1
  
            ELSE
              
              !! Increase S1 test counter of OK tests if no error occured
              IF(icycle==1) NumS1OK = NumS1OK + 1
              
              !! Increase the counter of detected LD if pretest significant
              IF(doAS4 .AND. PTAS4p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(1) = NumLD(1) + 1
              ENDIF
              IF(doAS1 .AND. PTAS1p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(2) = NumLD(2) + 1
              ENDIF
              IF(doDS4 .AND. PTDS4p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(3) = NumLD(3) + 1
              ENDIF
              IF(doDS1 .AND. PTDS1p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(4) = NumLD(4) + 1
              ENDIF
              IF(doPO4 .AND. PTPO4p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(5) = NumLD(5) + 1
              ENDIF
              IF(doPO1 .AND. PTPO1p <= level_S1) THEN
                rejected = .TRUE.
                IF(icycle==1) NumLD(6) = NumLD(6) + 1
              ENDIF
  
            ENDIF
  
          ENDIF
          
          !! IF THIS POINT REACHED, NO ERROR HAS OCCURRED
          !! Always indicate rejection in stage 1 when no pretest wanted 
          IF(skip_S1 .OR. reject_all_pretests) rejected = .TRUE.
          
          !! If S2 not wanted, skip it
          IF(only_S1) GOTO 25
          IF(.NOT.doS2 .AND. .NOT.doCS_all .AND. skip_S1) GOTO 24
          IF(.NOT.doS2 .AND. .NOT.doCS_all) GOTO 22
  
          !! If pretest didn't reject but score is still wanted by the user,
          !! indicate that it should be computed
          IF((.NOT.rejected .OR. .NOT.doS2) .AND. doCS_all) &
            only_CS = .TRUE.
            
          !! Call the score test if stage 1 rejected or when only score wanted
          IF(.NOT.rejected) THEN
            IF(only_CS) GOTO 21
            IF(.NOT.only_CS) GOTO 22
          ENDIF
          
          !! ****************************************************************** !!
          !! **                          STAGE  2                            ** !!
          !! ****************************************************************** !!
          !! **    IF LD DETECTED OR NO PRETEST PERFORMED DO SECOND TEST     ** !! 
          !! ****************************************************************** !!

          !! First check if there are not too many stage 2 tests already
          IF(nthreads==1 .AND. ntest_limit2 > 0 .AND. NumS2 > ntest_limit2) GOTO 24
  
          !! Increase the stage 2 test counters
          IF(icycle==1) NumS2 = NumS2 + 1
  
          21 CONTINUE
          
          !! ****************************************************************** !!
          !! Compute stage 2 statistics (or only score when only score wanted)
          !! ****************************************************************** !!
          IF(doAS .OR. doPO .OR. doCS) &
            CALL TestS2(sts=sts_all(1:nv), x=x_all(1:nv), y=y_all(1:nv), & 
                        nv=nv, n=nv, model=model_S2, CS=CS, CSp=CSp, &
                        ierrCS=ieCS, use_beta1=var_use_beta1, &
                        beta0=beta0, beta1=beta1, se_beta1=se_beta1, &
                        doAS4=doAS, doAS1=doAS, Debug=DebugAS, &
                        AS4=AS4, AS4p=AS4p, ierrAS4=ie2AS(1), AS1=AS1, AS1p=AS1p, & 
                        ierrAS1=ie2AS(2), min_ss=min_ss, &
                        var_bound=var_bound, var_cor_AS=var_cor_AS, &
                        ccco1=ccAS_co1, ccco2=ccAS_co2, &
                        ccca1=ccAS_ca1, ccca2=ccAS_ca2, Tn4=Tn4, &
                        V4=V4, V4_I=V4_I, V4_SqI=V4_SqI, Tn1=Tn1, &
                        V1=V1, V1_I=V1_I, V1_SqI=V1_SqI, Z1=Z1, &
                        varAS_sum=varAS_sum, varAS_sum_n=varAS_sum_n, &
                        WT1=WT11, cntrgrp=cntrgrp1, &
                        cov_grouped=cov_grouped, cov_divided=regres_divided, &
                        poststand=var_poststand, stand=ASstd, &
                        nonpar_probs=nonpar_probs)
                          
          !! Erase results that should not be reported
          IF(.NOT.doAS4) THEN
            AS4 = NAv
            AS4p = NAp
          ENDIF
          IF(.NOT.doAS1) THEN
            AS1 = NAv
            AS1p = NAp
          ENDIF
          IF(.NOT.doCS) THEN
            CS = NAv
            CSp = NAp
          ENDIF
          
          !! Get DS
          !print *,"C"
          IF(doDS) &
            CALL TestS2(sts2DS(1:nv2DS), x2DS(1:nv2DS), y2DS(1:nv2DS), nv2DS, &
                        nv2DS, model_S2, DS, DSp, ie2DS, var_use_beta1, &
                        beta0DS, beta1DS)
                            
          !! **************************************************************** !!
          !! CONTINUATION POINT IF ERROR OCCURRED => PARTS OF STAGE 2 SKIPPED !!
          !! **************************************************************** !!
          22 CONTINUE
  
          IF(icycle==1) THEN

            !! Determine rejection (if an error occured, pvalues are NAp)
            IF(MIN(AS4p,AS1p) <= level_S2) THEN
              NumEpi(1) = NumEpi(1) + 1
              NumEpi(2) = NumEpi(2) + 1
            ENDIF
            IF(DSp <= level_S2_nom) NumEpi(3:4) = NumEpi(3:4) + 1
            IF(CSp <= level_S2_nom) NumEpi(5:7) = NumEpi(5:7) + 1
  
          ENDIF
          
          23 CONTINUE
          
          !! Check for MAFs too big or small (too big has priority)
          IF(ANY((/ie1AS, ie1DS/)<=0) .AND. compute_maf) THEN
            IF(MIN(mafii, mafjj) < hundreth) THEN
              ie1AS = err_maf_sml
              ie1DS = err_maf_sml
            ENDIF
            IF(MAX(mafii, mafjj) >= half) THEN
              ie1AS = err_maf_big
              ie1DS = err_maf_big
            ENDIF
          ENDIF
  
          24 CONTINUE
          
          !! Check whether to incease counter of errors
          ie1 = (/ie1AS, ie1DS/)
          ie2 = (/ie2AS, ie2DS, ieCS/)
          IF(ANY(ie1>0) .OR. ANY(ie2>0)) NumErrs = NumErrs + 1
          
          25 CONTINUE
          
          !! ******************************************************** !!
          !! TEMPORARY CODE TO OUTPUT MULTIPLE TESTS ON THE SAME DATA !!
          !! ******************************************************** !!
          
        !! ****************************************************************** !!
          IF(icycle==1) THEN
            PT_co_X = PTAS4
            PTp_co_X = PTAS4p
            AS_co_X = AS4
            ASp_co_X = AS4p
            DS_co_X = DS
            DSp_co_X = DSp
          ELSEIF(icycle==2) THEN
            PT_co_Y = PTAS4
            PTp_co_Y = PTAS4p
            AS_co_Y = AS4
            ASp_co_Y = AS4p
            DS_co_Y = DS
            DSp_co_Y = DSp
          ELSEIF(icycle==3) THEN
            PT_co_0 = PTAS4
            PTp_co_0 = PTAS4p
            AS_co_0 = AS4
            ASp_co_0 = AS4p
            DS_co_0 = DS
            DSp_co_0 = DSp
          ELSEIF(icycle==4) THEN
            PT_cc = PTAS4
            PTp_cc = PTAS4p
            AS_cc = AS4
            ASp_cc = AS4p
            DS_cc = DS
            DSp_cc = DSp
          ENDIF
  
        !! End of do loop over icycle
        ENDDO
        !! ****************************************************************** !!

        !! Don't save results if no output file wanted
        IF(no_output) CYCLE

        !! Don't save results if no output file wanted
        IF(.NOT.report_all .AND. ALL((/ie1,ie2/)==err_not_done)) CYCLE
        
        !! Save the current results into the temporary result file only if they  
        !! have a chance of getting reported in the output file
        min_p1 = MIN(PTAS4p,PTAS1p,PTDS4p,PTDS1p,PTPO4p,PTPO1p)
        min_p2 = MIN(AS4p,AS1p,DSp,CSp)
        any_neg_err = ANY((/ie1AS, ie1AS, ie2AS, ie1DS, ie2DS, ieCS, ie1a4, ie1PO/)<=0)
        any_pos_err = ANY((/ie1AS, ie1AS, ie2AS, ie1DS, ie2DS, ieCS, ie1a4, ie1PO/)>0)
        out_it = .FALSE.
        IF(report_all .AND. any_neg_err)       out_it = .TRUE.
        IF(report_errs .AND. any_pos_err)      out_it = .TRUE.
        IF(min_p1 < olim1 .OR. min_p2 < olim2) out_it = .TRUE.

        IF(.NOT.out_it) CYCLE

        !write(70,'(A)') "G"

!$OMP CRITICAL
        nrec = nrec + 1
        k = nrec
!$OMP END CRITICAL

        !print *,"report_all=", report_all
        !print *,"report_errs=", report_errs
        !print *, "olim1=", olim1, "olim2=", olim2, "min_p1=", min_p1, "min_p2=", min_p2
        !print *,"ii=", ii, "jj=", jj, "k=", k
        
        !! Check if all ok with where a record will be stored
        IF(k > nrow) & 
          CALL PrntE("ARRAY BOUNDS EXCEEDED! (k="//i2cp(k)//", nrow="//&
                       i2cp(nrow)//")", Q=.TRUE.)

        IF(.NOT.tmp(k)%Empty) & 
          CALL PrntE("SOMETHING WENT WRONG AND I AM ABOUT TO OVEWRITE A"//&
                     " AN EXISTING RECORD! (k="//i2cp(k)//", nrow="//i2cp(nrow)//")", Q=.TRUE.)

        if(.false.) then
          print *,"----------------"
          print *,"ii=", ii, "jj=", jj, "PTAS1p=", PTAS1p, "AS4p=", AS4p, &
                "AS1p=", AS1p, "PTDS4p=", PTDS4p, "PTDS1p=", PTDS1p, "DSp=", DSp 
        endif

        !! Store the obtained results into shared variables
        tmp(k)%Empty    = .FALSE.
        
        tmp(k)%PTAS4    = PTAS4
        tmp(k)%PTAS4p   = PTAS4p
        tmp(k)%errPTAS4 = INT(ie1AS(1), 2)
        tmp(k)%PTAS1    = PTAS1
        tmp(k)%PTAS1p   = PTAS1p
        tmp(k)%errPTAS1 = INT(ie1AS(2), 2)
        tmp(k)%PTASNca  = PTASNca
        tmp(k)%PTASNco  = PTASNco
        
        tmp(k)%PTDS4    = PTDS4
        tmp(k)%PTDS4p   = PTDS4p
        tmp(k)%errPTDS4 = INT(ie1DS(1), 2)
        tmp(k)%PTDS1    = PTDS1
        tmp(k)%PTDS1p   = PTDS1p
        tmp(k)%errPTDS1 = INT(ie1DS(2), 2)
        tmp(k)%PTDSNca  = PTDSNca
        tmp(k)%PTDSNco  = PTDSNco
        
        tmp(k)%AS4      = AS4
        tmp(k)%AS4p     = AS4p
        tmp(k)%errAS4   = INT(ie2AS(1), 2)
        tmp(k)%AS1      = AS1
        tmp(k)%AS1p     = AS1p
        tmp(k)%errAS1   = INT(ie2AS(2), 2)
        tmp(k)%ASNca    = ASNca
        tmp(k)%ASNco    = ASNco
        
        tmp(k)%DS       = DS
        tmp(k)%DSp      = DSp
        tmp(k)%errDS    = INT(ie2DS, 2)
        tmp(k)%DSNca    = DSNca
        tmp(k)%DSNco    = DSNco
        
        tmp(k)%CS       = CS
        tmp(k)%CSp      = CSp
        tmp(k)%errCS    = INT(ieCS, 2)
        tmp(k)%CSNca    = CSNca
        tmp(k)%CSNco    = CSNco
        
        tmp(k)%PTPO4    = PTPO4
        tmp(k)%PTPO4p   = PTPO4p
        tmp(k)%errPTPO4 = INT(ie1PO(1), 2)
        tmp(k)%PTPO1    = PTPO1
        tmp(k)%PTPO1p   = PTPO1p
        tmp(k)%errPTPO1 = INT(ie1PO(2), 2)
        
        tmp(k)%beta0    = beta0
        tmp(k)%beta1    = beta1
        tmp(k)%se_beta1 = se_beta1
        tmp(k)%Tn4      = Tn4
        tmp(k)%Tn1      = Tn1
        tmp(k)%Pos      = (/ii, jj/) 
        tmp(k)%Debug    = DebugAS
        tmp(k)%PLevel   = (/level_S1, auto_level_optim/)
        tmp(k)%Lambda   = lambda
        tmp(k)%Slope    = slope
        tmp(k)%ASstd    = ASstd

        tmp(k)%PT4      = PT4
        tmp(k)%PT4p     = PT4p
        tmp(k)%errPT4   = INT(ie1a4, 2)
        tmp(k)%PT1      = PT1
        tmp(k)%PT1p     = PT1p
        tmp(k)%errPT1   = INT(ie1a1, 2)
        
        tmp(k)%AS_co_X  = AS_co_X
        tmp(k)%ASp_co_X = ASp_co_X
        tmp(k)%DS_co_X  = DS_co_X
        tmp(k)%DSp_co_X = DSp_co_X
        tmp(k)%PT_co_X  = PT_co_X
        tmp(k)%PTp_co_X = PTp_co_X
        tmp(k)%AS_co_Y  = AS_co_Y
        tmp(k)%ASp_co_Y = ASp_co_Y
        tmp(k)%DS_co_Y  = DS_co_Y
        tmp(k)%DSp_co_Y = DSp_co_Y
        tmp(k)%PT_co_Y  = PT_co_Y
        tmp(k)%PTp_co_Y = PTp_co_Y
        tmp(k)%AS_co_0  = AS_co_0
        tmp(k)%ASp_co_0 = ASp_co_0
        tmp(k)%DS_co_0  = DS_co_0
        tmp(k)%DSp_co_0 = DSp_co_0
        tmp(k)%PT_co_0  = PT_co_0
        tmp(k)%PTp_co_0 = PTp_co_0
        tmp(k)%AS_cc  = AS_cc
        tmp(k)%ASp_cc = ASp_cc
        tmp(k)%DS_cc  = DS_cc
        tmp(k)%DSp_cc = DSp_cc
        tmp(k)%PT_cc  = PT_cc
        tmp(k)%PTp_cc = PTp_cc

        !! Warn about zero marginal frequencies
        IF(ALL(ie1AS==0) .AND. zermarAS) THEN
          tmp(k)%ErrAS4 = err_zero_margin
          tmp(k)%ErrAS1 = err_zero_margin
          tmp(k)%ErrCS  = err_zero_margin
        ENDIF
        IF(ALL(ie1DS==0) .AND. zermarDS) &
          tmp(k)%ErrDS  = err_zero_margin

        !write(70,'(A)') "H"
      ENDDO ! jj

!$OMP END DO
!$OMP END PARALLEL
  
!!***********************!!
!!  END PARALLEL REGION  !!
!!***********************!! ,

      !! Skip saving results if no output file wanted
      IF(no_output) CYCLE

      !! If nrec is still small, which means there is enough space 
      !! left in the vectors, don't write anything yet
      !IF(nrow - nrec >= (nthreads+1)*n_incl_loci) CYCLE
      IF(nrow - nrec >= n_incl_loci) CYCLE

      !! Otherwise write to temporary output file if nrec is close to nrow
      k = nrec
      IF(k < 1) CYCLE
      CALL SaveTempFile(tmp(1:k), tmp_nl, tmp_cor, last=.FALSE.)
      
      !! Nulify result vectors and nrec
      nrec = 0
      tmp = NAtmp

    ENDDO   ! ii

    !! Check if the maximum number of tests has been reached
    IF(ij > maxntests) EXIT
      
  ENDDO     ! sample
  !! END TESTING LOOP

  !! Report the end of countdown
  CALL PrintCD(text, maxntests, maxntests, prcnt, NumErrs, testing=.TRUE., &
               show_pct=.FALSE.)

  !! Compute Total Running Time of the Computations
  CALL GetCurrentTime(time2f, time2)
  runtime = GetTimeDifference(time1, time2)
  
  !! Announce time of end of testing stage
  CALL Report(test_end=.TRUE., stoptime=time2f, runtime=runtime, nthreads=nthreads)
  
  !! If no tests performed and the reason is that all loci are on the same
  !! chromosome, print the warning
  IF(n_samechr_excl==ntested) CALL Report(all_on_same_chr=.TRUE.)
  
  !! Write to the output file what is left
  k = nrec
  IF(nrec > 0 .AND. .NOT.no_output) &
    CALL SaveTempFile(tmp(1:k), tmp_nl, tmp_cor, last=.TRUE.)

  !! Get multiple testing corrections
  NumTests = MAX(NumS1, NumS2) 

  !! No user value for MTC(3)
  DO i=1,6
    IF(MTC(i)<=zero) MTC(i) = NumLD(i)
  ENDDO
  IF(MTC(7)<=zero) MTC(7) = NumTests
  
  !! User value for correctS between 0 and 1
  DO i=1,SIZE(MTC)
    IF(MTC(i)>zero .AND. MTC(i)<one) MTC(i) = MTC(i) * NumTests
  ENDDO
  
  !! User value for MTC(1) between 0 and 1
  IF(MTC(1)>zero .AND. MTC(1)<one) MTC(1:6) = MTC(1:6) * NumTests
  
  !! Alter the multiple testing correction simultaneously for all tests 
  !! (typically make it larger by a factor 'MTC_factor')
  IF(MTC_factor/=one) MTC = MTC * MTC_factor
  
  !! Write the output file with corrected the p-values
  IF(.NOT.no_output) CALL WriteOutputFile(tmp_file, append=out_append)
          
  !! Print final report
  CALL Report(test_report=.TRUE., runtime=runtime)
  
  RETURN

END SUBROUTINE DoTests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_TEST
