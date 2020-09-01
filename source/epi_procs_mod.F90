MODULE EPI_REPORT

  USE EPI_MATH          !! MATHEMATICAL TOOLS
  USE EPI_INIT          !! DEFINITION AND INITIALIZATION OF VARIABLES

  IMPLICIT NONE
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Report(exec_start, exec_end, cmdline_call, used_args, input_opts, &
              data_summary, simul_input, test_only_neighbours, test_opts, &
              test_start, test_end, test_report, no_pretests, &
              all_pretests_rejected, cycle_length, no_auto_level, auto_level, &
              auto_level_value, delta, ntests_limit, all_excl, &
              all_but_one_excl, all_on_same_chr, partial_fam, &
              centering_problem, starttime, stoptime, runtime, &
              check_too_many_threads, check_no_parallelism, memory_hint, &
              nthreads)
!! Prints information on screen and into log file, uses global variables 
  IMPLICIT NONE
  INTEGER, PARAMETER                  :: lh = sect_width, rws = 55
  LOGICAL, INTENT(IN), OPTIONAL       :: exec_start, exec_end, &
                                         cmdline_call, used_args, &
                                         input_opts, data_summary, simul_input, &
                                         test_opts, test_start, test_end, &
                                         test_report, cycle_length, &
                                         ntests_limit, test_only_neighbours, &
                                         no_auto_level, auto_level, &
                                         auto_level_value, delta, &
                                         no_pretests, all_pretests_rejected, &
                                         partial_fam, all_but_one_excl, &
                                         all_excl, all_on_same_chr, &
                                         centering_problem, memory_hint, &
                                         check_too_many_threads, &
                                         check_no_parallelism
  CHARACTER(*), INTENT(IN), OPTIONAL  :: starttime, stoptime, runtime
  INTEGER, INTENT(IN), OPTIONAL       :: nthreads
  
  INTEGER(ikb)                        :: NLociPair, NErrOther 
  INTEGER                             :: i, j, k, A_N1, D_N1, A_N2, D_N2, &
                                         width, skip(mntl,2), rw, n_incl_loci, &
                                         n_excl_loci, n_tests, min_ntests 
  CHARACTER(mfl)                      :: arg
  CHARACTER(mntl)                     :: text, text1
  CHARACTER(mmtl), ALLOCATABLE        :: AT(:,:)
  CHARACTER(mltl)                     :: cmdline
  CHARACTER(lh)                       :: sstars
  CHARACTER(3)                        :: typ
  CHARACTER                           :: sep(2)
  
  INTEGER(ikb)                        :: NS2, NS2PostStand, &
                                         NS2OK, NS2NotOK, &
                                         NTest, NErrVar, NLowVar, &
                                         NNegVar, NLinDep, NS2Sex
  REAL(dpp)                           :: Correct, d, AvgSS

  DO i=1,LEN(sstars); sstars(i:i) = "*"; ENDDO

  !! Store add_commas in a shorter variable
  rw = report_width

  !! **************************************************************** !!
  !!             AUTOMATIC PRETEST LEVEL CALCULATION                  !! 
  !! **************************************************************** !!
  IF(PRESENT(no_auto_level)) THEN
    IF(no_auto_level) THEN
      CALL Prnt("No auto level calculation will be performed because"//&
                " of low number of tests", skip1=1, skip2=1)
    ENDIF
  ENDIF
   
  IF(PRESENT(auto_level)) THEN 
    IF(auto_level) THEN 
      IF(report_headers) THEN
        !CALL Prnt0(sstars, skip1=1)
        !CALL Prnt(text='AUTOMATIC S1 LEVEL CALCULATION', text1='>>>',text2='<<<', l=lh)
        !CALL Prnt0(sstars, skip2=1)
        CALL Prnt('>>> AUTOMATIC S1 LEVEL CALCULATION <<<', skip1=1, skip2=1)
      ENDIF

      !IF(WT1==T1co) THEN
        CALL Prnt("Assumptions for automatic S1 level calculation:", lead=3)
        IF(auto_level_DS) THEN
          CALL Prnt0("   - DS used in S2")
        ELSE
          CALL Prnt0("   - AS used in S2")
        ENDIF
        CALL Prnt0("   - all excess controls used for S1")
        CALL Prnt0("   - second step performed on the selected"//&
                        " nominal level that is Bonferroni corrected for"//&
                        " the expected number of tests (under full null of"//&
                        " no interactions and no background LD) in the"//&
                        " second phase", lead=5)
        IF(auto_level_maf) THEN
          CALL Prnt0("   - MAFs for each marker determined from the"//&
                          " control sample used in S1", lead=5)
        ELSE
          CALL Prnt0("   - MAFs for each marker determined from all"//&
                          " individuals in the entire sample", lead=5)
        ENDIF
        CALL Prnt0("   - interactions follow model '"//TRIM(GetModelName(ana_model))//&
                        "' with the following values of parameters", lead=4)
        CALL Prnt0("     - population prevalence      : "//r2c(auto_level_prev))
        CALL Prnt0("     - main effect 1 (odds ratio) : "//r2c(auto_level_OR1))
        CALL Prnt0("     - main effect 2 (odds ratio) : "//r2c(auto_level_OR2))
        CALL Prnt0("     - interaction   (odds ratio) : "//r2c(auto_level_OR))
        CALL Prnt0("     - assumed number of tests    : "//r2c(ntests_assumed))
        IF(auto_level_factor /= one) &
          CALL Prnt0("     - multiplication factor      : "//r2c(auto_level_factor))
        CALL Prnt0("")
      !ENDIF
      
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(auto_level_value)) THEN 
    IF(auto_level_value) THEN 
      CALL Prnt("Automatically determined S1 level is "//&
                TRIM(r2c(auto_level_optim/auto_level_factor)))
      IF(auto_level_factor /= one) &
        CALL Prnt("Due to the user given multiplication factor the"//&
                  " determined S1 level was changed to "//&
                  TRIM(r2c(auto_level_optim)))
    ENDIF 
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(delta)) THEN 
    IF(pss_nfiles>0) THEN
      CALL Prnt("AS&DS: S1 subsample determined from PSS file with ratio "//&
                TRIM(r2c(dAS)))
    ELSEIF(user_delta) THEN
      IF(doAS4 .OR. doAS1) &
        CALL Prnt("AS: S1 subsample determined randomly using the"//&
                  " user-given ratio "//TRIM(r2c(dAS)))
      IF(doDS4 .OR. doDS1) &
        CALL Prnt("DS: S1 subsample determined randomly using the"//&
                  " user-given ratio "//TRIM(r2c(dDS)))
    ELSEIF(auto_delta) THEN 
      CALL Prnt("AS&DS: S1 subsample determined randomly using the"//&
                " ratio "//TRIM(r2c(dAS)))
    ELSEIF(auto_d_incapable) THEN
      IF(doAS4 .OR. doAS1) &
        CALL Prnt("AS: S1 sample portion cannot be determined automatically,"//&
                  " using default ratio "//TRIM(r2c(dAS)))
      IF(doDS4 .OR. doDS1) &
        CALL Prnt("DS: S1 sample portion cannot be determined automatically,"//&
                  " using default ratio "//TRIM(r2c(dDS)))
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(ntests_limit)) THEN 
    IF(ntests_limit) & 
      CALL Prnt("Number of tests limited by user to "//i2c(ntest_limit)//" (approx.)")
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(cycle_length)) THEN 
    IF(cycle_length) THEN
      IF(temp_cycle > ped_nsamp * ntests) & 
        CALL Prnt("Results saved (approx.) every "//i2c(temp_cycle)//" tests")
    ENDIF
  ENDIF
  
  !! *********************************************************************** !!

  IF(PRESENT(test_only_neighbours)) THEN 
    IF(test_only_neighbours) THEN
      CALL PrntW("Only pairs on the same line in the SUBMAP file will be tested")
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(all_excl)) THEN 
    IF(all_excl) & 
      CALL PrntW("ALL LOCI WERE EXCLUDED FROM TESTING! NO TESTS WILL"//&
              " BE PERFORMED AND NO OUTPUT FILE PRODUCED!", skip1=1, &
              skip2=1)
  ENDIF
  
  !! *********************************************************************** !!

  IF(PRESENT(all_but_one_excl)) THEN 
    IF(all_but_one_excl) & 
      CALL PrntW("ALL EXCEPT ONE LOCI WERE EXCLUDED FROM TESTING!"//&
                 " NO TESTS CAN BE PERFORMED!", skip1=1, skip2=1)
  ENDIF
  
  !! *********************************************************************** !!

  IF(PRESENT(all_on_same_chr)) THEN 
    IF(all_on_same_chr) &
      CALL PrntW("All markers are located on the same chromosome, which"//&
                 " means that no S1 tests were performed. You can use flag"//&
                 " --pretest-same-chr to change this behavior", skip1=1, skip2=1)
      !CALL PrntW("ALL LOCI ARE LOCATED ON THE SAME CHROMOSOME! NO"//&
      !        " TESTS WERE PERFORMED AND NO OUTPUT FILE WAS PRODUCED!", &
      !        skip1=1, skip2=1)
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(no_pretests)) THEN 
    IF(no_pretests) & 
      CALL PrntW("S1 p-value threshold is set to 1.0, which makes all tested"//&
                 " pairs pass S1")
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(all_pretests_rejected)) THEN 
    IF(all_pretests_rejected) & 
      CALL Prnt("Flag --reject-all sets the value of the S1 level"//&
                " level to 1.0, which causes all S1 tests to be"//&
                " rejected irrespective of their p-values")
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(partial_fam)) THEN 
    IF(partial_fam) & 
      CALL PrntW("Only the first sample pedigree information saved into the"//&
                 " FAM file. For use of the data in another program the FAM"//&
                 " file will probably need to contain data for all samples")
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(centering_problem)) THEN 
    IF(centering_problem) & 
      CALL PrntW("Sample size used for centering in S1 is very low")
  ENDIF
  
  !! *********************************************************************** !!

#ifdef _OPENMP
  !! *********************************************************************** !!

  IF(PRESENT(memory_hint).AND.show_segfault_hint.AND.NUMBER_OF_THREADS>1) THEN
    IF(SIZE(X) > mem_hint_limit) &
      CALL Prnt("Hint: Genotype data arrays are large which may cause a"//&
             " segmentation fault (segfault) during OPENMP thread creation."//&
             " If this occurs try increasing the stack size. Under linux"//&
             " use 'ulimit -s' and/or 'export OMP_STACKSIZE=...' and/or"//&
             " 'export KMP_STACKSIZE=...'. Under Windows try using Visual"//&
             " Studio's 'editbin /STACK:<size> "//program_name//".exe. If"//&
             " the current call executes without crashing during OPENMP"//&
             " thread creation you can ignore this message.", skip1=1, skip2=1)
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(check_too_many_threads)) THEN 
    IF(check_too_many_threads .AND. NUMBER_OF_THREADS > NUMBER_OF_PROCESSORS) & 
      CALL PrntW("The selected number of threads ("//&
              i2c(NUMBER_OF_THREADS)//") to be used is higher than the"//&
              " number of available processors ("//i2c(NUMBER_OF_PROCESSORS)//&
              "), which may decrease the overall performance. Moreover, if"//&
              " the current analysis demands high amound of memory the high"//&
              " number of threads may cause segmentation faults.")
  ENDIF
  
  !! *********************************************************************** !!

  IF(PRESENT(check_no_parallelism)) THEN
    IF(NUMBER_OF_THREADS > 1 .AND. check_no_parallelism) THEN
      
      !! Check for too few tests relative to min_ntests_parallel and the 
      !! selected number of threads 
      min_ntests = MAX(min_ntests_parallel, NUMBER_OF_THREADS)
      IF(ntests < min_ntests) THEN
        NUMBER_OF_THREADS = 1
        CALL PrntW("The number of threads was changed to 1.", skip1=1)
        CALL Prnt("Reason: Parallelization is not needed because the total"//&
                  " number of tests is very low (less than "//&
                  i2c(min_ntests)//")")
      ENDIF
      
      !! Check for multi-sample optimized code
      IF(OMP_MULTISAMPLEOPTIMIZED) THEN
        NUMBER_OF_THREADS = 1
        CALL PrntW("The number of threads was changed to 1", skip1=1)
        CALL Prnt("Reason: The source code was compiled as multi-"//&
                  "sample, optimized which means that using multiple threads"//&
                  " does not improve the performance. If you want to improve"//&
                  " the performance by using multiple threads recompile the"//&
                  " code as single-sample optimized (set"//&
                  " OMP_MULTISAMPLEOPTIMIZED=.FALSE.)")
      ENDIF
    
    ENDIF
  ENDIF

  !! *********************************************************************** !!

#else

  !! Dummy code to avoid warning about unknown arguments when no openmp
  min_ntests = min_ntests_parallel
  IF(PRESENT(memory_hint)) THEN; ENDIF
  IF(PRESENT(check_too_many_threads)) THEN; ENDIF
  IF(PRESENT(check_no_parallelism) .AND. PRESENT(nthreads)) THEN; ENDIF

#endif

  !! *********************************************************************** !!

  IF(PRESENT(exec_start) .AND. PRESENT(starttime)) THEN 
    IF(exec_start) & 
      CALL Prnt("Execution started on "//TRIM(starttime(14:24))//" at "//starttime(1:12))
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(cmdline_call)) THEN 
    IF(cmdline_call) THEN
      cmdline = "" 
      DO k=1,COMMAND_ARGUMENT_COUNT() 
        CALL GET_COMMAND_ARGUMENT(k, arg)
        cmdline = TRIM(cmdline)//" "//TRIM(arg)
      ENDDO
      CALL Prnt("Call arguments:"//TRIM(cmdline), skip1=1)
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(used_args)) THEN
    IF(used_args) THEN    

      !! Announce unknown arguments first if any
      IF(nunknown_args>0) THEN
        CALL PrntW("THERE WERE SOME UNKNOWN ARGUMENTS!", skip1=1)
        CALL Prnt("UNKNOWN ARGUMENTS :", skip1=1, skip2=1)
        DO k=1,SIZE(args_used)
          IF(args_used(k)) CYCLE
          CALL Prnt0("  '"//TRIM(args_all(k))//"'", lead=4)
        ENDDO
      ENDIF

      !! Announce used arguments
      CALL Prnt("Effective arguments :", skip1=1, skip2=1)
      DO k=1,SIZE(args_selected,1)
        arg = "  "//TRIM(args_selected(k,1))//" "//TRIM(args_selected(k,2))
        CALL Prnt0(arg, lead=4)
      ENDDO

      CALL Prnt0("")
      
      !! Stop if any arguments not recognized and stop_on_unknown_args is true 
      IF(nunknown_args>0 .AND. stop_on_unknown_args) &
        CALL PrntE("UNKNOWN ARGUMENTS PRESENT!", Q=.TRUE.)
      
    ENDIF
  ENDIF
      
  !! *********************************************************************** !!

  IF(PRESENT(input_opts)) THEN
    IF(input_opts) THEN 

      !IF(report_headers) THEN
      !  CALL Prnt0(sstars, skip1=1)
      !  CALL Prnt0(text='DATA INPUT AND OUTPUT', text1='>>>', text2='<<<', l=lh)
      !  CALL Prnt0(sstars, skip2=1)
      !ENDIF
      !CALL Prnt0("")
      CALL Prnt0('>>> DATA <<<', skip1=1, skip2=0)
      
      IF(.NOT.simulate_data) THEN

        CALL Prnt("Input files:")
  
        !! Print filenames - PED
        IF(ped_nfiles>0) THEN
          DO k=1,SIZE(ped_file)
            IF(k==1) THEN
              CALL Prnt("PED <- ["//TRIM(ped_file(1))//"]", lead=8)
            ELSE
              CALL Prnt("       ["//TRIM(ped_file(k))//"]", lead=8)
            ENDIF
          ENDDO
        ENDIF
        
        !! Print filenames - BED
        IF(bed_nfiles>0) &
          CALL Prnt("BED <- ["//TRIM(bed_file(1))//"]", lead=8)
        
        !! Print filenames - MAP / BIM
        IF(map_nfiles>0) THEN
          IF(input_format==1) &
            CALL Prnt("BIM <- ["//TRIM(map_file(1))//"]", lead=8)
          IF(input_format==2) &
            CALL Prnt("MAP <- ["//TRIM(map_file(1))//"]", lead=8)
          DO k=2,SIZE(map_file)
            CALL Prnt0("       ["//TRIM(map_file(k))//"]", lead=8)
          ENDDO
        ELSE
          CALL Prnt("MAP <- No mapping file specified")
        ENDIF
  
        !! Print filenames - FAM (pedigree info)
        IF(fam_nfiles>0) &
          CALL Prnt("FAM <- ["//TRIM(fam_file)//"]", lead=8)
  
        !! Print filenames - SUBMAP (subset mapping file)
        IF(sub_nfiles>0) THEN
            CALL Prnt("SUB <- ["//TRIM(sub_file(1))//"]", lead=8)
          DO k=2,SIZE(sub_file)
            CALL Prnt0("       ["//TRIM(sub_file(k))//"]", lead=8)
          ENDDO
        ENDIF
  
        !! Print filenames - PSS (pretest selection status)
        IF(pss_nfiles>0) &
          CALL Prnt("PSS <- ["//TRIM(pss_file)//"]", lead=8)
  
        !CALL Prnt0("")
      ENDIF
      
      CALL Prnt("I/O files :")

      CALL Prnt("OUT -> ["//TRIM(out_file)//"]", lead=8)
      CALL Prnt("LOG -> ["//TRIM(log_file)//"]", lead=8)
      IF(do_out_maf) &
        CALL Prnt("MAF -> ["//TRIM(out_maf_file)//"]", lead=8)
      IF(do_out_pss) &
        CALL Prnt("PSS -> ["//TRIM(out_pss_file)//"]", lead=8)
      
      IF(save_input_data) THEN
        IF(save_binary) THEN
          CALL Prnt("BED <- ["//TRIM(save_ped_file)//"]", lead=8)
          CALL Prnt("BIM <- ["//TRIM(save_map_file)//"]", lead=8)
          CALL Prnt("FAM <- ["//TRIM(save_fam_file)//"]", lead=8)
        ELSE
          CALL Prnt("PED <- ["//TRIM(save_ped_file)//"]", lead=8)
          CALL Prnt("MAP <- ["//TRIM(save_map_file)//"]", lead=8)
        ENDIF
      ENDIF
      
      !CALL Prnt0("")

      !! Print file content info: PED FILES
      DO k=1,4
        IF(k==1) THEN
          IF(ped_nfiles==0) CYCLE
          typ = "PED"; sep = ped_cs
        ELSEIF(k==2) THEN
          IF(map_nfiles==0) CYCLE
          typ = "MAP"; sep = map_cs
        ELSEIF(k==3) THEN
          IF(sub_nfiles==0) CYCLE
          typ = "SUB"; sep = sub_cs
        ELSEIF(k==4) THEN
          IF(fam_nfiles==0) CYCLE
          typ = "FAM"; sep = fam_cs
        ENDIF
        
        text = typ//" file separator:"
        IF(ALL(sep == (/tab, space/))) THEN
          text = TRIM(text)//" [tab/space]"
        ELSEIF(sep(1) == tab) THEN
          text = TRIM(text)//" [tab]"
        ELSEIF(sep(1) == space) THEN
          text = TRIM(text)//" [space]"
        ELSE
          text = TRIM(text)//" ["//sep(1)
          IF(sep(1)/=sep(2)) text = TRIM(text)//"/"//sep(2)
          text = TRIM(text)//"]"
        ENDIF
        CALL Prnt(text, lead=3)
        
      ENDDO

      !! SEX CHROMOSOME CHECK
      IF(map_nfiles>0) THEN
        text = "Excluded chromosomes (if present):"
        IF(ANY(excluded_chr(special_chr))) THEN
          DO i=0,SIZE(excluded_chr)-1
            IF(excluded_chr(i)) text = TRIM(text)//" "//TRIM(chr_name(i))
          ENDDO
        ELSE
          text = TRIM(text)//" NONE"
        ENDIF
        CALL Prnt(text, lead=3)
  
        IF(ANY(.NOT.excluded_chr(sex_chr))) THEN
          IF(sex_chr_chck) THEN
            text = "Sex check active (males excluded for chromosome"//&
                   " X, females for chromosome Y)" 
          !ELSE
          !  text = " * Sex based exclusion for sex chromosomes is disabled"
          ENDIF
          CALL Prnt(text, lead=3)
        ENDIF
        
      ELSEIF(.NOT.simulate_data) THEN
        text = "No check for 'sex/special' chromosomes can be performed"//&
               " because no mapping information present" 
        CALL Prnt(text, lead=3)
      ENDIF

      !! SUBMAP FILE USED TO INCLUDE/EXCLUDE LOCI
      IF(sub_nfiles>0) THEN
        IF(sub_include) THEN
          CALL Prnt("Loci present in the submap file will be included")
        ELSE
          CALL Prnt("Loci present in the submap file will be removed")
        ENDIF
      ENDIF

      !! ASSUMING MALES OR FEMALES WHEN SEX MISSING
      IF(assume_male) &
        CALL Prnt("Assuming 'male' for individuals with missing sex"//&
                        " information (if any missing)", lead=3)
      IF(assume_fema) &
        CALL Prnt("Assuming 'female' for individuals with missing"//&
                        " sex information (if any missing)", lead=3)

      !! EXCLUSION BASED ON SEX
      IF(excl_male) &
        CALL Prnt("All males will be excluded", lead=3)
      IF(excl_fema) &
        CALL Prnt("All females will be excluded", lead=3)

      IF(excl_male .AND. excl_fema) &
        CALL PrntW("Both males and females were excluded from testing!")
      
      !! MINOR ALLELE FREQUENCY CHECKS
      IF(min_MAF>zero) THEN
        text = "Markers with MAF below "//TRIM(r2c(min_MAF))//" will be excluded"
        CALL Prnt(text, lead=3)
      ENDIF

      !CALL Prnt0("")
      
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(simul_input)) THEN
    IF(simul_input .AND. simulate_data) THEN 
      IF(reuse_data) THEN
        CALL Prnt("Using previously simulated data ...")
      ELSE
        CALL Prnt("INPUT DATA SIMULATION SETTINGS :", skip1=1, skip2=1)
        CALL Prnt0("  Total sample size          :   "//i2c(ped_ss))
        CALL Prnt0("  Number of cases            :   "//i2c(ped_nca))
        CALL Prnt0("  Number of controls         :   "//i2c(ped_nco))
        CALL Prnt0("  Number of samples          :   "//i2c(ped_nsamp))
        CALL Prnt0("  Total number of markers    :   "//i2c(ped_nloci))
        CALL Prnt0("  Causal markers             :   "//i2c(ncausalloci))
        CALL Prnt0("  Neutral markers            :   "//i2c(nneutralloci))
        CALL Prnt0("  Pairs of markers in LD     :   "//i2c(nLDpairs))
        CALL Prnt0("  Amount of LD               :   "//r2c(LD))
        CALL Prnt0("  * Causal pair simulation settings *")
        IF(allelefreqA == allelefreqB) THEN
          IF(fixed_MAF) THEN
            CALL Prnt0("  MAF (fixed value)          :   "//r2c(allelefreqA))
          ELSE
            CALL Prnt0("  MAF (lower bound)          :   "//r2c(allelefreqA))
          ENDIF
        ELSE
          IF(fixed_MAF) THEN
            CALL Prnt0("  MAF 1 (fixed value)        :   "//r2c(allelefreqA))
            CALL Prnt0("  MAF 2 (fixed value)        :   "//r2c(allelefreqB))
          ELSE
            CALL Prnt0("  MAF 1 (lower bound)        :   "//r2c(allelefreqA))
            CALL Prnt0("  MAF 2 (lower bound)        :   "//r2c(allelefreqB))
          ENDIF
        ENDIF
        CALL Prnt0("  Interaction model          :   "//TRIM(GetModelName(sim_model)))
        CALL Prnt0("  Population prevalence      :   "//r2c(prevalence))
        CALL Prnt0("  Main effect 1 (OR1)        :   "//r2c(OR1))
        CALL Prnt0("  Main effect 2 (OR2)        :   "//r2c(OR2))
        CALL Prnt0("  Interaction effect (OR)    :   "//r2c(OR))
        CALL Prnt0("")
      ENDIF
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(data_summary)) THEN
    IF(data_summary .AND. .NOT.no_data_summary) THEN 

      !CALL Prnt0(sstars, skip1=1)
      !CALL Prnt0(text='INPUT DATA SUMMARY', text1='>>>', text2='<<<', l=lh)
      !CALL Prnt0(sstars, skip2=1)
      CALL Prnt0('>>> INPUT DATA SUMMARY <<<', skip1=1, skip2=0)
      !CALL PrntF("INPUT DATA SUMMARY : ", skip1=1, skip2=1)

      CALL Prnt0(text1="Total sample size              :   ", &
                 text2=i2c(ped_nsamp*ped_ss), l=rws)
      CALL Prnt0(text1=" - valid status                :   ", &
                 text2=i2c(ped_nco+ped_nca), l=rws)
      CALL Prnt0(text1="   - controls                  :   ", &
                 text2=i2c(ped_nco), l=rws)
      CALL Prnt0(text1="   - cases                     :   ", &
                 text2=i2c(ped_nca), l=rws)
      CALL Prnt0(text1="   - males                     :   ", &
                 text2=i2c(ped_nmale), l=rws)
      CALL Prnt0(text1="   - females                   :   ", &
                 text2=i2c(ped_nfema), l=rws)
      CALL Prnt0(text1=" - missing status              :   ", &
                 text2=i2c(ped_nsamp*ped_ss-ped_nco-ped_nca), l=rws)
      CALL Prnt0(text1=" - missing sex                 :   ", &
                 text2=i2c(ped_nsamp*ped_ss-ped_nmale-ped_nfema), l=rws)
      
      NLociPair = INT(ped_nloci,ikb)*(INT(ped_nloci,ikb)-1)/2
      
      IF(ped_nsamp==1) THEN
        CALL Prnt0(text1="Total number of markers        :   ", &
                   text2=i2c(ped_nloci), l=rws)
        CALL Prnt0(text1="Total number of markers pairs  :   ", &
                   text2=i2c(NLociPair), l=rws)
      ELSE
        CALL Prnt0(text1="Number of samples              :   ", &                              
                   text2=i2c(ped_nsamp), l=rws)
        CALL Prnt0(text1="Number of markers per sample   :   ", &                              
                   text2=i2c(ped_nloci), l=rws)
        CALL Prnt0(text1="Number of pairs per sample     :   ", &
                   text2=i2c(NLociPair), l=rws)
      ENDIF
      
      !! Total number of included and excluded loci
      n_incl_loci = CountTrue(InclLoc)
      n_excl_loci = ped_nloci - n_incl_loci
      n_tests = n_incl_loci * (n_incl_loci-1) / 2  
  
      CALL Prnt0(text1="Number of included markers     :   ", &                              
                 text2=i2c(n_incl_loci), l=rws) 
      CALL Prnt0(text1="Number of excluded markers     :   ", &                              
                 text2=i2c(n_excl_loci), l=rws)
      IF(n_excl_loci>0) THEN
        CALL Prnt0(text1=" - excl. by user (chromosome)  :   ", &                              
                   text2=i2c(n_chr_excl), l=rws)
        CALL Prnt0(text1=" - excl. by user (negative bp) :   ", &                              
                   text2=i2c(n_negbp_excl), l=rws)
        CALL Prnt0(text1=" - excl. by user (submap file) :   ", &                              
                   text2=i2c(n_sub_excl), l=rws)
        CALL Prnt0(text1=" - excl. for low MAF           :   ", &                              
                   text2=i2c(n_maf_excl), l=rws)
        CALL Prnt0(text1=" - excl. for invalid RS number :   ", &                              
                   text2=i2c(n_rs_excl), l=rws)
                        
        IF(n_maxntests_excl>0) &
          CALL Prnt0(text1=" - excl. by --maxntests        :   ", &                              
                     text2=i2c(n_maxntests_excl), l=rws)
      ENDIF
      CALL Prnt0(text1="Number of included pairs       :   ", &
                 text2=i2c(ped_nsamp*n_tests), l=rws)

    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(test_opts)) THEN
    IF(test_opts) THEN    

      IF(report_headers) THEN
        !CALL Prnt0(sstars, skip1=1)
        !CALL Prnt0(text='TESTING OPTIONS', text1='>>>', text2='<<<', l=lh)
        !CALL Prnt0(sstars, skip2=1)
        CALL Prnt0('>>> ANALYSIS <<<', skip1=1, skip2=0)
      ENDIF

      !! **************************************************************** !!
      !!                        STAGE 1 SETTINGS                          !! 
      !! **************************************************************** !!
      IF(no_testing) THEN
        CALL Prnt("No tests will be performed", skip2=1)
      ELSE
      
        IF(.NOT.doS1) THEN
          CALL Prnt("No S1 tests will be performed")     
          !CALL Prnt("STAGE 1 : No S1 tests will be performed")     
        ELSE

          !CALL Prnt("STAGE 1 : ")

          !text = ""
          !IF(WT1 == T1cc) text = "Chisquare test (cases+controls)"
          !IF(WT1 == T1co) text = "Chisquare test (controls)"
          !IF(WT1 == T1ca) text = "Chisquare test (cases)"
          !IF(WT1 == T1po) text = "Chisquare test (pooled)"
          !IF(WT1 == T1sc) text = "Score test"
          !CALL Prnt("STAGE 1 : "//text, skip2=1)
          
          if(.false.) then
          text = ""
          IF(WT1 == T1cc) THEN
            IF(dAS==one) &
              text = " * AS : Tests of difference in dependence between cases"//&
                     " and controls using all individuals in the two samples"
            IF(dAS<one) &
              text = " * AS : Tests of difference in dependence between cases"//&
                     " and controls using only portion of both samples"
            IF(dDS==one) &
              text = " * DS : Tests of difference in dependence between cases"//&
                     " and controls using all individuals in the two samples"
            IF(dDS<one) &
              text = " * DS : Tests of difference in dependence between cases"//&
                     " and controls using only portion of both samples"
          ELSEIF(ANY(WT1 == (/T1co,T1ca/))) THEN
            IF(WT1 == T1co) text1 = "controls"
            IF(WT1 == T1ca) text1 = "cases"
            IF(dAS==one) &
              text = " * AS : Tests of independence using full sample of "//TRIM(text1)
            IF(dAS<one) &
              text = " * AS : Tests of independence using a subset of "//TRIM(text1)
            CALL Prnt0(text, lead=3)
            IF(dDS==one) &
              text = " * DS : Tests of independence using full sample of "//TRIM(text1)
            IF(dDS<one) &
              text = " * DS : Tests of independence using a subset of "//TRIM(text1)
            CALL Prnt0(text, lead=3)
          ELSEIF(WT1 == T1po) THEN
            text = " * AS&DS: Tests of independence in the pooled sample of"//&
                   " cases and controls using all individuals in the two"//&
                   " samples"
          ELSEIF(WT1 == T1sc) THEN
            text = " * AS&DS: Score tests for interactions between markers using only"//&
                   " subset of the sample of controls in S1 (analysis model '"//&
                   TRIM(GetModelName(ana_model))//"')"
            CALL Prnt0(text, lead=3)
          ENDIF
          endif
  
  
          text = "S1 tests will be performed using"
          IF(WT1==T1cc) text = TRIM(text)//" both controls and cases"
          IF(WT1==T1co) text = TRIM(text)//" only controls"
          IF(WT1==T1ca) text = TRIM(text)//" only cases"
          IF(WT1==T1cc) text = TRIM(text)//" both controls and cases (pooled sample)"
          CALL Prnt(text)
           
          !text = " * Controlled error rate       :  "//TRIM(ErrorRate1)
          !CALL Prnt0(text)
  
          text = "S1 p-value threshold is"
          IF(auto_level_use) THEN
            text = TRIM(text)//" determined automatically"
          ELSE
            text = TRIM(text)//" "//r2c(level_S1_nom)
          ENDIF
          CALL Prnt(text)
  
          text = "S1 subsample is determined"
          IF(pss_nfiles>0) THEN
            text = TRIM(text)//" by PSS file"
          ELSEIF(fix_subsamp) THEN
            text = TRIM(text)//" randomly once and fixed for all S1 tests"
          ELSE
            text = TRIM(text)//" randomly for each S1 test"
          ENDIF
          CALL Prnt(text)
  
          DO i=1,2
            IF(i==1) THEN
              IF(.NOT.doAS4 .AND. .NOT.doAS1) CYCLE
              text = "AS"
              d = dAS
            ELSE 
              IF(.NOT.doDS4 .AND. .NOT.doDS1) CYCLE
              text = "DS"
              d = dDS
            ENDIF
          
            text = "Sample ratio in S1 for "//text(1:2)//" tests is"
            IF(pss_nfiles>0) THEN
              text = TRIM(text)//" determined from PSS file (--delta ignored)"
            ELSE
              IF(auto_delta) THEN
                text = TRIM(text)//" determined automatically"
              ELSEIF(d>=zero .AND. d<=one) THEN
                text = TRIM(text)//" "//r2c(d)
              ELSE
                text = TRIM(text)//" unknown"
              ENDIF
            ENDIF
            CALL Prnt(text)
          ENDDO
        
          !! VARIANCE COMPUTED UNDER LE/GENERAL ASSUMPTIONS
          IF(var_indep) &
            CALL Prnt("Variance in S1 computed UNDER independence", lead=3)
          !ELSE
          !  CALL Prnt("Variance in S1 computed WITHOUT independence", lead=3)
          !ENDIF
  
        ENDIF

        !! **************************************************************** !!
        !!                        STAGE 2 SETTINGS                          !! 
        !! **************************************************************** !!
        IF(.NOT.doS2) THEN
          IF(doCS_all) THEN
            CALL Prnt0("In S2 only CS will be performed.")
          ELSE
            CALL Prnt0("No S2 tests will be performed.")
          ENDIF     
        ELSE

          !CALL Prnt("STAGE 2 :", skip1=1)

          !text = ""
          !IF(WT2 == T2cc) text = "Chisquare test (cases+controls)"
          !IF(WT2 == T2sc) text = "Score test"
          !CALL Prnt("STAGE 2 : "//text, skip1=1, skip2=1)
          
          !text = ""
          !IF(WT2==T2sc) THEN
          !  text = " * Test : Score tests for interactions between loci that are"//&     
          !         " performed independently of the S1 tests only for those pairs of"//&
          !         " loci that were rejected in S1 (analysis model '"//&
          !         TRIM(GetModelName(ana_model))//"')"
          !ELSEIF(WT2==T2cc) THEN
          !  text = " * Test : Chisquare test of difference in LD between cases and"//&     
          !         " controls will be performed independently of the S1 tests and only"//&
          !         " for those pairs loci that were rejected in S1"
          !ENDIF
          !CALL Prnt0(text, lead=3)
          
          !text = " * Indepedence of the two steps achieved by regression of the post-test"//&
          !       " statistic on the vector that generates the S1 statistic"
          !CALL Prnt0(text, lead=3)
          
          IF(var_poststand .AND. WT2==T2sc) THEN
            text = " * If the estimate of the variance of AS"//&
                   " is smaller than "//TRIM(r2c(var_bound))//", the statistic will be"//&
                   " standardized by the average variance of other markers with"//&
                   " sufficiently high variance estimate"
            CALL Prnt0(text, lead=3)
          ENDIF
          
          CALL Prnt("Significance level for S2 is "//r2c(level_S2_nom))
          CALL Prnt("Variance error bound for S2 is "//r2c(var_bound))
          
          IF((doAS4 .OR. doAS1) .AND. WT1>0) THEN
            IF(cntrgrp==0) THEN
              CALL Prnt("Pre-regression centering within AS in S2 is disabled", lead=3)
            ELSE
              text = "Pre-regression centering within AS in S2 will be done using"
              IF(WT1==T1cc) &
                text = TRIM(text)//" controls and cases"
              IF(WT1==T1co) THEN
                IF(cntrgrp==1) text = TRIM(text)//" controls"
                IF(cntrgrp==2) text = TRIM(text)//" cases"
                IF(cntrgrp==3) text = TRIM(text)//" controls and cases"
              ENDIF
              IF(WT1==T1ca) THEN
                IF(cntrgrp==1) text = TRIM(text)//" cases"
                IF(cntrgrp==2) text = TRIM(text)//" controls"
                IF(cntrgrp==3) text = TRIM(text)//" controls and cases"
              ENDIF
              IF(WT1==T1cc) &
                text = TRIM(text)//" CONTROLS AND CASES (POOLED)"
              CALL Prnt(text, lead=3)
            ENDIF
  
            IF(var_grouped) THEN
              text = "Variance of S1 generating vector is estimated using all"
              IF(WT1==T1co) text = TRIM(text)//" controls"
              IF(WT1==T1ca) text = TRIM(text)//" cases"
              IF(WT1==T1cc) text = TRIM(text)//" controls and cases"
              IF(WT1==T1po) text = TRIM(text)//" controls and cases (pooled)"
            ELSE
              text = "Variance of S1 generating vector is estimated using only"
              IF(WT1==T1co) text = TRIM(text)//" controls"
              IF(WT1==T1ca) text = TRIM(text)//" cases"
              IF(WT1==T1cc) text = TRIM(text)//" controls and cases"
              IF(WT1==T1po) text = TRIM(text)//" controls and cases (pooled)"
              text = TRIM(text)//" from S1"
            ENDIF
            CALL Prnt(text, lead=3)
    
          ENDIF
          
          if(.false.) then
          IF(doAS4 .OR. doAS1 .OR. doDS4 .OR. doDS1 .OR. doCS) &
            CALL Prnt0(" * Multiple testing correction :")

          DO k=1,7
            IF(k==1 .AND. .NOT.doAS4) CYCLE
            IF(k==2 .AND. .NOT.doAS1) CYCLE
            IF(k==3 .AND. .NOT.doDS4) CYCLE
            IF(k==4 .AND. .NOT.doDS1) CYCLE
            IF(k==5 .AND. .NOT.doPO4) CYCLE
            IF(k==6 .AND. .NOT.doPO1) CYCLE
            IF(k==7 .AND. .NOT.doCS) CYCLE
            IF(k==1) text = "   - AS4 :"
            IF(k==2) text = "   - AS1 :"
            IF(k==3) text = "   - DS4 :"
            IF(k==4) text = "   - DS1 :"
            IF(k==5) text = "   - PO4 :"
            IF(k==6) text = "   - PO1 :"
            IF(k==7) text = "   - CS :"
            IF(MTC(k)>=one) THEN
              text = TRIM(text)//" Bonferroni correction by user selected value "//&
                     TRIM(r2c(MTC(k)))
            ELSEIF(MTC(1)>zero .AND. MTC(k)<one) THEN
              text = TRIM(text)//" Bonferroni correction factor "//TRIM(r2c(MTC(k)))
            ELSEIF(k<=4) THEN
              text = TRIM(text)//" Bonferroni correction by the number of tests in S2"
            ELSE
              text = TRIM(text)//" Bonferroni correction by the total number of tests"
            ENDIF
            CALL Prnt0(text, lead=9)
          ENDDO
          endif

        ENDIF
  
        IF(doCS .AND. .NOT.doCS_all .AND. doS1 .AND. level_S1_nom<one) &
          CALL Prnt("CS will be computed only for the pairs rejected in S1", lead=3)

        IF(var_poststand .AND. WT2/=T2sc) THEN
          CALL Prnt("Post-standardization of the S2 statistic"//&
                     " is not neccessary for the selected S2 test"//&
                     " and it will not be performed", lead=3)
        ENDIF
          
        !! **************************************************************** !!
        !!                        OTHER SETTINGS                            !! 
        !! **************************************************************** !!
        !CALL Prnt("ADDITIONAL SETTINGS :", skip1=1, skip2=1)
  
        CALL Prnt("Minimal sample size to perform a test is "//TRIM(r2c(min_ss)))
        CALL Prnt("Minimal number of controls to perform a test is "//TRIM(r2c(min_co)))
        CALL Prnt("Minimal number of cases to perform a test is "//TRIM(r2c(min_ca)))
  
        !! TESTING OF LOCI PAIR ON THE SAME CHROMOSOME
        IF(map_nfiles==0) THEN
          CALL Prnt("No mapping file specified. Performing all possible tests", lead=3)
        ELSE
          text = ""
          IF(.NOT.test_same_chr_S1) text = " not"
          CALL Prnt("Pairs of markers on the same chromosome will"//&
                     TRIM(text)//" be tested", lead=3)
        ENDIF
  
        !! CELL AND MARGINAL COUNT CHECKS
        IF(min_cell_cnt>zero) THEN
          IF(cell_count_cor>zero) THEN
            text = " * Pairs with genotype cell frequency below "//&
                   TRIM(r2c(min_cell_cnt))//" will have all genotype counts"//&
                   " increased by "//TRIM(r2c(cell_count_cor))//"" 
          ELSE
            text = " * Pairs with genotype cell frequency below "//&
                   TRIM(r2c(min_cell_cnt))//" will NOT be tested" 
          ENDIF
          CALL Prnt0(text, lead=3)
        ENDIF
        IF(min_marg_cnt>zero) THEN
          text = " * Pairs of markers with genotype marginal frequency below "//&
                 TRIM(r2c(min_marg_cnt))//" will NOT be tested" 
          CALL Prnt0(text, lead=3)
        ENDIF
      
        !! REPORTING BOUNDS FOR OUTPUT FILE
        IF(olim3 > zero .AND. olim3 <= one) THEN
          text = "Reporting only results with S2 corrected"//&
                 " p-values below "//TRIM(r2c(olim3))
        ELSEIF(olim2 > zero .AND. olim2 <= one) THEN
          text = "Reporting only results with S2 raw p-values"//&
                 " below "//TRIM(r2c(olim2))
        ELSEIF(olim1 > zero .AND. olim1 <= one) THEN
          text = "Reporting only results with S1 raw p-values"//& 
                 " below "//TRIM(r2c(olim1))
        ELSE
          text = "All results will be reported"
        ENDIF
        CALL Prnt(text, lead=3)

        !! THE CALCULATIONS OF N1 AND N2 SHOULD BE CORRECTED !!!
        IF(WT1 == T1cc .OR. WT1 == T1sc) THEN
          A_N1 = round_int(dAS*ped_ss)
          D_N1 = round_int(dDS*ped_ss)
        ELSE
          A_N1 = round_int(dAS*ped_ss/two)
          D_N1 = round_int(dDS*ped_ss/two)
        ENDIF
        
        A_N2 = ped_ss
        D_N2 = ped_ss - round_int(dDS*ped_ss)
        
        IF((doAS4 .OR. doAS1) .AND. A_N1<min_ss) THEN
          text = "Size of the sample used in S1 for AS may be very small"//&
                 " which may cause stage 1 tests to be unreliable!" 
          CALL PrntW(text, skip1=1, skip2=1)
        ENDIF
        IF((doDS4 .OR. doDS1) .AND. D_N1<min_ss) THEN
          text = "Size of the sample used in S1 for DS may be very small"//&
                 " which may cause stage 1 tests to be unreliable!" 
          CALL PrntW(text, skip1=1, skip2=1)
        ENDIF
        IF((doAS4 .OR. doAS1) .AND. doS2 .AND. A_N2<min_ss) THEN
          CALL PrntW("Size of the sample used in S2 may be very"//&
                    " small which may cause stage 2 tests to be unreliable!", &
                    skip1=1)
          CALL Prnt("Possible solution: Increase the S2 sample"//&
                    " size for example by setting --delta to a smaller value")
        ENDIF
        IF((doDS4 .OR. doDS1) .AND. doS2 .AND. D_N2<min_ss) THEN
          CALL PrntW("Size of the sample used in S2 may be very"//&
                    " small which may cause stage 2 tests to be unreliable!", &
                    skip1=1)
          CALL Prnt("Possible solution: Increase the S2 sample"//&
                    " size for example by setting --delta to a smaller value")
        ENDIF
        IF(doS2 .AND. doAS .AND. A_N2-A_N1<min_ss) THEN
          text = "The data used in S1 and S2 are"
          IF(A_N2-A_N1>0) text = TRIM(text)//" almost"
          text = TRIM(text)//" equal, which may cause stage 2 tests to be unreliable!"
          CALL PrntW(text, skip1=1)
          text = "Possible solution: Increase the difference between stage 1"//&
                 " and stage 2 sample sizes by setting the parameter --delta"//&
                 " to a smaller value"
          CALL Prnt(text)
        ENDIF
  
        !! Check for too many tests to be reported
        IF(ntests > test_warn_limit .AND. report_all) THEN    
          text = "Reporting of all test results for such high number"//&
                 " of possible test results ("//i2c(ntests)//&
                 ") may produce a high number of very large output files!"
          CALL PrntW(text, skip1=1)
          IF(report_auto_switch) THEN
            olim2 = one
            text = "Automatically switching to reporting only tests where"//&
                   " S2 tests were performed (i.e. --report2 1)"
            CALL Prnt0(text, lead=2)
            text = "Reporting bounds can be set to any value between 0"//&
                   " and 1 by using flags --report1, --report2, --report3"
            CALL Prnt(text, lead=2)
          ENDIF
        ENDIF    

      !! PARALLELIZATION SETTINGS
#ifdef _OPENMP
        CALL Prnt("Parallelization: "//i2c(NUMBER_OF_THREADS)//&
                   " thread(s) will be used", lead=3)
#else
        CALL Prnt("PARALLELIZATION IS NOT AVAILABLE (Switch --nthreads is ignored)", lead=3, skip1=1)
        CALL Prnt("Hint: To enable parallel processing recompile the"//&
                   " source code with an OPENMP-compliant compiler such as"//&
                   " GNU Fortran or Intel Fortran using the necessary flags.", &
                   lead=3, skip1=1)
#endif

      ENDIF ! ELSE to IF(no_testing)
    ENDIF ! test_opts
  ENDIF ! PRESENT(test_opts)

  !! *********************************************************************** !!

  IF(PRESENT(test_start)) THEN
    IF(test_start) THEN 
      IF(report_headers) THEN
        !CALL Prnt0(sstars, skip1=1)
        !CALL Prnt0(text='TESTING', text1='>>>', text2='<<<', l=lh)
        !CALL Prnt0(sstars, skip2=1)
        !CALL Prnt0('>>> TESTING <<<', skip1=1, skip2=0)
      ENDIF

      IF(PRESENT(starttime)) &
        CALL Prnt("Testing started on "//TRIM(starttime(14:24))//" at "//starttime(1:12))
    ENDIF
  ENDIF
  
  !! *********************************************************************** !!

  IF(PRESENT(test_end)) THEN
    IF(test_end) THEN 
      IF(PRESENT(stoptime)) &
        CALL Prnt("Testing finished on "//TRIM(stoptime(14:24))//" at "//stoptime(1:12))
      IF(PRESENT(runtime)) &
        CALL Prnt("Testing runtime was "//TRIM(runtime))

#ifdef _OPENMP
      IF(PRESENT(nthreads)) &
        CALL Prnt("Number of parallel threads used during testing was "//i2c(nthreads))
#endif
    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(test_report)) THEN
    IF(test_report .AND. .NOT.no_testing .AND. .NOT.no_testing_report) THEN
      IF(report_headers) THEN
        !CALL Prnt0(sstars, skip1=1)
        !CALL Prnt0(text='TESTING REPORT', text1='>>>', text2='<<<', l=lh)
        !CALL Prnt0(sstars, skip2=1)
        CALL Prnt0('>>> TESTING REPORT <<<', skip1=1, skip2=1)
      ENDIF

      IF(PRESENT(runtime)) &
        CALL Prnt0("   * Testing runtime : "//runtime, skip2=1)

      CALL Prnt0(text1="   * Excluded pairs on same chromosome   : ", &
                 text2=i2c(n_samechr_excl), l=rw)
      CALL Prnt0(text1="   * Included pairs on same chromosome   : ", &
                 text2=i2c(n_samechr_incl), l=rw)
      CALL Prnt0(text1="   * Skipped pairs with low sample size  : ", &
                 text2=i2c(NumLowSS), l=rw, skip2=1)

      IF(WT1==0) THEN
        CALL Prnt0("  STAGE 1  >  No S1 tests performed!", skip2=1)
        CALL Prnt0(text1="   * Pairs with corrected count error    : ", &
                   text2=i2c(NumLowCellNumCor), l=rw)
        CALL Prnt0(text1="   * Pairs with count error              : ", &
                   text2=i2c(NumLowCell+NumLowMarg), &
                        l=rw)
        IF(NumLowCell+NumLowMarg>0) THEN
          CALL Prnt0(text1="      - Cell count errors                : ", &
                     text2=i2c(NumLowCell), l=rw)
          CALL Prnt0(text1="      - Marginal count errors            : ", &
                     text2=i2c(NumLowMarg), l=rw)
        ENDIF
        GOTO 130
      ENDIF

      text = "  STAGE 1  >  Selected test : "//i2c(WT1) 
      CALL Prnt0(text, skip2=1)

      CALL Prnt0(text1="   * Attempted tests                     : ", &
                 text2=i2c(NumS1), l=rw)
      IF(NumS1==0) GOTO 130
      CALL Prnt0(text1="   * Successful tests                    : ", &
                 text2=i2c(NumS1OK), l=rw)
      CALL Prnt0(text1="   * Erroneous tests                     : ", &
                 text2=i2c(NumS1NotOK), l=rw)
      IF(NumS1NotOK>0) THEN
        CALL Prnt0(text1="      - Cell count errors                : ", &
                   text2=i2c(NumLowCell), l=rw)
        CALL Prnt0(text1="      - Marginal count errors            : ", &
                   text2=i2c(NumLowMarg), l=rw)
        CALL Prnt0(text1="      - Low sample size stage 1 errors   : ", &
                   text2=i2c(NumLowSS1), l=rw)
        CALL Prnt0(text1="      - Low sample size stage 2 errors   : ", &
                   text2=i2c(NumLowSS2), l=rw)
        NErrOther = NumS1NotOK - NumLowCell - NumLowMarg &
                          - NumLowSS1 - NumLowSS2
        CALL Prnt0(text1="      - Other errors                     : ", &
                   text2=i2c(NErrOther), l=rw)
      ENDIF

      CALL Prnt0(text1="   * S1 sample size ratio                : ", &
                 text2=r2c(dAS), l=rw)                 

      IF(NumS1_sex==0) THEN
        CALL Prnt0(text1="   * Average sample size (error-free)    : ", &
                   text2=i2c(INT(MeanSS(0))), l=rw)                 
      ELSE
        CALL Prnt0(text1="   * Tests without sex chromosome        : ", &
                   text2=i2c(NumS1-NumS1_sex), l=rw)
        CALL Prnt0(text1="   * Tests with sex chromosome           : ", &
                   text2=i2c(NumS1_sex), l=rw)
      ENDIF

      IF(.NOT.auto_level_use .OR. auto_level_single) THEN
        CALL Prnt0(text1="   * Level of significance               : ", &
                   text2=r2c(level_S1), l=rw)
      ELSE
        CALL Prnt0(text1="   * Level of significance               : ", &
                   text2="automatic", l=rw)
      ENDIF
      CALL Prnt0(text1="   * Number of significant statistics    :", &
                 text2=i2c(MAXVAL(NumLD))//" *", l=rw+2)
      IF(NumS1OK>0) THEN
        CALL Prnt0(text1="   * Ratio of significance               :", &
                   text2=r2c(REAL(MAXVAL(NumLD), dpp)/NumS1OK), l=rw)
      ELSE
        CALL Prnt0(text1="   * Ratio of significance               :", &
                   text2="0", l=rw)
      ENDIF
      
      130 CONTINUE
      
      text = "  STAGE 2  >  Selected test : "//i2c(WT2) 
      CALL Prnt0(text, skip1=1)

      !!! ***************************************************************!!!
      !!! *************** ADJUSTED SCORE TEST STATISTIC *****************!!!
      !!! ***************************************************************!!!
      
      !! Allocate AT
      CALL ResizeVar(AT, mmtl, 3)
      
      !! Loop over adjusted, disjoint and classical score
      DO i=1,7
      
        IF(i==1 .AND. .NOT.doAS4) CYCLE
        IF(i==2 .AND. .NOT.doAS1) CYCLE
        IF(i==3 .AND. .NOT.doDS4) CYCLE
        IF(i==4 .AND. .NOT.doDS1) CYCLE
        IF(i==5 .AND. .NOT.doPO4) CYCLE
        IF(i==6 .AND. .NOT.doPO1) CYCLE
        IF(i==7 .AND. .NOT.doCS) CYCLE
        AT = ""
        skip = 0
        AvgSS = zero

        IF((i==4 .AND. doDS4) .OR. i==5 .OR. i==6) CYCLE

        !! Assign values if adjusted score
        IF(i==1) THEN
          NS2 = NumS2AS4 
          NS2OK = NumS2AS4OK 
          NS2NotOK = NumS2AS4NotOK
          NS2PostStand = NumS2AS4PostStand
          NTest = MAX(NumS2AS4, NumS1)
          NErrVar = NumVarErr(i)
          NLowVar = NumLowVar(i)
          NNegVar = NumNegVar(i)
          NLinDep = NumLinDep(i)
          AvgSS = MeanSS(i)
        ELSEIF(i==2) THEN
          NS2 = NumS2AS1 
          NS2OK = NumS2AS1OK 
          NS2NotOK = NumS2AS1NotOK
          NS2PostStand = NumS2AS1PostStand
          NTest = MAX(NumS2AS1, NumS1)
          NErrVar = NumVarErr(i)
          NLowVar = NumLowVar(i)
          NNegVar = NumNegVar(i)
          NLinDep = NumLinDep(i)
          AvgSS = MeanSS(i)
        !! Assign values if disjoint score
        ELSEIF(i==3 .OR. i==4) THEN 
          NS2 = NumS2DS
          NS2OK = NumS2DSOK 
          NS2NotOK = NumS2DSNotOK
          NTest = MAX(NumS2DS, NumS1)
          NErrVar = NumVarErr(3)
          NLowVar = NumLowVar(3)
          NNegVar = NumNegVar(3)
          NLinDep = NumLinDep(3)
          AvgSS = MeanSS(3)
        ELSEIF(i==5 .OR. i==6) THEN
        
          CALL PrntE("Missing code for pooled tests", Q=.TRUE.) 

          NErrVar = NumVarErr(4)
          NLowVar = NumLowVar(4)
          NNegVar = NumNegVar(4)
          NLinDep = NumLinDep(4)
          AvgSS = MeanSS(4)

        ELSEIF(i==7) THEN 
          NS2 = NumS2CS
          NS2OK = NumS2CSOK 
          NS2NotOK = NumS2CSNotOK
          NTest = NumS2CS
          NErrVar = NumVarErr(4)
          NLowVar = NumLowVar(4)
          NNegVar = NumNegVar(4)
          NLinDep = NumLinDep(4)
          AvgSS = MeanSS(4)
        ENDIF
        
        NS2Sex = NumS2_sex(i)
        Correct = MTC(i)
        
        NErrOther = NS2NotOK - NLinDep - NErrVar
         
        !! First the header for the current results
        j = 1
        IF(i==1) THEN
          AT(j,1) = "  ** ADJUSTED-4 SCORE TEST (AS4) **"
        ELSEIF(i==2) THEN
          AT(j,1) = "  ** ADJUSTED-1 SCORE TEST (AS1) **"
        ELSEIF(i==3) THEN
          AT(j,1) = "  ** DISJOINT-4 SCORE TEST (DS4) **"
        ELSEIF(i==4) THEN
          AT(j,1) = "  ** DISJOINT-1 SCORE TEST (DS1) **"
        ELSEIF(i==5) THEN
          AT(j,1) = "  ** POOLED-4 SCORE TEST (PO4) **"
        ELSEIF(i==6) THEN
          AT(j,1) = "  ** POOLED-1 SCORE TEST (PO1) **"
        ELSEIF(i==7) THEN
          AT(j,1) = "  ** CLASSICAL SCORE TEST (CS) **"
        ENDIF
  
        skip(j,1) = 1
        skip(j,2) = 1
        j = j+1
        AT(j,1) = "   * Attempted tests                     :"
        AT(j,2) = i2c(NS2)

        IF(NS2==0) GOTO 131

        j = j+1
        AT(j,1) = "   * Successful tests                    :"
        AT(j,2) = i2c(NS2OK)
        
        !! Only for ADJUSTED SCORE
        IF((i==1 .OR. i==2) .AND. var_poststand) THEN
          j = j+1
          AT(j,1) = "       - Normally standardized statistics :"
          AT(j,2) = i2c(NS2OK - NS2PostStand)
          j = j+1
          AT(j,1) = "       - Post-standardized statistics     :"
          AT(j,2) = i2c(NS2PostStand)
          j = j+1
          AT(j,1) = "   * Number of variance problems         :"
          AT(j,2) = i2c(NErrVar)
          IF(NLowVar>0) THEN
            j = j+1
            AT(j,1) = "      - Low (positive) variance errors   :"
            AT(j,2) = i2c(NLowVar)
            j = j+1
            AT(j,1) = "      - Negative variance errors         :"
            AT(j,2) = i2c(NNegVar)
            j = j+1
            AT(j,1) = "      - Resolved variance problems       :"
            AT(j,2) = i2c(NS2PostStand)
            j = j+1
            AT(j,1) = "      - Unresolved variance problems     :"
            AT(j,2) = i2c(NErrVar - NS2PostStand)
          ENDIF
        ENDIF

        j = j+1
        AT(j,1) = "   * Erroneous tests                     :"
        AT(j,2) = i2c(NS2NotOK)
        IF(NS2NotOK>0) THEN
          j = j+1
          AT(j,1) = "      - Linear dependency errors         :"
          AT(j,2) = i2c(NLinDep)
          IF(i==1 .AND. var_poststand) THEN
            j = j+1
            AT(j,1) = "      - Unresolved variance errors       :"
            AT(j,2) = i2c(NErrVar - NS2PostStand)
          ELSE
            j = j+1
            AT(j,1) = "      - Variance errors                  :"
            AT(j,2) = i2c(NErrVar) 
            j = j+1
            AT(j,1) = "        - Low (positive) variance errors :"
            AT(j,2) = i2c(NLowVar)
            j = j+1
            AT(j,1) = "        - Negative variance errors       :"
            AT(j,2) = i2c(NNegVar)
          ENDIF
          j = j+1
          AT(j,1) = "      - Other errors                     :"
          AT(j,2) = i2c(NErrOther)
        ENDIF

        IF(NS2Sex==0) THEN
          j = j+1
          AT(j,1) = "   * Average sample size (error-free)    :"
          AT(j,2) = i2c(INT(AvgSS))
        ELSE
          j = j+1
          AT(j,1) = "   * Tests without sex chromosome        :"
          AT(j,2) = i2c(NS2-NS2Sex)
          j = j+1
          AT(j,1) = "   * Tests with sex chromosome           :"
          AT(j,2) = i2c(NS2Sex)
        ENDIF

        j = j+1         
        AT(j,1) = "   * Level of significance               :"
        AT(j,2) = r2c(level_S2)

        j = j+1
        AT(j,1) = "   * Multiple testing correction         :"
        AT(j,2) = i2c(INT(Correct,ikb))

        j = j+1
        AT(j,1) = "   * Significant p-values (raw)          :"
        AT(j,2) = i2c(NumEpi(i))

        j = j+1
        AT(j,1) = "   * Ratio of significance (raw)         :"
        !AT(j,2) = r2c(REAL(NumEpi(i), dpp) * Invert(REAL(NS2OK, dpp)))
        AT(j,2) = r2c(Divide(NumEpi(i), NS2OK))

        j = j+1
        AT(j,1) = "   * Significant p-values (corrected)    :"
        AT(j,2) = i2c(NumEpiCorr(i))
        AT(j,3) = " *"
        j = j+1
        AT(j,1) = "   * Ratio of significance (corrected)   :"
        !AT(j,2) = r2c(REAL(NumEpiCorr(i), dpp) * Invert(REAL(NTest, dpp)))
        AT(j,2) = r2c(Divide(NumEpiCorr(i), NTest))

        131 CONTINUE
         
        DO k=1,j
          IF(ALL(AT(k,:)=="")) CYCLE
          width = rw + LEN_TRIM(AT(k,3))   
          CALL Prnt0(text1=AT(k,1), text2=TRIM(AT(k,2))//AT(k,3), & 
                          skip1=skip(k,1), skip2=skip(k,2), l=width)
        ENDDO
        
        IF(i==3 .AND. NS2 > 0 .AND. AvgSS < min_ss) &
          CALL PrntW("DS may be unreliable due to small sample size!", skip1=1, skip2=1)
      ENDDO
      
      CALL Prnt0(sstars, skip1=1)

    ENDIF
  ENDIF

  !! *********************************************************************** !!

  IF(PRESENT(exec_end)) THEN
    IF(exec_end) THEN
      
      CALL Prnt0("")
      
      !! Print error announcement
      IF(nerrors>0) THEN
        IF(nerrors==1) THEN
          CALL Prnt("IMPORTANT: There was 1 ERROR!")
        ELSEIF(nerrors>1) THEN
          CALL Prnt("IMPORTANT: There were "//i2c(nerrors)//" ERRORS!")
        ENDIF
        IF(nerrors>nerrors_limit) &
          CALL Prnt("Only "//i2c(nerrors_limit)//" errors were reported to"//&
                    " the standard output. Check the log file to see all"//&
                    " of them", log=.FALSE.)
      ENDIF
      
      !! Print warning announcement
      IF(nwarnings>0) THEN
        IF(nwarnings==1) THEN
          CALL Prnt("IMPORTANT: There was 1 warning!")
        ELSEIF(nwarnings>1) THEN
          CALL Prnt("IMPORTANT: There were "//i2c(nwarnings)//" warnings!")
        ENDIF
        IF(nwarnings>nwarnings_limit) &
          CALL Prnt("Only "//i2c(nwarnings_limit)//" warnings were reported"//&
                    " to standard output. Check the log file to see all of"//&
                    " them", log=.FALSE.)
      ENDIF
      
      !! Print runtime info
      IF(PRESENT(stoptime)) &
        CALL Prnt("Execution finished on "//TRIM(stoptime(14:24))//" at "//stoptime(1:12))
      IF(PRESENT(runtime)) &
        CALL Prnt("Total program runtime was "//TRIM(runtime))

      !! Print error announcement
      IF(nerrors>0) THEN
        CALL Prnt("LIST OF ISSUED ERRORS :", &
                       screen=.FALSE., skip1=1)
        DO i=1,nerrors
          CALL Prnt0("  "//errors(i), screen=.FALSE., lead=2)
        ENDDO
        CALL Prnt0("", screen=.FALSE.)
      ENDIF
      
      !! Print warning announcement
      IF(nwarnings>0) THEN
        CALL Prnt("LIST OF ISSUED WARNINGS :", &
                       screen=.FALSE., skip1=1)
        DO i=1,nwarnings
          CALL Prnt0("  "//warnings(i), screen=.FALSE., lead=2)
        ENDDO
        CALL Prnt0("", screen=.FALSE.)
      ENDIF
      
    ENDIF  
  ENDIF  

  !! *********************************************************************** !!

  RETURN
    
END SUBROUTINE Report

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_REPORT
