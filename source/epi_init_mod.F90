MODULE EPI_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE EPI_UTILS        !! VARIOUS TOOLS
  USE EPI_SIMUL        !! NEEDED TO CONVERT MODEL NAME INTO A NUMBER
                       !! AND BACK (FUNCTIONS GetModelName, GetModelNumber)
                       !! AND FOR OTHER SIMULATION AND ANALYSIS MODEL VARIABLES 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      THIS CODE NEEDS TO BE PREPROCESSED BEFORE COMPILING!       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  CHARACTER(mfl), ALLOCATABLE :: ped_file(:), &
                                 bed_file(:), &
                                 map_file(:), &
                                 tmp_file(:), &
                                 sub_file(:), &
                                 res_file(:),&
                                 prev_ped_file(:), &
                                 prev_bed_file(:), &
                                 prev_map_file(:), &
                                 prev_sub_file(:), &
                                 sort_file(:)
  
  CHARACTER(mcl), ALLOCATABLE :: cmdlines(:), &
                                 useed_str(:)
  
  CHARACTER(mal), ALLOCATABLE :: args_all(:), &
                                 args_selected(:,:), &
                                 args_unknown(:)
  
  LOGICAL, ALLOCATABLE        :: args_used(:)

  CHARACTER(mfl), SAVE        :: file_base = "", &
                                 def_file_base = "", &
                                 fam_file = "", &
                                 cmd_file = "", &
                                 out_file = "", &
                                 out_pss_file = "", &
                                 out_maf_file = "", &
                                 out_prev_file = "", &
                                 pss_file = "", &
                                 sts_file = "", &
                                 save_file = "", &
                                 save_ped_file = "", &
                                 save_map_file = "", &
                                 save_fam_file = ""
                                 
  CHARACTER(mcl), SAVE        :: cmdline = "", &
                                 user_cmdline = ""
                                 
  CHARACTER(30), SAVE         :: input = ""
   
  CHARACTER(mtsl), SAVE       :: starttime_full = "", &
                                 stoptime_full = "", &
                                 runtime_total = ""
  
  CHARACTER, SAVE             :: ped_cs(2) = "", &
                                 map_cs(2) = "", &
                                 sub_cs(2) = "", &
                                 fam_cs(2) = "", &
                                 out_cs = "", &
                                 NA = "", &
                                 dec = "", &
                                 delim = "", &
                                 NA_colsep = ""
  
  CHARACTER(2), SAVE          :: out_sort_sel

  LOGICAL, SAVE :: premat_halt, ask_halt, no_halt, stop_on_unknown_args, &
    run_analysis, cmdfile_present, check_fexist, keep_temp, &
    save_input_data, save_binary, do_out_pss, do_out_maf, do_out_actual_prev, &
    save_full_fam, case_sensitive_flags, append_ts, append_pid, no_output, &
    simulate_data, simulate_HWE, no_testing_report, no_data_summary, &
    doS2, no_testing, no_data_input, replication_mode, only_S1, &
    report_headers, fixed_MAF, do_all_pretests, tmp_to_out_only, out_append, &
    pss_file_subsample, out_available, out_legend, out_header_names, &
    out_header_comment, out_sort, out_sort_tmp, out_sort_multiple_files, &
    out_sort_second, out_discard_unsorted, out_epi_effect, &
    out_loc_info, bed_major, ped_minor, LE_population, &
    rs_valid_chck, fix_subsamp, var_indep, var_use_beta1, &
    correct_all_cells, user_delta, simple_delta, auto_delta, auto_d_incapable, &
    auto_d_optimize, auto_d_maf_fixed, auto_d_oracle, &
    auto_level_use, auto_level_report, auto_maf1, auto_maf2, auto_maf3, &
    auto_level_DS, auto_level_merge, auto_level_maf, &
    auto_level_single, auto_level_weight_by_counts, out_level_S1, &
    nonpar_probs, ntests_limit, add_settings_filename, test_same_chr, & 
    add_settings_filename_long, test_same_chr_S1, keep_all_chr, keep_all_bp, &
    sex_chr_chck, excluded_chr(0:(SIZE(chr_name)-1)), excl_male, &
    excl_fema, exclude_invalid_rs, assume_male, assume_fema, &
    assume_status, report_auto_switch, report_all, report_errs, &
    doS1, doAS, doAS4, doAS1, doDS4, doDS1, doCS, doPO4, doPO1, &
    doPTAS4, doPTAS1, doPTDS4, doPTDS1, doPO4CS, doPO1CS, doPO, &
    doPT4, doPT1, reject_all_pretests, &
    always_pretest, doCS_all, var_poststand, poststand_by_varbound, &
    file_opened, out_header, out_statistic, out_errcode, &
    out_errcode_glob, out_ss, out_pvalcor, out_mtc, &
    analyze_files_together, reuse_data, tmp_sorted, &
    out_testid, out_minimal, bed_data_SNPmajor, out_bed_data_SNPmajor, &
    out_debug_info, out_variance, out_alleles, out_S1_vector, &
    check_bed_filesize, sep_ascii, no_causal_pair, priority_user_cmdline, &
    var_grouped, out_zipped, cov_grouped, compute_maf, out_maf, &
    out_maf_coca, sub_read_all, sub_make_pairs, sub_include, sub_no_chr, &
    sim_single_chr, user_set_out, user_set_save, regres_divided, &
    fam_single_sample

  REAL(dpp), SAVE :: starttime_glob, starttime, stoptime, max_runtime, &
    current_time, ntests_assumed, olim1, olim2, olim3, &
    dAS, dDS, delta_bound, delta_start, delta_end, delta_step, &
    level_S1, auto_level_optim, level_S1_nom, level_S2, &
    level_S2_nom, var_bound, var_cor_AS, varAS_sum, varAS_sum_n, &
    allelefreqA, allelefreqB, allelefreqC, min_MAF, auto_maf, prevalence, OR, &
    OR_min, OR_max, OR_step, OR1, OR2, LD, auto_level_prev, auto_level_OR1, &
    auto_level_OR2, auto_level_OR, auto_level_min, auto_level_factor, &
    auto_level_max, auto_d_prev, auto_d_min, auto_d_max, &
    auto_d_minlevel, auto_d_maxlevel, auto_d_OR1, auto_d_OR2, &
    auto_d_OR, auto_d_maf1, auto_d_maf2, auto_d_COR, &
    min_cell_cnt, min_marg_cnt, min_ss, min_co, min_ca, &
    cell_count_cor, MeanSS(0:4), MTC(7), MTC_factor, auto_d_frac, res_prev

  INTEGER, SAVE :: seed_start, seed_const, nseed, max_reruns, nunknown_args, &
    ana_model, ana_model_S1, sim_model, ped_maxncol, &
    res_maxncol, ncmdlines, result_maxlines, temp_cycle, & !, max_samplesize, & 
    ped_nsamp, ped_nloci, ped_ss, ped_nco, ped_nca, &
    prev_ped_nco, prev_ped_nca, ped_nmale, ped_nfema, &
    prev_ped_nmale, prev_ped_nfema, prev_ped_nsamp, prev_ped_nloci, &
    prev_ped_ss, ped_nskip, input_format, &
    out_ncol, out_sort_col, map_ncols, sub_ncols, fam_ncols, &
    map_nlines, sub_nlines, map_nskip, sub_nskip, fam_nskip, &
    pss_nskip, sts_nskip, result_nskip, WT1, WT2, &
    T1_df, T1_df_co, T1_df_po, T2_df, &
    cntrgrp, min_map_ncols, map_chcol, map_rscol, map_dscol, &
    map_bpcol, map_a1col, map_a2col, sub_chcol, sub_rscol, fam_sexcol, &
    fam_fidcol, fam_iidcol, fam_pidcol, fam_midcol, fam_stscol, &
    ped_nrepeat, ped_nfiles, bed_nfiles, map_nfiles, fam_nfiles, &
    pss_nfiles, sts_nfiles, out_nfiles, tmp_nfiles, results_nfiles, sub_nfiles,&
    cT, cTpval, cAS, cASpval, def_cT, def_cTpval, def_cAS, &
    def_cASpval, ncausalloci, nneutralloci, nLDpairs, n_sub_excl, &
    n_chr_excl, n_rs_excl, n_negbp_excl, n_maf_excl, &
    n_chr_incl, n_maxntests_excl, n_samechr_incl, &
    n_samechr_excl, auto_d_n, auto_d_nlevel, auto_d_ntests, &
    ndig_stat, ndig_pval, ndig_mafs, ndelta, sim_chr_start, &
    sim_rs_start, sim_fid_start, sim_iid_start 

  INTEGER(ikb), SAVE :: max_out_size, NumS1, NumS1OK, NumS1NotOK, & !tmp_nlines, &
    NumS2, NumS2AS4, NumS2AS4OK, NumS2AS4NotOK, &
    NumS2AS4PostStand, NumS2AS1, NumS2AS1OK, NumS2AS1NotOK, &
    NumS2AS1PostStand, NumS2DS, NumS2DSOK, NumS2DSNotOK, &
    NumS2CS, NumS2CSOK, NumS2CSNotOK, NumS1_sex, &
    NumS2_sex(7), NumVarErr(4), NumLowVar(7), &
    NumNegVar(7), NumLinDep(4), NumLowCell, NumLowCellNumCor, &
    NumLowMarg, NumLowSS, NumLowSS1, NumLowSS2, &
    NumErrs, NumLD(6), NumEpi(7), NumEpiCorr(7), NumEpiC(2), &
    ntests, ntest_limit, ntest_limit2

  INTEGER(iks), ALLOCATABLE    :: X(:,:), &
                                  sts(:), &
                                  sex(:), &
                                  pss_AS(:), &
                                  pss_DS(:), &
                                  sts1(:), &
                                  sex1(:), &
                                  pss_AS1(:), &
                                  pss_DS1(:)

  INTEGER, ALLOCATABLE         :: seed(:), &
                                  useed_int(:), &
                                  useed_par(:), &
                                  Included_pairs(:,:) 
                                  !ped_ncols(:), prev_ped_ncols(:),

  INTEGER(ikb), ALLOCATABLE    :: MAPARS(:), &
                                  SUBMAPARS(:)                                   

  LOGICAL, ALLOCATABLE         :: InclLoc(:)

  REAL(dpp), ALLOCATABLE       :: MAF(:,:), MAFco(:,:), MAFca(:,:)

  CHARACTER(mml), ALLOCATABLE  :: MAPA(:,:), &
                                  SUBMAPA(:,:) 

  CHARACTER, ALLOCATABLE       :: FAM(:,:)

  !! Matrix MatOpt will store the values of the optimal S1 level (1),  
  !! the non-centrality parameter of the pretest (2) and the slope of the
  !! disjoint score test (3) for various combinations of maf1 and maf2 so that
  !! they have to be computed only once for each combination of maf1 and maf2
  !! of which there are only 100x100 pairs if maf computed with precision of 
  !! 2 digits. If OPENMP used, each thread keeps separate copy of MatOpt
  !! and fills it with its own values (no sharing between threads). The matrix
  !! should be initialized to a negative value.
  REAL(dpp) :: MatOpt(5, 10**opt_ndig_mafs, 10**opt_ndig_mafs)

  PUBLIC  :: InitGlobalVars

  CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE InitGlobalVars() 
!! This subroutine assigns default values to global variables  
  IMPLICIT NONE

#ifdef __INTEL_COMPILER
  screen_advance_always = .TRUE.
  log_advance_always = .TRUE.
  countdown = .TRUE.
  countdown_jump = 5
#else
#ifdef __GFORTRAN__
  screen_advance_always = .FALSE.
  log_advance_always = .TRUE.
  countdown = .TRUE.
  countdown_jump = 1
#endif
#endif

  !! GENERAL PARAMETERS
#ifdef _OPENMP
  NUMBER_OF_PROCESSORS = OMP_GET_NUM_PROCS()
#else
  NUMBER_OF_PROCESSORS = 1
#endif

  NUMBER_OF_THREADS = 1

  CALL ResizeVar(useed_int, 1, -1)
  CALL ResizeVar(useed_par, 1, -1)

  stop_on_unknown_args = .TRUE.
  comment = def_comment
  max_runtime = -one
  seed_start = 0
  seed_const = 0
  max_reruns = 0
  no_testing = .FALSE.
  no_data_input = .FALSE.
  no_log = .FALSE.
  silent = .FALSE.        ! When true, then no output to console
  add_commas = .TRUE.     ! When true then the numbers on output will be printed with thousand separators (commas)
             
  user_set_out = .FALSE.
  user_set_save = .FALSE.
  simulate_data = .FALSE.
  simulate_HWE = .TRUE.
  reuse_data = .FALSE.
  fix_subsamp = .TRUE.
  fixed_MAF = .FALSE.
  do_out_maf = .FALSE.
  ask_halt = .TRUE.
  no_halt = .FALSE.
  keep_temp = .FALSE.
  tmp_to_out_only = .FALSE.
  append_ts = .TRUE.
  append_pid = .TRUE.
  add_settings_filename = .FALSE.
  add_settings_filename_long = .FALSE.
  no_openfile_iostat = .FALSE.
  pause_run = .FALSE.
  out_minimal = .FALSE.
  out_debug_info = .FALSE.
  out_variance = .FALSE.
  out_alleles = .FALSE.
  out_testid = .TRUE.
  out_S1_vector = .FALSE.
  check_fexist = .TRUE.
  check_bed_filesize = .TRUE.
  compute_maf = .FALSE.
  out_maf = .TRUE.
  out_maf_coca = .TRUE.
  pss_file_subsample = .FALSE.
  tmp_sorted = .FALSE.
  out_statistic = .TRUE.
  out_errcode = .TRUE.
  out_errcode_glob = .FALSE.
  out_ss = .TRUE.
  out_pvalcor = .TRUE.
  out_mtc = .TRUE.
  do_all_pretests = .FALSE.
  out_zipped = .FALSE.
  regres_divided = .FALSE.
  fam_single_sample = .TRUE.
  
  priority_user_cmdline = .TRUE.
  report_all = .FALSE.                            ! y/n all test results are outputted
  report_errs = .FALSE.                           ! y/n erroneous test results are outputted
  report_auto_switch = .FALSE.
  out_available = .FALSE.
  out_append = .FALSE.
  out_legend = .TRUE.
  out_header = .TRUE.
  out_header_comment = .FALSE.                    ! y/n the output header line is commented
  out_header_names = .TRUE.                       ! y/n the output header line has column names
  out_sort = .FALSE.
  out_sort_tmp = .FALSE.
  out_sort_multiple_files = .FALSE.
  out_sort_second = .TRUE.
  out_discard_unsorted = .TRUE.
  out_sort_sel = def_sort_sel
  out_epi_effect = .FALSE.
  out_loc_info = .TRUE.

  case_sensitive_flags = .TRUE.                   ! y/n call flags are case sensitive
  no_output = .FALSE.                             ! y/n no output
  no_testing_report = .TRUE.
  no_data_summary = .FALSE.
  report_headers = .TRUE.
  bed_data_SNPmajor = .TRUE.
  out_bed_data_SNPmajor = .TRUE.
  sep_ascii = .FALSE.

  ndig_stat = def_ndig_stat                       ! num of digits of output columns 'statistics'
  ndig_pval = def_ndig_pval                       ! num of digits of output columns 'p-value'
  ndig_mafs = def_ndig_mafs                       ! num of digits of output columns 'maf' 
  ped_maxncol = 100000                            ! maximum input file line length 
  res_maxncol = 100                               ! maximum input result file length
  max_out_size = FLOOR(max_file_size * mfsf, ikb) ! maximum size of output/temp files
                                        
  temp_cycle = df_temp_cycle                      ! output cycle size

  def_cT = 2                                      ! position of column 'statistics' for results input
  def_cTpval = 3                                  ! position of column 'p-value' for results input
  def_cAS = 4                                     ! position of main statistic for results input
  def_cASpval = 5                                 ! position of main p-value for results input
  
  run_analysis  = .FALSE.
  analyze_files_together = .FALSE.
  
  !! GENERAL I/O PARAMETERS
  CALL get_environment_variable("DELIMITER", delim)
  IF(LEN_TRIM(delim)==0) delim = def_delim

  dec = "."                                       ! Decimal point on input (NOT IMPLEMENTED!)
  NA_colsep = ""
  NA = def_NA

  !! SIMULATION VARIABLES
  save_input_data = .FALSE.
  save_binary     = .TRUE.
  allelefreqA     = NAneg
  allelefreqB     = NAneg
  allelefreqC     = NAneg
  prevalence      = def_simul_prev
  no_causal_pair  = .FALSE.
  ncausalloci     = 2
  nneutralloci    = 0
  nLDpairs        = 0
  
  OR  = one
  OR1 = one
  OR2 = one
  LD  = zero
  OR_min = -one
  OR_max = -one
  OR_step = zero
  sim_single_chr = .TRUE.
  sim_chr_start = 0
  sim_rs_start = 1
  sim_fid_start = 1
  sim_iid_start = 1
  auto_level_factor = one
  auto_level_min = minimum_test_level
  auto_level_max = maximum_test_level
  auto_level_prev = def_auto_level_prev
  auto_level_OR1 = def_auto_level_OR1
  auto_level_OR2 = def_auto_level_OR2
  auto_level_OR = def_auto_level_OR
  auto_d_min = def_auto_d_min
  auto_d_max = def_auto_d_max
  auto_d_n = def_auto_d_n 
  auto_d_minlevel = def_auto_d_minlevel
  auto_d_maxlevel = def_auto_d_maxlevel
  auto_d_nlevel = def_auto_d_nlevel 
  auto_d_ntests = def_auto_d_ntests 
  auto_d_prev = def_auto_d_prev
  auto_d_OR1 = def_auto_d_OR1
  auto_d_OR2 = def_auto_d_OR2
  auto_d_OR = def_auto_d_OR
  auto_d_maf1 = def_auto_d_maf1
  auto_d_maf2 = def_auto_d_maf2
  auto_d_COR = def_auto_d_COR
  auto_d_frac = one
  ana_model = def_ana_model
  ana_model_S1 = -1
  sim_model = def_sim_model
  
  !! PEDDING FILE variables
  bed_major = .TRUE.
  ped_minor = .TRUE.
  ped_nsamp = 1
  ped_ss = 0
  ped_nco = 0
  ped_nca = 0
  ped_nmale = 0
  ped_nfema = 0 
  
  ped_cs = def_ped_cs           ! Column separator in the input ped file
  ped_nskip = 0                     ! Number of lines to skip in input ped file
  CharNAgen = def_CharNAgen      ! NA string for status
  
  CharNAsts = def_CharNAsts    ! NA string for status
  CharCa = def_CharCa            ! Status case string in input data file
  CharCo = def_CharCo      ! Status control string in input data file
  
  CharNAsex = def_CharNAsex          ! NA string for status
  CharMa = def_CharMa            ! Sex male string in input data file
  CharFe = def_CharFe        ! Sex female string in input data file

  !! MAPPING FILE variables
  map_nskip = 0                     ! Number of lines to skip in input mapping file
  map_cs = def_map_cs           ! Column separator in the input mapping file
  map_chcol = def_map_chcol
  map_rscol = def_map_rscol
  map_dscol = def_map_dscol
  map_bpcol = def_map_bpcol
  map_a1col = def_map_a1col
  map_a2col = def_map_a2col

  min_map_ncols = MAX(def_map_rscol, def_map_chcol, def_map_dscol, &
                      def_map_bpcol, def_map_a1col, def_map_a2col)
  map_ncols = -1
  input_format = 0

  sub_nskip = 0
  sub_cs = def_map_cs           ! Column separator in the input submapping file
  sub_chcol = def_sub_chcol
  sub_rscol = def_sub_rscol
  sub_ncols = MAX(4, sub_chcol, sub_rscol)
  sub_read_all = .FALSE.
  sub_make_pairs = .FALSE.
  sub_no_chr = .FALSE.
  
  rs_valid_chck = .FALSE.
  
  fam_nskip = 0
  fam_cs = def_fam_cs           ! Column separator in the input family info file
  fam_fidcol = def_fidcol
  fam_iidcol = def_iidcol
  fam_pidcol = def_pidcol
  fam_midcol = def_midcol
  fam_sexcol = def_sexcol
  fam_stscol = def_stscol
  fam_ncols = MAX(def_famncol,fam_stscol,fam_sexcol)

  out_cs = def_out_cs
  out_sort_col = -1
  
  pss_nskip = 0
  sts_nskip = 0
  result_nskip = 0
  result_maxlines = 0
  
  !! SETTING OF PARAMETERS FOR STATISTICAL TESTS

  doS1 = .TRUE.
  doAS4 = .TRUE.              ! y/n adjusted score (AS4) gets reported
  doAS1 = .TRUE.              ! y/n adjusted score (AS1) gets reported
  doDS4 = .TRUE.              ! y/n disjoint score (DS4) gets reported
  doDS1 = .TRUE.              ! y/n disjoint score (DS1) gets reported
  doCS = .TRUE.               ! y/n classical score (CS) gets reported
  doPTAS4 = .TRUE.
  doPTAS1 = .TRUE.
  doPTDS4 = .TRUE.
  doPTDS1 = .TRUE.
  doCS_all = .FALSE.          ! And whether it should be computed for all tests

  doPT4 = .FALSE.
  doPT1 = .FALSE.
  doPO4 = .TRUE.
  doPO1 = .TRUE.
  doPO4CS = .TRUE.
  doPO1CS = .TRUE.

  keep_all_chr = .FALSE.
  keep_all_bp = .FALSE.
  sex_chr_chck = .TRUE.           ! If true, males will be excluded if from analysis of sex chromosomes
  excl_male = .FALSE.
  excl_fema = .FALSE.
  exclude_invalid_rs = .FALSE.
  excluded_chr = .FALSE.
  !excluded_chr(special_chr) = .TRUE.
  assume_male = .FALSE.
  assume_fema = .FALSE.
  assume_status = .FALSE.
  var_poststand = .FALSE.
  poststand_by_varbound = .FALSE.
  !first_pair_signif = .FALSE.
  nonpar_probs = .FALSE.
  out_level_S1 = .TRUE.
  var_grouped = .FALSE.
  cov_grouped = .FALSE.
  save_full_fam = .FALSE.
  do_out_pss = .TRUE.
  ntests_limit = .FALSE.
  do_out_actual_prev = .FALSE.

  auto_level_use = .FALSE.
  auto_level_report = .FALSE.
  auto_level_DS = .FALSE.
  auto_level_merge = .TRUE.
  auto_level_maf = .TRUE.
  auto_level_single = .TRUE.
  auto_level_weight_by_counts = .TRUE.
  
  auto_delta = .FALSE.
  auto_d_incapable = .FALSE.
  auto_d_optimize = .TRUE.
  auto_d_oracle = .FALSE.
  auto_d_maf_fixed = .TRUE.
  user_delta = .FALSE.
  simple_delta = .FALSE.
  
  auto_maf1 = .FALSE.
  auto_maf2 = .FALSE.
  auto_maf3 = .FALSE.

  cntrgrp = 1

  ntest_limit = -1
  ntest_limit2 = -1
  nwarnings = 0
  nerrors = 0
  nwarnings_limit = 25                  ! If there are too many warnings, stop printing them to the standard output
  nerrors_limit = 1000                  ! If there are too many errors, stop printing them to the standard output

  min_ss = hundred              ! If any test works with less than this number of observations warning is printed
  min_co = fifty                  ! If any test works with less than this number of observations warning is printed
  min_ca = fifty                     ! If any test works with less than this number of observations warning is printed
  !max_samplesize = 200000               ! Specifies the upper bound for a sample size on input
  test_same_chr_S1 = .FALSE.            ! y/n or not loci on the same chromosome will be tested
  test_same_chr = .FALSE.
  doS2 = .TRUE.                         ! If false only pretests will be performed
  only_S1 = .FALSE.
  replication_mode = .FALSE.            ! If true, only pairs of neighboring loci are tested 1-2, 3-4, etc.
  reject_all_pretests = .FALSE.
  always_pretest = .TRUE.               ! WARNING: THIS SHOULD NOT BE CHANGED TO FALSE UNLESS 
                                        ! NECCESSARY MODIFICATIONS TO THE CODE IN epi_testing_mod ARE MADE
                                        ! IF THIS WAS FALSE AND PRETEST LEVEL WAS EQUAL TO 1 AS A CONSEQUENCE
                                        ! OF OPTIMAL PRETEST DETERMINATION, THE ADJUSTED SCORE WOULD
                                        ! WOULD HAVE TO BE REPLACED WITH CLASSICAL SCORE WHICH WOULD POSSIBLY
                                        ! HAVE TO BE CORRECTED BY THE TOTAL NUMBER OF TESTS, BECAUSE OTHERWISE
                                        ! IT WOULD NOT BE INDEPENDENT OF THE PRETEST LEVEL (WHICH IS BASED ON MAFS)
                                        ! IF PRETEST LEVEL IS 1 BECAUSE USER SET IT TO BE, THERE SHOULD BE NO PROBLEM
                                        ! AND SIMPLE REPLACEMENT OF AS WITH S WILL SUFFICE
  delta_bound = 0.99_dpp
  var_bound = ten**(-5)
  var_cor_AS = one
  WT1 = def_T1                     ! set default S1 test
  T1_df = 4                        ! default degrees of freedom for pretest chisquare distribution
  T1_df_co = 4                     ! default degrees of freedom for controls-4 chisquare test
  T1_df_po = 4                     ! default degrees of freedom for pool-4 chisquare test
  T2_df = 4
  WT2 = def_T2           ! determines which phase 2 test is used
  !level_S1 = NAneg
  level_S1 = def_level1
  auto_level_optim = NAneg
  level_S2 = def_level2           ! Default level of the (adjusted) score tests
  MTC = -one
  MTC_factor = one
  ntests_assumed = NAneg
  dAS = def_delta                   ! Determines splitting of controls into two subsamples for AS
  dDS = def_delta                   ! Determines splitting of controls into two subsamples for DS
  delta_start = NAv2
  delta_end = NAv2 
  delta_step = NAv2
  ndelta = 1 
  min_MAF = zero                        ! Determines the initial lower bound for minor allele frequency
  correct_all_cells = .FALSE.
  min_cell_cnt = zero                 ! Determines what the minimal cell count has to be to test that pair of loci
  cell_count_cor = zero
  min_marg_cnt = zero               ! Determines what the minimal marginal count has to be to test that pair of loci
  olim1 = two                           ! Report only chisquare test p-values smaller than this bound (must be smaller than NAp) 
  olim2 = two                           ! Report only score test p-values smaller than this bound (must be smaller than NAp) 
  olim3 = two                           ! Report only adjusted score test p-values smaller than this bound (must be smaller than NAp)
  var_use_beta1 = .FALSE.               ! Estimates of variance use either H0 or H1 estimates of beta
  var_indep = .FALSE.            ! y/n S1 variance will be computed assuming LE
  res_prev = -one
  
  !! Assign negative NA value to MatOpt
  MatOpt = NAneg
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  INITIALIZE DATA VARIABLES  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL ResizeVar(ped_file, 1)
  CALL ResizeVar(bed_file, 1)
  CALL ResizeVar(map_file, 1)
  CALL ResizeVar(sub_file, 1)
  CALL ResizeVar(tmp_file, 1)
  CALL ResizeVar(res_file, 1)
  CALL ResizeVar(sort_file, 1)

  ped_nfiles = 0
  map_nfiles = 0
  bed_nfiles = 0
  fam_nfiles = 0
  pss_nfiles = 0
  sts_nfiles = 0
  sub_nfiles = 0
  out_nfiles = 0
  tmp_nfiles = 0

  map_nlines = 0
  sub_nlines = 0

  !! INITIALIZE filenames to empty strings
  ped_file      = ""
  bed_file      = ""
  map_file      = ""
  fam_file      = ""
  pss_file      = ""
  sts_file      = ""
  file_base     = ""
  out_file      = ""
  sort_file     = ""
  tmp_file      = ""
  log_file      = ""
  sub_file      = ""
  res_file      = ""
  out_pss_file  = ""
  out_maf_file  = ""
  out_prev_file = ""
  save_file     = ""
  save_ped_file = ""
  save_map_file = ""
  save_fam_file = ""
  
  def_file_base = "out_Epi"

  RETURN

END SUBROUTINE InitGlobalVars

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE InitCounters()
  IMPLICIT NONE

  !tmp_nlines = -1  

  !! Initialize counters
  NumErrs               = 0
  NumS1             = 0
  NumS1OK           = 0
  NumS1NotOK        = 0
  NumS1_sex         = 0
  NumS2             = 0
  NumS2AS4          = 0
  NumS2AS4OK        = 0
  NumS2AS4NotOK     = 0
  NumS2AS4PostStand = 0
  NumS2AS1          = 0
  NumS2AS1OK        = 0
  NumS2AS1NotOK     = 0
  NumS2AS1PostStand = 0
  NumS2DS           = 0
  NumS2DSOK         = 0
  NumS2DSNotOK      = 0
  NumS2CS           = 0
  NumS2CSOK         = 0
  NumS2CSNotOK      = 0
  NumS2_sex         = 0
  NumLinDep          = 0
  NumVarErr        = 0
  NumLowVar        = 0
  NumNegVar        = 0
  NumLowCell       = 0
  NumLowCellNumCor = 0
  NumLowMarg     = 0
  NumLowSS      = 0
  NumLowSS1     = 0
  NumLowSS2     = 0
  NumLD                 = 0
  NumEpi                = 0
  NumEpiCorr            = 0

  !! Initialize exclusion counters
  n_chr_excl       = 0
  n_rs_excl        = 0
  n_negbp_excl     = 0
  n_sub_excl       = 0
  n_maf_excl       = 0
  n_maxntests_excl = 0
  n_samechr_incl   = 0
  n_samechr_excl   = 0 
  
  MeanSS = zero
  
  varAS_sum   = zero
  varAS_sum_n = zero
  
  RETURN

END SUBROUTINE InitCounters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AskForCmdline(cmdline, exit_code)
  IMPLICIT NONE
  CHARACTER(mcl), INTENT(OUT) :: cmdline
  INTEGER, INTENT(OUT)        :: exit_code
  INTEGER                     :: ios

  !! Ask for input from user
  WRITE(*,*) ""
  3010 WRITE(*,'(A)', ADVANCE='NO') &
    "> Enter additional command line arguments (0 to exit): "
  READ(*,'(A)', IOSTAT=ios) cmdline
  
  IF(LEN_TRIM(cmdline)==0) THEN
    GOTO 3010
  ELSEIF(TRIM(cmdline)=="0") THEN
      exit_code = -1
  ELSE
    IF(FindLastSubstr(cmdline, "-cmdfile")>0) THEN
      exit_code = 1
    ELSE
      exit_code = 2
    ENDIF
  ENDIF
  
  RETURN

END SUBROUTINE AskForCmdline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE CheckCmdfile(cmdfile, cmd_present, lines, nlines, comment, userline)
  IMPLICIT NONE
  CHARACTER, INTENT(INOUT)                  :: comment
  CHARACTER(mcl), INTENT(INOUT)             :: userline
  CHARACTER(mfl), INTENT(OUT)               :: cmdfile
  LOGICAL, INTENT(OUT)                      :: cmd_present
  CHARACTER(mcl), ALLOCATABLE, INTENT(OUT)  :: lines(:)
  INTEGER, INTENT(OUT)                      :: nlines
  CHARACTER(mfl), ALLOCATABLE               :: args(:), args1(:)
  CHARACTER(mcl)                            :: cmdline
  INTEGER                                   :: nargs, i, ios
  INTEGER                                   :: nwords

  10 CONTINUE
  
  cmdfile = ""
  cmd_present = .FALSE.
  nlines = 0
  
  !! Count the number of arguments from command line
  nargs = COMMAND_ARGUMENT_COUNT()
  IF(nargs==0) THEN
    !! Print error announcement
    CALL PrntE("No command line arguments specified.", skip1=1)
    CALL Prnt("Use --help to see the list of possible arguments.", &
              Q=.TRUE., premature=.FALSE.)
  ENDIF

  !! Get arguments from command line
  CALL ResizeVar(args, nargs) 
  DO i=1,nargs 
    CALL GET_COMMAND_ARGUMENT(i, args(i))
  ENDDO

  !! Add the command line from userline
  IF(LEN_TRIM(userline)>0) THEN
    CALL SplitString(userline, char(32), args1, nwords)
    CALL ResizeVar(args, nargs+nwords)
    args(nargs+1:nargs+nwords) = args1
  ENDIF
  
  !! Total number of arguments
  nargs = SIZE(args)
  
  IF(nargs<=1) RETURN
  
  !! Check for presence
  DO i=1,nargs 
    SELECT CASE (args(i))
      CASE ("--no-openfile-iostat")
        no_openfile_iostat = .TRUE.
    END SELECT
  ENDDO
  
  !! Check for values
  DO i=2,nargs 
    SELECT CASE (args(i-1))
      CASE ("--ignore-unknown-args")
        stop_on_unknown_args = .FALSE.
      CASE ("--comment")
        comment = TRIM(ADJUSTL(args(i)))
      CASE ("--cmd", "--cmdfile", "--cmd-file", "--commandfile", "--command-file")
        cmdfile = TRIM(ADJUSTL(args(i)))
        cmd_present = .TRUE.
    END SELECT
  ENDDO
  DEALLOCATE(args)
  
  !! If --cmdfile found then read the entire cmd file
  IF(cmd_present) THEN
    CALL CheckFileExistence(cmdfile)
    CALL OpenFile(UNIT=ucmd, FILE=cmdfile, ACTION='R', STATUS='O', &
                  FORM='F', POSITION="R", ACCESS='S', IOSTAT=ios)

    DO WHILE (ios == 0)
      cmdline = ""
      !READ(cmdfile_unit,'(A)',ADVANCE='NO',SIZE=ls,IOSTAT=ios) cmdline
      READ(ucmd,'(A)', ADVANCE='NO', IOSTAT=ios) cmdline

      !! If end-of-line error, remove the error
      IF(ios==-2) ios = 0
      !IF(ios==-2 .AND. LEN_TRIM(cmdline)>0) ios = 0

      !! Exit if error
      IF(ios/=0) EXIT

      IF(cmdline(1:1)==comment) CYCLE
      IF(LEN_TRIM(cmdline)>0) THEN
        nlines = nlines+1
        CALL ResizeVar(lines, nlines)
        lines(nlines) = cmdline
      ENDIF

    ENDDO

    CLOSE(ucmd)
  ENDIF

  !! If --cmdfile found but empty then ask to terminate
  IF(cmd_present .AND. nlines==0) THEN
    CALL PrntE("Command file ["//TRIM(cmdfile)//"] is empty or"//&
                    " inaccessible.")
    CALL Terminate(ask=ask_halt, premature=.TRUE.)
    !! If still running, try again
    GOTO 10
  ENDIF

  RETURN

END SUBROUTINE CheckCmdfile
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetCmdLineArgs(cmdline, nrun, nrun_delta, nrun_OR, exit_code)
  IMPLICIT NONE
  
  CHARACTER(mcl), INTENT(IN)  :: cmdline
  INTEGER, INTENT(IN)         :: nrun
  INTEGER, INTENT(INOUT)      :: nrun_delta, nrun_OR
  INTEGER, INTENT(OUT)        :: exit_code
  CHARACTER(mal), ALLOCATABLE :: args(:)
  CHARACTER(mal)              :: arg, parg, carg, narg
  CHARACTER(mstl)             :: ped_cs1, map_cs1, sub_cs1, fam_cs1, out_cs1
  CHARACTER(mfl)              :: pedfile_list, mapfile_list, res_filelist, &
                                 subfile_list
  REAL(dpp)                   :: level1, level2, delta1, delta2, &
                                 olim1n, olim2n, olim3n, &
                                 tempcycle1, minmf, ssfac1, mf1, mf2, mf3
  LOGICAL                     :: show_help, carg_found, parg_found, cond1, &
                                 cond2, olim1_auto, olim2_auto, olim3_auto
  INTEGER                     :: i, j, nargs, nargsfound, nargsunknown, &
                                 nwords, ios, pedmaxncols1, WT11, WT21, &
                                 nthreads1, nloci1, & !nneutralloci1, & 
                                 nco1, nca1, ss1, T1_df1
  CHARACTER(mltl)             :: text

  exit_code = 0

  !! Initialize starting values of local option variables
  mapfile_list = ""
  pedfile_list = ""
  subfile_list = ""
  res_filelist = ""
  ped_cs1  = ""
  map_cs1  = ""
  sub_cs1  = ""
  fam_cs1  = ""
  out_cs1  = ""
  
  show_help    = .FALSE.
  nthreads1    = -1
  WT11         = -1
  T1_df1   = -1
  WT21         = -1
  tempcycle1   = -1
  pedmaxncols1 = -1
  nloci1       = -1
  ss1          = -1
  nco1         = -1
  nca1         = -1
  !nneutralloci1 = 0
  ssfac1 = -one
  level1 = -one
  level2 = -one
  delta1 = -one
  delta2 = -one
  olim1n = -one
  olim2n = -one
  olim3n = -one
  minmf  = -one
  mf1    = -one
  mf2    = -one
  mf3    = -one
  olim1_auto = .TRUE.
  olim2_auto = .FALSE.
  olim3_auto = .FALSE.

  !! Analysis: columns in the file
  cT      = 0
  cTpval  = 0
  cAS     = 0
  cASpval = 0

  ios = 0
     
  !! Read command line arguments     
  !nargs = iargc()
  nargs = COMMAND_ARGUMENT_COUNT()
  
  !! Allocate an array to store command line arguments
  CALL ResizeVar(args_all, nargs)
  
  !! Read arguments from the command line
  DO i=1,nargs 
    CALL GET_COMMAND_ARGUMENT(i, args_all(i))
  ENDDO

  !! Add the command line from cmd_file and those added by user 
  !! in the second and latter run
  IF(LEN_TRIM(cmdline)>0) THEN
    CALL SplitString(cmdline, char(32), args, nwords)
    CALL ResizeVar(args_all, nargs+nwords)
    IF(priority_user_cmdline) THEN
      args_all(nwords+1:nwords+nargs) = args_all(1:nargs)
      args_all(1:nwords) = args(1:nwords)
    ELSE
      args_all(1:nargs) = args_all(1:nargs)
      args_all(nargs+1:nwords+nargs) = args(1:nwords)
    ENDIF
    nargs = nargs + nwords
  ENDIF
  
  !! Used and unknown arguments
  nargsfound = 0
  nargsunknown = 0
  nunknown_args = nargs
  CALL ResizeVar(args_used, nargs, newvalue = .FALSE.)
  CALL ResizeVar(args_selected, nargs, 2, newvalue = "")
  CALL ResizeVar(args_unknown, nargs, newvalue = "")
  
  !! Search for known arguments
  parg = ""
  DO i=1,nargs
    arg = ADJUSTL(args_all(i))
    carg = arg
    narg = " "
    IF(i<nargs) narg = ADJUSTL(args_all(i+1))
    
    !! Turn switchers to lower case
    IF(.NOT.case_sensitive_flags) THEN
      CALL lowcase(carg)    ! currently read argument
      CALL lowcase(parg)    ! previously read argument
      CALL lowcase(narg)    ! next argument in line
    ENDIF
    
    !! CHECK FOR PRESENCE ONLY
    carg_found = .TRUE.
    SELECT CASE (carg)
      CASE ("--help", "--h")
        show_help = .TRUE.
      CASE ("--ignore-unknown-args")
        stop_on_unknown_args = .FALSE.
      CASE ("--no-halt")
        no_halt = .TRUE.
      CASE ("--halt")
        no_halt = .FALSE.
        ask_halt = .FALSE.
      CASE ("--silent")
        silent = .TRUE.
      CASE ("--overwrite")
        overwrite_files = .TRUE.
      CASE ("--no-overwrite")
        overwrite_files = .FALSE.
      CASE ("--keeptemp", "--keep-temp")
        keep_temp = .TRUE.
      CASE ("--deletetemp", "--delete-temp")
        keep_temp = .FALSE.
      CASE ("--nocountdown", "--no-countdown")
        countdown = .FALSE.
      CASE ("--noprogress", "--no-progress")
        countdown = .FALSE.
      CASE ("--nooutput", "--no-output")
        no_output = .TRUE.
      CASE ("--no-data-summary")
        no_data_summary = .TRUE.
      CASE ("--noreport", "--no-report")
        no_testing_report = .TRUE.
      CASE ("--show-report")
        no_testing_report = .FALSE.
      CASE ("--simulateinput", "--simulate-input", "--simulate", "--simulate-data")
        simulate_data = .TRUE.
        test_same_chr_S1 = .TRUE.
        test_same_chr = .TRUE.
        do_out_pss = .FALSE.
        input_format = 1
        out_loc_info = .FALSE.
      CASE ("--simulate-HWE", "--HWE")
        simulate_HWE = .TRUE.
      CASE ("--simulate-no-HWE", "--no-HWE", "--simulate-HWD", "--HWD")
        simulate_HWE = .FALSE.
      CASE ("--sim-single-chr")
        sim_single_chr = .TRUE.
      CASE ("--sim-multiple-chr")
        sim_single_chr = .FALSE.
      CASE ("--reusedata", "--reuse-data")
        reuse_data = .TRUE.
      CASE ("--save-input", "--saveinput", "--save-data", "--savedata", "--save")
        save_input_data = .TRUE.
      CASE ("--save-binary", "--savebinary", "--save-bed", "--make-bed")
        save_binary = .TRUE.
        save_input_data = .TRUE.
      CASE ("--save-plain", "--save-ped", "--make-ped")
        save_binary = .FALSE.
        save_input_data = .TRUE.
        save_full_fam = .TRUE.
      CASE ("--save-full-fam")
        save_full_fam = .TRUE.
      CASE ("--save-snp-major")
        out_bed_data_SNPmajor = .TRUE.
        save_input_data = .TRUE.
      CASE ("--save-indiv-major")
        out_bed_data_SNPmajor = .FALSE.
        save_input_data = .TRUE.
      CASE ("--simul-fixed-maf", "--fixed-maf", "--maf-fixed", "--fixedmaf")
        fixed_MAF = .TRUE.
      CASE ("--save-maf", "--freq")
        do_out_maf = .TRUE.
        compute_maf = .TRUE.
      CASE ("--out-pss", "--save-pss")
        do_out_pss = .TRUE.
      CASE ("--no-out-pss", "--no-save-pss")
        do_out_pss = .FALSE.
      CASE ("--out-prev", "--save-prev")
        do_out_actual_prev = .TRUE.
      CASE ("--do-test", "--do-testing")
        no_testing = .FALSE.
      CASE ("--no-test", "--no-testing", "--no-tests", "--notest")
        no_testing = .TRUE.
      CASE ("--no-pretest", "--no-S1", "--noS1")
        doS1 = .FALSE.
        doPTAS4 = .FALSE.
        doPTAS1 = .FALSE.
        doPTDS4 = .FALSE.
        doPTDS1 = .FALSE.
        doPO4 = .FALSE.
        doPO1 = .FALSE.
        doPT4 = .FALSE.
        doPT1 = .FALSE.
      CASE ("--dopretest", "--do-pretest", "--do-S1", "--doS1")
        doS1 = .TRUE.
      CASE ("--regression-divided")
        regres_divided = .TRUE.
      CASE ("--regression-not-divided")
        regres_divided = .FALSE.
      CASE ("--replication-mode")
        replication_mode = .TRUE.
        test_same_chr_S1 = .TRUE.
        report_all = .TRUE.
        !rs_valid_chck = .TRUE.
      CASE ("--random-subsample")
        fix_subsamp = .FALSE.
      CASE ("--fixed-subsample")
        fix_subsamp = .TRUE.
      CASE ("--reject-all", "--reject-all-pretests", "--rejectall")
        reject_all_pretests = .TRUE.
      CASE ("--only-pretest", "--only-S1")
        only_S1 = .TRUE.
      CASE ("--no-maintest", "--no-S2")
        doS2 = .FALSE.
        WT21 = 0
      CASE ("--do-maintest")
        doS2 = .TRUE.
      CASE ("--auto-delta-oracle")
        auto_d_oracle = .TRUE.
      CASE ("--auto-level")
        auto_level_use = .TRUE.
        auto_level_report = .TRUE.
      CASE ("--report-auto-level")
        auto_level_report = .TRUE.
      CASE ("--no-auto-level")
        auto_level_report = .FALSE.
        auto_level_use = .FALSE.
      CASE ("--auto-level-single")
        auto_level_single = .TRUE.
      CASE ("--auto-level-each")
        auto_level_single = .FALSE.
      CASE ("--auto-level-weight")
        auto_level_weight_by_counts = .TRUE.
      CASE ("--auto-level-noweight")
        auto_level_weight_by_counts = .FALSE.
      CASE ("--auto-level-adjusted")
        auto_level_DS = .FALSE.
      CASE ("--auto-level-disjoint")
        auto_level_DS = .TRUE.
      CASE ("--auto-level-merged")
        auto_level_merge = .TRUE.
      CASE ("--auto-level-split")
        auto_level_merge = .FALSE.
      CASE ("--auto-level-pretest-maf")
        auto_level_maf = .TRUE.
      CASE ("--auto-level-all-data-maf")
        auto_level_maf = .FALSE.
      CASE ("--analyze")
        run_analysis = .TRUE.
      CASE ("--analyze-as-one-file")
        analyze_files_together = .TRUE.
      CASE ("--out-append")
        out_append = .TRUE.
      CASE ("--outnolegend", "--out-no-legend")
        out_legend = .FALSE.
      CASE ("--outnoheader", "--out-no-header")
        out_header = .FALSE.
      CASE ("--outheadercommented", "--out-header-commented")
        out_header_comment = .TRUE.
      CASE ("--outheadernamed", "--out-header-named")
        out_header_names = .TRUE.
      CASE ("--out-debug")
        out_debug_info = .TRUE.
      CASE ("--out-variance")
        out_variance = .TRUE.
      CASE ("--out-alleles")
        out_alleles = .TRUE.
      CASE ("--out-zip", "--out-zipped")
        out_zipped = .TRUE.
      CASE ("--out-sort", "--out-sorted")
        out_sort = .TRUE.
        out_sort_multiple_files = .TRUE.
      CASE ("--out-no-sort", "--out-unsorted")
        out_sort = .FALSE.
      CASE ("--out-sort-only")
        out_sort = .TRUE.
        out_sort_tmp = .FALSE.
        out_sort_multiple_files = .TRUE.
        out_discard_unsorted = .FALSE.
        out_available = .TRUE.
        no_data_input = .TRUE.
        no_testing = .TRUE.
      CASE ("--tmp-sort-only")
        out_sort = .TRUE.
        out_sort_tmp = .TRUE.
        out_sort_multiple_files = .TRUE.
        out_available = .TRUE.
        no_data_input = .TRUE.
        no_testing = .TRUE.
        keep_temp = .TRUE.
      CASE ("--out-sort-keep-old")
        out_discard_unsorted = .FALSE.
      CASE ("--out-interaction", "--out-interactions")
        out_epi_effect = .TRUE.
      CASE ("--out-no-interaction", "--out-no-interactions")
        out_epi_effect = .FALSE.
      CASE ("--out-loc","--out-location-info")
        out_loc_info = .TRUE.
      CASE ("--out-maf")
        out_maf = .TRUE.
      CASE ("--out-maf-coca")
        out_maf_coca = .TRUE.
      CASE ("--out-maf-all")
        out_maf = .TRUE.
        out_maf_coca = .TRUE.
      CASE ("--no-out-maf")
        out_maf = .FALSE.
      CASE ("--no-out-maf-coca")
        out_maf_coca = .FALSE.

      CASE ("--out-debug-all", "--out-all-info")
        out_epi_effect = .TRUE.
        out_loc_info = .TRUE.
        out_maf = .TRUE.
        out_alleles = .TRUE.
        out_debug_info = .TRUE.
        out_S1_vector = .TRUE.
        out_testid = .TRUE.
        out_statistic = .TRUE.
        out_errcode = .TRUE.
        out_ss = .TRUE.
        out_pvalcor = .TRUE.
        out_mtc = .TRUE.

      CASE ("--out-basic")
        out_testid = .FALSE.
        out_statistic = .FALSE.
        out_pvalcor = .FALSE.
        out_mtc = .FALSE.
        out_debug_info = .FALSE.
        out_epi_effect = .FALSE.
        out_loc_info = .FALSE.
        out_alleles = .FALSE.
        out_level_S1 = .FALSE.
        out_debug_info = .FALSE.
        out_S1_vector = .FALSE.
        doPTDS4 = .FALSE.
        doPTDS1 = .FALSE.
        doPO4CS = .FALSE.
        doPO1CS = .FALSE.
        doPT4 = .FALSE.
        doPT1 = .FALSE.

      CASE ("--out-minimal", "--out-minimalistic")
        out_minimal = .TRUE.
        out_testid = .FALSE.
        out_legend = .FALSE.
        out_statistic = .FALSE.
        out_errcode = .FALSE.
        out_errcode_glob = .TRUE.
        out_ss = .FALSE.
        out_pvalcor = .FALSE.
        out_mtc = .FALSE.
        out_debug_info = .FALSE.
        out_epi_effect = .FALSE.
        !out_loc_info = .FALSE.
        out_maf = .FALSE.
        out_maf_coca = .FALSE.
        out_alleles = .FALSE.
        out_level_S1 = .FALSE.
        out_S1_vector = .FALSE.
        ndig_pval = 2
        ndig_mafs = 2
        doPO4CS = .FALSE.
        doPO1CS = .FALSE.
        doPT4 = .FALSE.
        doPT1 = .FALSE.
        
      CASE ("--out-statistic")
        out_statistic = .TRUE.
      CASE ("--out-errcode", "--out-error-code")
        out_errcode = .TRUE.
      CASE ("--out-global-errcode", "--out-global-error-code")
        out_errcode = .TRUE.
        out_errcode_glob = .TRUE.
      CASE ("--out-indiv-errcode")
        out_errcode = .TRUE.
        out_errcode_glob = .FALSE.

      CASE ("--testsamechr", "--test-same-chr", "--test-same-chromosomes", &
            "--same-chr", "--pretestsamechr", "--pretest-same-chr", &
            "--pretest-same-chromosomes")
        test_same_chr_S1 = .TRUE.
        test_same_chr = .TRUE.
      CASE ("--skip-pretest-same-chr")
        test_same_chr_S1 = .FALSE.
      CASE ("--skip-all-same-chr", "--skip-same-chr")
        test_same_chr_S1 = .FALSE.
        test_same_chr = .FALSE.
      CASE ("--nolog", "--no-log")
        no_log = .TRUE.
        
      CASE ("--do-AS", "--adjusted", "--reportadjusted", "--report-adjusted", "--get-adjusted", "--report-AS")
        doAS4 = .TRUE.
        doAS1 = .TRUE.
      CASE ("--do-AS4", "--adjusted4", "--reportadjusted4", "--report-adjusted4", "--get-adjusted4", "--report-AS4")
        doAS4 = .TRUE.
      CASE ("--do-AS1", "--adjusted1", "--reportadjusted1", "--report-adjusted1", "--get-adjusted1", "--report-AS1")
        doAS1 = .TRUE.
      CASE ("--do-DS", "--disjoint", "--reportdisjoint", "--report-disjoint", "--get-disjoint", "--report-DS")
        doDS4 = .TRUE.
        doDS1 = .TRUE.
      CASE ("--do-DS4", "--disjoint4", "--reportdisjoint4", "--report-disjoint4", "--get-disjoint4", "--report-DS4")
        doDS4 = .TRUE.
      CASE ("--do-DS1", "--disjoint1", "--reportdisjoint1", "--report-disjoint1", "--get-disjoint1", "--report-DS1")
        doDS1 = .TRUE.
      CASE ("--do-PO", "--report-pooled")
        doPO4 = .TRUE.
        doPO1 = .TRUE.
      CASE ("--do-CS", "--reportscore", "--report-score", "--out-score", "--reportclassical", "--report-classical", "--report-CS")
        doCS = .TRUE.
      CASE ("--no-adjusted", "--no-AS", "--noAS")
        doAS4 = .FALSE.
        doAS1 = .FALSE.
      CASE ("--no-adjusted4", "--no-AS4", "--noAS4")
        doAS4 = .FALSE.
      CASE ("--no-adjusted1", "--no-AS1", "--noAS1")
        doAS1 = .FALSE.
      CASE ("--no-disjoint", "--no-DS", "--noDS")
        doDS4 = .FALSE.
        doDS1 = .FALSE.
      CASE ("--no-disjoint4", "--no-DS4", "--noDS4")
        doDS4 = .FALSE.
      CASE ("--no-disjoint1", "--no-DS1", "--noDS1")
        doDS1 = .FALSE.
      CASE ("--no-pooled", "--no-PO", "--noPO")
        doPO4 = .FALSE.
        doPO1 = .FALSE.
      CASE ("--no-pooled4", "--no-PO4", "--noPO4")
        doPO4 = .FALSE.
      CASE ("--no-pooled1", "--no-PO1", "--noPO1")
        doPO1 = .FALSE.
      CASE ("--no-score", "--no-CS", "--noCS", "--noscore", "--no-classical", "--noclassical")
        doCS = .FALSE.
      CASE ("--allscore", "--all-score", "--alwaysscore", "--all-classical", "--always-score")
        doCS_all = .TRUE.
        doCS = .TRUE.
      CASE ("--do-all-two-stage")
        doAS4 = .TRUE.
        doAS1 = .TRUE.
        doDS4 = .TRUE.
        doDS1 = .TRUE.
        doPO4 = .TRUE.
        doPO1 = .TRUE.
      CASE ("--do-all")
        doAS4 = .TRUE.
        doAS1 = .TRUE.
        doDS4 = .TRUE.
        doDS1 = .TRUE.
        doPO4 = .TRUE.
        doPO1 = .TRUE.
        doCS_all = .TRUE.
        doCS = .TRUE.
        
      CASE ("--nocentering", "--no-centering", "--center-by-none")
        cntrgrp = 0
      CASE ("--centerbyself", "--center-by-self")
        cntrgrp = 1
      CASE ("--centerbyother", "--center-by-other")
        cntrgrp = 2
      CASE ("--centerbyall", "--center-by-all")
        cntrgrp = 3
      CASE ("--gvar", "--groupvar", "--group-var", "--groupedvariance", "--grouped-variance")
        var_grouped = .TRUE.
      CASE ("--ngvar", "--nongroupvar", "--non-group-var", "--nongroupedvariance", "--non-grouped-variance")
        var_grouped = .FALSE.
      CASE ("--gcov", "--groupcov", "--group-cov", "--groupedcovariance", "--grouped-covariance")
        cov_grouped = .TRUE.
      CASE ("--ngcov", "--nongroupcov", "--non-group-cov", "--nongroupedcovariance", "--non-grouped-covariance")
        cov_grouped = .FALSE.
      CASE ("--nonparametric")
        nonpar_probs = .TRUE.
      CASE ("--reportall", "--report-all", "--outputall", "--output-all", "--out-all")
        report_all = .TRUE.
        report_errs = .TRUE.
      CASE ("--reportallresults", "--out-all-results")
        report_all = .TRUE.
      CASE ("--reporterrors", "--report-errors", "--out-errors")
        report_errs = .TRUE.
      CASE ("--hideerrors", "--hide-errors")
        report_errs = .FALSE.
      CASE ("--sexchrcheck", "--sex-chr-check")
        sex_chr_chck = .TRUE.
      CASE ("--nosexchrcheck", "--no-sex-chr-check")
        sex_chr_chck = .FALSE.
      CASE ("--includespecialchr", "--include-special-chr")
        excluded_chr(special_chr) = .FALSE.
      CASE ("--excludespecialchr", "--exclude-special-chr")
        excluded_chr(special_chr) = .TRUE.
      CASE ("--keepallchr", "--keep-all-chr")
        keep_all_chr = .TRUE.
      CASE ("--keepallbp", "--keep-all-bp")
        keep_all_bp = .TRUE.
      CASE ("--includesexchr", "--include-sex-chr")
        excluded_chr(sex_chr) = .FALSE.
      CASE ("--excludesexchr", "--exclude-sex-chr")
        excluded_chr(sex_chr) = .TRUE.
      CASE ("--includeunplaced", "--include-unplaced")
        excluded_chr(0) = .FALSE.
      CASE ("--excludeunplaced", "--exclude-unplaced")
        excluded_chr(0) = .TRUE.
      CASE ("--excludemales", "--exclude-males")
        excl_male = .TRUE.
      CASE ("--excludefemales", "--exclude-females")
        excl_fema = .TRUE.
      CASE ("--excludeinvalidrs", "--exclude-invalid-rs")
        exclude_invalid_rs = .TRUE.
      CASE ("--submap-all-rs", "--submap-all-names", "--submap-no-rs-check", "--submap-read-all")
        sub_read_all = .TRUE.
      CASE ("--submap-read-one")
        sub_read_all = .FALSE.
      CASE ("--submap-pairs", "--submap-make-pairs")
        sub_make_pairs = .TRUE.
        report_all = .TRUE.
        report_errs = .TRUE.
      CASE ("--submap-no-chr")
        sub_no_chr = .TRUE.
        sub_chcol = 0
        sub_rscol = 1
      CASE ("--rs-check")
        rs_valid_chck = .TRUE.
      CASE ("--no-rs-check")
        rs_valid_chck = .FALSE.
      CASE ("--assume-males")     
        assume_male = .TRUE.
        assume_fema = .FALSE.
      CASE ("--assume-females")
        assume_fema = .TRUE.
        assume_male = .FALSE.
      CASE ("--assume-status")
        assume_status = .TRUE.
      CASE ("--var-use-beta1", "--var-beta1")
        var_use_beta1 = .TRUE.
      CASE ("--var-use-beta0", "--var-beta0")
        var_use_beta1 = .FALSE.
      CASE ("--var-assume-indep")
        var_indep = .TRUE.
      CASE ("--no-var-assume-indep", "--var-assume-general")
        var_indep = .FALSE.
      CASE ("--poststandardize", "--post-standardize")
        var_poststand = .TRUE.
      CASE ("--poststandbybound", "--post-standbybound")
        poststand_by_varbound = .TRUE.
      CASE ("--poststandbymean", "--post-stand-by-mean")
        poststand_by_varbound = .FALSE.
      CASE ("--nocausalpair", "--no-causal-pair")
        no_causal_pair = .TRUE.
      CASE ("--nocommas", "--no-commas")
        add_commas = .FALSE.
      CASE ("--no-openfile-iostat")
        no_openfile_iostat = .TRUE.
      CASE ("--pause")
        pause_run = .TRUE.
      CASE ("--adddetails", "--add-details")
        add_settings_filename = .TRUE.
      CASE ("--addlongdetails", "--add-long-details")
        add_settings_filename = .TRUE.
        add_settings_filename_long = .TRUE.
      CASE ("--addtimestamp", "--add-time-stamp")
        append_ts = .TRUE.
      CASE ("--notimestamp", "--no-timestamp", "--no-ts")
        append_ts = .FALSE.
      CASE ("--addpid", "--add-pid")
        append_pid = .TRUE.
      CASE ("--nopid", "--no-pid")
        append_pid = .FALSE.
      CASE ("--allpretests", "--all-pretests")
        do_all_pretests = .TRUE.
      CASE ("--compress-pval")
        ndig_pval = -1
      CASE ("--tmp-to-out-only")
        tmp_to_out_only = .TRUE.
        keep_temp = .TRUE.
      CASE ("--fam-all-samples")
        fam_single_sample = .FALSE.
      CASE ("--fam-single-sample")
        fam_single_sample = .TRUE.

      CASE DEFAULT
        carg_found = .FALSE.
    END SELECT
    
    !! Save the argument if successfully identified
    IF(carg_found) THEN
      nargsfound = nargsfound + 1
      nunknown_args = nunknown_args - 1
      args_selected(nargsfound,1) = TRIM(arg)
      args_used(i) = .TRUE.
      !GOTO 111
    ENDIF

    !! ****************************** !!
    !!  CHECK FOR PRESENCE AND VALUE  !!
    !! ****************************** !!
    parg_found = .TRUE.
    SELECT CASE (parg)
      CASE ("--cmd", "--cmdfile", "--cmd-file", "--commandfile", "--command-file")
        !! Processed elsewhere, this time only register the command so that it
        !! does not show up as unknown argument
      CASE ("--randomseed", "--seed", "--rseed")
        CALL SplitString(arg, ":", useed_str, nseed)
        CALL ResizeVar(useed_int, nseed)
        DO j=1,nseed; useed_int(j) = c2i(useed_str(j)); ENDDO
      CASE ("--randomseed-omp", "--seed-omp", "--rseed-omp")
        CALL SplitString(arg, ":", useed_str, nseed)
        CALL ResizeVar(useed_par, nseed)
        DO j=1,nseed; useed_par(j) = c2i(useed_str(j)); ENDDO
      CASE ("--randomseedstart", "--seedstart", "--seed-start")
        seed_start = MAX(10, c2i(arg))
      CASE ("--randomseedconst", "--seedconst", "--seed-const")
        seed_const = MAX(10, c2i(arg))
      !! Input and output files
      CASE ("--max-runtime", "--max-time", "--maxruntime")
        max_runtime = MAX(zero, c2r(arg))
      CASE ("--max-runtime-hours", "--max-time-hours", "--max-runtime-h", "--max-time-h")
        max_runtime = MAX(zero, c2r(arg))
      CASE ("--max-runtime-minutes", "--max-time-minutes", "--max-runtime-m", "--max-time-m")
        max_runtime = MAX(zero, c2r(arg)) / 60
      CASE ("--max-runtime-seconds", "--max-time-seconds", "--max-runtime-s", "--max-time-s")
        max_runtime = MAX(zero, c2r(arg)) / 3600
      CASE ("--comment")
        comment = arg(1:1)
      CASE ("--nrepeat", "--nrepeats", "--numrepeat", "--numrepeats")
        IF(arg=="unlimited" .OR. arg=="max" .OR. arg=="auto") THEN
          max_reruns = MAX_REPEATS
        ELSE
          max_reruns = MAX(1,MIN(MAX_REPEATS, c2i(arg)))
        ENDIF
      CASE ("--nthreads", "--nthread")
        nthreads1 = c2i(arg)
        
      !! Input/Output filenames
      CASE ("--file")
        input_format = 2
        ped_file = TRIM(arg)
        map_file = TRIM(arg)
        fam_file = TRIM(arg)
        ped_nfiles = 1
        map_nfiles = 1
        fam_nfiles = 1
      CASE ("--ped")
        input_format = 2
        ped_file = TRIM(arg)
        ped_nfiles = 1
      CASE ("--ped-list")
        input_format = 2
        pedfile_list = TRIM(arg)
      CASE ("--map")
        input_format = 2
        map_file = TRIM(arg)
        map_nfiles = 1
      CASE ("--map-list")
        input_format = 2
        mapfile_list = TRIM(arg)

      CASE ("--bfile")
        input_format = 1
        bed_file = TRIM(arg)
        map_file = TRIM(arg)
        fam_file = TRIM(arg)
        bed_nfiles = 1
        map_nfiles = 1
        fam_nfiles = 1
      CASE ("--bed")
        input_format = 1
        bed_file = TRIM(arg)
        bed_nfiles = 1
      CASE ("--fam")
        fam_file = TRIM(arg)
        fam_nfiles = 1
      CASE ("--bim")
        map_file = TRIM(arg)
        map_nfiles = 1
        input_format = 1
      CASE ("--pss", "--in-pss")
        pss_file = TRIM(arg)
        pss_file_subsample = .TRUE.
        fix_subsamp = .TRUE.
        pss_nfiles = 1
        delta1 = half
      CASE ("--sts", "--status")
        sts_file = TRIM(arg)
        sts_nfiles = 1
      CASE ("--filter", "--extract", "--submap", "--submap-include")
        sub_file = TRIM(arg)
        sub_include = .TRUE.
        sub_nfiles = 1
      CASE ("--filter-out", "--exclude", "--submap-exclude")
        sub_file = TRIM(arg)
        sub_include = .FALSE.
        sub_nfiles = 1
      CASE ("--out", "--output")
        file_base = TRIM(arg)
        save_file = TRIM(arg)
        user_set_out = .TRUE.
      CASE ("--resultfile", "--resultsfile", "--result-file", "--results-file")
        res_file = TRIM(arg)
      CASE ("--resultfilelist", "--result-file-list")
        res_filelist = TRIM(arg)
      CASE ("--resultfilemaxlines", "--result-file-max-lines")
        result_maxlines = MAX(c2i(arg),0)
      CASE ("--savefile", "--sfile", "--save-file")
        save_file = TRIM(arg)
        save_input_data = .TRUE.
        user_set_save = .TRUE.
      !! Input file settings
      CASE ("--max-tests", "--max-ntests", "--maxntests", "--maxntest")
        IF(IsInteger(arg)) ntest_limit = lc2i(arg)
      CASE ("--max-tests-phase2", "--max-ntests-phase2", "--maxntests-phase2", &
            "--maxntest-phase2", "--maxntests2", "--maxntest2")
        IF(IsInteger(arg)) ntest_limit2 = lc2i(arg)
      CASE ("--max-out-size", "--max-output-file-size")
        max_out_size = MAX(1000, lc2i(arg))
      CASE ("--mapncols", "--bimncols" , "--map-ncols", "--bim-ncols")
        map_ncols = MAX(1, c2i(arg))
      CASE ("--submapncols", "--submap-ncols")
        sub_ncols = MAX(1, c2i(arg))
      CASE ("--submaprscol", "--submap-rs-col")
        sub_rscol = MAX(1, c2i(arg))
      CASE ("--submapchcol", "--submap-ch-col", "--submapchrcol", "--submap-chr-col")
        sub_chcol = MAX(1, c2i(arg))
      CASE ("--submap-list")
        subfile_list = TRIM(arg)
      CASE ("--famncols", "--fam-ncols")
        fam_ncols = MAX(1, c2i(arg))
      CASE ("--famfidcol", "--fam-fid-col")
        fam_fidcol = MAX(0, c2i(arg))
      CASE ("--famiidcol", "--fam-iid-col")
        fam_iidcol = MAX(0, c2i(arg))
      CASE ("--fampidcol", "--fam-pid-col")
        fam_pidcol = MAX(0, c2i(arg))
      CASE ("--fammidcol", "--fam-mid-col")
        fam_midcol = MAX(0, c2i(arg))
      CASE ("--famsexcol", "--fam-sex-col")
        fam_sexcol = MAX(0, c2i(arg))
      CASE ("--famstatuscol", "--fam-status-col", "--fam-sts-col")
        fam_stscol = MAX(1,c2i(arg))
      CASE ("--sep")
        ped_cs1 = TRIM(arg)
        map_cs1 = TRIM(arg)
        sub_cs1 = TRIM(arg)
        fam_cs1 = TRIM(arg)
      CASE ("--pedsep", "--ped-sep")
        ped_cs1 = TRIM(arg)
      CASE ("--mapsep", "--map-sep")
        map_cs1 = TRIM(arg)
      CASE ("--subsep", "--sub-sep", "--submap-sep")
        sub_cs1 = TRIM(arg)
      CASE ("--famsep", "--fam-sep")
        fam_cs1 = TRIM(arg)
      CASE ("--outsep", "--out-sep")
        out_cs1 = TRIM(arg)

      CASE ("--out-sort-by", "--out-sort-col")
        IF(arg=="0") THEN
          out_sort_col = 0        
          out_sort_sel = "0"
        ELSEIF(IsInteger(arg)) THEN
          out_sort_col = c2i(arg)
          out_sort_sel = ""
        ELSE
          out_sort_sel = upcasef(arg(1:1))
        ENDIF

      CASE ("--missing-genotype", "--missing-geno")
        CharNAgen = TRIM(arg)
      CASE ("--missing-sex")
        CharNAsex = TRIM(arg)
      CASE ("--missing-status")
        CharNAsts = TRIM(arg)
      CASE ("--statuscase")
        CharCa = TRIM(arg)
      CASE ("--statuscontrol")
        CharCo = TRIM(arg)
      CASE ("--sexmale")
        CharMa = TRIM(arg)
      CASE ("--sexfemale")
        CharFe = TRIM(arg)
      CASE ("--delim")
        delim = TRIM(arg)
      CASE ("--pedskip")
        ped_nskip = MAX(0,c2i(arg))
      CASE ("--mapskip")
        map_nskip = MAX(0,c2i(arg))
      CASE ("--submapskip")
        sub_nskip = MAX(0,c2i(arg))
      CASE ("--famskip")
        fam_nskip = MAX(0,c2i(arg))
      CASE ("--pssskip")
        pss_nskip = MAX(0,c2i(arg))
      CASE ("--stsskip")
        sts_nskip = MAX(0,c2i(arg))
      CASE ("--resultfileskip")
        result_nskip = MAX(0,c2i(arg))
      CASE ("--pedmaxncols")
        pedmaxncols1 = c2i(arg)

      CASE ("--samplesize")
        ss1 = MAX(0,c2i(arg))
      CASE ("--ncontrols", "--ncontrol", "--nco")
        nco1 = MAX(0,c2i(arg))
      CASE ("--ncases", "--ncase", "--nca")
        nca1 = MAX(0,c2i(arg))
      CASE ("--ssfactor")
        ssfac1 = MAX(zero,c2r(arg))
      CASE ("--min-ss", "--minss", "--minsamplesize", "--min-samplesize")
        min_ss = MAX(zero,c2r(arg))
      CASE ("--mincontrols", "--min-controls")
        min_co = MAX(zero,c2r(arg))
      CASE ("--mincases", "--min-cases")
        min_ca = MAX(zero,c2r(arg))
      CASE ("--nsamples")
        ped_nsamp = MAX(1,c2i(arg))
      CASE ("--tempcycle", "--temp-cycle", "--cycle")
        tempcycle1 = c2r(arg)

      CASE ("--ndig-stat", "--ndigit-stat", "--ndigits-stat")
        ndig_stat = MAX(0,c2i(arg))
      CASE ("--ndig-pval", "--ndigit-pval", "--ndigits-pval")
        ndig_pval =  MAX(0,c2i(arg))
      CASE ("--ndig-mafs", "--ndigit-mafs", "--ndigits-mafs", "--ndig-maf", "--ndigit-maf", "--ndigits-maf")
        ndig_mafs =  MAX(0,c2i(arg))

      !! Testing parameters
      CASE ("--S1-sample", "--S1-data")
        SELECT CASE (lowcasef(arg))
        CASE ("ca", "case", "cases", "affected")
          WT11 = T1ca
        CASE ("co", "control", "controls", "unaffected")
          WT11 = T1co
        CASE ("caco", "coca", "cases+controls", "both")
          WT11 = T1cc
        END SELECT
      CASE ("--pretest", "--pre-test", "-test1")
        WT11 = c2i(arg)
        
      CASE ("--method")
        SELECT CASE (arg)
        CASE ("all")
          doAS4 = .TRUE.
          doAS1 = .TRUE.
          doDS4 = .TRUE.
          doDS1 = .TRUE.
          doPO4 = .TRUE.
          doPO1 = .TRUE.
          doCS = .TRUE.
        CASE ("AS4")
          doAS4 = .TRUE.
          doAS1 = .FALSE.
          doDS4 = .FALSE.
          doDS1 = .FALSE.
          doPO4 = .FALSE.
          doPO1 = .FALSE.
          doCS = .FALSE.
        CASE ("AS1")
          doAS4 = .FALSE.
          doAS1 = .TRUE.
          doDS1 = .FALSE.
          doDS4 = .FALSE.
          doPO4 = .FALSE.
          doPO1 = .FALSE.
          doCS = .FALSE.
        CASE ("DS4")
          doAS4 = .FALSE.
          doAS1 = .FALSE.
          doDS1 = .FALSE.
          doDS4 = .TRUE.
          doPO4 = .FALSE.
          doPO1 = .FALSE.
          doCS = .FALSE.
        CASE ("DS1")
          doAS4 = .FALSE.
          doAS1 = .FALSE.
          doDS1 = .TRUE.
          doDS4 = .FALSE.
          doPO4 = .FALSE.
          doPO1 = .FALSE.
          doCS = .FALSE.
        CASE ("PO4")
          doAS4 = .FALSE.
          doAS1 = .FALSE.
          doDS1 = .FALSE.
          doDS4 = .FALSE.
          doPO4 = .TRUE.
          doPO1 = .FALSE.
          doCS = .FALSE.
        CASE ("PO1")
          doAS4 = .FALSE.
          doAS1 = .FALSE.
          doDS1 = .FALSE.
          doDS4 = .FALSE.
          doPO4 = .FALSE.
          doPO1 = .TRUE.
          doCS = .FALSE.
        CASE ("CS")
          doAS4 = .FALSE.
          doAS1 = .FALSE.
          doDS1 = .FALSE.
          doDS4 = .FALSE.
          doPO4 = .FALSE.
          doPO1 = .FALSE.
          doCS = .TRUE.
        END SELECT

      CASE ("--pretestdf", "--test1df")
        T1_df1 = c2i(arg)
      CASE ("--level1", "--pretestlevel", "--plevel", "--alpha1", "--a1")
        IF(TRIM(arg)=="auto") THEN
          auto_level_use = .TRUE.
          auto_level_report = .TRUE.
        ELSE
          level1 = c2r(arg)
          auto_level_use = .FALSE.
        ENDIF
      CASE ("--test", "--test2", "--maintest", "--main-test", "--posttest", "--post-test")
        WT21 = c2i(arg)
      CASE ("--level2")
        level2 = c2r(arg)
      CASE ("--assumedntests", "--assumentests", "--assumetests")
        ntests_assumed = c2r(arg)
      CASE ("--mtc1", "--correct1")
        MTC(1:3) = c2r(arg)
      CASE ("--mtc-AS")
        MTC(1:2) = c2r(arg)
      CASE ("--mtc-DS")
        MTC(3:4) = c2r(arg)
      CASE ("--mtc-PO")
        MTC(5:6) = c2r(arg)
      CASE ("--mtc-CS")
        MTC(7) = c2r(arg)
      CASE ("--mtc-factor")
        MTC_factor = c2r(arg)
      CASE ("--outreport1", "--out-report1", "--report1", "--output1")
        IF(TRIM(arg)=="auto") THEN
          olim1_auto = .TRUE.
          olim1n = -one
        ELSE
          olim1_auto = .FALSE.
          olim1n = c2r(arg)
        ENDIF
      CASE ("--outreport2", "--out-report2", "--report2", "--output2")
        IF(TRIM(arg)=="auto") THEN
          olim2_auto = .TRUE.
          olim2n = -one
        ELSE
          olim2_auto = .FALSE.
          olim2n = c2r(arg)
        ENDIF
      CASE ("--outreport3", "--out-report3", "--report3", "--output3")
        IF(TRIM(arg)=="auto") THEN
          olim3_auto = .TRUE.
          olim3n = -one
        ELSE
          olim3_auto = .FALSE.
          olim3n = c2r(arg)
        ENDIF
      CASE ("--delta")
        IF(TRIM(arg)=="auto-simple") THEN
          auto_delta = .TRUE.
          user_delta = .FALSE.
          simple_delta = .TRUE.
        ELSE IF(TRIM(arg)=="auto" .OR. TRIM(arg)=="auto-optim") THEN
          auto_delta = .TRUE.
          user_delta = .FALSE.
          simple_delta = .FALSE.
        ELSE
          delta1 = c2r(arg)
        ENDIF
      CASE ("--delta-simple-frac")
        auto_d_frac = c2r(arg)
      CASE ("--delta-DS", "--deltaDS")
        delta2 = c2r(arg)
      CASE ("--delta-start", "--delta-min")
        delta_start = MIN(one, MAX(zero, c2r(arg)))
      CASE ("--delta-end", "--delta-max")
        delta_end = MIN(one, MAX(zero, c2r(arg)))
      CASE ("--delta-increase", "--delta-incr", "--delta-add", "--delta-step")
        delta_step = MIN(one, MAX(-one, c2r(arg)))
        !delta_step = c2r(arg)

      CASE ("--model")
        IF(IsInteger(arg)) THEN
          ana_model = c2i(arg)
          sim_model = c2i(arg)
        ELSE
          ana_model = GetModelNumber(arg)
          sim_model = GetModelNumber(arg)
        ENDIF                                
      CASE ("--analysis-model", "--analysismodel", "--amodel")
        IF(IsInteger(arg)) THEN
          ana_model = c2i(arg)
        ELSE
          ana_model = GetModelNumber(arg)
        ENDIF
      CASE ("--analysis-model-1", "--analysismodel1", "--amodel1", "--amodel-1")
        IF(IsInteger(arg)) THEN
          ana_model_S1 = c2i(arg)
        ELSE
          ana_model_S1 = GetModelNumber(arg)
        ENDIF
      CASE ("--minmaf", "--min-maf")
        minmf = c2r(arg)
      CASE ("--center-group", "--group-center")
        cntrgrp = c2i(arg)
      CASE ("--maf")
        IF(TRIM(arg)=="auto") THEN
          auto_maf1 = .TRUE.
          auto_maf2 = .TRUE.
          auto_maf3 = .TRUE.
        ELSE
          mf1 = c2r(arg)
          mf2 = c2r(arg)
          mf3 = c2r(arg)
        ENDIF
      CASE ("--maf1")
        IF(TRIM(arg)=="auto") THEN
          auto_maf1 = .TRUE.
        ELSE
          mf1 = c2r(arg)
        ENDIF
      CASE ("--maf2")
        IF(TRIM(arg)=="auto") THEN
          auto_maf2 = .TRUE.
        ELSE
          mf2 = c2r(arg)
        ENDIF
      CASE ("--maf3")
        IF(TRIM(arg)=="auto") THEN
          auto_maf3 = .TRUE.
        ELSE
          mf3 = c2r(arg)
        ENDIF

      CASE ("--auto-delta-ntests")
        auto_d_ntests = c2i(arg)
      CASE ("--auto-delta-correction")
        auto_d_COR = c2r(arg)
      CASE ("--auto-delta-OR")
        auto_d_OR = c2r(arg)
      CASE ("--auto-delta-OR1")
        auto_d_OR1 = c2r(arg)
      CASE ("--auto-delta-OR2")
        auto_d_OR2 = c2r(arg)
      CASE ("--auto-delta-prev")
        auto_d_prev = c2r(arg)
      CASE ("--auto-delta-maf1")
        auto_d_maf1 = c2r(arg)
      CASE ("--auto-delta-maf2")
        auto_d_maf2 = c2r(arg)
      CASE ("--auto-delta-min")
        auto_d_min = c2r(arg)
      CASE ("--auto-delta-max")
        auto_d_max = c2r(arg)
      CASE ("--auto-delta-n")
        auto_d_n = c2i(arg)
      CASE ("--auto-delta-minlevel")
        auto_d_min = c2r(arg)
      CASE ("--auto-delta-maxlevel")
        auto_d_max = c2r(arg)
      CASE ("--auto-delta-nlevel")
        auto_d_n = c2i(arg)

      CASE ("--autolevelfactor")
        auto_level_factor = c2r(arg)
      CASE ("--autolevelmin")
        auto_level_min = MIN(MAX(c2r(arg), minimum_test_level), &
                                 one - minimum_test_level)
      CASE ("--autolevelmax")
        auto_level_max = MIN(MAX(c2r(arg), minimum_test_level), &
                                 maximum_test_level)
      CASE ("--autolevelprevalence")
        auto_level_prev = c2r(arg)
      CASE ("--autolevelOR1")
        auto_level_OR1 = MAX(c2r(arg), min_OR)
      CASE ("--autolevelOR2")
        auto_level_OR2 = MAX(c2r(arg), min_OR)
      CASE ("--autolevelOR")
        auto_level_OR = MAX(c2r(arg), min_OR)

      CASE ("--simulation-model", "--simulationmodel", "--smodel")
        IF(IsInteger(arg)) THEN
          sim_model = c2i(arg)
        ELSE
          sim_model = GetModelNumber(arg)
        ENDIF
      CASE ("--prevalence", "--prev")
        prevalence = c2r(arg)
      CASE ("--OR", "--RR")
        OR = MAX(min_OR, c2r(arg))
      CASE ("--OR-min", "--RR-min", "--OR-start")
        OR_min = MAX(min_OR, c2r(arg))
      CASE ("--OR-max", "--RR-max", "--OR-end")
        OR_max = MAX(zero, c2r(arg))
      CASE ("--OR-increase", "--RR-increase", "--OR-incr", "--RR-incr", &
            "--OR-add", "--RR-add", "--OR-step", "--RR-step")
        OR_step = c2r(arg)
      CASE ("--OR1", "--RR1")
        OR1 = MAX(min_OR,c2r(arg))
      CASE ("--OR2", "--RR2")
        OR2 = MAX(min_OR,c2r(arg))
      CASE ("--beta")
        OR = EXP(c2r(arg))
      CASE ("--beta1")
        OR1 = EXP(c2r(arg))
      CASE ("--beta2")
        OR2 = EXP(c2r(arg))
      CASE ("--sim-chr-start")
        sim_chr_start = MIN(MAX(c2i(arg), 0), 26)
      CASE ("--sim-rs-start")
        sim_rs_start = MAX(c2i(arg), 0)
      CASE ("--sim-id-start")
        sim_fid_start = MAX(c2i(arg), 0)
        sim_iid_start = MAX(c2i(arg), 0)
      CASE ("--sim-fid-start")
        sim_fid_start = MAX(c2i(arg), 0)
      CASE ("--sim-iid-start")
        sim_iid_start = MAX(c2i(arg), 0)
      CASE ("--LD")                                                    
        LD = c2r(arg)
        simulate_HWE = .FALSE.
      CASE ("--nloci", "--nmarkers")
        nloci1 = MAX(0,c2i(arg))
      !CASE ("--nneutralloci")
      !  nneutralloci1 = MAX(0,c2i(arg))
      CASE ("--nLDpairs")
        nLDpairs = MAX(0,c2i(arg))
      CASE ("--minmargincount", "--min-margin-count")
        min_marg_cnt = MAX(zero,c2r(arg))
      CASE ("--mincellcount", "--min-cell-count")
        min_cell_cnt = MAX(zero,c2r(arg))
      CASE ("--cellcorrection", "--cell-correction")
        cell_count_cor = MAX(zero,c2r(arg))
      CASE ("--min-cell-correct-all")
        correct_all_cells = .TRUE.
      CASE ("--varbound", "--var-bound", "--min-var")
        var_bound = MAX(epstol,c2r(arg))
      CASE ("--varcorrection", "--var-correction", "--varcorr", "--var-corr")
        var_cor_AS = MAX(one,c2r(arg))
      CASE ("--colt1")
        cT = c2i(arg)
      CASE ("--colp1")
        cTpval = c2i(arg)
      CASE ("--colt2")
        cAS = c2i(arg)
      CASE ("--colp2")
        cASpval = c2i(arg)
      CASE ("--out-prev-res")
        res_prev = c2r(arg)
      CASE DEFAULT
        parg_found = .FALSE.
    END SELECT

    !! Save the argument if successfully identified
    IF(parg_found) THEN
      nargsfound = nargsfound + 1
      nunknown_args = nunknown_args - 2
      args_selected(nargsfound,1) = TRIM(parg)
      args_selected(nargsfound,2) = TRIM(arg)
      args_used(i-1) = .TRUE.
      args_used(i) = .TRUE.
    ENDIF
    
    !! Save the arguments if NOT successfully identified
    IF(.NOT.carg_found .AND. .NOT.parg_found) THEN
      nargsunknown = nargsunknown + 1
      args_unknown(nargsunknown) = TRIM(parg)//" "//TRIM(carg)
    ENDIF

    !111 CONTINUE 
    parg = arg
    
  ENDDO
  
  !! PROCESS SELECTED OPTIONS
#ifdef _OPENMP
  IF(nthreads1>0) NUMBER_OF_THREADS = nthreads1
#endif
  
  IF(silent) THEN
    countdown = .FALSE.
    append_ts = .TRUE.
    append_pid = .TRUE.
  ENDIF

  !! Remove unnecessary rows in args_selected 
  IF(nargsfound>0) THEN
    CALL ResizeVar(args_selected, nargsfound, 2)
  ELSE
    CALL PrntE("No valid command line arguments found!", skip1=1)
    CALL Prnt("Use --help to see the list of possible arguments.", &
                   Q=.TRUE., premature=.FALSE.)
  ENDIF
  
  !! Remove unnecessary space in args_unknown
  IF(nargsunknown>0) CALL ResizeVar(args_unknown, nargsunknown)

  !! Remove duplicate arguments and keep the last ones of the duplicates    
  !CALL RemoveDuplicates(args_selected)
  
  !! TEMPORARILLY, WARNINGS ARE TURNED OFF!
  IF(nunknown_args/=0 .AND. stop_on_unknown_args .AND. .FALSE.) THEN
    IF(nunknown_args > 0) PRINT '(A)',lt//"Command line error: Too few arguments."
    IF(nunknown_args < 0) PRINT '(A)',lt//"Command line error: Too many arguments."
    STOP
  ENDIF

  !! SHOW HELP IF SELECTED AND STOP
  IF(show_help) THEN
    no_log = .TRUE.
    CALL PrintHelpScreen()
    STOP
  ENDIF
  
  !! ======================================================================= !!
  !!                           DATA INPUT SETTINGS                           !!
  !! ======================================================================= !!

  !! ----------------------------------------------------------------------- !!
  !!                               INPUT FILES                               !!  
  !! ----------------------------------------------------------------------- !!

  !! Get filenames from user input
  CALL GetFilenames(ped_nfiles, uped, ped_file, pedfile_list, delim, comment)
  CALL GetFilenames(map_nfiles, umap, map_file, mapfile_list, delim, comment)
  CALL GetFilenames(sub_nfiles, umap, sub_file, subfile_list, delim, comment)

  !! User specified both ped and bed files -> QUIT
  IF(ped_nfiles>0 .AND. bed_nfiles>0) &
    CALL PrntE("Pedding file(s) (--ped, --pedlist) and binary ped"// &
               " file(s) (--bed) cannot be inputted simultaneously.", &
                   log=silent, Q=.TRUE., premature=.FALSE.)
  
  !! Warn about missing map and fam files
  IF(bed_nfiles>0 .AND. (map_nfiles==0 .OR. fam_nfiles==0)) THEN
    IF(map_nfiles==0) & 
      CALL PrntE("MAP file cannot be missing when BED file specified.", skip1=1)
    IF(fam_nfiles==0) & 
      CALL PrntE("PED file cannot be missing when BED file specified.", skip1=1)
    CALL Prnt0("", Q=.TRUE., premature=.FALSE.)
  ENDIF

  !! ----------------------------------------------------------------------- !!
  !!                          INPUT FILE SETTINGS                            !!  
  !! ----------------------------------------------------------------------- !!
  !! Make sure map_ncols or map_rscol are not too small
  IF(map_ncols<0) THEN
    IF(input_format==1) map_ncols = def_bim_ncols
    IF(input_format==2) map_ncols = def_map_ncols
  ENDIF
  
  IF(map_rscol<1) &
    CALL PrntE("RS column of a mapping file must be positive!", Q=.TRUE., &
               premature=.FALSE.)
  
  map_ncols = MAX(map_ncols, map_rscol, map_chcol)

  !! Make sure sub_ncols or sub_rscol are not too small
  IF(sub_rscol==1 .AND. sub_chcol==sub_rscol) sub_chcol = 0 
  IF(sub_rscol<1) &
    CALL PrntE("RS column of a submap file must be positive!", Q=.TRUE., &
               premature=.FALSE.)
  
  sub_ncols = MAX(sub_ncols, sub_rscol, sub_chcol)

  !! Make sure fam_fidcol, fam_iidcol, fam_sexcol or fam_stscol are not too large
  IF(fam_ncols == 1 .AND. fam_sexcol == def_sexcol) fam_sexcol = 0  
  text = ""
  IF(fam_ncols < fam_fidcol) text = TRIM(text)//" FID"
  IF(fam_ncols < fam_iidcol) text = TRIM(text)//" IID"
  IF(fam_ncols < fam_pidcol) text = TRIM(text)//" PID"
  IF(fam_ncols < fam_midcol) text = TRIM(text)//" MID"
  IF(fam_ncols < fam_sexcol) text = TRIM(text)//" SEX"
  IF(fam_ncols < fam_stscol) text = TRIM(text)//" STATUS"

  IF(LEN_TRIM(text)>0) &
    CALL PrntE("Given position of"//TRIM(text)//" column(s) is larger than"//&
               " the total number of columns.", Q=.TRUE., premature=.FALSE.)

  !! Column separators
  IF(sep_ascii) THEN
    IF(IsInteger(ped_cs1)) ped_cs = CHAR(c2i(ped_cs1))
    IF(IsInteger(map_cs1)) map_cs = CHAR(c2i(map_cs1))
    IF(IsInteger(sub_cs1)) sub_cs = CHAR(c2i(sub_cs1))
    IF(IsInteger(fam_cs1)) fam_cs = CHAR(c2i(fam_cs1))
  ENDIF

  IF(ped_cs1=="\t")    ped_cs = tab
  IF(ped_cs1=="t")     ped_cs = tab
  IF(ped_cs1=="tab")   ped_cs = tab
  IF(ped_cs1=="\s")    ped_cs = space
  IF(ped_cs1=="s")     ped_cs = space
  IF(ped_cs1=="space") ped_cs = space
  IF(map_cs1=="\t")    map_cs = tab
  IF(map_cs1=="t")     map_cs = tab
  IF(map_cs1=="tab")   map_cs = tab
  IF(map_cs1=="\s")    map_cs = space
  IF(map_cs1=="s")     map_cs = space
  IF(map_cs1=="space") map_cs = space
  IF(fam_cs1=="\t")    fam_cs = tab
  IF(fam_cs1=="t")     fam_cs = tab
  IF(fam_cs1=="tab")   fam_cs = tab
  IF(fam_cs1=="\s")    fam_cs = space
  IF(fam_cs1=="s")     fam_cs = space
  IF(fam_cs1=="space") fam_cs = space
  IF(sub_cs1=="\t")    sub_cs = tab
  IF(sub_cs1=="t")     sub_cs = tab
  IF(sub_cs1=="tab")   sub_cs = tab
  IF(sub_cs1=="\s")    sub_cs = space
  IF(sub_cs1=="s")     sub_cs = space
  IF(sub_cs1=="space") sub_cs = space
  
  !! Make sure there is not a conflict between 'subset mapping file being paired
  !! up' and whether 'markers listed in it are included or excluded'
  IF(sub_make_pairs) sub_include = .TRUE.

  !! Filesize variables
  IF(pedmaxncols1>0) ped_maxncol = pedmaxncols1

  !! Cycle for storing data into temporary files 
  IF(INT(tempcycle1)>0) temp_cycle = INT(tempcycle1)
  
  !! Save user given sample sizes into global variables
  IF(ss1 > 0) ped_ss = ss1
  IF(nco1 > 0)  ped_nco = nco1
  IF(nca1 > 0)  ped_nca = nca1

  !! Sample and subsample sizes
  IF(ped_nco>0 .AND. ped_nca>0) &
    ped_ss = ped_nco + ped_nca
  IF(ped_ss>0 .AND. ped_nco>0 .AND. ped_nca<0) &
    ped_nca = ped_ss - ped_nco
  IF(ped_ss>0 .AND. ped_nca>0 .AND. ped_nco<0) &
    ped_nco = ped_ss - ped_nca
  
  IF(ssfac1 > zero .AND. nrun>1) THEN
    reuse_data = .FALSE.
    !! ssfac1 is an integer that is to be added to both controls and cases
    IF(INT(ssfac1) == ssfac1 .AND. ssfac1 > 100*one) THEN
      ped_nco  = INT((nrun-1)*ssfac1 + ped_nco)  
      ped_nca     = INT((nrun-1)*ssfac1 + ped_nca)
    !! ssfac1 is a real quotient by which sample size is to be increased
    ELSE
      ped_nco  = INT((nrun-1)*ssfac1*ped_nco + ped_nco)  
      ped_nca     = INT((nrun-1)*ssfac1*ped_nca + ped_nca)
    ENDIF
    !! Adjust the total sample size
    ped_ss = ped_nco + ped_nca
  ENDIF
    
  !! Check for input MAF
  IF(minmf >= zero .AND. minmf <= half) min_MAF = minmf

  !! Adjust things when new data are to be simulated
  IF(simulate_data .AND. .NOT.reuse_data) THEN
    
    IF(no_causal_pair) ncausalloci = 0

    !! Check for missing or non-sensical values in nloci1
    IF(nloci1 < 1) &
      CALL PrntE("Total number of loci for simulation must be specified"//&
                      " and cannot be smaller than 1.", Q=.TRUE., &
                      premature=.FALSE.)
    IF(nloci1 < ncausalloci+2*nLDpairs) &
      CALL PrntE("Specified loci counts are inconsistent. Check the"//&
                      " values of --nloci, --nLDpairs and the presence of"//&
                      " --nocausalpair.", Q=.TRUE., premature=.FALSE.)

    !! Determine the number of neutral loci from nloci1 and nLDpairs
    IF(nloci1 < 2) THEN
      nneutralloci = nloci1
      no_causal_pair = .TRUE.
    ELSEIF(no_causal_pair) THEN
      nneutralloci = nloci1 - 2*nLDpairs
    ELSE
      nneutralloci = nloci1 - 2*nLDpairs - 2
    ENDIF
  
    !!! ALLELE FREQUENCIES SETTINGS
    !! Set allele frequencies for simulation
    IF(mf1>=min_MAF .AND. mf1<=half) allelefreqA = mf1
    IF(mf2>=min_MAF .AND. mf2<=half) allelefreqB = mf2
    IF(mf3>=min_MAF .AND. mf3<=half) allelefreqC = mf3
    
    !! If any not given, determine automatically
    IF(allelefreqA <= zero) auto_maf1 = .TRUE. 
    IF(allelefreqB <= zero) auto_maf2 = .TRUE. 
    IF(allelefreqC <= zero) auto_maf3 = .TRUE. 
  
    IF(allelefreqA > zero) auto_maf1 = .FALSE. 
    IF(allelefreqB > zero) auto_maf2 = .FALSE. 
    IF(allelefreqC > zero .OR. nneutralloci==0) auto_maf3 = .FALSE.
    
    !! If maf1 specified and maf3 is to be automatic, just put it equal to maf1
    IF(auto_maf3 .AND. .NOT.auto_maf1) THEN
      allelefreqC = allelefreqA
      auto_maf3 = .FALSE.
    ENDIF
  
    !! Automatically determine MAF based on sample size
    IF(auto_maf1 .OR. auto_maf2 .OR. auto_maf3) THEN
      auto_maf = (five/MAX(ped_ss/2, INT(min_ss)))**quarter
      auto_maf = round(MAX(min_MAF, MIN(one-min_MAF, auto_maf)), 2)
      IF(auto_maf1) allelefreqA = auto_maf
      IF(auto_maf2) allelefreqB = auto_maf
      IF(auto_maf3) allelefreqC = auto_maf
    ENDIF
  
  ENDIF
  
  !! ----------------------------------------------------------------------- !!
  !!                        TESTS, LEVELS, DELTA, ETC.                       !!  
  !! ----------------------------------------------------------------------- !!
  !! Test levels
  IF(delta1>=zero .AND. delta1<=one) THEN
    dAS = delta1
    dDS = delta1
    auto_delta = .FALSE.
    user_delta = .TRUE.
  ENDIF
  IF(level1>zero .AND. level1<=one) level_S1 = level1
  IF(level2>zero .AND. level2<=one) level_S2 = level2
  
  IF(delta2>=zero .AND. delta2<=one) dDS = delta2

  !! Check for non-sense delta_step, delta_start, delta_end input
  IF(delta_step == zero) &
    CALL PrntE("Value given by --delta-step must be non-zero!", Q=.TRUE.)
  IF(delta_start < delta_end .AND. delta_step <= zero) &
    CALL PrntE("With --delta-start smaller than --delta-end the value"//&
                    " given by --delta-step must be positive!", Q=.TRUE.)
  IF(delta_start > delta_end .AND. delta_step >= zero) &
    CALL PrntE("With --delta-start larger than --delta-end the value"//&
                    " given by --delta-step must be negative!", Q=.TRUE.)

  !! If not all non-zero, then none non-zero
  IF(delta_step == NAv2) THEN
    delta_start = NAv2
    delta_end = NAv2
  ENDIF
  
  !! Set automatic value of reruns if required and not set by user
  IF(ALL((/delta_step,delta_start,delta_end/)/=NAv2) .AND. max_reruns==0) &
    max_reruns = MAX_REPEATS
    
  !! --delta-min preceeds (or can act as) --delta
  IF(delta_start /= NAv2) THEN
    dAS = delta_start
    dDS = delta_start
  ENDIF
  IF(delta_start==NAv2 .AND. dAS>=zero) delta_start = dAS
  
  !! Adjuste OR step if OR max and min are the same
  IF(ABS(OR_max - OR_min) < rerun_eps) OR_step = zero

  !! --OR-min preceeds (or can act as) --OR
  IF(OR_min >= zero) OR = OR_min

  !! Check for no reruns needed at all
  IF(OR_step==zero .AND. delta_step==NAv2 .AND. nrun>1 .AND. &
    max_reruns == MAX_REPEATS) &
    CALL Prnt("Stopping because OR and delta steps are 0 and the number of"//&
              " repeats is automatically determined.", Q=.TRUE., skip1=1, &
              skip2=1, premature=.FALSE., dellog=.TRUE.)

  !! If multiple runs and delta_step or OR_step given, modify delta and/or OR
  IF(nrun>1 .AND. (delta_step/=NAv2 .OR. (simulate_data .AND. OR_step/=zero))) THEN
  
    nrun_delta = nrun_delta + 1
    dAS = dAS + nrun_delta * delta_step
    dDS = dDS + nrun_delta * delta_step
    
    IF(delta_step>zero) THEN
      cond1 = dAS <= MIN(delta_end + rerun_eps, one) 
      cond2 = dAS > MIN(delta_end + rerun_eps, one) .AND. &
              dAS - delta_step < delta_end - rerun_eps
    ELSE
      cond1 = dAS >= MIN(delta_end - rerun_eps, one) 
      cond2 = dAS < MAX(delta_end - rerun_eps, -one) .AND. &
              dAS - delta_step > delta_end + rerun_eps
    ENDIF
    
    IF(delta_step/=NAv2 .AND. (cond1 .OR. cond2)) THEN
      
      reuse_data = .TRUE.
      
      IF(cond2) THEN
        dAS = delta_end
        dDS = delta_end
      ENDIF
    
    ELSE
    
      nrun_delta = 0
      dAS = delta_start
      dDS = delta_start
      nrun_OR = nrun_OR + 1
      reuse_data = .FALSE.
    
    ENDIF
    
    OR = OR + nrun_OR * OR_step

    !! Check for stopping conditions
    text = ""
    IF(.NOT.simulate_data .AND. delta_step/=NAv2 .AND. nrun_delta==0) THEN
      text = "Stopping because all delta values repeated and automatic"//&
             " number of repeats is active." 
    ELSEIF(OR_step==zero .AND. nrun_OR>0 .AND. max_reruns == MAX_REPEATS) THEN
      text = "Stopping because OR step is 0 and automatic number of"//&
             " repeats is active."
    ELSEIF(OR_max>zero .AND. OR>OR_max+rerun_eps) THEN
      text = "User given maximum OR exceeded."             
    ELSEIF(OR_min>zero .AND. OR<OR_min-rerun_eps) THEN
      text = "User given minimum OR exceeded."
    ENDIF
    IF(LEN_TRIM(text)>0) &
      CALL Prnt(text, Q=.TRUE., skip1=1, skip2=1, premature=.FALSE., &
                     dellog=.TRUE.)
  ENDIF

  !! Check if user wants delta to range over multiple values
  IF(delta_step /= NAv2) THEN
    auto_delta = .FALSE.
    user_delta = .TRUE.
    fix_subsamp = .FALSE.
  ENDIF 
  
  !! If pretest level is NA, then try to get automatically
  auto_level_use = level_S1 == NAneg
  
  !! If all pretests are to be rejected, set level_S1 to 1
  IF(reject_all_pretests) THEN
    level_S1 = one
    auto_level_use = .FALSE.
  ENDIF
  
  !! Remember nominal levels
  level_S1_nom  = level_S1
  level_S2_nom = level_S2

  !! Select pretest
  !! Controls only with splitting and chisquare test
  IF(ANY(WT11==(/0,T1co,T1ca,T1cc,T1po,T1sc/))) WT1 = WT11

  !! Disable the pretests
  IF(WT1==0) doS1 = .FALSE.

  !! Select main test
  IF(ANY(WT21==(/0,T2sc,T2cc/))) WT2 = WT21

  !! Disable main tests
  doS2 = WT2 > 0

  IF(WT2==T2sc .AND. WT1==T1sc) THEN
    doAS4 = .FALSE.
    doAS1 = .FALSE.
  ENDIF
  
  !! Change parameters if dAS and/or dDS too low or high
  IF(.NOT.auto_delta .AND. dAS==zero .AND. dDS==zero) doS1 = .FALSE.
    
  IF(dAS > delta_bound) THEN
    IF(cntrgrp==1)  cntrgrp = 2
    IF(WT1 == T1cc) cntrgrp = 0
  ENDIF
  
  !! Disable reporting of adjusted and disjoint pretests if not sensible
  IF(WT1 == T1po) THEN
    doAS4 = .FALSE.
    doAS1 = .FALSE.
    doDS4 = .FALSE.
    doDS1 = .FALSE.
    doPO4 = .TRUE.
    doPO1 = .TRUE.
    doPT4 = .FALSE.
    doPT1 = .FALSE.
    doCS = .TRUE.
    cntrgrp = 0
  ENDIF
  
  IF(dDS > delta_bound) THEN
    doDS4 = .FALSE.
    doDS1 = .FALSE.
  ENDIF

  !! Assign default S1 degrees of freedom
  IF(WT1 == T1co) T1_df = T1co_df
  IF(WT1 == T1ca) T1_df = T1ca_df
  IF(WT1 == T1cc) T1_df = T1cc_df
  IF(WT1 == T1cc) T1_df = T1po_df
  
  !! Change the pretest degrees of freedom to user input
  IF(T1_df1>0 .AND. T1_df1==INT(T1_df1)) T1_df=INT(T1_df1)
  
  !! If random pretest subsample, no pss file saving
  IF(.NOT.fix_subsamp) THEN
    do_out_pss = .FALSE.
    IF(pss_file_subsample) &
      CALL PrntE("Pretest sample cannot be random when S1 selection status"//&
                 " (PSS) file specified.", Q=.TRUE., premature=.FALSE.)
  ENDIF
  
  !! pss is inputted, adjust some variables
  IF(pss_nfiles>0) THEN
    auto_delta = .FALSE.
    user_delta = .FALSE.
  ENDIF

  !! If none of the tests is wanted, don't do any testing
  IF(.NOT.ANY((/doAS4, doAS1, doDS4, doDS1, doPO4, doPO1, doCS/))) &
    no_testing = .TRUE.
   
  !! NO TESTS AT ALL WILL BE PERFORMED !!  
  IF(no_testing) THEN
    add_settings_filename = .FALSE.
    do_out_pss = .FALSE.
    ana_model = 1 
  ENDIF

  !! NO PRETEST WILL BE PERFORMED !!  
  IF(.NOT.doS1) THEN
    WT1  = 0
    dAS = zero
    dDS = zero
    auto_delta = .FALSE.
    auto_level_use = .FALSE.
    level_S1 = one
    level_S1_nom = one
    doAS4 = .FALSE.
    doAS1 = .FALSE.
    doDS4 = .FALSE.
    doDS1 = .FALSE.
    doPT4 = .FALSE.
    doPT1 = .FALSE.
    doPO4 = .FALSE.
    doPO1 = .FALSE.
    do_out_pss = .FALSE.
  ENDIF

  !! NO POST-TEST WILL BE PERFORMED !!  
  IF(.NOT.doS2) THEN
    WT2  = 0
    !dAS = one
    !dDS = one
    !auto_delta = .FALSE.
    cntrgrp = 0
    doCS = .FALSE.
  ENDIF
  
  !! Disable pretests if no two-stage method selected
  IF(.NOT.ANY((/doAS4,doAS1,doDS4,doDS1,doPO4,doPO1/))) doS1 = .FALSE.
  
  !! Make sure the model selection fits
  IF(ana_model<1 .OR. ana_model > SIZE(model_names)) &
    ana_model = def_ana_model
  IF(sim_model<1 .OR. sim_model > SIZE(model_names)) &
    sim_model = def_sim_model
    
  IF(ana_model_S1<1) ana_model_S1 = ana_model
    
  !! Make sure missing sim_model does not cause error when not relevant 
  IF(.NOT.simulate_data) sim_model = 1
    
  !! If user made no model (analysis and/or simulation) selection and default 
  !! model numbers are negative (set in module PARAMETERS), then stop
  IF(ana_model<1 .OR. sim_model<1) THEN
    IF(ana_model<1) THEN
      CALL PrntE("Analysis model not specified!", Q=.FALSE.)
      CALL Prnt("Use --model or --amodel to set it.")
    ENDIF
    IF(sim_model<1) THEN
      CALL PrntE("Simulation model not specified!")
      CALL Prnt("Use --model or --smodel to set it.")
    ENDIF
    CALL Prnt0("", Q=.TRUE., premature=.FALSE.)
  ENDIF

  !! ----------------------------------------------------------------------- !!
  !!                         AUTO LEVEL DETERMINATION                        !!  
  !! ----------------------------------------------------------------------- !!
  !! Make sure that minor allele frequencies get computed for auto level
  IF(.NOT.no_testing) THEN
    IF(out_maf) compute_maf = .TRUE.
    IF(auto_level_report .AND. auto_level_single) compute_maf = .TRUE.
    IF(auto_level_report .AND. .NOT.auto_level_maf) compute_maf = .TRUE.
  ENDIF 
  IF(min_MAF>zero) compute_maf = .TRUE. 

  !! Make sure the odds ratios for automatic determination of OR are within 
  !! the allowed range (must be positive!)
  auto_level_OR1 = MIN(MAX(auto_level_OR1, min_OR), max_OR)
  auto_level_OR2 = MIN(MAX(auto_level_OR2, min_OR), max_OR)
  auto_level_OR = MIN(MAX(auto_level_OR, min_OR), max_OR)
  
  !! Check for non-sense values in auto_level_min and auto_level_max
  IF(auto_level_max <= auto_level_min) &
    CALL PrntE("The value given by --autolevelmin must be smaller than"//&
                    " the value given by --autolevelmax!", Q=.TRUE., &
                    premature=.FALSE.)
                    
  !! Don't check for existence of output files if data will be simulated
  IF(simulate_data) check_fexist = .FALSE. 

  !! ======================================================================= !!
  !!                            RESULTS REPORTING                            !!  
  !! ======================================================================= !!

  !! Check for minimalistic output and change certain settings
  IF(out_minimal) THEN
    IF(var_poststand) &
      CALL PrntE("Post-standardization of AS is not allowed when only"//&
                 " minimalistic output is requested.")
    out_debug_info = .FALSE.

    IF(doAS4) doPTDS4 = .FALSE.
    IF(doAS1) doPTDS1 = .FALSE.
    IF(doDS4) doDS1 = .FALSE.
    
  ENDIF
  
  !! First, nulify the variables if report bounds will be set by user
  IF(olim1_auto) olim1 = level_S1
  IF(olim2_auto) olim2 = one
  IF(olim3_auto) olim3 = level_S2
  
  IF(MAX(olim1n, olim2n, olim3n)>zero .OR. report_all) THEN
    olim1 = two
    olim2 = two
    olim3 = two
  ENDIF
  
  !! Then, set the value given by user (limit by 2 since NAp must be bigger)
  IF(.NOT.report_all) THEN
    IF(olim1n > zero) olim1 = MIN(olim1n, two)
    IF(olim2n > zero) olim2 = MIN(olim2n, two)
    IF(olim3n > zero) olim3 = MIN(olim3n, two)
  ENDIF
  
  !! If all bounds are larger than 1, everything will be reported
  report_all = MIN(olim1, olim2, olim3) > one
    
  IF(COUNT((/olim1, olim2, olim3/) <= one)>1) &
    CALL PrntE("Conflicting reporting bounds.")

  IF(ALL(out_sort_sel/=sort_selections) .AND. out_sort_sel/="") THEN
    out_sort_sel = def_sort_sel
    CALL PrntE("Unknown output file sort selection. Using default value.")
  ENDIF
  
  !! Decide about possible sort column
  IF(out_sort_sel=="A" .AND. .NOT.(doAS4 .OR. doAS1)) THEN
    IF(doDS4 .OR. doDS1) THEN
      out_sort_sel = "D"
    ELSEIF(doCS) THEN
      out_sort_sel = "C"
    ELSEIF(doS1) THEN
      out_sort_sel = "P"
    ELSEIF(doPT4) THEN
      out_sort_sel = "PA"
    ELSEIF(doPT1) THEN
      out_sort_sel = "PB"
    ELSEIF(doPO4) THEN
      out_sort_sel = "PC"
    ELSEIF(doPO1) THEN
      out_sort_sel = "PD"
    ELSE
      out_sort_sel = "0"
    ENDIF
  ENDIF
  IF(out_sort_sel=="D" .AND. .NOT.(doDS4 .OR. doDS1)) THEN
    IF(doAS4 .OR. doAS1) THEN
      out_sort_sel = "A"
    ELSEIF(doCS) THEN
      out_sort_sel = "C"
    ELSEIF(doS1) THEN
      out_sort_sel = "P"
    ELSEIF(doPT4) THEN
      out_sort_sel = "PA"
    ELSEIF(doPT1) THEN
      out_sort_sel = "PB"
    ELSEIF(doPO4) THEN
      out_sort_sel = "PC"
    ELSEIF(doPO1) THEN
      out_sort_sel = "PD"
    ELSE
      out_sort_sel = "0"
    ENDIF
  ENDIF
  IF(out_sort_sel=="C" .AND. .NOT.doCS) THEN
    IF(doAS4 .OR. doAS1) THEN
      out_sort_sel = "A"
    ELSEIF(doDS4 .OR. doDS1) THEN
      out_sort_sel = "D"
    ELSEIF(doS1) THEN
      out_sort_sel = "P"
    ELSEIF(doPT4) THEN
      out_sort_sel = "PA"
    ELSEIF(doPT1) THEN
      out_sort_sel = "PB"
    ELSEIF(doPO4) THEN
      out_sort_sel = "PC"
    ELSEIF(doPO1) THEN
      out_sort_sel = "PD"
    ELSE
      out_sort_sel = "0"
    ENDIF
  ENDIF
  IF(out_sort_sel=="P" .AND. .NOT.doS1) THEN
    IF(doPO4) THEN
      out_sort_sel = "PB"
    ELSEIF(doAS4 .OR. doAS1) THEN
      out_sort_sel = "A"
    ELSEIF(doDS4 .OR. doDS1) THEN
      out_sort_sel = "D"
    ELSEIF(doCS) THEN
      out_sort_sel = "C"
    ELSEIF(doPT4) THEN
      out_sort_sel = "PA"
    ELSEIF(doPT1) THEN
      out_sort_sel = "PB"
    ELSEIF(doPO4) THEN
      out_sort_sel = "PC"
    ELSEIF(doPO1) THEN
      out_sort_sel = "PD"
    ELSE
      out_sort_sel = "0"
    ENDIF
  ENDIF
  
  !! ======================================================================= !!
  !!                         RESULT FILE ANALYSIS                            !!  
  !! ======================================================================= !!
  IF(run_analysis) THEN
    !! Read the names of output files to analyze and verify existence
    CALL GetFilenames(results_nfiles, ures, res_file, res_filelist, delim, &
                      comment)
    RETURN
  ELSE
    IF(ALLOCATED(res_file)) DEALLOCATE(res_file)
  ENDIF
  
  RETURN
     
END SUBROUTINE GetCmdLineArgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CheckMissingFiles()
  CHARACTER(mfl) :: fname, fkind
  LOGICAL        :: file_exists
  INTEGER        :: i

  !! Check for missing input files
  DO i=1,7
  
    IF(simulate_data .AND. i>=1 .AND. i<=4) CYCLE
    
    fname = ""
    IF(i==1 .AND. ped_nfiles>0) THEN
      fname = ped_file(1)
      fkind = "PED"
    ELSEIF(i==2 .AND. bed_nfiles>0) THEN
      fname = bed_file(1)
      fkind = "BED"
    ELSEIF(i==3 .AND. map_nfiles>0) THEN
      fname = map_file(1)
      fkind = "MAP"
    ELSEIF(i==4 .AND. fam_nfiles>0) THEN
      fname = fam_file
      fkind = "FAM"
    ELSEIF(i==5 .AND. sub_nfiles>0) THEN
      fname = sub_file(1)
      fkind = "SUBMAP"
    ELSEIF(i==6 .AND. pss_nfiles>0) THEN
      fname = pss_file
      fkind = "PSS"
    ELSEIF(i==7 .AND. sts_nfiles>0) THEN
      fname = sts_file
      fkind = "STS"
    ENDIF
    
    IF(LEN_TRIM(fname)==0) CYCLE 
    
    INQUIRE(FILE=fname, EXIST=file_exists)
    IF(.NOT.file_exists) &
      CALL PrntE("File ["//TRIM(fname)//"] ("//TRIM(fkind)//") does not exist.", Q=.TRUE.)
  
  ENDDO

END SUBROUTINE CheckMissingFiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetFilenames()
  IMPLICIT NONE
#ifdef __INTEL_COMPILER
  INTEGER, EXTERNAL :: GETPID
#endif
  CHARACTER(mfl) :: add, details, tlev1, tlev2, tdelta, tmaf
  !INTEGER        :: i

  details = ""

  !! Add process ID to filenames
  IF(append_pid) &
    details = "pid"//TRIM(i2cp(GETPID()))//TRIM(details) 
    
  !! Add time stamp to filenames
  IF(append_ts) &
    details = TRIM(timestamp)//TRIM(details)

  !! Add separator
  IF(append_pid .OR. append_ts) & 
    details = "_"//TRIM(details)
    
  !! Define testing settings that will be attached to filenames
  IF(.NOT.run_analysis .AND. add_settings_filename) THEN

    IF(doS1) &
      details = "_cg="//TRIM(i2cp(cntrgrp))//TRIM(details)

    !! Add more detailed info into file names
    IF(add_settings_filename_long) THEN

      IF(var_use_beta1) THEN
        details = "_varb1=y"//TRIM(details)
      ELSE
        details = "_varb1=n"//TRIM(details)
      ENDIF

      !! Level info
      tlev1 = r2c(level_S1)
      tlev2 = r2c(level_S2)
      IF(auto_level_use) tlev1 = "auto"
      details = "_a1="//TRIM(tlev1)//"_a2="//TRIM(tlev2)//TRIM(details)
    
      !! Test selection info
      IF(WT2>0)  details = "_t2="//i2cp(WT2)//TRIM(details)
      IF(WT2==0) details = "_t2=n"//TRIM(details)
      IF(WT1>0)  details = "_t1="//i2cp(WT1)//TRIM(details)
      IF(WT1==0) details = "_t1=n"//TRIM(details)

    ENDIF

    !! Add extra info if data is simulated
    IF(simulate_data) THEN

      IF(var_use_beta1 .AND. .NOT.add_settings_filename_long) &
        details = "_varb1=y"//TRIM(details)

      IF(simulate_HWE) THEN
        details = "_hwe=y"//TRIM(details)
      ELSE
        details = "_hwe=n"//TRIM(details)
      ENDIF
      
      IF(allelefreqA > zero) THEN
        tmaf = "_f1="//TRIM(r2c(round(allelefreqA,5)))
      ELSE
        tmaf = "_f1=na"
      ENDIF
      
      IF(allelefreqB > zero) THEN
        tmaf = TRIM(tmaf)//"_f2="//TRIM(r2c(round(allelefreqB,5)))
      ELSE
        tmaf = TRIM(tmaf)//"_f2=na"
      ENDIF
      
      IF(nneutralloci>0) THEN
        IF(allelefreqC > zero) THEN
          tmaf = TRIM(tmaf)//"_f3="//TRIM(r2c(round(allelefreqC,5)))
        ELSE
          tmaf = TRIM(tmaf)//"_f3=na"
        ENDIF
      ENDIF
      
      IF(fixed_MAF) THEN
        tmaf = TRIM(tmaf)//"_fix"
      ELSE
        tmaf = TRIM(tmaf)//"_lbnd"
      ENDIF
      
      details = TRIM(tmaf)//TRIM(details)
      IF(LD /= zero) details = "_LD="//TRIM(r2c(LD))//TRIM(details)

      details = "_ca="//TRIM(i2cp(ped_nca))//TRIM(details)
      details = "_co="//TRIM(i2cp(ped_nco))//TRIM(details)

      details = "_OR="//TRIM(r2c(OR, mindigs=1))//&
                "_OR1="//TRIM(r2c(OR1, mindigs=1))//&
                "_OR2="//TRIM(r2c(OR2, mindigs=1))//TRIM(details)

    ENDIF
    
    !! Pretest sample size ratio info
    IF(doS1) THEN
    
      IF(pss_nfiles>0 .OR. auto_delta) THEN 
        IF(pss_nfiles>0) THEN
          tdelta = "=pssfile"
        ELSE IF(auto_delta) THEN
          IF(simple_delta) THEN
            tdelta = "=auto-s"
          ELSE
            tdelta = "=auto-o"
          ENDIF
        ENDIF
      ELSE
        IF(dAS==dDS) THEN
          tdelta = "="//TRIM(r2c(dAS))
        ELSE
          tdelta = "AS="//TRIM(r2c(dAS))//"_dDS="//TRIM(r2c(dDS))
        ENDIF
      ENDIF
      details = "_d"//TRIM(tdelta)//TRIM(details) 
      
    ENDIF
    
    !! Analysis model info
    details = "_am="//TRIM(GetModelName(ana_model))//TRIM(details)
    !! Analysis model info
    details = "_sm="//TRIM(GetModelName(sim_model))//TRIM(details)

  ENDIF
  
  !! Remove extension from input files (if any present)
  IF(input_format>0) THEN
    IF(input_format==1) THEN
      CALL RemoveSubstrEnd(bed_file, bed_ext)
      CALL RemoveSubstrEnd(map_file, bim_ext)
      bed_file = TrimStr(bed_file)//bed_ext
      map_file = TrimStr(map_file)//bim_ext
    ELSEIF(input_format==2) THEN 
      CALL RemoveSubstrEnd(ped_file, ped_ext)
      CALL RemoveSubstrEnd(map_file, map_ext)
      ped_file = TrimStr(ped_file)//ped_ext
      map_file = TrimStr(map_file)//map_ext
    ENDIF
    CALL RemoveSubstrEnd(fam_file, fam_ext)
    fam_file = TrimStr(fam_file)//fam_ext
  ENDIF
  
  !! Remove extension from input pss file (if any present)
  IF(pss_nfiles>0) THEN
    CALL RemoveSubstrEnd(pss_file, pss_ext)
    pss_file = TrimStr(pss_file)//pss_ext
  ENDIF
  
  !! Remove space from arrays that could have spaces not removed by TrimStr 
  CALL DeleteAllSubstr(ped_file, space)
  CALL DeleteAllSubstr(bed_file, space)
  CALL DeleteAllSubstr(map_file, space)
  CALL DeleteAllSubstr(pss_file, space)
  CALL DeleteAllSubstr(sub_file, space)

  !! Set/Get LOG FILE name
  IF(.NOT.user_set_out) file_base = TRIM(def_file_base)//details

  !! Append extension to output file names
  log_file = TRIM(file_base)//log_ext
  out_file = TRIM(file_base)//out_ext
  tmp_file = TRIM(file_base)//tmp_ext
  sort_file = TRIM(file_base)//"_sorted"//out_ext
  out_pss_file = TRIM(file_base)//pss_ext 
  out_maf_file = TRIM(file_base)//maf_ext
  out_prev_file = TRIM(file_base)//prev_ext
  
  !! Remove any possible spaces
  CALL ReplaceAllSubstr(log_file, space, "_")
  CALL ReplaceAllSubstr(out_file, space, "_")
  CALL ReplaceAllSubstr(tmp_file, space, "_")
  CALL ReplaceAllSubstr(out_pss_file, space, "_")
  CALL ReplaceAllSubstr(out_maf_file, space, "_")
  CALL ReplaceAllSubstr(out_prev_file, space, "_")

  !! Set output data file names
  IF(save_input_data) THEN
    
    !! If not set by user, define filename to save data into
    IF(.NOT.user_set_out .AND. .NOT.user_set_save) THEN     
      
      IF(simulate_data) THEN
        save_file = "Epi_simdata_sm="//&
                    TRIM(GetModelName(sim_model))//&
                    "_nLE="//i2cp(nneutralloci)//&
                    "_nLD="//i2cp(nLDpairs)//&
                    "_n="//i2cp(ped_ss)//&
                    "_OR="//TRIM(r2c(OR, mindigs=1))// &
                    "_OR1="//TRIM(r2c(OR1, mindigs=1))//&
                    "_OR2="//TRIM(r2c(OR2, mindigs=1))//&
                    "_LD="//TRIM(r2c(LD))//&
                    "_Prev="//TRIM(r2c(round(prevalence,2)))
      ELSE
        save_file = TRIM(program_name)//"_outdata"
      ENDIF
    
    ENDIF

    !! Add simulation details
    IF(.FALSE.) THEN
      IF(simulate_data) THEN
  
        !! HWE used in simulation or not
        IF(simulate_HWE) THEN
          details = TRIM(details)//"_hwe=y"
        ELSE
          details = TRIM(details)//"_hwe=n"
        ENDIF
  
        !! Add allele info
        IF(allelefreqA > zero) THEN
          details = TRIM(details)//"_f1="//TRIM(r2c(round(allelefreqA,2)))
        ELSE
          details = TRIM(details)//"_f1=auto"
        ENDIF
        IF(allelefreqB > zero) THEN
          details = TRIM(details)//"_f2="//TRIM(r2c(round(allelefreqB,2)))
        ELSE
          details = TRIM(details)//"_f2=auto"
        ENDIF
        IF(nneutralloci>0) THEN
          IF(allelefreqC > zero) THEN
            details = TRIM(details)//"_f3="//TRIM(r2c(round(allelefreqC,2)))
          ELSE
            details = TRIM(details)//"_f3=auto"
          ENDIF
        ENDIF
        IF(fixed_maf) THEN
          details = TRIM(details)//"_fix"
        ELSE
          details = TRIM(details)//"_lbnd"
        ENDIF
        
      ENDIF
                  
      !!! Add time stamp to filenames
      !IF(append_ts) &
      !  save_file = TRIM(save_file)//"_"//TRIM(timestamp)
      !!! Add process ID to filenames
      !IF(append_pid) &
      !  save_file = TRIM(save_file)//"_pid="//TRIM(i2cp(GETPID()))
        
    ENDIF
    
    !! Append the details
    save_file = TRIM(save_file)//TRIM(details)
    
    !! Remove all possible spaces
    CALL ReplaceAllSubstr(save_file, space, "_")
    
    !! Define names
    IF(save_binary) THEN
      save_ped_file = TRIM(save_file)//bed_ext
      save_map_file = TRIM(save_file)//bim_ext 
      save_fam_file = TRIM(save_file)//fam_ext
    ELSE 
      save_ped_file = TRIM(save_file)//ped_ext
      save_map_file = TRIM(save_file)//map_ext
    ENDIF

  ENDIF
  
  !! Modify the output and log file names if repeated run
  IF(ped_nrepeat > 0 .AND. .NOT.append_ts) THEN
    add = "_"//i2cp(ped_nrepeat)
    CALL ChangeFilename(log_file, log_ext, add)
    CALL ChangeFilename(out_file, out_ext, add)
    CALL ChangeFilename(out_pss_file, pss_ext, add)
    CALL ChangeFilename(out_maf_file, maf_ext, add)
    CALL ChangeFilename(out_prev_file, prev_ext, add)
    IF(save_binary) THEN
      CALL ChangeFilename(save_ped_file, bed_ext, add)
      CALL ChangeFilename(save_map_file, bim_ext, add)
      CALL ChangeFilename(save_fam_file, fam_ext, add)
    ELSE
      CALL ChangeFilename(save_ped_file, ped_ext, add)
      CALL ChangeFilename(save_map_file, map_ext, add)
    ENDIF
  ENDIF

  RETURN
     
END SUBROUTINE SetFilenames

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrintHelpScreen()
  IMPLICIT NONE
  INTEGER, PARAMETER           :: max_nflag = 200, nlead = 1
  INTEGER, PARAMETER           :: lenempty = 1
  CHARACTER(mstl)              :: switch(max_nflag)
  CHARACTER(mltl)              :: label(max_nflag), text
  CHARACTER(mstl)              :: alias(max_nflag), emptysp
  CHARACTER(mstl), ALLOCATABLE :: words(:)
  INTEGER                      :: i, j, n, swlen
  
  emptysp = ""
  alias = ""
  n=0
  
  !! Read in the help information
  n=n+1
  switch(n) = "--help or --h"
  label(n) = "Prints this help screen."
  n=n+1
  switch(n) = "--cmdfile"
  label(n) = "Follow by file name of a file that contains setting switchers." //&
             " Each line of this file is processed separately as a single" //&
             " call of the program."
  n=n+1
  switch(n) = "--nrepeat"
  label(n) = "Follow by how many times the analysis should be rerun (mostly usefull in simulation)."
  n=n+1
  switch(n) = "--silent"
  label(n) = "Supresses any output to the console."
  n=n+1
  switch(n) = "--no-termination"
  label(n) = "If present the program does not exit after performing "//&
    "analysis and asks user for additional command line input."
  n=n+1
  switch(n) = "--terminate"
  label(n) = "Negates the --no-termination switch, i.e. forces "//&
    "termination after analysis."
  n=n+1
  switch(n) = "--overwrite"
  label(n) = "Disables file existence check and overwrites existing"//&
    " output and log files."
  n=n+1
  switch(n) = "--nocountdown"
  label(n) = "Supresses testing countdown from being printed."
  n=n+1
  switch(n) = "--nooutput"
  label(n) = "No output file with test results will be produced."
  n=n+1
  switch(n) = "--analyze"
  label(n) = "Runs the analysis of output files of previous runs."
  n=n+1
  switch(n) = "--simulateinput"
  label(n) = "If present the input data will be simulated. Requires"//&
    " --samplesize with positive value and makes --ped and"//&
    " --map switchers redundant. User can specify the"//&
    " name of the log file by using --log to avoid generic"//&
    " name of the log file."
  n=n+1
  switch(n) = "--saveinput"
  label(n) = "If present the simulated input data will be save into"//&
    " a file."
  n=n+1
  switch(n) = "--reusedata"
  label(n) = "If present the input data from previous run will be used,"//&
    " provided the previous run included --no-termination or the runs are"//&
    " specified through --cmdfile. Applies both to simulated data and data"//&
    " read from a file."
  n=n+1
  switch(n) = "--analyzeasonefile"
  label(n) = "Analysis mode only. If present all files on input of the"//&
    " analysis will be treated as one big file and analyzed together."//&
    " (Use only with --analyze)"
  n=n+1
  switch(n) = "--outnoheader"
  label(n) = "Output file: No header will be printed."
  n=n+1
  switch(n) = "--outheadercommented"
  label(n) = "Output file: Header will be printed but commented."
  n=n+1
  switch(n) = "--outheadernamed"
  label(n) = "Output file: Header will have a name for each column."
  n=n+1
  switch(n) = "--firstsignificant"
  label(n) = "If present the first pair in a sample is considered"//&
    " truly interacting and reports of successfull"//&
    " detection are printed."
  n=n+1
  switch(n) = "--testsamechr"
  label(n) = "Presence indicates that pairs of loci with both loci"//&
    " on the same chromosome will also be tested."
  n=n+1
  switch(n) = "--skipsamechr"
  label(n) = "Presence indicates that pairs of loci with both loci"//&
    " on the same chromosome will not be tested. Default behaviour."
  n=n+1
  switch(n) = "--nolog"
  label(n) = "Supresses printing of the log file."
  n=n+1
  switch(n) = "--reportall"
  label(n) = "All test results including erroneous are printed in the"//&
    " output file regardless of the values specified by --report1,"//&
    " --report2, --report3."
  n=n+1
  switch(n) = "--reportallresults"
  label(n) = "All non-erroneous test results are printed in the"//&
    " output file regardless of the values specified by --report1,"//&
    " --report2, --report3."
  n=n+1
  switch(n) = "--reporterrors"
  label(n) = "Defaultly, erroneous tests are not printed in the output"//&
    " file. Adding this switch will print them (unless --reportall comes later)."
  n=n+1
  switch(n) = "--hideerrors"
  label(n) = "This switch causes erroneous tests are not printed in"//&
    " the output file. Default behavior."
  n=n+1
  switch(n) = "--reportscore"
  label(n) = "If present the program also reports classical score"//&
    " test results for the phase 2 tests that are performed."
  n=n+1
  switch(n) = "--noscore"
  label(n) = "If present the program will not report the results of"//&
    " classical score tests."
  n=n+1
  switch(n) = "--allscore"
  label(n) = "If present the program computes classical score test"//&
    " for all possible pairs, not only pairs that are rejected in the pretest."
  n=n+1
  switch(n) = "--nocentering"
  label(n) = "If present the pretest vector will NOT be centered"//&
    " for regression during the main test. Default behavior is to center it."
  n=n+1
  switch(n) = "--centerbyself"
  label(n) = "If present the pretest vector will be centered by the"//&
    " individuals from the same group (cases or controls) for the purposes"//&
    " of regression during the main test. Otherwise, the other group is used"//&
    " to center the pretest vector. Default behavior."
  n=n+1
  switch(n) = "--centerbyother"
  label(n) = "If present the pretest vector will be centered by the"//&
    " individuals from the other group (cases or controls) for the purposes"//&
    " of regression during the main test. Otherwise, the same group is used"//&
    " to center the pretest vector. Applicable only for pretest 3 and 4."
  n=n+1
  switch(n) = "--nonparametric"
  label(n) = "If present the probabilities of being case or control"//&
    " (conditioned on a genotype) will be computed"//&
    " nonparametrically."
  n=n+1

  !! Input and output files
  switch(n) = "--uniquefilenames"
  label(n) = "If present the output and log files will contain input"//&
    " parameter information in them. (Default behavior)."
  n=n+1
  switch(n) = "--comment"
  label(n) = "Allows to specify the comment character (length 1)."//&
    " Default value "//TRIM(comment)//"."
  n=n+1
  !switch(n) = "--NAoutput"
  !label(n) = "Allows to specify the NA character for output files, not"//&
  !  " the input files (length 1). Default value "//TRIM(NA)//"."
  !n=n+1
  switch(n) = "--nthreads"
#ifdef _OPENMP
  label(n) = "Allows to specify the number of parallel threads to run (integer)."
#else
  label(n) = "Not available in this version. The source code must be compiled"//&
    " with an OpenMP-compliant compiler such as Gfortran (flag --openmp)."
#endif      
  n=n+1
  switch(n) = "--ped"
  label(n) = "Specifies the name of an input ped file (relative or absolute path)."
  n=n+1
  switch(n) = "--pedlist"
  label(n) = "Specifies the name of a file with a list of ped files The path"//&
    " to this file will be added to filenames in it (one file per row)."
  n=n+1
  switch(n) = "--bed"
  label(n) = "Specifies the name of an input binary ped file (relative or absolute path)."
  n=n+1
  switch(n) = "--map"
  label(n) = "Specifies the name of an input mapping file (relative or absolute path)."
  n=n+1
  switch(n) = "--mapncols"
  label(n) = "Specifies the number of columns in the given map files."
  n=n+1
  switch(n) = "--maplist"
  label(n) = "Specifies the name of a file with a list of mapping files. The path"//&
    " to this file will be added to filenames in it (one file per row)."
  n=n+1
  switch(n) = "--filter, --submap"
  label(n) = "Specifies the name of a file with a list of loci (subset of specified mapping"//&
    " file) for which to perform the tests (relative or absolute path). If missing"//&
    " considered equal to the specified mapping file."
  n=n+1
  switch(n) = "--submaplist"
  label(n) = "Specifies the name of a file with a list of mapping subset files. The path"//&
    " to this file will be added to filenames in it (one file per row)."
  n=n+1
  switch(n) = "--submapncols"
  label(n) = "Specifies the number of columns in the submap file(s)."
  n=n+1
  switch(n) = "--submaprscol"
  label(n) = "Specifies the position of a RS column within the"//&
    " submap file(s)."
  n=n+1
  switch(n) = "--fam"
  label(n) = "Specifies the name of an input pedigree (family) info and"//&
             " status file (relative or absolute path)."
  n=n+1
  switch(n) = "--famncols"
  label(n) = "Specifies the number of columns in a pedigree (fam) file(s) or"//&
    " number of non-genetic columns in a ped file"//&
    " (positive integer). In a standard linkage format there"//&
    " are 6 non-genetic columns (pedigree ID, individual ID,"//&
    " father's ID, mother's ID, sex, status). Default value "//&
    i2cp(fam_ncols)//"."
  n=n+1
  switch(n) = "--famsexcol"
  label(n) = "Specifies the position of a sex status column within"//&
    " the fam file(s). Put 0 if sex column not present."
  n=n+1
  switch(n) = "--famstatuscol"
  label(n) = "Specifies the position of a status column within"//&
    " the fam file(s)."
  n=n+1
  switch(n) = "--pss"
  label(n) = "Specifies the name of pretest selection status file (relative or absolute path)."
  n=n+1
  switch(n) = "--output"
  label(n) = "Specifies the name of output file to print the test results into."
  n=n+1
  switch(n) = "--savefile"
  label(n) = "Specifies the root name of file to which the input data"//&
    " will be save saved if --saveinput also present. Plain text ped"//&
    " file with extension .ped (or .bed if --savebinary present),"//&
    " mapping file with extension .map (or .bim) and fam file with"//&
    " extension .fam (generated only if --savebinary present)."
  n=n+1
  switch(n) = "--log"
  label(n) = "Specifies the name of log file (relative or absolute path)."
  n=n+1
  switch(n) = "--resultfile"
  label(n) = "Specifies the name of a file with results of previous"//&
    " program runs to be analyzed. (Use only with --analyze)"
  n=n+1
  switch(n) = "--resultfilelist"
  label(n) = "Specifies the name of a file with a list of result"//&
    " files of previous program runs to be analyzed."//&
    " (Use only with --analyze)"
  n=n+1
  switch(n) = "--resultfilemaxlines"
  label(n) = "Specifies the maximum number of lines read from input"//&
    " file(s) (0 means no limit). Default value: 0."//&
    " (Use only with --analyze)"
  n=n+1

  !! Input file parameters
  switch(n) = "--sepascii"
  label(n) = "If present the specified numeric values of separators will"//&
    " be converted to characters according to ASCII table. Default behavior."
  n=n+1
  switch(n) = "--sep"
  label(n) = "Specifies the separator or its ASCII code of the"//&
    " ped/map/submap/fam files (character or integer)."//&
    " If ASCII code given, then --sepascii needs to be also present."//&
    " Default value tabulator."
  n=n+1
  switch(n) = "--pedsep"
  label(n) = "Specifies the separator or its ASCII code of the ped"//&
    " files (character or integer)."//&
    " If ASCII code given, then --sepascii needs to be also present."//&
    " Default value tabulator."
  n=n+1
  switch(n) = "--mapsep"
  label(n) = "Specifies the separator or its ASCII code of the mapping"//&
    " files (character or integer)."//&
    " If ASCII code given, then --sepascii needs to be also present."//&
    " Default value tabulator."
  n=n+1
  switch(n) = "--submapsep"
  label(n) = "Specifies the separator or its ASCII code of the mapping"//&
    " subset files (character or integer)."//&
    " If ASCII code given, then --sepascii needs to be also present."//&
    " Default value tabulator."
  n=n+1
  switch(n) = "--famsep"
  label(n) = "Specifies the separator or its ASCII code of the fam"//&
    " files (character or integer)."//&
    " If ASCII code given, then --sepascii needs to be also present."//&
    " Default value tabulator."
  n=n+1
  switch(n) = "--statuscase"
  label(n) = "Specifies the character for case status in the input file(s)."//&
    " Default value "//def_CharCa//"."
  n=n+1
  switch(n) = "--statuscontrol"
  label(n) = "Specifies the character for control status in the input file(s)."//&
    " Default value "//def_CharCo//"."
  n=n+1
  switch(n) = "--missing-genotype"
  label(n) = "Specifies the character for NA genotype in the input file(s)."//&
    " Default value "//def_CharNAgen//"."
  n=n+1
  switch(n) = "--missing-status"
  label(n) = "Specifies the character for NA status in the input file(s)."//&
    " Default value "//def_CharNAsts//"."
  n=n+1
  switch(n) = "--missing-sex"
  label(n) = "Specifies the character for NA sex in the input file(s)."//&
    " Default value "//def_CharNAsex//"."
  n=n+1
  switch(n) = "--sexmale"
  label(n) = "Specifies the character for male sex in the input file(s)."//&
    " Default value "//def_CharMa//"."
  n=n+1
  switch(n) = "--sexfemale"
  label(n) = "Specifies the character for female sex in the input file(s)."//&
    " Default value "//def_CharFe//"."
  n=n+1
  switch(n) = "--delim"
  label(n) = "Specifies the file path delimiter of given paths (e.g. forwardslash or backslash)."
  n=n+1
  switch(n) = "--pedskip"
  label(n) = "Specifies the number of lines to be skipped in each ped file excluding"//&
    " commented lines (integer). Default value 0."
  n=n+1
  switch(n) = "--mapskip"
  label(n) = "Specifies the number of lines to be skipped in each mapping file excluding"//&
    " commented lines (integer). Default value 0."
  n=n+1
  switch(n) = "--submapskip"
  label(n) = "Specifies the number of lines to be skipped in each subset mapping file excluding"//&
    " commented lines (integer). Default value 0."
  n=n+1
  switch(n) = "--famskip"
  label(n) = "Specifies the number of lines to be skipped in each fam file excluding"//&
    " commented lines (integer). Default value 0."
  n=n+1
  switch(n) = "--pssskip"
  label(n) = "Specifies the number of lines to be skipped in the pretest "//&
    " selection status file excluding commented lines (integer). Default value 0."
  n=n+1
  switch(n) = "--resultfileskip"
  label(n) = "Specifies the number of lines to be skipped in each input file excluding"//&
    " commented lines (integer). Default value 0. (Use only with --analyze)"
  n=n+1
  switch(n) = "--pedmaxncols"
  label(n) = "Specifies the maximum number of columns in a ped file (integer)."//&
    " Default value "//i2cp(ped_maxncol)//"."
  n=n+1
  switch(n) = "--minsamplesize"
  label(n) = "Specifies the lower bound for the sample size that needs"//&
    " to be met in order for a test to be performed."//&
    " (non-negative integer)."
  n=n+1
  switch(n) = "--mincontrols"
  label(n) = "Specifies the lower bound for the size of sample of "//&
    " controls that needs to be met in order for a test to be performed."//&
    " (non-negative integer)."
  n=n+1
  switch(n) = "--mincases"
  label(n) = "Specifies the lower bound for the size of sample of "//&
    " cases that needs to be met in order for a test to be performed."//&
    " (non-negative integer)."
  n=n+1
  switch(n) = "--samplesize"
  label(n) = "Specifies the sample size if multiple samples are"//&
    " contained within ped file(s) or when input data"//&
    " is to be simulated specified by --simulateinput"//&
    " (positive integer)."
  n=n+1
  switch(n) = "--ncontrols"
  label(n) = "Specifies the number of controls in the sample."
  n=n+1
  switch(n) = "--ncases"
  label(n) = "Specifies the number of cases in the sample."
  n=n+1
  switch(n) = "--nsamples"
  label(n) = "Specifies the number of input data samples that should"//&
    " be simulated. Use only with --simulateinput."//&
    " (positive integer)."
  n=n+1

  !! Testing parameters
  switch(n) = "--disjoint"
  label(n) = "Determines whether for pretest and second test disjoint subsamples of"//&
    " the original sample will be used."//&
    " Default behavior is to use entire sample in both tests."
  n=n+1
  switch(n) = "--reportdisjoint"
  label(n) = "If present test statistic for the second test computed from"//&
    " the part of the sample not used in the first step will also be"//&
    " computed."
  n=n+1
  switch(n) = "--nopretest"
  label(n) = "If present no pretests will be performed."
  n=n+1
  switch(n) = "--dopretest"
  label(n) = "If present pretests will be performed (Default behavior)."
  n=n+1
  switch(n) = "--nomaintest"
  label(n) = "If present only pretests and no second step testing will be performed."
  n=n+1
  switch(n) = "--domaintest"
  label(n) = "If present second step testing will be performed (Default behavior)."
  n=n+1
  switch(n) = "--pretest"
  label(n) = "Selects the pretest. Possible values: 1,2,3,4."//&
    " Default value: "//i2cp(def_T1)//"."
  n=n+1
  switch(n) = "--test"
  label(n) = "Selects the main test performed to detect interactions. Possible values: 1,2."//&
    " Default value: "//i2cp(def_T2)//"."
  n=n+1
  switch(n) = "--test2"
  label(n) = "Selects the second main test to be performed to detect interactions."//&
    " Possible values: 1,2. Default value: 0 (no second main test)."
  n=n+1
  switch(n) = "--errorrate1"
  label(n) = "Selects the type of error to be controlled during pretest (FWER, FDR, pFDR)."//&
    " Default value FWER (currently not supported)."
  n=n+1
  switch(n) = "--errorrate2"
  label(n) = "Selects the type of error to be controlled during pretest (FWER, FDR, pFDR)."//&
    " Default value FWER (currently not supported)."
  n=n+1
  switch(n) = "--level1"
  label(n) = "Determines the level of significance of the pretest."//&
    " Possible values: real number between 0 and 1 or string auto."//&
    " Default value: "//TRIM(r2c(level_S1))//"."//&
    " If value 'auto' the pretest level (level1) is determined"//&
    " from the main test level (level2) and the total"//&
    " number of tests (N) as level1=SQRT(level2/N)."
  n=n+1
  switch(n) = "--level2"
  label(n) = "Determines the level of significance of the main test."//&
    " Possible values: real numbers between 0 and 1."//&
    " Default value: "//TRIM(r2c(level_S2))//"."
  n=n+1
  switch(n) = "--pretestdf"
  label(n) = "Sets the degrees of freedom for chisquare distribution"//&
    " quantiles in pretest. Default value "//i2cp(T1_df)//"."
  n=n+1

  switch(n) = "--correct1"
  label(n) = "Determines the correction factor for multiple testing for the adjusted"//&
    " or disjoint second step statistic. If integer than p-values multiplied"// &
    " by the given value. If real number between 0 and 1 than the correction"//&
    " factor is taken as the percentage of the total number of tests or the"//&
    " percentage of the value given by --correct2."
  n=n+1
  switch(n) = "--correct2"
  label(n) = "Determines the correction factor for multiple testing for the classical"//&
    " score second step statistic. If integer than p-values multiplied"// &
    " by the given value. If real number between 0 and 1 than the correction"//&
    " factor is taken as the percentage of the total number of tests."
  n=n+1
  switch(n) = "--report1"
  label(n) = "Determines the report bound for raw p-values of the pretest"//&
    " (real value between 0 and 1). If all results are to be printed select"// &
    " any value larger than 1. Default behavior: report all."
  n=n+1
  switch(n) = "--report2"
  label(n) = "Determines the report bound for raw p-values of the main test (real value between 0 and 1). "// &
    "If all results are to be printed select any value larger than 1. Superseeds --report1."//&
    " Default behavior: report all."
  n=n+1
  switch(n) = "--report3"
  label(n) = "Determines the report bound for corrected p-values of the main test."//&
    " Possible values: real number between 0 and 1."// &
    " If all results are to be printed select any value larger than 1."//&
    " Superseeds both --report1 and --report2. Default behavior: report all."
  n=n+1
  switch(n) = "--reportlogiterr"
  label(n) = "Determines whether errors encountered during logistic regression"//&
    " will be reported to the console. Default behavior: errors are"//&
    " not reported to console, only to the output and log files."
  n=n+1
  switch(n) = "--delta"
  label(n) = "Determines the portion of controls to be used during pretest"//&
    " by adjusted score and disjoint score. Possible values: real numbers"//&
    " between 0 and 1 or 'auto'."
  n=n+1
  switch(n) = "--deltaAS"
  label(n) = "Determines the portion of controls to be used during pretest"//&
    " by adjusted score. Possible values: real numbers between 0 and 1 or"//&
    " 'auto'. Default value "//TRIM(r2c(dAS))//"."
  n=n+1
  switch(n) = "--deltaDS"
  label(n) = "Determines the portion of controls to be used during pretest"//&
    " by disjoint score. Possible values: real numbers between 0 and 1 or"//&
    " 'auto'. Default value "//TRIM(r2c(dDS))//"."
  n=n+1
  switch(n) = "--analysismodel"
  label(n) = "Specifies the penetrance model used during analysis."
  alias(n) = "--amodel, --model"
  n=n+1
  switch(n) = "--minmaf"
  label(n) = "Specifies a lower bound for the minor allele frequency (MAF)."//&
    " If MAF for given pair of loci is lower than this bound the pair"//&
    " will not be tested. When used with --simulateinput it"//&
    " specifies the lower bound or the fixed MAF (when "//&
    " --fixedmaf present) Possible values: real numbers between 0 and 1."
  n=n+1
  switch(n) = "--mincellcount"
  label(n) = "Specifies a lower bound for the genotype frequency. If"//&
    " for a given pair of loci any cell in the genotype contingency"//&
    " table has frequency lower than this bound the pair will not be tested."//&
    " Possible values: non-negative integers. Default value "//&
    TRIM(r2c(min_cell_cnt))
  IF(min_cell_cnt==zero) label(n) = TRIM(label(n))//" (no exclusion)"
  label(n) = TRIM(label(n))//"."
  n=n+1
  switch(n) = "--minmargincount"
  label(n) = "Specifies a lower bound for the marginal genotype frequency. If"//&
    " for given pair of loci any marginal frequency in the genotype contingency"//&
    " table has frequency lower than this bound the pair will not be tested."//&
    " Possible values: non-negative integers. Default value "//&
    TRIM(r2c(min_marg_cnt))
  IF(min_marg_cnt==zero) label(n) = TRIM(label(n))//" (no exclusion)"
  label(n) = TRIM(label(n))//"."
  n=n+1
  switch(n) = "--cellcorrection"
  label(n) = "Sets the count by which the cells of the pretest contingency"//&
    " table will be increased if any of the cells is below --mincellcount."//&
    " (non-negative integer). Default value 0."
  n=n+1
  switch(n) = "--varbound"
  label(n) = "Sets the lower bound for a value of an estimate for the variance of adjusted"//&
    " score statistic (non-negative real). Default value 0."
  n=n+1
  switch(n) = "--sexchrcheck"
  label(n) = "If present males will be excluded from analysis when"//&
                   " 23rd chromosome analyzed. (Default behavior)"
  n=n+1
  switch(n) = "--nosexchrcheck"
  label(n) = "If present males will be included into analysis when"//&
                   " 23rd chromosome analyzed."
  n=n+1
  switch(n) = "--excludemales"
  label(n) = "If present all males will be excluded from analysis."
  n=n+1
  switch(n) = "--excludefemales"
  label(n) = "If present all females will be excluded from analysis."
  n=n+1
  switch(n) = "--assumemales"
  label(n) = "If present all individuals with missing sex status will be assumed to be males."
  n=n+1
  switch(n) = "--assumefemales"
  label(n) = "If present all individuals with missing sex status will be assumed to be females."
  n=n+1
  switch(n) = "--poststandardize"
  label(n) = "If present those pairs where the variance of adjusted score was estimated as"//&
    " too small (lower than the bound given by --varbound) is standardized by the"//&
    " average variance of all other pairs or by the variance bound if "//&
    " --standardbybound is present."
  n=n+1
  switch(n) = "--poststandbymean"
  label(n) = "If present and --poststandardize is also present the standardization will be"//&
    " done with the average variance of those pairs of loci with variance above"//&
    " the variance bound. Default behavior."
  n=n+1
  switch(n) = "--poststandbybound"
  label(n) = "If present and --poststandardize is also present the standardization will be"//&
    " done with the variance low bound and not with the average variance of other"//&
    " pairs of loci."
  n=n+1
  switch(n) = "--tempcycle"
  label(n) = "Specifies how often results will be written into temporary output file (integer)."//&
    " In some cases lowering this value may help resolve segmentation faults."//&
    " Default value "//i2cp(temp_cycle)//"."
  n=n+1
  switch(n) = "--maxfilesize"
  label(n) = "Specifies the maximum size of one (temporary) output file (integer)."//&
    " Default value "//i2cp(max_out_size)//"."
  n=n+1
  switch(n) = "--addtimestamp"
  label(n) = "If present the names of the log and output files will contain time stamp."
  n=n+1
  switch(n) = "--notimestamp"
  label(n) = "If present the names of the log and output files will not contain time stamp."
  n=n+1
  switch(n) = "--colT1"
  label(n) = "Analysis only: Determines the position (column) of the pretest"//&
    " statistic. Default value: "//i2cp(def_cT)//"."
  n=n+1
  switch(n) = "--colP1"
  label(n) = "Analysis only: Determines the position (column) of the pretest"//&
    " p-value. Default value: "//i2cp(def_cTpval)//"."
  n=n+1
  switch(n) = "--colT2"
  label(n) = "Analysis only: Determines the position (column) of the main test"//&
    " statistic. Default value: "//i2cp(def_cAS)//"."
  n=n+1
  switch(n) = "--colP2"
  label(n) = "Analysis only: Determines the position (column) of the main test"//&
    " p-value. Default value: "//i2cp(def_cASpval)//"."
  switch(n) = "--prevalence"
  label(n) = "Simulation only. Specifies the disease prevalence."//&
    " Possible values: real number between 0 and 1."
  n=n+1
  switch(n) = "--simulationmodel"
  label(n) = "Simulation only. Specifies the interaction model used."
  alias(n) = "--smodel"
  n=n+1
  switch(n) = "--OR"
  label(n) = "Simulation only. Specifies the interaction odds ratio."//&
    " Possible values: positive real number."
  n=n+1
  switch(n) = "--OR1"
  label(n) = "Simulation only. Specifies the first main effect odds ratio."//&
    " Possible values: positive real number."
  n=n+1
  switch(n) = "--OR2"
  label(n) = "Simulation only. Specifies the second main effect odds ratio."//&
    " Possible values: positive real number."
  n=n+1
  switch(n) = "--LD"
  label(n) = "Simulation only. Specifies the amound of LD between"//&
    " simulated pairs of loci (their number is specified by --nLDpairs."
  n=n+1
  switch(n) = "--nocausalpair"
  label(n) = "Simulation only. If present all simulated loci will be"//&
                   " neutral."
  n=n+1
  switch(n) = "--nloci"
  label(n) = "Simulation only. Specifies the total number of causal"//&
                   " and neutral loci."
  n=n+1
  !switch(n) = "--nneutralloci"
  !label(n) = "Simulation only. Specifies the number of neutral loci."
  !n=n+1
  switch(n) = "--nLDpairs"
  label(n) = "Simulation only. Specifies the number of pairs of"//&
    " loci that are in LD specified by --LD."
  n=n+1
  switch(n) = "--maf"
  label(n) = "Simulation only. Specifies the minor allele frequency"//&
    " (MAF) for data simulation mode (--simulateinput)."
  n=n+1
  switch(n) = "--fixedmaf"
  label(n) = "Simulation only. When present the value given by"//&
    " --maf is the fixed MAF for simulation. Use only with --simulateinput"
  switch(n) = "--no-openfile-iostat"
  label(n) = "Debugging: If present files will be open without IOSTAT."//&
    " Helps when debugging file open errors."
    
  
  CALL PrintProgramHeader(usto)
  CALL Prnt("List of available arguments:", skip2=1)

  !! Find the longest switch
  swlen = 0
  DO i=1,n
    swlen = MAX(swlen, LEN_TRIM(switch(i)))
  ENDDO
  
  !! Print available switchers and help
  j = 0
  DO i=1,n
    j = j+1
    IF(j>help_limit) THEN
      WRITE(*,'(A)') ""
      WRITE(*,'(A)', ADVANCE='NO') " Press enter to continue ..."
      READ(*,'(A)')      
      WRITE(*,'(A)') ""
      j = 0
    ENDIF

    !! Assign what to print
    text = switch(i)
    text = emptysp(1:nlead)//text(1:swlen)//emptysp(1:lenempty)//label(i)

    !! Print help item
    CALL Prnt0(text, lead=nlead+swlen+lenempty)

    !! If the switch has aliases, print them
    IF(LEN_TRIM(alias(i))>0) THEN
      text = emptysp(1:nlead+swlen+lenempty)//"Aliases: "//TRIM(alias(i))
      CALL Prnt0(text, lead=nlead+swlen+lenempty)
    ENDIF

  ENDDO

  !! Print signature
  CALL Prnt("Author: Jakub Pecanka, VU University Amsterdam", skip1=1)
  CALL Prnt("Please send comments and report bugs to "//&
                 TRIM(contact_email), skip2=1)

  IF(ALLOCATED(words)) DEALLOCATE(words)
  
  RETURN
    
END SUBROUTINE PrintHelpScreen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE EPI_INIT
