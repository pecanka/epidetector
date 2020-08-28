MODULE EPI_PARAMS
!! This file is intended for inclusion whenever the below variables and 
!! constants are needed. The idea is to have them defined in one place only
!! for convenience of changing their values.

  !USE PREPROCESSED

!! def_ accessibility is PUBLIC  
  PUBLIC

  !! PROGRAM NAME AND VERSION
  CHARACTER(11), PARAMETER    :: program_name = "EpiDetector"
  CHARACTER(8), PARAMETER     :: program_version = "0.1 beta"
  CHARACTER(50), PARAMETER    :: contact_email = "EpiDetector@pecanka.net"
  
  !! PROGRAM FULL NAME
  INTEGER, PARAMETER          :: fnl = LEN(program_name)+LEN(program_version)+1
  CHARACTER(fnl), PARAMETER   :: program_fullname = program_name//" "//&
                                                    program_version
  
  !! Define a constant which identifies the system (taken from 
  !! http://sourceforge.net/p/predef/wiki/OperatingSystems/ or
  !! http://nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system
  !! List of predefined macros can be obtained by calling: gcc -dM -E - < /dev/null

#if defined(__linux__)
  CHARACTER(15), PARAMETER    :: sys_type = "GNU Linux"
  CHARACTER(10), PARAMETER    :: sys_bit = "32/64bit"
  CHARACTER(15), PARAMETER    :: sys_envir = ""
#elif defined(__MINGW64__)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "64bit"
  CHARACTER(15), PARAMETER    :: sys_envir = "with MinGW"
#elif defined(__MINGW32__)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "32bit"
  CHARACTER(15), PARAMETER    :: sys_envir = "with MinGW"
#elif defined(__CYGWIN64__)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "64bit"
  CHARACTER(15), PARAMETER    :: sys_envir = "with Cygwin"
#elif defined(__CYGWIN32__)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "32bit"
  CHARACTER(15), PARAMETER    :: sys_envir = "with Cygwin"
#elif defined(__CYGWIN__)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "32bit"
  CHARACTER(15), PARAMETER    :: sys_envir = "with Cygwin"
#elif defined(_WIN64)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "64bit"
  CHARACTER(15), PARAMETER    :: sys_envir = ""
#elif defined(_WIN32)
  CHARACTER(15), PARAMETER    :: sys_type = "Windows"
  CHARACTER(10), PARAMETER    :: sys_bit = "32bit"
  CHARACTER(15), PARAMETER    :: sys_envir = ""
#else                   
  CHARACTER(15), PARAMETER    :: sys_type = "unknown system"
  CHARACTER(10), PARAMETER    :: sys_bit = ""
  CHARACTER(15), PARAMETER    :: sys_envir = ""
#endif

  !! Default type for REALs
  INTEGER, PARAMETER :: &
    dpp = SELECTED_REAL_KIND(12,60), &
    !dpp = SELECTED_REAL_KIND(16, 60)
    !dpp = SELECTED_REAL_KIND(13)     ! returns 8
    !dpp = KIND(1.0D0)     ! returns 8
    dpd = SELECTED_REAL_KIND(12,60)     ! returns 8
    !dpd = SELECTED_REAL_KIND(13)     ! returns 8

  !! Type of REALs used by cdflib (must be the same as in biomath_constants_mod)
  !INTEGER, PARAMETER :: dpcdf = KIND(0.0D0)
  !INTEGER, PARAMETER :: dpcdf = SELECTED_REAL_KIND(13)
  
! Note: 'dpp' does not have to be the same as 'dp' set by module LSQ (called by  
!       Logistic). To give at least 12 decimal digit representation of floating  
!       point numbers, set it to SELECTED_REAL_KIND(12, 60). This should be  
!       adequate for most problems except the fitting of polynomials. dpp and dp
!       are being set this way so that the same code can be run on PCs and Unix 
!       systems, which will usually represent floating-point numbers in 'double
!       precision', and other systems with larger word lengths which will
!       give similar accuracy in `single precision'.
  
  INTEGER, PARAMETER :: &
    iks = 1, &
    ikn = 4, &
    !ikl = 16, &
    !ik = SELECTED_INT_KIND(1), &
    !ikn = SELECTED_INT_KIND(4), &
    !ikl = SELECTED_INT_KIND(8), &
    !ikb = SELECTED_INT_KIND(16), &
    ikb = 8, &
    iks_len = 5, &
    ikn_len = 7, &
    ikb_len = 21, &
    !ik = ikn, &
    dmF = 4, &                     ! Dimension of beta, Fisher, score, ...
    dmV = 9, &                     ! N of two-locus genotypes (determines dimensions of many variables) 
    dmS = 3, &                     ! Number genotypes (determines dimensions of many variables)
    sDebug = 16, &                 ! Size of debug info vector (must be at least 7!!!)
    MAX_REPEATS = 1000000              

  !! ERROR NUMBERS
  !! If an event should not stop the execution of the program, give it a negative 
  !! error number. Then all statistics will get computed and the error will be
  !! considered more of a warning in the output file (such as zero marginals).
  INTEGER, PARAMETER :: &

    err_not_done = -1, &           ! Error code for test not performed
    err_low_cell_cor = -4, &       ! Error code for low contingency table cell count after correction (TestPhase1)
    err_zero_margin = -5, &        ! Error code for zero contingency table cell count (TestPhase1)
    err_no_data = -100, &
    err_maf_big = -311, &
    err_maf_sml = -312, &

    err_low_margin = 5, &          ! Error code for low contingency table marginal count (TestPhase1)    
    err_lin_dep = 6, &             ! Error code for logistic regression linear dependence problem
    err_low_cell = 4, &            ! Error code for low contingency table cell count (TestPhase1)
    err_low_ss_S1 = 101, &         ! Error code for low sample size in S1   
    err_low_ss_S2 = 102, &         ! Error code for low sample size in S1
    err_no_co = 201, &
    err_no_ca = 202, &
    err_no_cc = 203, &
    err_low_co = 211, &
    err_low_ca = 212, &
    err_low_cc = 213, &
    
    err_all_zeros = 888, & 
    err_low_var_AS = 801, & 
    err_neg_var_AS = 901, &
    err_low_var_DS = 802, &
    err_neg_var_DS = 902, &
    err_low_var_CS = 803, &
    err_neg_var_CS = 903, &
      
    errs_low_ss(8) = (/ &          ! A vector of all sample size related error codes (for simple checking)
       err_low_ss_S1, &       
       err_low_ss_S2, &
       err_no_cc, &
       err_no_co, &
       err_no_ca, &
       err_low_cc, &
       err_low_co, &
       err_low_ca/)
                               
  !! I/O ERROR NUMBERS
  INTEGER, PARAMETER :: &
    err_eof(3) = (/-1, 3, 5002/), &
    err_peof(5) = (/36, 39, 213, 5001, 5003/), &
    err_eol(1) = (/-2/), &
    err_nonexistfile(1) = (/2/), &
    err_openfile(1) = (/13/), &
    err_unclosedfile(2) = (/5002, 104/)
            !! ERROR NUMBERS NOTES
            !! eof: 5002 in gfortran
            !! past eof: 5003 gfortran, 213 g95, 36, 39 ifort
            !! reading non-existant record: 36
            !! error during read: 39

!       Gfortran Error Codes
!       
!       When using READ in fortran error codes can returned in the value IOSTAT. 
!       These values are
!       
!       -3    LIBERROR_FIRST = -3,
!       -2    LIBERROR_EOR = -2,
!       -1    LIBERROR_END = -1,
!       0     LIBERROR_OK = 0, 
!       5000  LIBERROR_OS = 5000, 
!       5001  LIBERROR_OPTION_CONFLICT
!       5002  LIBERROR_BAD_OPTION
!       5003  LIBERROR_MISSING_OPTION
!       5004  LIBERROR_ALREADY_OPEN
!       5005  LIBERROR_BAD_UNIT
!       5006  LIBERROR_FORMAT
!       5007  LIBERROR_BAD_ACTION
!       5008  LIBERROR_ENDFILE
!       5009  LIBERROR_BAD_US
!       5010  LIBERROR_READ_VALUE
!       5011  LIBERROR_READ_OVERFLOW
!       5012  LIBERROR_INTERNAL
!       5013  LIBERROR_INTERNAL_UNIT
!       5014  LIBERROR_ALLOCATION
!       5015  LIBERROR_DIRECT_EOR
!       5016  LIBERROR_SHORT_RECORD
!       5017  LIBERROR_CORRUPT_FILE
!       5018  LIBERROR_INQUIRE_INTERNAL_UNIT
!       
!       as defined in gcc/fortran/libgfortran.h. (Taken from gfortran 4.6)

  !! PROGRAM SETTINGS 
  LOGICAL, PARAMETER :: &
    show_segfault_hint = .TRUE., &  
    OMP_MULTISAMPLEOPTIMIZED = .FALSE.
       ! If true a warning will be printed about possible stack overflow with
       ! large arrays and parallel computation telling user how to fix it
       
  INTEGER(ikb), PARAMETER :: &
    max_file_size = HUGE(1)   ! File size cannot be bigger than this
  
  REAL(dpp), PARAMETER :: &
    mfsf = 0.9_dpp           ! max file size factor (must be between 0 and 1)

  REAL(dpd), PARAMETER :: &
    zero_d = 0.0_dpd, &
    one_d = 1.0_dpd, &
    two_d = 2.0_dpd, &
    ten_d = 10.0_dpd
  
  REAL(dpp), PARAMETER :: &
    pi = 3.141592653589793_dpp, &
    zero = 0.0_dpp, &
    quarter = 0.25_dpp, &
    half = 0.5_dpp,&
    one = 1.0_dpp, & 
    two = 2.0_dpp, &
    three = 3.0_dpp, &
    four = 4.0_dpp, &
    five = 5.0_dpp, &
    six = 6.0_dpp, &
    seven = 7.0_dpp, &
    eight = 8.0_dpp, &
    nine = 9.0_dpp, &
    ten = 10.0_dpp, &
    fifty = 50.0_dpp, &
    hundred = 100_dpp, &
    tenth = 0.1_dpp, &
    hundreth = 0.01_dpp, &
    thousandth = 0.001_dpp, &
    one3(1:3) = (/one, one, one/), &
    !epstol = ten**7*EPSILON(zero_d)       
    epstol = EPSILON(1.0)       
       !! Store negligable number for checking (almost) equality
       !! Note: If not EPSILON() multiplied by 100 pseudoinverse 
       !! sometimes complains about negative eigenvalues 

  !! UNITS AND VARIABLE LENGTHS
  INTEGER, PARAMETER :: &
    system_filename_limit = 255, &
    test_warn_limit = 10**8, &
    mem_warn_limit = 10**9, &
    mem_hint_limit = 10**7, &
    df_temp_cycle = 500000, &
    min_ntests_parallel = 100, &
    
    !! Various string lengths
    mfl = 1500, &        ! max file name length
    mal = 1000, &        ! max cmd line argument length
    mcl = 10000, &       ! max cmd line length
    mml = 50, &          ! max MAPA item length
    mmll= 200, &         ! max MAPA linelength
    mfml = 50, &         ! max FAM item length     
    mfmll = 500, &       ! max FAM file line length
    mifll = 50000, &     ! max input file line length
    mhtl = 10000, &      ! max huge text length
    mltl = 2000, &       ! max long text length
    mntl = 500, &        ! max neutral text length
    mmtl = 200, &        ! max medium text length
    mstl = 50, &         ! max short text length
    mttl = 20, &         ! max tiny text length
    mel = 10000, &       ! max empty space length
    mtsl = 24, &         ! max time string length
    mextl = 10, &        ! max file extension length
    moil = 20, &         ! max output item length
    mllc = 500, &        ! max log line count 
    
    !! Output units
    usto = 6, &
    ucmd = 10, &
    uped = 20, &
    ubed = 21, &
    umap = 30, &
    ufam = 31, &
    umaf = 33, &
    uout = 40, &
    utmp = 41, &
    ures = 42, &
    ulog = 50, &
    ufre = 11, &
    udel = 99, &
    
    !! Screen and log file widths
    log_width = 79, &    ! Sets the width of text in a log file
    head_width = 75, &   ! Sets the width of the program header (For standard console this should not exceed 79)
    sect_width = 40, &   ! Sets the width of the section header (For standard console this should not exceed 79)
    stdo_width = 999, &  ! Sets the width of text on screen (For standard console this should not exceed 79)
    report_width = 64    ! Sets the width of report tables (screen and log) (Should be at least 63!)

  !! PARAMETERS ETC.
  
  CHARACTER, PARAMETER :: &
    tab = char(9), &
    space = char(32), &
    colon = ":", &
    semicolon = ";", &
    comma = ",", & 
    def_comment = "#", &
    def_NA = "*", &
    def_ped_cs(2) = (/tab, space/), &
    def_map_cs(2) = (/tab, space/), &
    def_fam_cs(2) = (/tab, space/), &
    def_out_cs = tab, &
    def_delim = "\"

  CHARACTER(2), PARAMETER :: &
    lt = "> "            ! lt stores the text that is printed at the 
                         ! beginning of log file and screen lines by Prnt

  CHARACTER(mextl), PARAMETER :: &
    log_ext = ".log", &
    out_ext = ".epi", &
    ped_ext = ".ped", &
    bed_ext = ".bed", &
    map_ext = ".map", &
    bim_ext = ".bim", &
    fam_ext = ".fam", &
    maf_ext = ".frq", &
    prev_ext = ".prev", &
    pss_ext = ".pss", &
    res_ext = ".res", &
    tmp_ext = ".tmp"

  !! TEST NUMBERS AND MODELS
  INTEGER, PARAMETER :: &
    T1sc = 1, &
    T1cc = 2, &
    T1co = 3, &
    T1ca = 4, &
    T1po = 5, &
    T1cc_df = 4, &
    T1co_df = 4, &
    T1ca_df = 4, &
    T1po_df = 4, &
    T2sc = 1, &
    T2cc = 2, &
    def_T1 = T1ca, &
    def_T2 = T2sc, &
    def_ana_model = 1, &
    def_sim_model = 1, &
    def_map_ncols = 4, &
    def_bim_ncols = 6, &
    def_map_chcol = 1, &
    def_map_rscol = 2, &
    def_map_dscol = 3, &
    def_map_bpcol = 4, &
    def_map_a1col = 5, &
    def_map_a2col = 6, &
    def_sub_chcol = 1, &
    def_sub_rscol = 2, &
    def_auto_d_n = 25, &
    def_auto_d_nlevel = 100, &
    def_auto_d_ntests = 1000, &
    def_auto_d_maxnreps = 5
    
  CHARACTER, PARAMETER :: &
    def_sort_sel = "A", &      ! Default statistics by whose p-values the output file is ordered
    sort_selections(5) = (/"A","D","C","P","0"/) 
       ! Possible values for sort selection of pvalues: 
       !   A (adjusted), D (disjoint), C (classical), P (pooled), 0 (no sort)
  
  INTEGER :: imm
    
  REAL(dpp), PARAMETER :: & 
    NAp = nine, &                         ! NA stand-in for p-value (must be bigger than 2!)
    NAv = HUGE(dpp), &                    ! NA stand-in for statistic (0 seems safest for this)
    NAv2 = -123456.0_dpp, &               ! Alternative NA stand-in for statistic (in case 0 does not fit)
    MINvalue = 1e-99_dpp, &               ! The value of 0 is safest for this
    NAneg = -nine, &                      ! Must be negative
    NApos = nine, &                       ! Must be positive (and should be bigger than 1)
    rerun_eps = 1e-6_dpp, &               ! "Epsilon" for delta
    def_delta = half, &                   ! Default S1 portion (no S1)
    def_level1 = 0.01_dpp, &              ! Default S1 level
    def_level2 = 0.05_dpp, &              ! Default main test level
    minimum_test_level = 1.0E-10_dpp, &   ! Minimum allowed level of significance for a test
    maximum_test_level = 0.999_dpp, &     ! Maximum allowed level of significance for a test  
    min_OR = hundreth**2, &               ! Minimum allowed odds ratio in simulation
    max_OR = hundred, &                   ! Maximum allowed odds ratio in simulation
    def_simul_prev = hundreth, &
    def_auto_d_min = zero, &
    def_auto_d_max = 0.95_dpp, &
    def_auto_d_minlevel = 1e-7_dpp, &
    def_auto_d_maxlevel = one, &
    def_auto_d_COR = 1e8_dpp, &
    def_auto_d_prev = hundreth, &
    def_auto_d_OR1 = one, &
    def_auto_d_OR2 = one, &
    def_auto_d_OR = 1.4_dpp, &
    def_auto_d_maf1 = 0.25_dpp, &
    def_auto_d_maf2 = 0.25_dpp, &
    def_auto_level_prev = hundreth, &
    def_auto_level_OR1 = one, &
    def_auto_level_OR2 = one, &
    def_auto_level_OR = 1.5_dpp
       
  INTEGER, PARAMETER :: &
    help_limit = 20, &
    def_fidcol = 1, &
    def_iidcol = 2, &
    def_pidcol = 3, &
    def_midcol = 4, &
    def_sexcol = 5, &
    def_stscol = 6, &
    def_famncol = 6, &
    delta_precision = 4, &
    sex_chr(2) = (/23, 24/), &
    autosomal_chr = 25, &
    mitochondrial_chr = 26, &                   
    !special_chr(5) = (/0, &
    special_chr(4) = (/sex_chr, &
                       autosomal_chr, &
                       mitochondrial_chr/)
                       
  CHARACTER(2), PARAMETER :: &
    chr_name(0:26) = (/"0 ", "1 ", "2 ", "3 ", &
                       "4 ", "5 ", "6 ", "7 ", &
                       "8 ", "9 ", "10", "11", &
                       "12", "13", "14", "15", &
                       "16", "17", "18", "19", &
                       "20", "21", "22", "X ", &
                       "Y ", "XY", "MT"/)
  
  !! INPUT DATA PARAMETERS
  INTEGER(iks), PARAMETER :: &
    NumNA = -1_iks, & 
    NumExcl = -2_iks, &
    NumEmpty = -3_iks, &
    NumCo = 0_iks, &
    NumCa = 1_iks, &
    NumMa = 1_iks, &
    NumFe = 2_iks

  CHARACTER, PARAMETER :: &
    def_CharNAgen = "0", &
    def_CharNAsts = "9", &
    def_CharNAsex = "9", &   
    def_CharCa = "2", &   
    def_CharCo = "1", &   
    def_CharMa = "1", &   
    def_CharFe = "2"

  !! FORMAT FOR STATISTICS AND PVALUES
  INTEGER, PARAMETER :: &
    def_ndig_stat = 3, &
    def_ndig_pval = 3, &
    def_ndig_mafs = 3, &
    opt_ndig_mafs = 2 
  
  !! Recommended relationship between field width 'W' and the number of  
  !! fractional digits 'D' in a format descriptor is 'W>=D+3'.
  CHARACTER(10), PARAMETER :: &
    fm_maf = 'F6.3', &
    fm_lvl = 'ES12.4'
    
  !! Define type for a temporary file record    
  TYPE TMPTYPE

    REAL(dpp)  :: PTAS4          ! Pretest statistic AS-PT-4
    REAL(dpp)  :: PTAS4p         ! Pretest p-value AS-PT-4                      
    REAL(dpp)  :: PTAS1          ! Pretest statistic AS-PT-1
    REAL(dpp)  :: PTAS1p         ! Pretest p-value AS-PT-1                      
    INTEGER    :: PTASNca        ! Counts of cases used by PTAS
    INTEGER    :: PTASNco        ! Counts of controls used by PTAS
    INTEGER(2) :: errPTAS4       ! Error indicator for AS-PT-4
    INTEGER(2) :: errPTAS1       ! Error indicator for AS-PT-1
    REAL(dpp)  :: PTDS4          ! Pretest statistic DS-PT-4
    REAL(dpp)  :: PTDS4p         ! Pretest p-value DS-PT-4                     
    REAL(dpp)  :: PTDS1          ! Pretest statistic DS-PT-1
    REAL(dpp)  :: PTDS1p         ! Pretest p-value DS-PT-1                      
    INTEGER    :: PTDSNca        ! Counts of cases used by PTDS
    INTEGER    :: PTDSNco        ! Counts of controls used by PTDS
    INTEGER(2) :: errPTDS4       ! Error indicator for DS-PT-4
    INTEGER(2) :: errPTDS1       ! Error indicator for DS-PT-1
    REAL(dpp)  :: PT4            ! Pretest control-4 statistic
    REAL(dpp)  :: PT4p           ! Pretest control-4 p-value                      
    REAL(dpp)  :: PT1            ! Pretest control-1 statistic
    REAL(dpp)  :: PT1p           ! Pretest control-1 p-value                      
    INTEGER(2) :: errPT4         ! Error indicator for PT4
    INTEGER(2) :: errPT1         ! Error indicator for PT1
    REAL(dpp)  :: PTPO4          ! Pretest pool-4 statistic
    REAL(dpp)  :: PTPO4p         ! Pretest pool-4 p-value                      
    REAL(dpp)  :: PTPO1          ! Pretest pool-1 statistic
    REAL(dpp)  :: PTPO1p         ! Pretest pool-1 p-value                      
    INTEGER(2) :: errPTPO4       ! Error indicator for PTPO4
    INTEGER(2) :: errPTPO1       ! Error indicator for PTPO1
    REAL(dpp)  :: AS4            ! Adjusted score statistic for pretest-4
    REAL(dpp)  :: AS4p           ! Adjusted score p-value for pretest-4 
    REAL(dpp)  :: AS1            ! Adjusted score statistic for pretest-1
    REAL(dpp)  :: AS1p           ! Adjusted score p-value for pretest-1
    INTEGER    :: ASNca          ! Counts of cases used by AS
    INTEGER    :: ASNco          ! Counts of controls used by AS
    INTEGER(2) :: errAS4         ! Error indicator for AS-4
    INTEGER(2) :: errAS1         ! Error indicator for AS-1
    REAL(dpp)  :: DS             ! Disjoint score statistic
    REAL(dpp)  :: DSp            ! Disjoint score p-value
    INTEGER    :: DSNca          ! Counts of cases used by DS
    INTEGER    :: DSNco          ! Counts of controls used by DS
    REAL(dpp)  :: CS             ! Classical score statistic
    REAL(dpp)  :: CSp            ! Classical score p-value
    INTEGER    :: CSNca          ! Counts of cases used by CS
    INTEGER    :: CSNco          ! Counts of controls used by CS
    INTEGER(2) :: errDS          ! Error indicator for DS
    INTEGER(2) :: errCS          ! Error indicator for CS
    REAL(dpp)  :: Tn4(dmV)       ! Pretest vector for chisquare-4
    REAL(dpp)  :: Tn1(dmV)       ! Pretest vector for chisquare-1
    REAL(dpp)  :: beta0(0:2)     ! Estimate of beta (under H_0)
    !REAL(dpp)  :: se_beta0(0:2)  ! Estimate of std. error of beta0
    REAL(dpp)  :: beta1(0:3)     ! Estimate of beta (general)
    REAL(dpp)  :: se_beta1(0:3)  ! Estimate of std. error of beta1
    !REAL(dpp)  :: MAF(4)         ! Two pairs of minor allele frequencies
    REAL(dpp)  :: PLevel(2)      ! Used and optimal pretest level
    REAL(dpp)  :: Lambda         ! Estimated non-centrality parameter of pretest statistic
    REAL(dpp)  :: Slope          ! Estimated slope of score statistics 
    REAL(dpp)  :: Debug(sDebug)  ! Debugging information (variance, etc.)
    INTEGER    :: Pos(2)         ! Position of the two loci in MAPA
    !INTEGER    :: Chr(2)         ! Chromosomes
    !INTEGER    :: Rs(2)          ! RS numbers
    LOGICAL    :: ASstd          ! Indicator of AS already standardized
    !LOGICAL    :: Signif         ! Indicator whether the current pair is known to be significant
    LOGICAL    :: Empty          ! Indicator of NA record
                 
    REAL(dpp)  :: AS_co_X, ASp_co_X, DS_co_X, DSp_co_X, PT_co_X, PTp_co_X, &
                  AS_co_Y, ASp_co_Y, DS_co_Y, DSp_co_Y, PT_co_Y, PTp_co_Y, &
                  AS_co_0, ASp_co_0, DS_co_0, DSp_co_0, PT_co_0, PTp_co_0, &
                  AS_cc, ASp_cc, DS_cc, DSp_cc, PT_cc, PTp_cc   

  END TYPE TMPTYPE
  
  !! Define type for a temporary file record in minimalistic format    
  TYPE TMPTYPEmin

    REAL(dpp)  :: PTAS4p         ! Pretest p-value AS-PT-4                      
    REAL(dpp)  :: PTAS1p         ! Pretest p-value AS-PT-1                      
    REAL(dpp)  :: PTDS4p         ! Pretest p-value DS-PT-4                     
    REAL(dpp)  :: PTDS1p         ! Pretest p-value DS-PT-1                      
    REAL(dpp)  :: PT4p           ! Pretest control-4 p-value                      
    REAL(dpp)  :: PT1p           ! Pretest control-1 p-value                      
    REAL(dpp)  :: PTPO4p         ! Pretest pool-4 p-value                      
    REAL(dpp)  :: PTPO1p         ! Pretest pool-1 p-value                      
    REAL(dpp)  :: AS4p           ! Adjusted score p-value for pretest-4 
    REAL(dpp)  :: AS1p           ! Adjusted score p-value for pretest-1
    REAL(dpp)  :: DSp            ! Disjoint score p-value
    REAL(dpp)  :: CSp            ! Classical score p-value
    INTEGER(2) :: errPTAS4       ! Error indicator for AS-PT-4
    INTEGER(2) :: errPTAS1       ! Error indicator for AS-PT-1
    INTEGER(2) :: errPTDS4       ! Error indicator for DS-PT-4
    INTEGER(2) :: errPTDS1       ! Error indicator for DS-PT-1
    INTEGER(2) :: errPT4         ! Error indicator for PTPO4
    INTEGER(2) :: errPT1         ! Error indicator for PTPO1
    INTEGER(2) :: errPTPO4       ! Error indicator for PTPO4
    INTEGER(2) :: errPTPO1       ! Error indicator for PTPO1
    INTEGER(2) :: errAS4         ! Error indicator for AS-4
    INTEGER(2) :: errAS1         ! Error indicator for AS-1
    INTEGER(2) :: errDS          ! Error indicator for DS
    INTEGER(2) :: errCS          ! Error indicator for CS
    REAL(dpp)  :: PLevel(2)      ! Used and optimal pretest level
    INTEGER    :: Pos(2)           ! Position of the first and second loci in MAPA

  END TYPE TMPTYPEmin

  !! Define empty record parameter
  TYPE(TMPTYPE), PARAMETER :: NAtmp = TMPTYPE(  &
    NAv, &     ! AS-PT4        
    NAp, &      ! AS-PT4p   
    NAv, &     ! AS-PT1        
    NAp, &      ! AS-PT1p   
    -1, &          ! Nca PTAS
    -1, &          ! Nco PTAS
    -1, &          ! Error AS-PT4
    -1, &          ! Error AS-PT1
    NAv, &     ! DS-PT4        
    NAp, &      ! DS-PT4p   
    NAv, &     ! DS-PT1        
    NAp, &      ! DS-PT1p   
    -1, &          ! Nca PTDS
    -1, &          ! Nco PTDS
    -1, &          ! Error DS-PT4
    -1, &          ! Error DS-PT1
    NAv, &     ! PT4        
    NAp, &      ! PT4p
    NAv, &     ! PT1        
    NAp, &      ! PT1p
    -1, &          ! Error PT4
    -1, &          ! Error PT1
    NAv, &     ! PTPO4        
    NAp, &      ! PTPO4p
    NAv, &     ! PTPO1        
    NAp, &      ! PTPO1p
    -1, &          ! Error PTPO4
    -1, &          ! Error PTPO1
    NAv, &     ! AS4   
    NAp, &      ! AS4p  
    NAv, &     ! AS1   
    NAp, &      ! AS1p  
    -1, &          ! Nca AS
    -1, &          ! Nco AS
    -1, &          ! Error AS4
    -1, &          ! Error AS1
    NAv, &     ! DS   
    NAp, &      ! DSp  
    -1, &          ! Nca DS
    -1, &          ! Nco DS
    NAv, &     ! CS
    NAp, &      ! CSp
    -1, &          ! Nca CS
    -1, &          ! Nco CS
    -1, &          ! Error DS
    -1, &          ! Error CS
    NAv, &     ! Tn4
    NAv, &     ! Tn1
    NAv, &     ! beta0
    !-one, &        ! se_beta0
    NAv, &     ! beta1
    -one, &        ! se_beta1
    !NAneg, &       ! MAF
    NAneg, &       ! PLevel
    NAneg, &       ! Lambda
    NAneg, &       ! Slope
    NAv, &     ! Debug
    -1, &          ! Pos
    !-1, &          ! Chr
    !-1, &          ! Rs
    !.FALSE., &     ! Signif
    .FALSE., &     ! ASstd
    .TRUE., &      ! Empty
    !
    NAv, &     ! AS_co_X        
    NAp, &      ! etc.   
    NAv, &             
    NAp, &          
    NAv, &             
    NAp, &         
    !
    NAv, &     ! AS_co_Y        
    NAp, &         
    NAv, &             
    NAp, &          
    NAv, &             
    NAp, &         
    !
    NAv, &     ! AS_co_0        
    NAp, &         
    NAv, &             
    NAp, &          
    NAv, &             
    NAp, &         
    !
    NAv, &     ! AS_cc        
    NAp, &         
    NAv, &             
    NAp, &          
    NAv, &             
    NAp )        

  !! Define empty record parameter
  TYPE(TMPTYPEmin), PARAMETER :: NAtmpmin = TMPTYPEmin(  &
    NAp, &      ! AS-PT4p   
    NAp, &      ! AS-PT1p   
    NAp, &      ! DS-PT4p   
    NAp, &      ! DS-PT1p   
    NAp, &      ! PT4p
    NAp, &      ! PT1p
    NAp, &      ! PTPO4p
    NAp, &      ! PTPO1p
    NAp, &      ! AS4p  
    NAp, &      ! AS1p  
    NAp, &      ! DSp  
    NAp, &      ! CSp
    -1, &          ! Error AS-PT4
    -1, &          ! Error AS-PT1
    -1, &          ! Error DS-PT4
    -1, &          ! Error DS-PT1
    -1, &          ! Error PT4
    -1, &          ! Error PT1
    -1, &          ! Error PTPO4
    -1, &          ! Error PTPO1
    -1, &          ! Error AS4
    -1, &          ! Error AS1
    -1, &          ! Error DS
    -1, &          ! Error CS
    NAneg, &       ! PLevel
    -1 &           ! Pos 
    )        

END MODULE EPI_PARAMS