!! ====================================================================== !!
!! ===                                                                === !!
!! ====                          EpiDetector                         ==== !!
!! ====                              by                              ==== !!
!! ====                         Jakub Pecanka                        ==== !!
!! ===                                                                === !!
!! ====================================================================== !!

PROGRAM EpiDetector 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! INITIALIZATION OF VARIABLES AND MODULES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE EPI_RESULTS                       !! ROUTINES FOR ANALYZING OUTPUT
  USE EPI_TEST                          !! DATA TESTING ROUTINES
  !USE EPI_NONTESTING_PROCEDURES        !! NON-TESTING ROUTINES
  !USE EPI_DATA_OUTPUT_INPUT            !! INPUT/OUTPUT ROUTINES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PARAMETERS AND VARIABLES INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  
  INTEGER :: ncmd, exit_code, nrun, nrun_OR, nrun_delta
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PROGRAM BODY START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !! Get the global start time (for purposes of checking maximum runtime)
  CALL GetCurrentTime(time_h=starttime_glob)

  !! Initialize the counter of reruns
  nrun = 0
  nrun_delta = 0
  nrun_OR = 0
  
  !! Initialize random seed
  CALL SetSeed(seed=seed, useed_int=useed_int, useed_par=useed_par, &
               start=seed_start, const=seed_const, init=.TRUE., get=.TRUE.)
               
  !! Assign intial values for global variables
  CALL InitGlobalVars()
  
  comment = def_comment  ! variable with the comment sign for input files
  user_cmdline = ""
  
  !! Point of return to process cmd file or if no termination is true and
  !! user_cmdline contains --cmdfile argument
  1000 CONTINUE   !! Point of return: run again WITH re-reading cmdfile
  
  premat_halt = .FALSE.
  overwrite_files = .FALSE.
  ncmd = 0

  !! Check whether --cmdfile is present among command line arguments
  CALL CheckCmdfile(cmd_file, cmdfile_present, cmdlines, ncmdlines, comment, user_cmdline)

  !! Point of return to process cmd file or if no termination is true
  2000 CONTINUE   !! Point of return: run again WITHOUT re-reading cmdfile
  
  nrun = nrun + 1
  !! If rerun, announce additional analysis, then switch rerun indicator
  IF(nrun > 1) THEN
    IF(pause_run) &
      CALL Prnt("Execution paused. Press enter to continue ...", skip1=1, &
                advance='NO', wait=.TRUE.)
    CALL Prnt("Starting analysis number "//TRIM(i2cp(nrun))//" ...", & 
              log=.FALSE., skip1=1, skip2=1)
  ENDIF

  premat_halt = .FALSE.

  !! Get Execution start time and timestamp
  CALL GetCurrentTime(starttime_full, starttime)
  CALL GetTimeStamp(timestamp, starttime_full)

  log_file = ""

  !********************************************************************!
  !!      REMEMBER VALUES OF DATA VARIABLES FROM PREVIOUS RUNS        !!
  !********************************************************************!
  IF(simulate_data .OR. ALLOCATED(ped_file) .OR. ALLOCATED(bed_file)) &
    CALL RememberData()

  !! Read the next line from the command file (if any left)
  cmdline = ""
  IF(cmdfile_present .AND. ncmd < ncmdlines) THEN
    ncmd = ncmd + 1
    cmdline = cmdlines(ncmd)
  ENDIF

  !! Merge user command line and command line read from cmdfile
  IF(priority_user_cmdline) THEN
    cmdline = TRIM(cmdline)//" "//TRIM(user_cmdline)
  ELSE
    cmdline = TRIM(user_cmdline)//" "//TRIM(cmdline)
  ENDIF
  
  !! Assign initial values for global variables
  CALL InitGlobalVars()

  exit_code = 0

  !! Get command line arguments and set various options and settings
  CALL GetCmdLineArgs(cmdline, nrun, nrun_delta, nrun_OR, exit_code)
  
  !! If run_analysis is true call subroutine to analyze given output files
  IF(run_analysis) GOTO 4000

  !********************************************************************!
  !!                SET FILE NAMES AND CREATE LOG FILE                !!
  !********************************************************************!

  !! Check current input (binary) ped files and if the names are the same
  !! and if they are, keep the old data
  IF(.NOT.no_data_input) CALL CheckForSameInputFiles()

  !! Set the output and log file names (it is done here because these names may
  !! contain information that has just been recalled from previous runs
  CALL SetFilenames()

  !! Create log file (normal creation, this is where it should be created)
  CALL CreateLogFile(log_file, check_fexist, starttime_full, timestamp)

  !! Print execution start time
  CALL Report(exec_start=.TRUE., starttime=starttime_full)
  !! Print all submitted command line arguments
  CALL Report(cmdline_call=.TRUE.)
  !! Print used arguments
  CALL Report(used_args=.TRUE.)

  !! Re-initialize and announce random seed in case user specified it
  CALL SetSeed(seed=seed, useed_int=useed_int, useed_par=useed_par, &
               start=seed_start, const=seed_const, init=.TRUE., get=.TRUE., &
               nthreads=NUMBER_OF_THREADS, announce=.TRUE.)
               
  !! If no_data_input is true skip all data input and testing
  IF(no_data_input) GOTO 2500

  !! If a conversion of old tmp file to out file is requested, do it and quit
  IF(tmp_to_out_only) THEN
  
    CALL Prnt("Converting temporary result file ["//TRIM(tmp_file(1))//&
              "] to output file ["//TRIM(out_file)//"] ...")
    
    !! Make sure the file exists (was given by user properly)
    CALL CheckFileExistence(tmp_file(1))
    
    !! Do the conversion
    CALL WriteOutputFile(tmp_file)
    GOTO 2500
    
  ENDIF
      
  !********************************************************************!
  !!              PREPARE VARIABLES THAT CONCERN DATA INPUT           !!
  !********************************************************************!

  !! Make sure the input files exist
  CALL CheckMissingFiles()
  
  !! Make sure data variable sizes are ready
  CALL PrepareDataVariables(exit_code)
  IF(exit_code /= 0) GOTO 3000

  !! Print input options
  CALL Report(input_opts=.TRUE.)
  
  CALL FlushOutput((/usto, ulog/))
 
  !! Read-in or simulate genetical data 
  CALL GetData()
  
  !! Process the read-in data for loci exclusion based on MAF or other factors
  !! (Includes the calculation of MAF) 
  CALL ProcessData(exit_code)
  IF(exit_code /= 0) GOTO 3000
  
  !********************************************************************!
  !     SAVE THE INPUT DATA INTO A FILE BEFORE RUNNING THE ANALYSIS    !
  !********************************************************************!

  !! Save input data of maf information (only if new)
  IF(.NOT.reuse_data .OR. nrun==1) CALL SaveData()

  !! Print input data report    
  CALL Report(data_summary=.TRUE.)

  !********************************************************************!
  !  THE FOLLOWING PROCEDURE PERFORMS (TWO-STEP) GENOME WIDE ANALYSIS  ! 
  !********************************************************************!

  CALL FlushOutput((/usto, ulog/))
  
  IF(.NOT.no_testing) THEN
  
    !! Check for large arrays and announce a hint (OPENMP only)
    CALL Report(memory_hint=.TRUE.)
    
    CALL ResizeVar(sts1, ped_ss)
    CALL ResizeVar(sex1, ped_ss)
    CALL ResizeVar(pss_AS1, ped_ss)
    CALL ResizeVar(pss_DS1, ped_ss)
    sts1 = sts(1:ped_ss)
    sex1 = sex(1:ped_ss)
    pss_AS1 = pss_AS(1:ped_ss)
    pss_DS1 = pss_AS(1:ped_ss)
    
#ifdef _OPENMP
    !! Set the number of threads for parallel computing
    CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)
#endif

    !! Perform the actual testing
    CALL DoTests(X, sts1, sex1, ped_ss, ped_nco, ped_nca, ped_nsamp, dAS, dDS, &
                 pss_AS1, pss_DS1, temp_cycle, MAPA, MAF, InclLoc, &
                 ana_model_S1, ana_model)
                 !ana_model_S1, ana_model, seed_user, seed_start, seed_const)
    
  ENDIF
                      
  !********************************************************************!
  !!  ALL TESTS ARE PERFORMED AT THIS POINT AND TESTING IS FINISHED   !!
  !********************************************************************!
  
  2500 CONTINUE
  
  !! Sort the output file (if no more than 1 output file produced)
  IF(out_sort .AND. out_available .AND. .NOT.out_zipped .AND. (out_nfiles==1 .OR. &
  out_sort_multiple_files)) THEN
    CALL SortFile(sort_file, temp=out_sort_tmp, announce=.TRUE., skip1=.TRUE.)
    IF(.NOT.out_sort_tmp) CALL WriteOutputFile(sort_file)
  ENDIF

  !! Compute Total Running Time of the Computations
  CALL GetCurrentTime(stoptime_full, stoptime)
  runtime_total = GetTimeDifference(starttime, stoptime)

  !! Print execution stop time
  CALL Report(exec_end=.TRUE., starttime=starttime_full, stoptime=stoptime_full, &
                   runtime=runtime_total)
                   
  !********************************************************************!
  !!        DECIDE WHETHER THE PROGRAM SHOULD START FROM TOP          !!
  !********************************************************************!

  !! Get current total runtime and check global runtime against max_runtime
  CALL GetCurrentTime(time_h=current_time)
  IF(max_runtime>zero .AND. current_time - starttime_glob >= max_runtime) THEN
    CALL PrntW("User set maximum runtime of "//TRIM(r2c(max_runtime))//" hours exceeded!", &
               skip1=1, skip2=1)
    nrun = HUGE(nrun)
  ENDIF

  !! If the current parameter setting should be repeated 'max_reruns' times
  !! start from top with ncmd decreased by 1 so that the current cmd line gets 
  !! repeated
  IF(nrun < MIN(max_reruns, MAX_REPEATS)) THEN
    ncmd = MAX(0, ncmd - 1)
    GOTO 2000   !! 2000: Run again without re-reading cmdfile
  ENDIF

  !! Otherwise process the next line in cmd_line and start from top
  IF(cmdfile_present .AND. ncmd < ncmdlines) GOTO 2000
  
  !! No repetition selected, which means that the program will terminate (if
  !! no_halt is false) or asks the user for more commands (if no_halt is true)
  3000 CONTINUE 

  !! If no_halt is true, keep running and ask for new parameters
  IF(no_halt) THEN
    CALL AskForCmdline(user_cmdline, exit_code)
    !! Evaluate exit codes
    IF(exit_code == 1)  GOTO 1000   !! 1000: Run again including reading cmdfile
    IF(exit_code == 2)  GOTO 2000   !! 2000: Run again without reading cmdfile
    IF(exit_code == -1) ask_halt = .FALSE.  !! No new user input -> quit
  ELSE
    ask_halt = .FALSE.
  ENDIF 
  
  !! Print execution stop time
  CALL Terminate(ask=ask_halt, premature=premat_halt)

  !! If Terminate did not stop the program, start from the top
  GOTO 1000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                              ANALYSIS MODE                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  4000 CONTINUE

  !! First set output file names
  CALL SetFilenames()

  !! Check for log and output filenames
  IF(.NOT.no_log .AND. log_file=="") THEN
    IF(results_nfiles>0) THEN
      log_file=TRIM(res_file(1))//"_analysis"//log_ext
    ELSE
      log_file="EpiDetector_"//TRIM(timestamp)//"_analysis"//log_ext
    ENDIF
    !CALL CreateLogFile(log_file)
  ENDIF
  
  !! Create log file
  CALL CreateLogFile(log_file, check_fexist, starttime_full, timestamp)

  !! Throw an error if no files on input
  IF(results_nfiles<=0) &
    CALL PrntE("No input files to analyze.", skip1=1, premature=.FALSE., Q=.TRUE.)

  !! Print execution start time
  CALL Report(exec_start=.TRUE., starttime=starttime_full)

  !! Re-Initialize random seed in case user specified the random seed
  CALL SetSeed(seed=seed, start=seed_start, const=seed_const, init=.TRUE., &
               get=.TRUE., announce=.TRUE.)
               
  !! Run Analysis (the only time the program gets here)
  CALL AnalyzeResults(res_file, out_header)

  !! Print execution stop time
  CALL GetCurrentTime(stoptime_full, stoptime)
  runtime_total = GetTimeDifference(starttime, stoptime)
  CALL Report(exec_end=.TRUE., stoptime=stoptime_full)
  CALL Dealloc(res_file)
  CALL FlushOutput((/usto, ulog/))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PROGRAM BODY END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM EpiDetector
