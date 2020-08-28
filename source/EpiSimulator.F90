PROGRAM EpiSimulator

  USE UTILITIES
  USE SIMULATION

! Generate a case-control sample and writes the output into a pedfile and a 
! mapfile in PLINK or other format format assuming:
!
!    1-biallelic loci
!    2-random mating (i.e. Hardy and Weinberg equilibrium)
!    3-user-specified amount of LD between the first two loci
!    4-user specified interaction effect between the first two loci for 
!      determining the probability that the individual is case vs. control
!    5-choice between:
!       1-a uniform allele frequency distribution with a lower boundary for MAF
!       2-a fixed minor allel frequency (same for both "causal" loci)
! 
!   Note: If the user would wish to allow the major alleles to be the "causal" 
!         alleles than uncomment lines ?? and ?? in the code ... !!!
! 
! In addition to the first two loci that may be in LD (in both cases and 
! controls) and may contribute to the susceptibility for the disease.
! The user can specify additional markers to be simulated which are completely 
! random, i.e. unlinked and no effect on the phenotype
! 
! Also note that some defaults and limitations are built in for the minimal 
! sample size and allele frequency such that the EXPECTED number of individuals 
! with the least frequent two-locus genotype should be at least five. The 
! actual number may differ!


  IMPLICIT NONE

  CHARACTER(12), PARAMETER    :: program_name = "EpiSimulator"
  CHARACTER(8), PARAMETER     :: program_version = "0.0.1"
  REAL(dpp), PARAMETER        :: min_freq = 0.01_dpp
  INTEGER, PARAMETER          :: min_ss = 80
  INTEGER, ALLOCATABLE        :: seed(:)
  CHARACTER, ALLOCATABLE      :: cases(:,:), controls(:,:)
  INTEGER(ik), ALLOCATABLE    :: X(:,:)

  CHARACTER       :: model
  CHARACTER(mfl)  :: pedfile, mapfile, filename, text
  INTEGER         :: nsamples, nneutralloci, nLDpairs, ss, &
                     isample, prev_percent, nca, nco, nloci   
  REAL(dpp)       :: allelefreq, RR, RR1, RR2, prevalence, LD, rnum, & 
                     starttime, stoptime
  CHARACTER(24)   :: starttime_full, stoptime_full, runtime_total
  LOGICAL         :: appendfile, onlygenes, fixedallelefreq, nomapfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PROGRAM BODY                                                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  no_log = .TRUE.
  !! Set and store random seed
  CALL SetSeed(seed, initialize=.TRUE., get=.TRUE.)
  
  !! Print header
  CALL Prnt("", empty=.TRUE.)
  CALL Prnt("***********************************", empty=.TRUE.)
  CALL Prnt(">>> "//TRIM(program_name)//", version "//&
                   TRIM(program_version)//" <<<", empty=.TRUE.)
  CALL Prnt("***********************************", empty=.TRUE.)
  CALL Prnt("", empty=.TRUE.)

  !! Get start time
  CALL GetCurrentTime(starttime_full, starttime)

  !! Get parameters from cmd line
  CALL GetCmdLineArgs()

#ifdef _OPENMP
  !! Set the number of threads for parallel computing
  CALL OMP_SET_NUM_THREADS(NUMBER_OF_THREADS)
#endif

  !! Get parameters from cmd line
  CALL Options()

  !! Allocate data variables
  IF(ALLOCATED(cases)) DEALLOCATE(cases)
  IF(ALLOCATED(controls)) DEALLOCATE(controls)
  IF(ALLOCATED(X)) DEALLOCATE(X)
  ALLOCATE(cases(nca,2*nloci))
  ALLOCATE(controls(nco,2*nloci))
  ALLOCATE(X(nca+nco,nloci))
  
  !! Start simulation
  prev_percent=-1
  text = ""
 	CALL Prnt("Simulation started ...") 
 	IF(nsamples<=1) report_countdown = .FALSE.
  IF(report_countdown) THEN
    CALL Prnt("Progress: ", ADVANCE='NO')
    IF(.NOT.silent) CALL PrintCountdown(text, 0, nsamples, prev_percent)
  ENDIF
  
  !! Run simulation loop
  DO isample = 1, nsamples
    CALL GenerateCasesAndControls(nca, nco, fixedallelefreq, & 
                  allelefreq, RR, RR1, RR2, LD, model, prevalence, &
                  nneutralloci, nLDpairs, X, .FALSE., cases, controls)
    CALL WriteSimulationOutput(cases, controls, nca, nco, pedfile, &
                  mapfile, nomapfile, onlygenes)
    IF(report_countdown) CALL PrintCountdown(text, isample, nsamples, &
                                             prev_percent)
  ENDDO

  !! Get end time
  CALL GetCurrentTime(stoptime_full, stoptime)
  !! Get runtime
  runtime_total = GetTimeDifference(starttime, stoptime)
  
 	CALL Prnt("Simulation finished successfully.") 
 	CALL Prnt("Duration: "//runtime_total) 
 	CALL Prnt("Program terminated.") 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINITIONS OF SUBROUTINES                                              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GetCmdLineArgs()	
!! Process command line options
  IMPLICIT NONE
  INTEGER                       :: nargs, i
  CHARACTER(mfl), ALLOCATABLE   :: args(:)
  CHARACTER(mfl)                :: parg, carg

	!! Set default values
  NUMBER_OF_THREADS = 1
  nsamples = 1
	nneutralloci = 0
	nLDpairs = 0
	ss = 0
	nca = -1
  nco = -1
  allelefreq = -one
	prevalence = min_freq
	RR = 1
	RR1 = 1
	RR2 = 1
	LD = 0
	model = "A"
	fixedallelefreq = .FALSE.
	nomapfile = .TRUE.
	appendfile = .TRUE.
	onlygenes = .TRUE.
	silent = .FALSE.
			
  nargs = COMMAND_ARGUMENT_COUNT()
  IF(nargs==0) THEN
    CALL Prnt("NOTE: No command line arguments specified. "//&
                   "Using default values.")
    allelefreq = MAX(min_freq,MIN(one-min_freq,(five/ss)**(one/four)))
    RETURN
  ENDIF

  IF(ALLOCATED(args)) DEALLOCATE(args)
  ALLOCATE(args(nargs))
  DO i=1,nargs 
    CALL GET_COMMAND_ARGUMENT(i, args(i))
  ENDDO

  carg = ""
  DO i=1,nargs
    parg = carg
    carg = args(i)
    
    !! CHECK FOR PRESENCE ONLY
    SELECT CASE (carg)
      CASE ("--fixedallelefreq", "-fixedallelefreq")
        fixedallelefreq = .TRUE.
      CASE ("--nomapfile", "-nomapfile")
        nomapfile = .TRUE.
      CASE ("--mapfile", "-mapfile")
        nomapfile = .FALSE.
      CASE ("--overwrite", "-overwrite")
        appendfile = .FALSE.
      CASE ("--append", "-append")
        appendfile = .TRUE.
      CASE ("--onlygenes", "-onlygenes")
        onlygenes = .TRUE.
      CASE ("--plinkformat", "-plinkformat")
        onlygenes = .FALSE.
      CASE ("--silent", "-silent")
        silent = .TRUE.
    END SELECT
      
    !! CHECK FOR PRESENCE AND VALUE
    SELECT CASE (parg)
      CASE ("--nthreads", "-nthreads")
        NUMBER_OF_THREADS = MAX(char2int(carg),1)
      CASE ("--nsamples", "-nsamples")
        nsamples = MAX(char2int(carg),1)
      CASE ("--nneutralloci", "-nneutralloci")
        nneutralloci = MAX(char2int(carg),0)
      CASE ("--nLDpairs", "-nLDpairs")
        nLDpairs = MAX(char2int(carg),0)
      CASE ("--ss", "-ss")
        ss = MAX(char2int(carg),0)
      CASE ("--nca", "-nca")
        nca = MAX(char2int(carg),0)
      CASE ("--nco", "-nco")
        nco = MAX(char2int(carg),0)
      
      CASE ("--allelefreq", "-allelefreq")
        allelefreq = MIN(MAX(char2real(carg),min_freq),one-min_freq)
      CASE ("--prevalence", "-prevalence")
        prevalence = MIN(MAX(char2real(carg),min_freq),one)
      CASE ("--RR", "-RR")
        RR = MAX(char2real(carg),zero)
      CASE ("--RR1", "-RR1")
        RR2 = MAX(char2real(carg),zero)
      CASE ("--RR2", "-RR2")
        RR2 = MAX(char2real(carg),zero)
      CASE ("--LD", "-LD")
        LD = char2real(carg)
    END SELECT

  ENDDO
  DEALLOCATE(args)

  !! Assign proper values to nca, nco, ss based on input
	IF(nca>-1 .OR. nco>-1) THEN
	  nca = MAX(nca,0)
    nco = MAX(nco,0)
    ss = nca + nco
  ENDIF

  IF(ss < min_ss) &
    CALL Prnt("Warning: Specified sample size is very low.")
  
  IF(nca == -1) nca = ss/2
  IF(nco == -1) nco = ss-ss/2

  !! Determine the allele frequency from sample size
  IF(allelefreq<zero) &
    allelefreq = MAX(min_freq,MIN(one-min_freq,&
                     (five/MAX(ss/2, min_ss))**(one/four)))
  
  !! Total number of loci
  nloci = nneutralloci+2*(nLDpairs+1)

  RETURN

END SUBROUTINE GetCmdLineArgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Options()
!! Set options and prints announcements
  IMPLICIT NONE

#ifdef _OPENMP
  CALL Prnt("NUMBER OF PARALLEL THREADS: "//&
                 TRIM(int2char(NUMBER_OF_THREADS)))  
#endif

  filename = "EPIC.simulated.model="//model//".LEloci="//&
             TRIM(int2char(nneutralloci))//".LDloci="//&
             TRIM(int2char(nLDpairs))//".n="//&
             TRIM(int2char(ss))//".RR="//TRIM(real2char(RR))//".RR1="//&
             TRIM(real2char(RR1))//".RR2="//TRIM(real2char(RR2))//".LD="//&
             TRIM(real2char(LD))//".AlleleF="//&
             TRIM(real2char(round(allelefreq,2)))//".Prev="//&
             TRIM(real2char(round(prevalence,2)))
             
  pedfile = TRIM(filename)//".ped" 
  mapfile = TRIM(filename)//".map" 

	CALL Prnt("Output ped file: [ "//TRIM(pedfile)//" ]")
  IF(nomapfile) THEN
   	CALL Prnt("Output map file: none")
  ELSE
    CALL Prnt("Output map file: [ "//TRIM(mapfile)//" ]") 
  ENDIF 
  IF(appendfile) THEN 
    CALL Prnt("Note: Previously existing output files with the same"//&
                   " names will be appended.")
  ELSE
    CALL Prnt("Note: Previously existing output files with the same"//&
                   " names will be removed.")
  ENDIF
  
  !! Create the output file if appendfile is .false.
  IF(.NOT.appendfile) THEN ! -e $filename.".ped"){
    CALL OpenFile(UNIT=uped, FILE=pedfile, ACTION='W', STATUS='U', &
                  FORM='F', POSITION="R", ACCESS='S')
   	CLOSE(uped)
    IF(.NOT.nomapfile) THEN
      CALL OpenFile(UNIT=uped, FILE=pedfile, ACTION='W', STATUS='U', &
                    FORM='F', POSITION="R", ACCESS='S')
   	  CLOSE(umap)
 	  ENDIF
  ENDIF
  
  CALL Prnt("Simulation parameters : ")
  CALL Prnt("", empty=.TRUE.)
  CALL Prnt("   Total sample size       : "//TRIM(int2char(ss)),&
                 empty=.TRUE.)
  CALL Prnt("   Number of cases         : "//TRIM(int2char(nca)),&
                 empty=.TRUE.)
  CALL Prnt("   Number of controls      : "//TRIM(int2char(nco)),&
                 empty=.TRUE.)
  CALL Prnt("   Number of samples       : "//TRIM(int2char(nsamples)),&
                 empty=.TRUE.)
  CALL Prnt("   Total number of loci    : "//TRIM(int2char(nloci)),&
                 empty=.TRUE.)
  CALL Prnt("   Pairs of loci in LD     : "//TRIM(int2char(nLDpairs)),&
                 empty=.TRUE.)
  CALL Prnt("   Amount of LD            : "//TRIM(real2char(LD)),&
                 empty=.TRUE.)
  CALL Prnt("   Neutral loci            : "//TRIM(int2char(nneutralloci)),&
                 empty=.TRUE.)
  CALL Prnt("   Minor allele frequency  : "//TRIM(real2char(allelefreq)),&
                 empty=.TRUE.)
  CALL Prnt("   Prevalence              : "//TRIM(real2char(prevalence)),&
                 empty=.TRUE.)
  CALL Prnt("   Interaction model       : "//model, empty=.TRUE.)
  CALL Prnt("   Interaction effect (RR) : "//TRIM(real2char(RR)),&
                 empty=.TRUE.)
  CALL Prnt("   Main effect 1 (RR1)     : "//TRIM(real2char(RR1)),&
                 empty=.TRUE.)
  CALL Prnt("   Main effect 2 (RR2)     : "//TRIM(real2char(RR2)),&
                 empty=.TRUE.)
  CALL Prnt("", empty=.TRUE.)

  RETURN

END SUBROUTINE Options

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END PROGRAM EpiSimulator
