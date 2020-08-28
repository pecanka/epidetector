MODULE EPI_RESULTS

  USE EPI_INIT          !! DEFINITION AND INITIALIZATION OF VARIABLES

  IMPLICIT NONE
  
  CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE AnalyzeResults(filename, noheader)
  IMPLICIT NONE
  INTEGER(ikb), PARAMETER       :: ten_i = 10_ikb, five_i = 5_ikb
  CHARACTER(*), INTENT(INOUT)   :: filename(:)
  LOGICAL, INTENT(IN)           :: noheader
  CHARACTER(res_maxncol*50)     :: line
  !CHARACTER(50)                 :: line(res_maxncol)
  CHARACTER(50)                 :: header(res_maxncol), text
  CHARACTER(50), ALLOCATABLE    :: words(:) 
  CHARACTER(mfl)                :: filename_res
  CHARACTER(3)                  :: advance
  REAL(dpp), ALLOCATABLE        :: output(:,:)
  REAL(dpp), ALLOCATABLE        :: alpha1(:), alpha1b(:), alpha2(:)
  REAL(dpp)                     :: oline(res_maxncol), AS, ASpval, T, Tpval
  INTEGER                       :: iostatus, iostatus2, ifile, i1, i2,& 
                                   nalpha1, nalpha1b, nalpha2, ncorrAS, &
                                   nanalyzed, icor, ilevel, j, ncol, &
                                   prev_ncol, pncols, nlines(SIZE(filename))
  INTEGER(ikb)                  :: i, prev_percent, total_nl, &
                                   file_nl, read_nl  
  INTEGER(ikb), ALLOCATABLE     :: corrAS(:)
  
    !! Allocate the vectors with levels and correction factors
    nalpha1 = 7
    nalpha2 = 7
    IF(MTC(1)>=one) nalpha1 = MAX(INT(LOG10(MTC(1))),1)
    IF(MTC(3)>=one) nalpha2 = MAX(INT(LOG10(MTC(3))),1)

    IF(ALLOCATED(alpha1))  DEALLOCATE(alpha1)
    IF(ALLOCATED(alpha1b)) DEALLOCATE(alpha1b)
    IF(ALLOCATED(alpha2))  DEALLOCATE(alpha2)
    ALLOCATE(alpha1(nalpha1))
    ALLOCATE(alpha1b(nalpha1))
    ALLOCATE(alpha2(nalpha2))
    
    !! Set the levels for pretest and main test
    alpha1 = (/ (level_S1 * ten**(-i), i=0,nalpha1-1) /)
    alpha2 = (/ (level_S2 * ten**(-i), i=0,nalpha2-1) /)

    !! Set the correction factors
    ncorrAS = 2*6
    IF(ALLOCATED(corrAS))  DEALLOCATE(corrAS)
    ALLOCATE(corrAS(ncorrAS))
    corrAS = INT( (/ (ten_i**i, five_i*ten_i**i, i=6,11) /), KIND=ikb)

    !! Number of columns in output before corrAS kicks in
    pncols = 8

    !! Allocate the output array    
    IF(ALLOCATED(output)) DEALLOCATE(output)
    ALLOCATE(output(nalpha1*nalpha2,ncorrAS+pncols))
    
    !! Print analysis mode announcement
    CALL Prnt(TRIM(program_name)//" is running in the output file"//&
                   " analysis mode.")

    !! Print log file and output file information
    CALL Prnt("LOG FILE    : ["//TRIM(log_file)//"]")
    DO ifile=1,SIZE(filename)
      IF(ifile==1) &
        CALL Prnt("INPUT FILES : ["//TRIM(filename(ifile))//"]")
      IF(ifile>1) &
        CALL Prnt0("                ["//TRIM(filename(ifile))//"]")
    ENDDO

    !! If columns not specified
    CALL Prnt("Input file column information :", skip1=1, skip2=1)
    IF(cT==0) THEN
      cT = def_cT
      CALL Prnt0("  * Pretest statistics column not specified."//&
                      " Using default value "//TRIM(i2c(def_cT))//".")
    ELSE
      CALL Prnt0("  * Pretest statistics column   : "//TRIM(i2c(cT)))
    ENDIF
    IF(cTpval==0) THEN
      cTpval = def_cTpval
      CALL Prnt0("  * Pretest p-value column not specified. Using"//&
                      " default value "//TRIM(i2c(def_cTpval))//".")
    ELSE                
      CALL Prnt0("  * Pretest p-value column      : "//TRIM(i2c(cTpval)))
    ENDIF
    IF(cAS==0) THEN
      cAS = def_cAS
      CALL Prnt0("  * Main test statistics column not specified."//&
                      " Using default value "//TRIM(i2c(def_cAS))//".")
    ELSE
      CALL Prnt0("  * Main test statistics column : "//TRIM(i2c(cAS)))
    ENDIF
    IF(cASpval==0) THEN
      cASpval = def_cASpval
      CALL Prnt0("  * Main test p-value column not specified."//&
                      " Using default value "//TRIM(i2c(def_cASpval))//".")
    ELSE
      CALL Prnt0("  * Main test p-value column    : "//TRIM(i2c(cASpval)))
    ENDIF
    
    !! Inform about skipping lines
    IF(result_nskip > 0) &
      CALL Prnt("Skipping first "//TRIM(i2c(result_nskip))//&
                     " of each input file.", skip1=1)
    
    !! Announce start of analysis
    CALL Prnt("Performing analysis of the input file(s) ... ", skip1=1, &
                   advance='NO')

    !! Initialize variables
    header = ""
    nanalyzed = 0
    total_nl = 0
    ilevel = 0

    !! Get number of lines in each file
    CALL GetFileNumLines(filename, nlines, terminate=.FALSE.)

    !! Announce analysis end
    CALL Prnt0("Finished!")

    !! Nulify the file counter
    ifile = 0
    
    !! POINT OF RETURN FOR READING OF THE NEXT FILE
    100 CONTINUE
    !! Next file number
    ifile = ifile + 1
    !! If all files read then go to the end
    IF(ifile>SIZE(filename)) GOTO 200
    !! If file empty, skip it
    IF(nlines(ifile)==0) GOTO 100

    !! Announce analysis    
    CALL Prnt("Analyzing file ["//TRIM(filename(ifile))//"] ... ", &
                   advance='NO')

    !! Open the output file
    CALL OpenFile(UNIT=uout, FILE=filename(ifile), ACTION='R', &
                  STATUS="O", FORM='F', POSITION="R", ACCESS='S')
    
    nanalyzed = nanalyzed + 1

    !! Get the correction factor if not given by user
    IF(MTC(1)<one) THEN
      IF(analyze_files_together)      MTC(1) = SUM(nlines)
      IF(.NOT.analyze_files_together) MTC(1) = nlines(ifile)
    ENDIF

    !! Initialize variables
    advance = 'NO'
    read_nl = -1
    prev_ncol = -1
    prev_percent = -1
    IF(noheader) read_nl = 0
    file_nl = 0
    iostatus = 0
    IF(ifile==1 .OR. .NOT. analyze_files_together) output = zero
    text = ""

    DO WHILE(iostatus==0)
    
      !! Print countdown report
      IF(countdown) &
        CALL PrintCD(text, file_nl, INT(nlines(ifile),ikb), prev_percent, finish=.TRUE.)
                            
      !! Get the pretest levels
      alpha1b = alpha1
      nalpha1b = nalpha1
   
      !! If reading the last line advance=yes so that eof does not give an error
      !IF(file_nl == nlines(ifile) - 1) advance = "YES"
      !! Read the following line
      line = ""
      READ(uout, '(A)', ADVANCE=advance, IOSTAT=iostatus) line
      !! Check for end-of-line error
      IF(iostatus == -2) iostatus = 0
      !! If other error then end-of-file then exit
      IF(iostatus/=0 .AND. ALL(iostatus/=err_eof)) EXIT
      !! If end-of-file error and nothing was read, then exit
      IF(ANY(iostatus==err_eof) .AND. LEN_TRIM(line)==0) EXIT
      !! Increase counter of lines
      file_nl = file_nl + 1
      !! Skip given number of lines
      IF(file_nl<result_nskip) CYCLE
      !! Impose user selected limit on maximum number of lines read
      IF(file_nl>result_maxlines .AND. result_maxlines>0) EXIT
      
      !! Check for comment character
      IF(LEN_TRIM(line)>0) THEN
        IF(line(1:1)==comment) CYCLE
      ENDIF
      
      !! Remove leading spaces
      !line = ADJUSTL(line)

      IF(LEN_TRIM(line)==0) THEN
        CALL Prnt("Warning: Skipping line "//TRIM(i2c(file_nl))//&
                       ". It appears to be empty.", skip1=1, warning=.TRUE.)
        CYCLE
      ENDIF
      
      !! Devide the line into "words" according to tabs
      CALL SplitString(line, tab, words, ncol)

      !! Check for improper number of columns
      IF(ncol/=prev_ncol .AND. prev_ncol/=-1) THEN
        CALL Prnt("ERROR: Line "//TRIM(i2c(file_nl))//&
                       " appears to have "//TRIM(i2c(ncol))//&
                       " instead of "//TRIM(i2c(prev_ncol))//&
                       " columns. Line skipped.", warning=.TRUE.)
        CYCLE
      ENDIF
      
      !! Increase valid line counter
      read_nl = read_nl + 1
      
      !! Check for space in header and oline
      IF(ncol>SIZE(header) .OR. ncol>SIZE(oline)) &
        CALL Prnt("ERROR: Not enough space to store the result line.", &
                       skip1=1, Q=.TRUE.)

      !! Read line into header 
      IF(read_nl==0) THEN
        header(1:ncol) = words(1:ncol)
        CYCLE
      ENDIF

      !! Remember the number of columns
      prev_ncol = ncol
      ncol = 0

      !! Assign NA values
      oline = NAv
      oline(cTpval)  = NAp
      oline(cASpval) = NAp
      !! Read line into oline 
      DO i = 1, prev_ncol
        IF(LEN_TRIM(words(i))==0) CYCLE
        ncol = ncol+1
        IF(TRIM(words(i))==NA) CYCLE
        oline(i) = char2real(words(i), iostatus2)
      ENDDO

      !! Check whether NA in pretest level doesn't just mean "no pretest" and if
      !! so, then modify the stored pretest p-value and the pretest levels
      !! The "no pretest" situation is implied by non-NA values in ASpval 
      IF(oline(cTpval)==NAp .AND. oline(cASpval)/=NAp) THEN
        oline(cTpval) = zero
        alpha1b = one
        nalpha1b = 1
      ENDIF
      
      AS = oline(cAS)
      ASpval = oline(cASpval)
      T = oline(cT)
      Tpval = oline(cTpval)
      
      !! Check for pretest p-values out of range
      IF(Tpval<zero .OR. (Tpval>one .AND. Tpval/=NAp)) THEN
        CALL Prnt("Warning: Nonsensical pretest p-value "//&
                TRIM(r2c(Tpval))//" found on line "//&
                TRIM(i2c(file_nl))//"!", skip1=1, warning=.TRUE.)
        CALL Prnt("Suggestion: Make sure the position of pretest test"//&
                       " p-value column is correctly specified.")
        IF(.NOT.silent) CALL Terminate(ask=.TRUE., announce=.TRUE.)
      ENDIF

      !! Check for main test p-values out of range
      IF(ASpval<zero .OR. (ASpval>one .AND. ASpval/=NAp)) THEN
        CALL Prnt("Warning: Nonsensical main test p-value "//&
                TRIM(r2c(ASpval))//" found on line "//&
                TRIM(i2c(file_nl))//"!", skip1=1, warning=.TRUE.)
        CALL Prnt("Suggestion: Make sure the position of main test"//&
                       " p-value column is correctly specified.")
        IF(.NOT.silent) CALL Terminate(ask=.TRUE., announce=.TRUE.)
      ENDIF
      
      ilevel = 0
      DO i1 = 1, nalpha1b
        DO i2 = 1, nalpha2
          ilevel = ilevel + 1
          !! 1: Total number of tests
          output(ilevel,1) = output(ilevel,1) + 1
          !! 2: Level pretest
          output(ilevel,2) = alpha1b(i1)
          !! 3: Level main test
          output(ilevel,3) = alpha2(i2)

          !! If pretest did not reject, then skip
          IF(.NOT.(Tpval<alpha1b(i1))) CYCLE

          !! 4: Significant pretests
          output(ilevel,4) = output(ilevel,4) + 1
          !! 5: Power of pretest 
          !! (computed later as output(:,4) / output(:,1) )
          !! 6: Significant main tests
          IF(ASpval<alpha2(i2)) output(ilevel,6) = output(ilevel,6) + 1
          !! 7: Power of main test 
          !! (computed later as output(:,7) = output(:,6) / output(:,4) )
          !! 8: Significant main tests after correction for alpha*N
          IF(ASpval*MTC(1)*alpha1b(i1) < alpha2(i2)) & 
            output(ilevel,8) = output(ilevel,8) + 1
          !! 9:? Significant main tests after corrections from corrAS
          DO icor = pncols+1, MIN(SIZE(output,2),SIZE(corrAS)+pncols)
            IF(ASpval*corrAS(icor-pncols)*alpha1b(i1) < alpha2(i2)) & 
              output(ilevel,icor) = output(ilevel,icor) + 1
          ENDDO
        ENDDO
      ENDDO
      
    ENDDO
    
    CLOSE(uout)

    !! Number of lines read over all files
    total_nl = total_nl + read_nl

    !! No lines found to be analyzed
    IF(read_nl <= 0) &
      CALL Prnt("Warning: File ["//TRIM(filename(ifile))//"] appears "//&
                     "to have no data to analyze. No analysis performed on "//&
                     "this file.", skip1=1, warning=.TRUE.)
  
    !! Write output only if files are to be analyzed separately or when last
    !! file has been analyzed
    IF(.NOT.analyze_files_together .OR. ifile==SIZE(filename)) THEN
    
      !! Ad 5: Empirical level/power of pretest
      output(:,5) = output(:,4) / output(:,1)
      !! AS 7: Empirical level/power of main test
      output(:,7) = output(:,6) / MAX(output(:,4),one)
      
      IF(analyze_files_together) THEN
        filename_res = out_file
      ELSE
        filename_res = TRIM(filename(ifile))//"_analysis"//res_ext
      ENDIF

      !! Open file for writing output of analysis
      CALL OpenFile(UNIT=ures, FILE=filename_res, ACTION='W', &
                    STATUS="U", FORM='F', POSITION="R", ACCESS='S', &
                    overwrite=overwrite_files, chck_exist=check_fexist)
  
      WRITE(ures,'(15A)', ADVANCE='NO') "Correction", tab, "NTests1",&
          tab, "Level1", tab, "Level2", tab, "NSignif1", tab, &
          "EmpLevel1", tab, "NSignif2", tab, "EmpLevel2"
      WRITE(ures,'(2A)',ADVANCE='NO') tab, "'C=Level1*Correction'"
      DO icor=1,ncorrAS
        WRITE(ures,'(2A)',ADVANCE='NO') & 
          tab, "'C="//TRIM(i2c(corrAS(icor)))//"*Level1'"
      ENDDO
      WRITE(ures,'(A)') ""
  
      !! Write the number of rejections
      DO i=1,ilevel
        WRITE(ures,'(A)', ADVANCE='NO') TRIM(r2c(MTC(1)))
        DO j=1,SIZE(output,2)
          WRITE(ures,'(2A)',ADVANCE='NO') &
            tab, TRIM(r2c(output(i,j)))
        ENDDO
        WRITE(ures,'(A)') ""
      ENDDO
  
      CLOSE(ures)
    
    ENDIF

    !! Announce end of current file analysis (only the log file if on screen
    !! "finished" will be printed by PrintCD)
    CALL Prnt0("Finished!", screen=.NOT.countdown)

    !! Read the next file
    GOTO 100

    !! Skip here if all files were read
    200 IF(ALLOCATED(output)) DEALLOCATE(output)

    !! No lines found to be analyzed
    IF(nanalyzed <= 0 .OR. total_nl<=0) THEN
      CALL Prnt("ERROR: NO file was analyzed and NO output was produced!.")
    ELSE
      CALL Prnt("Analysis of all files finished!")
    ENDIF
    
    RETURN
  
END SUBROUTINE AnalyzeResults

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_RESULTS
