MODULE EPI_UTILS

  USE EPI_PARAMS
  !USE LOGISTICREG

#ifdef _OPENMP
  USE OMP_LIB
  !USE PAR_ZIG_MOD 
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PARAMETERS AND INTERFACES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CHARACTER(mfl), SAVE         :: log_file = ""
  CHARACTER(24), SAVE          :: timestamp = ""
  LOGICAL, SAVE                :: silent, &
                                  no_log, &
                                  overwrite_files, &
                                  countdown, &
                                  log_advance_always, & 
                                  screen_advance_always, &
                                  no_openfile_iostat, &
                                  log_opened, pause_run, &
                                  add_commas = .TRUE., &
                                  add_spaces = .FALSE.
                                  
  INTEGER, SAVE                :: NUMBER_OF_THREADS, &
                                  NUMBER_OF_PROCESSORS, &
                                  nwarnings, &
                                  nerrors, &
                                  nwarnings_limit, &
                                  nerrors_limit, &
                                  countdown_jump
                                  
  CHARACTER, SAVE              :: comment = "#", &
                                  CharNAsts, &
                                  CharNAsex, &
                                  CharNAgen, &
                                  CharCa, &
                                  CharCo, &
                                  CharMa, &
                                  CharFe
                                 
  CHARACTER(mltl), ALLOCATABLE :: warnings(:), errors(:) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERFACES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE TrimStr
    MODULE PROCEDURE TrimStr1, TrimStrN
  END INTERFACE

  INTERFACE upcase
    MODULE PROCEDURE upcase1, upcaseN
  END INTERFACE

  INTERFACE upcasef
    MODULE PROCEDURE upcase1f, upcaseNf
  END INTERFACE

  INTERFACE lowcase
    MODULE PROCEDURE lowcase1, lowcaseN
  END INTERFACE

  INTERFACE lowcasef
    MODULE PROCEDURE lowcase1f, lowcaseNf
  END INTERFACE
  
  INTERFACE DelAllSubstr
    MODULE PROCEDURE DelAllSubstr1, DelAllSubstrN
  END INTERFACE
  
  INTERFACE DeleteAllSubstr
    MODULE PROCEDURE DelAllSubstr1, DelAllSubstrN
  END INTERFACE
  
  INTERFACE ReplaceAllSubstr
    MODULE PROCEDURE ReplaceAllSubstr1, ReplaceAllSubstrN
  END INTERFACE

  INTERFACE SplitString
    MODULE PROCEDURE SplitString1, SplitStringN
  END INTERFACE

  INTERFACE RemoveSubstrEnd
    MODULE PROCEDURE RemoveSubstrEnd1, RemoveSubstrEndN
  END INTERFACE

  INTERFACE DeleteFile
    MODULE PROCEDURE DeleteFile1, DeleteFileN
  END INTERFACE 

  !! Change the size of allocatable arrays

  INTERFACE Dealloc
    MODULE PROCEDURE Dealloc1, Dealloc2
  END INTERFACE

  INTERFACE ResizeVar
    MODULE PROCEDURE Resize1Real, Resize1Char
    MODULE PROCEDURE Resize1Int_ikn, Resize1Int_iks, Resize1Int_ikb !, Resize1IntL
    MODULE PROCEDURE Resize1Logi, Resize1TEMP, Resize1TEMPmin

    MODULE PROCEDURE Resize2Real, Resize2Char
    MODULE PROCEDURE Resize2Int_ikn, Resize2Int_iks, Resize2Int_ikb !, Resize2IntL

    MODULE PROCEDURE Resize3Int_ikn, Resize3Int_iks, Resize3Int_ikb !, Resize3IntL
    MODULE PROCEDURE Resize3Char, Resize3Real
  END INTERFACE

  !! Convert integer to string with possible formating 
  INTERFACE int2char
    MODULE PROCEDURE si2c, ni2c, li2c
    MODULE PROCEDURE si2cN, ni2cN, li2cN
  END INTERFACE

  !! Allias for int2char
  INTERFACE i2c
    MODULE PROCEDURE si2c, ni2c, li2c
    MODULE PROCEDURE si2cN, ni2cN, li2cN
  END INTERFACE
  
  !! Convert integer to string without any formating
  INTERFACE i2cp
    MODULE PROCEDURE si2cp, ni2cp, bi2cp!, li2cp
  END INTERFACE
  
  !! Convert string to integer (normal integer type)
  INTERFACE c2i
    MODULE PROCEDURE nchar2int, nchar2intN
  END INTERFACE

  !! Allias for char2int (normal integer type)
  INTERFACE c2i8
    MODULE PROCEDURE nchar2int8, nchar2int8N
  END INTERFACE

  !! Convert string to integer (long integer type)
  INTERFACE lc2i
    MODULE PROCEDURE lchar2int, lchar2intN
  END INTERFACE

  !! Convert real to string
  INTERFACE r2c
    MODULE PROCEDURE real2char_sp, real2char_dp
  END INTERFACE

  !! Convert string to real
  INTERFACE c2r
    MODULE PROCEDURE char2real1, char2realN
  END INTERFACE

  !! Get the number of lines in a file
  INTERFACE GetFileNumLines
    MODULE PROCEDURE GetFileNumLines1, GetFileNumLinesN
  END INTERFACE

  !! Erase text specified by giving the text or the length of the text
  INTERFACE EraseText
    MODULE PROCEDURE EraseTextL, EraseTextT
  END INTERFACE

  !! Procedure that replace '=' when there are many assignments of the same 
  !! value
  INTERFACE SetVal
    MODULE PROCEDURE SetValReal1, SetValRealN, SetValRealNM
    MODULE PROCEDURE SetValInt1, SetValIntN
    MODULE PROCEDURE SetValBool1, SetValBoolN
  END INTERFACE
  
  INTERFACE PrintCD
    MODULE PROCEDURE PrintCD1, PrintCD2
  END INTERFACE

 CONTAINS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! SUBROUTINES AND FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PURE FUNCTION IND(bool1, bool2) RESULT(num)
!! Function which converts a logical input into numerical representation. If at
!! least one of bool1 and bool2 is true it returns 1 otherwise it returns 0.
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: bool1, bool2
  INTEGER             :: num
  
  num = 0
  
  IF(bool1 .OR. bool2) num = 1

  RETURN
  
END FUNCTION IND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!**************************************************************************!!  
!!**************          FILE RELATED ROUTINES         ********************!! 
!!**************************************************************************!!  

RECURSIVE SUBROUTINE OpenFile(UNIT, FILE, ACTION, FORM, STATUS, POSITION, &
      ACCESS, RECL, IOSTAT, chck_exist, chck_miss, overwrite, ask, &
      announce, islog)
  IMPLICIT NONE
  INTEGER, INTENT(IN)              :: UNIT
  CHARACTER(*), INTENT(INOUT)      :: FILE
  CHARACTER, INTENT(IN)            :: ACTION, FORM
  CHARACTER, INTENT(IN), OPTIONAL  :: STATUS, POSITION, ACCESS
  INTEGER, INTENT(IN), OPTIONAL    :: RECL
  INTEGER, INTENT(OUT), OPTIONAL   :: IOSTAT
  LOGICAL, INTENT(IN), OPTIONAL    :: chck_exist, chck_miss, ask, &
                                      announce, islog
  LOGICAL, INTENT(INOUT), OPTIONAL :: overwrite
  CHARACTER(LEN(FILE))             :: filename
  CHARACTER(20)                    :: action1, form1, status1, position1, &
                                      access1
  CHARACTER                        :: response1
  INTEGER                          :: ios1, iname
  LOGICAL                          :: file_exists, chck_exist1, ask1, &
                                      announce1, chck_miss1, islog1
  
  !! Assign default values
  action1 = 'READWRITE'
  status1 = 'UNKNOWN'
  access1 = 'SEQUENTIAL'
  position1 = 'REWIND'
  form1 = ''

  !! Read in user input values
  IF(upcasef(ACTION)=="R")     action1 = 'READ'
  IF(upcasef(ACTION)=="W")     action1 = 'WRITE'
  IF(upcasef(FORM)=="F")       form1   = 'FORMATTED'
  IF(upcasef(FORM)=="U")       form1   = 'UNFORMATTED'
  IF(PRESENT(STATUS)) THEN
    IF(upcasef(STATUS)=="N")   status1 = 'NEW'
    IF(upcasef(STATUS)=="O")   status1 = 'OLD'
    IF(upcasef(STATUS)=="S")   status1 = 'SCRATCH'
    IF(upcasef(STATUS)=="U")   status1 = 'UNKNOWN'
    IF(upcasef(STATUS)=="R")   status1 = 'REPLACE'
  ENDIF
  IF(PRESENT(POSITION)) THEN
    IF(upcasef(POSITION)=="R") position1 = 'REWIND'
    IF(upcasef(POSITION)=="A") position1 = 'APPEND'
  ENDIF
  IF(PRESENT(ACCESS)) THEN
    IF(upcasef(ACCESS)=="S")   access1 = 'SEQUENTIAL'
    IF(upcasef(ACCESS)=="R")   access1 = 'STREAM'
    IF(upcasef(ACCESS)=="D")   access1 = 'DIRECT'
  ENDIF

  filename = FILE
  iname = 1
  ios1 = 0
  
  !! Assign intial values for variables
  ask1 = .FALSE.                ! asking for new filename
  announce1 = .TRUE.            ! announcing premature termination
  chck_exist1 = .FALSE.    ! check file existence before opening it
  chck_miss1 = .FALSE.      ! check whether file missing before opening it
  islog1 = .FALSE.              ! file to be opened is a log file
  
  !! Change the values based on input
  IF(PRESENT(chck_exist)) chck_exist1 = chck_exist
  IF(PRESENT(chck_miss))  chck_miss1 = chck_miss
  IF(PRESENT(islog))      islog1 = islog
  IF(PRESENT(ask))        ask1 = ask
  IF(PRESENT(announce))   announce1 = announce

  !! If proper conditions are met then skip existence check
  IF(action1=="READ") chck_exist1 = .FALSE.  
  IF(PRESENT(overwrite)) THEN
    IF(overwrite) chck_exist1 = .FALSE.
  ENDIF

  10 CONTINUE

  !! Save new file name into variable FILE
  FILE = filename

  !! If existence or missingness check not wanted, go directly to opening it
  IF(.NOT.chck_exist1 .AND. .NOT.chck_miss1) GOTO 20

  !! Check for existence of file
  INQUIRE(FILE=filename, EXIST=file_exists)

  !! Problem with file: FILE MISSING
  IF(chck_miss1 .AND. .NOT.file_exists) THEN
    
    !! Announce missing file
    CALL PrntE("Cannot find file ["//TRIM(filename)//"].", log=.NOT.islog1)
    
    !! If silent, then don't ask and try to resolve the name conflict
    IF(silent) CALL Terminate(log=.NOT.islog1) 

    !! Ask for new filename is wanted, otherwise ask for termination
    IF(ask1) THEN
      CALL AskForNewFilename(filename)
    ELSE
      CALL AskUser(response1, retry=.TRUE.)
      IF(response1=='N') CALL Terminate(ask=.TRUE., log=.NOT.islog1)
    ENDIF
    !! If still running go to the beggining
    GOTO 10
  ENDIF

  !! Problem with file: FILE ALREADY EXISTS
  IF(chck_exist1 .AND. file_exists) THEN
  
    !! If silent, then don't ask and try resolve the name conflict
    IF(silent) THEN
      CALL ChangeFilename(filename,"", "_"//TRIM(i2cp(iname)), findext=.TRUE.)
      iname = iname + 1
      GOTO 10
    ENDIF
    
    !! Otherwise, ask whether to rewrite
    CALL Prnt("File ["//TRIM(filename)//"] already exists.", skip1=1)

    IF(PRESENT(overwrite)) THEN
      CALL AskUser(response1, overwriteall=.TRUE.)
      IF(response1=="A") overwrite = .TRUE.
    ELSE
      CALL AskUser(response1, overwrite=.TRUE.)
    ENDIF
    
    IF(response1=="N") THEN
      !! Ask for new filename is wanted, otherwise ask for termination
      IF(ask1) THEN
        CALL AskForNewFilename(filename)
      ELSE
        CALL Terminate(announce=.NOT.islog1, ask=.TRUE.)
      ENDIF
      !! If still running go to the beggining
      GOTO 10
    ENDIF
    
  ENDIF

  !! Try opening the file
  20 CONTINUE
  
  !! Use DIRECT ACCESS to open the file
  IF(access1=='DIRECT') THEN
    IF(.NOT.PRESENT(RECL)) & 
      CALL PrntE("Direct file access requires RECL!", Q=.TRUE.)
    IF(RECL<=0) & 
      CALL PrntE("Record length (RECL) must be positive!", Q=.TRUE.)
    
    !! If RECL ok, open the file
    IF(no_openfile_iostat) THEN
      OPEN(UNIT=unit, FILE=filename, STATUS=status1, FORM=form1, &
           ACCESS=access1, RECL=RECL)
    ELSE
      OPEN(UNIT=unit, FILE=filename, STATUS=status1, FORM=form1, &
           ACCESS=access1, RECL=RECL, IOSTAT=ios1)
    ENDIF
  !! Use SEQUENTIAL ACCESS to open the file
  ELSE
    IF(no_openfile_iostat) THEN
      OPEN(UNIT=unit, FILE=filename, ACTION=action1, STATUS=status1, &
           FORM=form1, POSITION=position1)
    ELSE
      OPEN(UNIT=unit, FILE=filename, ACTION=action1, STATUS=status1, &
           FORM=form1, POSITION=position1, IOSTAT=ios1)
    ENDIF
  ENDIF
           
  !! Announce problems
  IF(ios1/=0) THEN
    IF(announce1) THEN
      !! Announce error
      CALL PrntE("Cannot open file ["//TRIM(filename)//"] (IOSTAT "//&
                 TRIM(i2cp(ios1))//").", log=.NOT.islog1)
      !! Try to give a hint on what type of error
      IF(ANY(ios1==err_unclosedfile)) &
        CALL Prnt("NOTE: The file appears to be already opened and"//&
                  " blocked by a process. It is possible the file is"//&
                  " blocked by this process, in which case there is"//&
                  " probably a bug in the code. Please report this to"//&
                  " the author with a detailed description of when this"//&
                  " error occurred. Thank you.", log=.NOT.islog1, lead=2)
      IF(ANY(ios1==err_openfile)) &
        CALL Prnt("NOTE: The file appears to be already opened and"//&
                  " blocked by a process. Unlock it and try again.", &
                  log=.NOT.islog1, lead=2)
      IF(ANY(ios1==err_nonexistfile)) THEN
        CALL Prnt("REASON: File does not exist.", log=.NOT.islog1, lead=2)
        IF(LEN_TRIM(filename)>system_filename_limit) &
          CALL Prnt("Hint: The 'non-existent file' error can occur"//&
                    " also if the file name (including path) is excessively"//&
                    " long. If you think this is the reason for the current"//&
                    " error try changing the file name or moving the file into"//&
                    " a directory with a shorter path.", log=.NOT.islog1, lead=2)
      ENDIF
    ENDIF
  ENDIF
    
  !! Return the ios
  IF(PRESENT(IOSTAT)) THEN
    IOSTAT = ios1
    RETURN
  ENDIF
  
  !! If still here, either terminate or ask user for confirmation
  IF(ios1/=0) THEN
    !! If silent, just terminate
    IF(silent) CALL Terminate(announce=.NOT.islog1, ask=.FALSE.)

    !! Otherwise ask user to retry
    CALL AskUser(response1, retry=.TRUE.)
    IF(response1=="Y") GOTO 10

    !! Ask for a new file name or terminate
    IF(ask1) THEN
      CALL AskForNewFilename(filename)
    ELSE
      CALL Terminate(ask=.TRUE., log=.NOT.islog1)
    ENDIF
    
    !! If not terminated by now, try again
    GOTO 10
  ENDIF
  
  RETURN

END SUBROUTINE OpenFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ExcessFileSize(fn, funit, filesize, max_size, nl, ext, &
                          check_size, form, access, rec_size)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: fn(:)
  INTEGER, INTENT(IN)                      :: funit
  INTEGER(ikb), INTENT(OUT)                :: filesize
  INTEGER(ikb), INTENT(IN)                 :: max_size
  INTEGER(ikb), INTENT(INOUT)              :: nl
  CHARACTER(*), INTENT(IN)                 :: ext
  LOGICAL, INTENT(IN), OPTIONAL            :: check_size
  CHARACTER, INTENT(IN), OPTIONAL          :: form, access
  INTEGER, INTENT(IN), OPTIONAL            :: rec_size
  INTEGER                                  :: nf, ios
  CHARACTER(mfl)                           :: fbase
  CHARACTER                                :: form1, access1
  CHARACTER(10)                            :: apx
  LOGICAL                                  :: is_opened, check_size1
  
  !! Value in form should be one of 'U' or 'F'
  form1 = 'U'
  IF(PRESENT(form)) form1 = form

  !! Value in access should be one of 'D' or 'S'
  access1 = 'D'
  IF(PRESENT(access)) access1 = access
  
  !! Select whether file size check is performed
  check_size1 = .TRUE.
  IF(PRESENT(check_size)) check_size1 = check_size
  
  !! Get the current number of files
  nf = SIZE(fn)

  !! Unless user disabled it, check for filesize exceeding the maximum size
  IF(check_size1) THEN
    
    IF(PRESENT(rec_size)) THEN
      
      filesize = nl * rec_size
      IF(filesize < max_size) RETURN
    
    ELSE
    
      !! Close the file first so that its size can be obtained
      INQUIRE(FILE=fn(nf), OPENED=is_opened)      
      IF(is_opened) CLOSE(funit)
    
      !! Obtain the size
      filesize = GetFileSize(fn(nf))
    
      !! If small enough, then reopen it and return
      IF(filesize < max_size) THEN
        IF(PRESENT(rec_size)) THEN
          CALL OpenFile(UNIT=funit, FILE=fn(nf), ACTION='W', STATUS="O", &
                        FORM=form1, ACCESS=access1, POSITION='A', RECL=rec_size)
        ELSE
          CALL OpenFile(UNIT=funit, FILE=fn(nf), ACTION='W', STATUS="O", &
                        FORM=form1, ACCESS=access1, POSITION='A')
        ENDIF
        RETURN
      ENDIF
    
    ENDIF
  
  ENDIF
  
  !! If the file is too big (approaching maximum filesize), open new one
  IF(filesize >= max_size) THEN
    
    !! Close the old file
    INQUIRE(FILE=fn(nf), OPENED=is_opened)      
    IF(is_opened) CLOSE(funit)
    
    !! Rename the first file to match the other files
    !IF(.FALSE.) THEN
    IF(nf==1) THEN
      fbase = fn(1)
      CALL RemoveSubstrEnd(fn(1), ext)
      fn(1) = TRIM(fn(1))//"_001"//TRIM(ext)
      CALL RenameFile(fbase, fn(1), ios)
      IF(ios/=0) & 
        CALL PrntE("Error occurred during renaming of file ["//TRIM(fbase)//&
                     "] to ["//TRIM(fn(1))//"] (IOSTAT "//i2cp(ios)//").", Q=.TRUE.)
    ENDIF
    !ENDIF 
    
    !! Reset the current file line counter 
    nl = 1
    
    !! Increase the number of output files
    nf = nf + 1
    CALL ResizeVar(fn, nf)
    fbase = fn(1)

    !! Drop the extention
    CALL RemoveSubstrEnd(fbase, ext)

    !! Drop the '_00?' counter
    fbase = fbase(1:LEN_TRIM(fbase)-4)    
    fn(nf) = fbase
    
    !! Change the fn for the new file
    CALL RemoveSubstrEnd(fn(nf), ext)
    WRITE(apx, '(A,I0.3)') "_", nf
    fn(nf) = TRIM(fn(nf))//TRIM(apx)//TRIM(ext)
    
    !! Open a new output file
    IF(PRESENT(rec_size)) THEN
      CALL OpenFile(UNIT=funit, FILE=fn(nf), ACTION='W', STATUS="U", &
                    FORM=form1, ACCESS=access, POSITION="R", RECL=rec_size)
    ELSE
      CALL OpenFile(UNIT=funit, FILE=fn(nf), ACTION='W', STATUS="U", &
                    FORM=form1, ACCESS=access, POSITION="R")
    ENDIF
  
  ENDIF

  RETURN
  
END SUBROUTINE ExcessFileSize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DeleteFile1(filename, confirm, ios)
!! Deletes a file by opening it and closing it with status delete
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: filename
  LOGICAL, INTENT(IN), OPTIONAL  :: confirm
  INTEGER, INTENT(OUT), OPTIONAL :: ios
  INTEGER                        :: ios1
  LOGICAL                        :: file_exists
  CHARACTER                      :: response

  !! Check if file exists
  INQUIRE(FILE=filename, EXIST=file_exists)
  
  !! If it doesn't exist just return
  IF(.NOT.file_exists) THEN
    IF(PRESENT(ios)) ios = -9
    RETURN
  ENDIF

  !! Ask for confirmation of delete
  IF(PRESENT(confirm)) THEN
    IF(file_exists .AND. confirm) THEN
      CALL AskUser(response, confirm_delete=.TRUE.)
      IF(response=='N') THEN
        CALL Prnt("File ["//TRIM(filename)//"] was not deleted.")
        RETURN
      ENDIF
    ENDIF
  ENDIF
  
  !! Delete the file
  CALL OpenFile(UNIT=udel, FILE=filename, ACTION='W', STATUS="O", &
                FORM='U', POSITION="R", ACCESS='S', IOSTAT=ios1)
                
  IF(ios1==0) CLOSE(udel, STATUS='DELETE', IOSTAT=ios1)
  
  IF(PRESENT(ios)) ios = ios1
  
  RETURN

END SUBROUTINE DeleteFile1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DeleteFileN(filename, confirm, ios)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: filename(:)
  LOGICAL, INTENT(IN), OPTIONAL  :: confirm
  INTEGER, INTENT(OUT), OPTIONAL :: ios
  INTEGER                        :: i
  
  !! Loop over given file names and delete them
  DO i=1,SIZE(filename)
    CALL DeleteFile1(filename(i), confirm, ios)
  ENDDO
  
  RETURN

END SUBROUTINE DeleteFileN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RenameFile(old, new, ios, confirm, announce)
!! Renames a file (currently works only on GNU Fortran)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: old, new
  INTEGER, INTENT(OUT), OPTIONAL :: ios
  LOGICAL, INTENT(IN), OPTIONAL  :: confirm, announce
  LOGICAL                        :: exists, announce1
  CHARACTER                      :: response
  
  announce1 = .FALSE.
  IF(PRESENT(announce)) announce1 = announce

  IF(announce1) CALL Prnt("Renaming file ["//TRIM(old)//"] to ["//TRIM(new)//"] ...")

  IF(PRESENT(ios)) ios = 0 
  
  !! Check if the old file exists
  INQUIRE(FILE=old, EXIST=exists)
  IF(.NOT.exists) CALL PrntE("File ["//TRIM(old)//"] does not exist.", Q=.TRUE.)

  !! Ask for confirmation of rename if new file already exists
  INQUIRE(FILE=new, EXIST=exists)
  IF(PRESENT(confirm)) THEN
    IF(exists .AND. confirm) THEN
      CALL Prnt("File ["//TRIM(new)//"] already exists.")
      CALL AskUser(response, overwrite=.TRUE.)
      IF(response/='Y') THEN
        CALL Prnt("File ["//TRIM(new)//"] was NOT overwritten.")
        IF(PRESENT(ios)) ios = -8
        RETURN
      ENDIF
    ENDIF
  ENDIF
  
  !!! If the "new" file already exists, delete it before renaming the old one to it
  !IF(exists) CALL DeleteFile(new)

  !! Rename the file
  CALL RENAME(old, new)
  
  RETURN
  
END SUBROUTINE RenameFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CheckFileExistence(fn, fnlist, filepath, filelistpath, msg, &
                              delim, exists, quit, skip1, skip2)
!! Check existence of files given in fn (scalar) and fnlist (array)
  IMPLICIT NONE
  CHARACTER(*), OPTIONAL, INTENT(IN) :: fn, fnlist(:), filepath, &
                                        filelistpath, msg, delim
  LOGICAL, INTENT(OUT), OPTIONAL     :: exists
  LOGICAL, INTENT(IN), OPTIONAL      :: quit
  INTEGER, INTENT(IN), OPTIONAL      :: skip1, skip2
  CHARACTER                          :: slash
  CHARACTER(mfl)                     :: fname
  CHARACTER(mmtl)                    :: msg1
  LOGICAL                            :: file_exists, quit1
  INTEGER                            :: i, length, skip11, skip22
  
  msg1 = ""
  quit1 = .TRUE.
  IF(PRESENT(quit)) quit1 = quit
  IF(PRESENT(msg)) msg1 = " ("//TRIM(msg)//")"
  
  !! Get system delimiter
  IF(PRESENT(delim)) THEN
    IF(LEN_TRIM(delim)==0) THEN
      CALL get_environment_variable("DELIMITER", slash)
    ELSE 
      slash = delim
    ENDIF
  ELSE
    CALL get_environment_variable("DELIMITER",slash)
  ENDIF
  
  !! Set how much skips will be made when printing
  skip11 = 1
  skip22 = 1
  IF(PRESENT(skip1)) skip11 = MAX(0, skip1)
  IF(PRESENT(skip2)) skip22 = MAX(0, skip2)

  !! Check for existence of file with the name filename
  IF(PRESENT(fn)) THEN
    IF(PRESENT(filepath)) THEN
      length = LEN_TRIM(filepath)
      IF(filepath(length:length)==slash) slash=""
      fname = TRIM(filepath)//TRIM(slash)//TRIM(fn)
    ELSE
      fname = fn
    ENDIF
    INQUIRE(FILE=fname, EXIST=file_exists)
    IF(PRESENT(exists)) exists=file_exists
    IF(.NOT.file_exists) &
      CALL Prnt("File ["//TRIM(fname)//"]"//TRIM(msg1)//" does not exist.", &
                     skip1=skip11, skip2=skip22, Q=quit1)
  ENDIF
  
  !! Check for existence of files listed in file with name filenamelist
  IF(PRESENT(fnlist)) THEN
    IF(SIZE(fnlist)>0) THEN
      DO i=1,SIZE(fnlist)
        IF(PRESENT(filelistpath)) THEN
          length = LEN_TRIM(filelistpath)
          IF(filelistpath(length:length)==slash) slash=""
          fname = TRIM(filelistpath)//TRIM(slash)//TRIM(fnlist(i))
        ELSE
          fname = fnlist(i)
        ENDIF  
        INQUIRE(FILE=fname, EXIST=file_exists)
        IF(.NOT.file_exists) &
          CALL Prnt("File ["//TRIM(fname)//"] does not exist.", &
                         skip1=skip11, skip2=skip22, Q=quit1)
      ENDDO
    ENDIF
  ENDIF
  
  RETURN

END SUBROUTINE CheckFileExistence

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ChangeFilename(filename, extension, add, findext)
!! Changes filename by adding add at the end of filename before extension 
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: filename
  CHARACTER(*), INTENT(IN)      :: extension
  CHARACTER(*), INTENT(IN)      :: add
  LOGICAL, INTENT(IN), OPTIONAL :: findext
  CHARACTER(mextl)              :: extension1
  INTEGER                       :: length, pos

  extension1 = extension
  length = LEN_TRIM(filename)
  
  !! Find out the extension
  IF(LEN_TRIM(extension1)==0) THEN
    pos = 0
    IF(PRESENT(findext)) THEN
      IF(findext) THEN
        extension1 = ""
        pos = FindLastSubstr(filename,".")
        IF(pos>0) extension1 = filename(pos:)
      ENDIF
    ENDIF
  ELSE  
    pos = FindLastSubstr(filename,extension1)
  ENDIF
  
  !! Change the file name
  IF(pos>0) THEN
    filename = filename(1:MAX(1,pos-1))//TRIM(add)//TRIM(extension1)
  ELSE
    filename = TRIM(filename)//TRIM(add)
  ENDIF
  
  RETURN
    
END SUBROUTINE ChangeFilename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION GetFileSize(filename) RESULT(fs)
!! Returns size of a given file
!! WARNING: If compiled with Intel Fortran (version at least 8) the size of one
!! record is 4 bytes, unless '-assume byterecl' added to ifort compilation line  
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: filename
  INTEGER(ikb)                :: fs, maxlen, pos1, pos2, i
  INTEGER                     :: ios
  LOGICAL                     :: file_exists
  !CHARACTER                   :: one_byte
  INTEGER(1)                  :: one_byte
  
  !! Set fs to a value that indicates non-existent file
  fs = -1
  
  !! Check for file existence
  INQUIRE(FILE=filename, EXIST=file_exists)
  IF(.NOT.file_exists) RETURN

  !!! Try inquiring the size intrinsically (INQUIRE works in G95)
  !INQUIRE(FILE=filename, SIZE=fs)

  !! If inquire returns other value than -1, exit
  !! Otherwise, try getting the size by "splitting interval"
  IF(fs/=-1) RETURN

  !! Maximum filesize is given by the type of maxlen
  maxlen = HUGE(maxlen)
  
  !! Open with record length 1:
  CALL OpenFile(UNIT=ufre, FILE=filename, ACTION='R', STATUS='O', &
                FORM='U', POSITION="R", ACCESS='D', RECL=1, ask=.FALSE.)

  !! Binary search for file end.
  !! If we read fine, we're below file end,
  !! if we get an error, we're above file end.
  pos1 = 0
  pos2 = maxlen

  10 CONTINUE 
  
  !! Found the position of the last byte of the file, save it and exit
  IF(pos2 - pos1 <= 1) THEN
    fs = pos1
    CLOSE(ufre)
    RETURN
  ENDIF

  !! Last byte was not yet found, try again
  i = (pos1 + pos2)/2
  READ(UNIT=ufre, REC=i, IOSTAT=ios) one_byte

  !! We read without error, which means file size >= I
  IF(ios == 0) THEN
  
    pos1 = i
    GOTO 10
  
  !! End-of-file or Past-End-of-file errors encountered, thus file size < i
  ELSEIF(.TRUE. .OR. ANY(ios==(/err_eof, err_peof/))) THEN
  
    pos2 = i
    GOTO 10
  
  !! Encounted other error than EOF:
  ELSE    
  
    fs = -1
    CLOSE(ufre)
    CALL PrntE("The size of file ["//TRIM(filename)//"] could not be"//&
                 " determined. Last READ returned IOSTAT value "//&
                 TRIM(i2cp(ios))//".", Q=.TRUE.)
    RETURN
  
  ENDIF

END FUNCTION GetFileSize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetFileNumLines1(filename, nlines, quit)
!! Changes filename by adding add at the end of filename before extension 
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: filename
  INTEGER, INTENT(OUT)          :: nlines
  LOGICAL, INTENT(IN), OPTIONAL :: quit
  INTEGER                       :: ios
  LOGICAL                       :: quit1, file_exists

  quit1 = .TRUE.
  IF(PRESENT(quit)) quit1 = quit

  nlines = 0
  
  !! Check file existence
  CALL CheckFileExistence(filename, exists=file_exists, quit=quit1)
  
  !! If the program didn't quit but the file does not exist, return
  IF(.NOT.file_exists) RETURN

  !! Open the file
  CALL OpenFile(UNIT=ufre, FILE=filename, ACTION='R', STATUS="O", &
                FORM='F', POSITION="R", ACCESS='S')
  
  !! Read lines until ios not zero
  DO
    READ(ufre, *, IOSTAT=ios)
    IF(ios/=0) EXIT
    nlines = nlines + 1
  ENDDO
  
  CLOSE(ufre)

  RETURN
    
END SUBROUTINE GetFileNumLines1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetFileNumLinesN(filename, nlines, terminate)
  CHARACTER(*), INTENT(INOUT)   :: filename(:)
  INTEGER, INTENT(OUT)          :: nlines(:)
  LOGICAL, INTENT(IN), OPTIONAL :: terminate
  INTEGER                       :: i
  
  IF(SIZE(filename)/=SIZE(nlines) .OR. SIZE(filename)<=0) &
    CALL PrntE("Array sizes do not match (GetFileNumLinesN)", &
                   Q=.TRUE., skip1=1)
                     
  DO i=1,SIZE(filename)
    CALL GetFileNumLines1(filename(i), nlines(i), terminate)
  ENDDO
  
  RETURN

END SUBROUTINE GetFileNumLinesN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadTable(filename, CONTENT)
  INTEGER, PARAMETER                       :: funit = 62
  CHARACTER, PARAMETER                     :: LF = CHAR(10), CR = CHAR(13)
  CHARACTER(*), INTENT(INOUT)              :: filename
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: CONTENT(:,:) 
  CHARACTER, ALLOCATABLE                   :: TEXT(:) 
  INTEGER, ALLOCATABLE                     :: EOL(:)
  INTEGER(ikb)                             :: filesize
  INTEGER                                  :: ios, n, i, j, k, l, nr, nc
  LOGICAL                                  :: is_win
  
  !! Allocate the output array
  CALL ResizeVar(CONTENT, 1, 1)

  !! Get the total file size
  filesize = GetFileSize(filename)

  !! Stop if file too big for fortran
  IF(filesize > HUGE(INT(1))) CALL PrntE("File is too big!")
  
  !! If the file is empty, just leave
  IF(filesize==0) RETURN

  !! Otherwise read the contents    
  n = INT(filesize) 
  CALL ResizeVar(TEXT, n)
  
  !! Read the whole file at once
  OPEN(UNIT=funit, FILE=filename, STATUS="OLD", ACCESS="STREAM")
  !CALL OpenFile(UNIT=funit, FILE=filename, ACTION='R', STATUS="O", &
  !              FORM='F', POSITION="R", ACCESS='R')
  READ(funit, POS=1, IOSTAT=ios) TEXT
  
  !! Locate end-of-line symbols. On windows, EOL consists of CR+LF, while on 
  !! Linux EOL is only LF, where CR stands for "carriage-return" (CHAR(13)) 
  !! and LF stands for "line feed" (CHAR(10)). 
  CALL ResizeVar(EOL, 1000)
  k = 0
  is_win = .FALSE.
  DO i=1,n
  
    !! Look for the LF symbol
    IF(TEXT(i)==LF) THEN
    
      !! Increase EOL counter
      k = k+1

      !! Check for windows encoding
      IF(i>1 .and. TEXT(i-1)==CR) is_win = .TRUE.
      
      !! Make sure EOL is large enough
      IF(k>SIZE(EOL)) CALL ResizeVar(EOL, SIZE(EOL)+1000)
      
      !! Store the position of the first/only EOL character
      IF(is_win) THEN
        EOL(k) = i-1
      ELSE
        EOL(k) = i
      ENDIF
      
    ENDIF
    
  ENDDO
  
  !! Check for the file not ending with eol symbol
  !print *,"n=", n
  !print *,"k=", k
  !write (2,'(I0)') eol
  IF(is_win .AND. EOL(k)<n-1) THEN
    k = k+1
    IF(k>SIZE(EOL)) CALL ResizeVar(EOL, SIZE(EOL)+1)
    EOL(k) = n
  ELSEIF(.NOT.is_win .AND. EOL(k)<n) THEN
    k = k+1
    IF(k>SIZE(EOL)) CALL ResizeVar(EOL, SIZE(EOL)+1)
    EOL(k) = n+1
  ENDIF
  CALL ResizeVar(EOL, k)
  
  !! Determine the number of columns and rows
  nr = k                  
  IF(k==1) THEN
    nc = EOL(1)-1
  ELSEIF(is_win) THEN
    nc = MAX(EOL(1), MAXVAL(EOL(2:)-EOL(1:k-1)-1)) - 1 
  ELSE
    nc = MAX(EOL(1), MAXVAL(EOL(2:)-EOL(1:k-1))) - 1 
  ENDIF
  
  !print *,"is_win=", is_win
  !write (3,'(I0)') eol
  !print *,"nr=", nr
  !print *,"nc=", nc
    
  !! Allocate the array
  CALL ResizeVar(CONTENT, nr, MAX(1,nc), "")
  
  IF(nc==0) RETURN
    
  !! Fill the array up with values from TEXT
  l = 0
  k = 0
  DO i=1,nr
    DO j=1,nc
      l = l+1
      IF(TEXT(l)==LF) EXIT
      CONTENT(i,j) = TEXT(l)
    ENDDO
    k = k+1
    l = EOL(k)
    IF(is_win) l = l+1
  ENDDO
  
  RETURN
  
END SUBROUTINE ReadTable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReadFile(fn, fu, CONTENT, read_nlines, file_nlines, ncols, colsep, &
                    nskip, skip_ncols, cmt, skip_cmted, attach, max_nl, &
                    only_cmted, close_file, ios, fixed_ncols, ignore_ncols, &
                    learn_ncols, read_ncols, progress, warn_empty)
  IMPLICIT NONE
  INTEGER, PARAMETER                       :: chunk = 2500
  CHARACTER(*), INTENT(INOUT)              :: fn
  INTEGER, INTENT(IN)                      :: fu
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: CONTENT(:,:)
  INTEGER, INTENT(INOUT)                   :: read_nlines, file_nlines, ncols
  CHARACTER, INTENT(INOUT), OPTIONAL       :: colsep(:)
  CHARACTER, INTENT(IN), OPTIONAL          :: cmt
  LOGICAL, INTENT(IN), OPTIONAL            :: only_cmted, skip_cmted, &
                                              attach, close_file, fixed_ncols, &
                                              ignore_ncols, learn_ncols, &
                                              progress, warn_empty, read_ncols
  INTEGER, INTENT(IN), OPTIONAL            :: nskip, skip_ncols, max_nl
  INTEGER, INTENT(OUT), OPTIONAL           :: ios
  CHARACTER(mifll)                         :: line
  CHARACTER(mml), ALLOCATABLE              :: words(:)
  INTEGER                                  :: nlines, clines, max_nl1, &
                                              nskipped, ncolsfound, ios1, &
                                              CONTENT_size, nseps, skip_ncols1
  INTEGER(ikb)                             :: percent
  CHARACTER                                :: firstchar, cmt1, response
  CHARACTER, ALLOCATABLE                   :: colsep1(:)
  LOGICAL                                  :: is_opened, close_file1, &
                                              attach1, skip_cmted1, &
                                              only_cmted1, ignore_ncols1, &
                                              read_ncols1, fixed_ncols1, &
                                              learn_ncols1, progress1, warn_empty1
  CHARACTER(mstl)                          :: text

  !! Get number of lines in each file
  CALL GetFileNumLines(fn, max_nl1, quit=.FALSE.)

  !! If max_nl is not present read as many lines as possible
  IF(PRESENT(max_nl)) max_nl1 = max_nl 

  !! Initial values for variables
  close_file1   = .TRUE.
  cmt1          = def_comment
  attach1       = .TRUE.
  skip_cmted1   = .TRUE.
  only_cmted1   = .FALSE.
  fixed_ncols1  = .FALSE.
  ignore_ncols1 = .FALSE.
  learn_ncols1  = .FALSE.
  read_ncols1   = .FALSE.
  progress1     = .FALSE.
  warn_empty1   = .TRUE.
  skip_ncols1   = 0
  
  IF(PRESENT(colsep)) THEN
    CALL ResizeVar(colsep1, 1, "")
  ELSE
    CALL ResizeVar(colsep1, SIZE(colsep))
  ENDIF
  
  !! Change the value according to input
  IF(PRESENT(close_file))   close_file1 = close_file
  IF(PRESENT(cmt))          cmt1 = cmt
  IF(PRESENT(colsep))       colsep1 = colsep
  IF(PRESENT(attach))       attach1 = attach
  IF(PRESENT(skip_cmted))   skip_cmted1 = skip_cmted
  IF(PRESENT(only_cmted))   only_cmted1 = only_cmted
  IF(PRESENT(fixed_ncols))  fixed_ncols1 = fixed_ncols
  IF(PRESENT(ignore_ncols)) ignore_ncols1 = ignore_ncols
  IF(PRESENT(learn_ncols))  learn_ncols1 = learn_ncols
  IF(PRESENT(read_ncols))   read_ncols1 = read_ncols
  IF(PRESENT(progress))     progress1 = progress
  IF(PRESENT(warn_empty))   warn_empty1 = warn_empty
  IF(PRESENT(skip_ncols))   skip_ncols1 = skip_ncols
  
  !! Don't warn about an empty file if only commented lines are to be read
  IF(only_cmted1) warn_empty1 = .FALSE.

  !! Warn about wrong column count
  IF(ncols<=0) THEN
    IF(fixed_ncols1) THEN
      CALL PrntE("Value in 'ncols' must be positive (ReadFile)!", Q=.TRUE.)
    ELSE
      ncols = 1
      learn_ncols1 = .TRUE.
    ENDIF
  ENDIF

  !! Open the output file if not open and exists)
  INQUIRE(FILE=fn, OPENED=is_opened)      
  IF(.NOT.is_opened) &
    CALL OpenFile(UNIT=fu, FILE=fn, ACTION='R', STATUS="O", FORM='F',&
                  POSITION="R", ACCESS='S')
                  
  !! Set counters
  nskipped = 0
  clines = 0
  nlines = 0

  !! If CONTENT should not be overwritten and it has data in it, store its size
  IF(attach1 .AND. ALLOCATED(CONTENT)) nlines = SIZE(CONTENT,1)
  
  !! Make sure the array is not unallocated at the end of the run  
  IF(.NOT.ALLOCATED(CONTENT)) &
    CALL ResizeVar(CONTENT, MAX(1, max_nl1), ncols, "")
    !CALL ResizeVar(CONTENT, MAX(1, MIN(max_nl1, chunk)), ncols, "")
    
  !! Initialize CONTENT_size and progress variables
  CONTENT_size = SIZE(CONTENT, 1)   
  text = ""
  percent = -1
  ios1 = 0
  
  !! Read file
  DO WHILE(ios1 == 0 .AND. clines < max_nl1)

    !! Increase the read line counter
    clines = clines + 1
    
    !! Print progress
    IF(progress1) CALL PrintCD(text, clines, max_nl1, percent, finish=.TRUE.)

    !! Read the next line of the input file
    line = ""
    READ(fu, '(A)', IOSTAT=ios1, ADVANCE='NO') line
    file_nlines = file_nlines + 1
    
    !! Make sure the content inside line is not streatching its length to almost
    !! maximum, which could mean loss data 
    IF(LEN_TRIM(line)>LEN(line)-10) &
      CALL PrntW("Input file line with length very near the maximum"//&
                 " supported length of "//TRIM(i2cp(mifll))//" found!"//&
                 " Make sure this is OK, some data might have been"//&
                 " not read in! Hint: Change input to a binary format"//&
                 " or increase the value in 'mifll' and recompile.", &
                 skip1=1, skip2=1)
    
    !! Remove end-of-line status
    IF(ANY(ios1==err_eol)) ios1 = 0

    !! If error quit
    IF(ios1/=0) EXIT
    
    !! Save the first character
    firstchar = line(1:1)

    !! If commented line should be skipped, skip this one if commented
    IF(skip_cmted1 .AND. firstchar==cmt1) CYCLE
    
    !! If only commented lines should be read and this one does not start with
    !! a comment, skip it
    IF(only_cmted1 .AND. firstchar/=cmt1) CYCLE

    !! Before reading, skip specified number of lines
    IF(PRESENT(nskip)) THEN
      IF(nskipped < nskip) THEN
        nskipped = nskipped + 1
        CYCLE
      ENDIF 
    ENDIF
  
    !! Empty line, announce
    IF(LEN_TRIM(line)==0) &
      CALL PrntE("Line "//i2c(file_nlines)//" of the file ["//TRIM(fn)//&
                 "] appears to be empty.", Q=.TRUE.)
    
    !! Increase total line counter
    read_nlines = read_nlines + 1

    !! Increase line counter 
    nlines = nlines + 1
    
    !! Make sure CONTENT still has space to store the next line
    IF(nlines > CONTENT_size) THEN
      CONTENT_size = CONTENT_size + chunk
      CALL ResizeVar(CONTENT, CONTENT_size, ncols, "")
    ENDIF

    !! If only 1 columns wanted, save the line into CONTENT and go to next line
    IF(ncols == 1 .AND. .NOT.learn_ncols1) THEN
      CONTENT(nlines,:) = TRIM(line)
      CYCLE
    ENDIF

    !! Make sure the given separator is indeed the separator of the first line
    !! If it doesn't appear to be, change it to one of the candidates
    IF(nlines==1 .AND. .NOT.learn_ncols1 .AND. SIZE(colsep1)==1) THEN
      nseps = ncols - 1 + skip_ncols1
      CALL CheckSeparator(colsep1(1), nseps, TRIM(line), findsep=.TRUE.)
      IF(PRESENT(colsep)) colsep = colsep1
      IF(learn_ncols1) THEN
        ncols = nseps + 1 - skip_ncols1
        CALL ResizeVar(CONTENT, CONTENT_size, ncols, "")
      ENDIF
    ENDIF
    
    !! Split line according to colsep
    CALL SplitString(TRIM(line), colsep1, words, ncolsfound, ios1)
    ncolsfound = ncolsfound - skip_ncols1
    
    IF(nlines==1 .AND. learn_ncols1) ncols = ncolsfound
    
    !! If other error occured, quit and print error message
    IF(ios1/=0) &
      CALL PrntE("An error occured when splitting line "//i2c(file_nlines)//&
                 " of the file ["//TRIM(fn)//"]! (IOSTAT="//TRIM(i2cp(ios1))//&
                 ")", Q=.TRUE.)
         
    !! Check for wrong number of columns
    IF(ncolsfound > ncols .AND. read_ncols1) THEN
      ncolsfound = ncols
    ELSEIF(ncolsfound /= ncols) THEN
      IF(fixed_ncols1) THEN
        CALL PrntE("Line "//i2c(file_nlines)//" of the file ["//TRIM(fn)//&
                   "] contains "//i2c(ncolsfound + skip_ncols1)//" columns"//&
                   " when "//i2c(ncols+skip_ncols1)//" columns expected.", &
                   Q=.TRUE.)
      ELSE
        IF(ignore_ncols1) THEN
          !! If this is the first line in CONTENT which is all empty, change
          !! ncols to ncolsfound, otherwise make sure ncols is not too small
          IF(SIZE(CONTENT,1)==1 .AND. ALL(LEN_TRIM(CONTENT(1,:))==0)) THEN
            ncols = ncolsfound
          ELSE
            ncols = MAX(ncols, ncolsfound)
          ENDIF 
        ELSE
          CALL PrntW("Line "//i2c(file_nlines)//" of the file ["//&
                  TRIM(fn)//"] contains "//i2c(ncolsfound+skip_ncols1)//&
                  " columns when "//i2c(ncols+skip_ncols1)//" columns"//&
                  " expected.", skip1=1)
          IF(PRESENT(ios)) THEN
            CALL AskUser(response, cont=.TRUE.)
            IF(response=="Y") THEN
              ios1 = 0
              ncols = ncolsfound
              CALL PrntW("Expected number of columns was changed to"//&
                                i2c(ncols+skip_ncols1)//".")
            ELSE
              ios1 = 1005
              EXIT
            ENDIF
          ENDIF
          CALL Prnt0("", Q=.TRUE.)
        ENDIF
      ENDIF
    ENDIF
    
    !! Make sure CONTENT is big enough to fit the next line
    IF(SIZE(CONTENT,2) < ncolsfound) &
      CALL ResizeVar(CONTENT, CONTENT_size, ncolsfound, "")
      
    !! Save the current line into CONTENT
    CONTENT(nlines,1:ncolsfound) = words(1+skip_ncols1:ncolsfound+skip_ncols1)
    
  ENDDO
  
  !! Remove extra space for CONTENT
  CALL ResizeVar(CONTENT, MAX(1,nlines), ncols)
  
  !! If no lines were read return an empty array of size 1
  IF(nlines==0) CONTENT = ""

  !! If an error or an end of file occurred, close the file
  IF(ios1/=0) close_file1 = .TRUE.

  !! Print end of reading
  !IF(progress1) CALL PrintCD(text, max_nl1, max_nl1, percent, finish=.TRUE.)
      
  !! If close_file is true, close the file and perform checks for warnings
  IF(close_file1) THEN
    INQUIRE(FILE=fn, OPENED=is_opened)      
    IF(is_opened) CLOSE(fu)

    !! Warn about an empty file
    IF(warn_empty1 .AND. read_nlines == 0) THEN
      CALL PrntW("No data was read from file ["//TRIM(fn)//"]!")

      !! Explanation
      IF(ios1/=0) THEN
        IF(ANY(ios1==err_eof)) THEN
          CALL Prnt("REASON: End-of-file appeared before any lines were"//&
                         " read (IOSTAT "//TRIM(i2cp(ios1))//").")
          CALL Prnt("NOTE: If the file is not empty or all lines are"//&
                         " commented make sure that the end-of-line symbols"//&
                         " are consistent with the current operating system.")
        ELSE
          CALL Prnt("Unknown error (IOSTAT "//TRIM(i2cp(ios1))//").")
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  
  !! Remove end-of-file status
  IF(ANY(ios1==err_eof)) ios1 = 0

  !! Report status if ios present
  IF(PRESENT(ios)) ios = ios1

  RETURN

END SUBROUTINE ReadFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE WriteFile(fn, fu, CONTENT, sep, attach, progress, ios)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)     :: fn
  INTEGER, INTENT(IN)             :: fu
  CHARACTER(*), INTENT(IN)        :: CONTENT(:,:)
  LOGICAL, INTENT(IN), OPTIONAL   :: attach, progress
  CHARACTER, INTENT(IN), OPTIONAL :: sep
  INTEGER, INTENT(OUT), OPTIONAL  :: ios
  LOGICAL                         :: attach1, progress1
  CHARACTER                       :: position, status, sep1
  CHARACTER(mstl)                 :: text
  !CHARACTER(mltl)                 :: line
  INTEGER                         :: i, j, ncol, nrow, ios1
  INTEGER(ikb)                    :: percent
  
  !! Decide whether to rewrite or attach an already existing file
  attach1 = .FALSE.
  IF(PRESENT(attach)) attach1 = attach
  status = "U"
  IF(attach1) THEN
    position = "A"
  ELSE
    position = "R"
    status = "R"
  ENDIF
  
  !! Set the progress counter and separator
  progress1 = .FALSE.
  sep1      = tab
  IF(PRESENT(progress)) progress1 = progress
  IF(PRESENT(sep))      sep1 = sep
  
  !! Open the file
  CALL OpenFile(UNIT=fu, FILE=fn, ACTION='W', STATUS=status, FORM='F',&
                POSITION=position, ACCESS='S')
                
  !! Set variables
  text = ""
  percent = -1
  nrow = SIZE(CONTENT,1)
  ncol = SIZE(CONTENT,2)

  !! Write CONTENT to file
  DO i=1,nrow

    !! Print progress
    IF(progress1) &
      CALL PrintCD(text, INT(i,ikb), INT(nrow,ikb), percent, finish=.TRUE.)

    !! Write the line into file
    WRITE(fu, '(A)', IOSTAT=ios1, ADVANCE='NO') TRIM(CONTENT(i,1))
    DO j=2,ncol
      WRITE(fu, '(A)', IOSTAT=ios1, ADVANCE='NO') sep1//TRIM(CONTENT(i,j))
      IF(ios1/=0) GOTO 100
    ENDDO
    WRITE(fu,'(A)', IOSTAT=ios1) 
    
    !!! Put together the line columns and separator
    !line = TRIM(CONTENT(i,1))
    !DO j=2,ncol
    !  line = TRIM(line)//sep1//CONTENT(i,j)
    !ENDDO
    !
    !!! Write the line into file
    !WRITE(fu, '(A)', IOSTAT=ios1) TRIM(line)
    !IF(ios1/=0) GOTO 100
    
  ENDDO
  
  !! Close the file and exit
  CLOSE(fu)  
  RETURN
  
  !! An error occurred 
  100 CONTINUE

  !! Close the file
  CLOSE(fu)

  !! Return ios
  IF(PRESENT(ios)) THEN
    ios = ios1
    RETURN
  ENDIF
  
  !! Announce error
  IF(ios1/=0) &
    CALL PrntE("Unknown error occurred during writing to file ["//&
               TRIM(fn)//"] (IOSTAT "//i2cp(ios1)//").", Q=.TRUE.)

END SUBROUTINE WriteFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetFilenames(nfiles, fu, file, filelist, delim, comment, quit)
!! Extracts filenames from lines of the file in filelist
  IMPLICIT NONE
  INTEGER, INTENT(OUT)                      :: nfiles
  INTEGER, INTENT(IN)                       :: fu
  CHARACTER(*), ALLOCATABLE, INTENT(INOUT)  :: file(:)
  CHARACTER(*), INTENT(INOUT)               :: filelist
  CHARACTER, INTENT(IN)                     :: delim, comment
  LOGICAL, INTENT(IN), OPTIONAL             :: quit
  CHARACTER(mfl)                            :: filelist_path, line
  INTEGER                                   :: pathsplit, ios
  LOGICAL                                   :: file_exists, quit1
      
  quit1 = .TRUE.
  IF(PRESENT(quit)) quit1 = quit

  nfiles = 0        

  IF(SIZE(file)>0 .AND. LEN_TRIM(file(1))>0) nfiles = nfiles + 1

  IF(LEN_TRIM(filelist)>0) THEN

    filelist_path = ""
    pathsplit = FindLastSubstr(filelist,delim)
    IF(pathsplit>0) filelist_path = filelist(1:pathsplit)

    !! Check whether the file exists
    CALL CheckFileExistence(filelist, exists=file_exists, quit=quit1)

    !! Open it
    CALL OpenFile(UNIT=fu, FILE=filelist, ACTION='R', STATUS='O', FORM='F')
    ios = 0
    DO WHILE(ios==0)
      READ(fu, *, IOSTAT=ios) line
      IF(ios==0) THEN
        !! If line commented or empty skip to the next one 
        IF(line(1:1)==comment .OR. LEN_TRIM(line)==0) CYCLE
        nfiles=nfiles+1
        IF(nfiles>SIZE(file)) CALL ResizeVar(file, nfiles)
        file(nfiles) = TRIM(filelist_path) // TRIM(line)
      ENDIF
    ENDDO
    CLOSE(fu)
  ENDIF

  RETURN
      
END SUBROUTINE GetFilenames

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE AskForNewFilename(filename, log)
  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL :: log
  CHARACTER(*), INTENT(INOUT)   :: filename
  CHARACTER(LEN(filename))      :: filename1
  CHARACTER                     :: response
  LOGICAL                       :: log1
  
  log1 = .TRUE.
  IF(PRESENT(log)) log1 = log
  
  !! Ask user whether to enter new file name
  CALL AskUser(response, newfile=.TRUE.)

  !! If response is NO, then terminate
  IF(response=="N") CALL Terminate(ask=.TRUE., log=log1)
  
  !! If response is NO, then ask for the new filename
  IF(response=="Y") THEN
    12 CONTINUE 
    WRITE (usto,'(A)') &
      "> NOTE: When entering new filename you can use '*' at the beggining"//&
      " of your input as a mask for the previous file name."
    WRITE (usto,'(A)', ADVANCE='NO') "> Enter new filename: "
    READ(*,'(A)') filename1
    IF(LEN_TRIM(filename1)==0) THEN
      CALL Terminate(ask=.TRUE.)
      GOTO 12
    ENDIF
    IF(filename1(1:1)=="*") THEN
      IF(LEN_TRIM(filename1)>1) THEN
        filename1 = filename1(2:)
      ELSE
        filename1 = ""
      ENDIF
      filename1 = TRIM(filename)//TRIM(filename1)
    ENDIF
    filename = filename1
    RETURN
  ENDIF

  RETURN

END SUBROUTINE AskForNewFilename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE CreateLogFile(log_file, chck_exist, time_full, stamp)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)        :: log_file
  LOGICAL, INTENT(IN), OPTIONAL      :: chck_exist
  CHARACTER(*), INTENT(IN), OPTIONAL :: time_full, stamp
  CHARACTER(24)                      :: time_full1, stamp1
  REAL(dpp)                          :: time

  !! Get creation time
  IF(PRESENT(time_full)) THEN
    time_full1 = time_full
  ELSE
    CALL GetCurrentTime(time_full1, time)
  ENDIF

  !! If no log file given determine randomly
  IF(LEN_TRIM(log_file)==0) THEN
    IF(PRESENT(stamp)) THEN
      stamp1 = stamp
    ELSE
      CALL GetTimeStamp(stamp1, time_full1)
    ENDIF
    log_file = program_name//"_"//TRIM(stamp1)//log_ext
  ENDIF
  
  !! Create log file
  CALL OpenFile(UNIT=ulog, FILE=log_file, ACTION='W', STATUS="U", &
                FORM='F', POSITION="R", ACCESS='S', ask=.TRUE., &
                chck_exist = chck_exist, overwrite=overwrite_files, &
                islog=.TRUE.)
  CLOSE(ulog)

  !! Print log file header
  CALL Prnt(TRIM(program_fullname)//" log file created at "//time_full1, &
                 screen=.FALSE.)
  
  !! Print header
  CALL PrintProgramHeader()

  RETURN

END SUBROUTINE CreateLogFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION TrimStr1(str) RESULT(tstr)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: str
  CHARACTER(LEN_TRIM(str)) :: tstr
  
  tstr = TRIM(str)
  RETURN
  
END FUNCTION TrimStr1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION TrimStrN(str) RESULT(tstr)
!! This function DOES NOT WORK PROPERLY. I cannot figure out a way of returning
!! variable lengths to each of the components of str
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)   :: str(:)
  CHARACTER(MaxLenTrim(str)) :: tstr(SIZE(str))
  INTEGER                    :: i
  
  DO i=1,SIZE(str)
    tstr(i) = TRIM(str(i))
  ENDDO
  RETURN
  
END FUNCTION TrimStrN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE RemoveSubstrEnd1(str, substr, trimmed)
!! Removes substr from the end of scaler str (if the end fits)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: str
  CHARACTER(*), INTENT(IN)       :: substr
  LOGICAL, INTENT(OUT), OPTIONAL :: trimmed
  INTEGER                        :: l, sl

  !! Assume no removal will be done
  IF(PRESENT(trimmed)) trimmed = .FALSE.
  
  !! Quit if string sizes do not allow any removal
  IF(LEN_TRIM(substr)==0 .OR. LEN_TRIM(str)==0) RETURN
  IF(LEN_TRIM(substr) > LEN_TRIM(str)) RETURN
  
  !! Check for presence of substr at the end of str
  l = LEN_TRIM(str)
  sl = LEN_TRIM(substr)
  
  !! Presence found
  IF(str(l-sl+1:) == substr) THEN
    str = str(:l-sl)
    IF(PRESENT(trimmed)) trimmed = .TRUE.
  ENDIF 
  
  RETURN
  
END SUBROUTINE RemoveSubstrEnd1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE RemoveSubstrEndN(str, substr, trimmed)
!! Removes substr from the end of all components of vector str (if the end fits) 
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: str(:)
  CHARACTER(*), INTENT(IN)       :: substr
  LOGICAL, INTENT(OUT), OPTIONAL :: trimmed
  INTEGER                        :: i

  DO i=1,SIZE(str)
    CALL RemoveSubstrEnd1(str(i), substr, trimmed)
  ENDDO
    
  RETURN
  
END SUBROUTINE RemoveSubstrEndN           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

PURE FUNCTION MaxLenTrim(str) RESULT(l)
!! Returns the maximum of trimmed lengths of components of vector str 
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: str(:)
  INTEGER                  :: i, l

  l = 0
  DO i=1,SIZE(str)
    IF(LEN_TRIM(str(i))>l) l = LEN_TRIM(str(i)) 
  ENDDO
    
  RETURN
  
END FUNCTION MaxLenTrim           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!**************************************************************************!!  
!!**************          TIME RELATED ROUTINES         ********************!! 
!!**************************************************************************!!  

SUBROUTINE GetCurrentTime(time_full, time_s, time_h)
!! Returns current time and date time_full (format HH:MM:SS.sss MON DAY YEAR),
!! in seconds in time_s and/or in hours in time_h. If all optional variables
!! are missing, then the user friendly formatted time gets printed to stdout 
  IMPLICIT NONE
#ifdef __INTEL_COMPILER
  INTEGER, EXTERNAL :: TIME
#endif
  CHARACTER(*), INTENT(OUT), OPTIONAL :: time_full
  REAL(dpp), INTENT(OUT), OPTIONAL    :: time_s, time_h
  CHARACTER(24)                       :: time_full1
  CHARACTER(8)                        :: dat
  CHARACTER(10)                       :: tim
  REAL(dpp)                           :: time_num, time_ms, time_in_s
  !INTEGER                             :: TIME
  !CHARACTER(3)                        :: months(12)

  !! Store month names
  !months = (/"Jan","Feb","Mar","Apr","May","Jun", &
  !           "Jul","Aug","Sep","Oct","Nov","Dec" /)
   
  !! Get the current time, date and time in seconds
  CALL DATE_AND_TIME(dat, tim)
  time_in_s = REAL(TIME(), dpp)

  !! Store the time numberically and compute miliseconds   
  time_num = char2real(tim)
  time_ms = REAL(time_num - INT(time_num), dpp)

  !! Save time and date in user friendly format
  time_full1 = tim(1:2)//colon//tim(3:4)//colon//tim(5:)//space//&
               dat(1:4)//"-"//dat(5:6)//"-"//dat(7:8)
               !dat(1:8)
               !months(c2i(dat(5:6)))//space//dat(7:8)//space//dat(1:4)

  !! If wanted, return time and date in user friendly format
  IF(PRESENT(time_full)) time_full = time_full1

  !! If wanted, return time also in the number of seconds and miliseconds
  IF(PRESENT(time_s)) time_s = time_in_s + time_ms

  !! If wanted, return time also in hours
  IF(PRESENT(time_h)) time_h = ( time_in_s + time_ms ) / 3600 

  !! Print the time in user friendly format               
  IF(.NOT.ANY((/PRESENT(time_full), PRESENT(time_s), PRESENT(time_h)/))) &
    WRITE(*,'(A)') time_full1
              
  RETURN

END SUBROUTINE GetCurrentTime
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION GetTimeDifference(timeA, timeB) RESULT(time_dif)
!! Computes a suitably formatted absolute value of the difference of times  
!! timeB - timeA (input in seconds). If both timeA and timeB are missing, it 
!! returns formatted zero time and if one of these is missing, it returns a 
!! formatted version of the one present. 
  IMPLICIT NONE
  REAL(dpp), INTENT(IN), OPTIONAL :: timeA, timeB
  REAL(dpp)                       :: time, time1, time2
  CHARACTER(24)                   :: time_dif
  CHARACTER(5)                    :: hours
  INTEGER                         :: difh, difm, difs, difss

    !! Assign values
    time1 = zero
    time2 = zero
    IF(PRESENT(timeA)) time1 = timeA
    IF(PRESENT(timeB)) time2 = timeB

    !! Get the time difference is seconds
    time = ABS(time2 - time1)
    
    !! Get the time differences in hours, minutes
    difh  = INT(time/3600)
    difm  = INT((time-difh*3600)/60)    
    difs  = INT(time-difh*3600-difm*60)
    difss = INT((time-difh*3600-difm*60-difs)*1000 + half)
    
    WRITE (hours,'(I5.2)') difh
    CALL DelAllSubstr(hours,space)

    WRITE(time_dif,'(2A,I2.2,A,I2.2,A,I3.3)') TRIM(hours),":",difm,":",difs,".",difss

END FUNCTION GetTimeDifference
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE GetTimeStamp(stamp, time_full)
  IMPLICIT NONE
  CHARACTER(*), INTENT(OUT)          :: stamp
  CHARACTER(*), INTENT(IN), OPTIONAL :: time_full
  CHARACTER(24)                      :: time_full1
  REAL(dpp)                          :: time_num 
  
  !! If time_full is present and long enough, don't get it
  IF(PRESENT(time_full)) THEN
    IF(LEN(time_full)>=14) THEN
      time_full1 = time_full
      GOTO 10
    ENDIF
  ENDIF
  
  !! Otherwise get time_full
  CALL GetCurrentTime(time_full1, time_num)

  10 CONTINUE

  !! Create time stamp
  stamp = lowcasef(time_full1(14:)//time_full1(:8)//time_full1(10:12))
  
  !! Remove spaces and colons from stamp
  CALL ReplaceAllSubstr(stamp, space, "_")
  CALL ReplaceAllSubstr(stamp, colon, "_")
  CALL ReplaceAllSubstr(stamp, "_", "")
  
  RETURN

END SUBROUTINE GetTimeStamp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!**************************************************************************!!  
!!**************        FILE DEALLOCATION ROUTINES      ********************!! 
!!**************************************************************************!!  

SUBROUTINE Dealloc1(A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12)
!! Deallocates up to 13 arrays
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: A0(:)
  CHARACTER(*), INTENT(INOUT), OPTIONAL, ALLOCATABLE :: A1(:), A2(:), A3(:), &
      A4(:), A5(:), A6(:), A7(:), A8(:), A9(:), A10(:), A11(:), A12(:)
                                        
#include "epi_inc_deallocate.F90"

  RETURN
END SUBROUTINE Dealloc1       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Dealloc2(A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12)
!! Deallocates up to 12 arrays
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT), ALLOCATABLE :: A0(:,:)
  CHARACTER(*), INTENT(INOUT), OPTIONAL, ALLOCATABLE :: A1(:,:), A2(:,:), &
      A3(:,:), A4(:,:), A5(:,:), A6(:,:), A7(:,:), A8(:,:), A9(:,:), &
      A10(:,:), A11(:,:), A12(:,:)
                                        
#include "epi_inc_deallocate.F90"

  RETURN
END SUBROUTINE Dealloc2       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!                         ONE DIMENSIONAL OBJECTS                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Real(vector, newlength, newvalue, allvalue)
  !! Changes size of a real vector to newlength  
  IMPLICIT NONE
#define _Resize1Real
  CHARACTER(11), PARAMETER                  :: routine_name = "Resize1Real"
  REAL(dpp), ALLOCATABLE, INTENT(INOUT)     :: vector(:)
  REAL(dpp), INTENT(IN), OPTIONAL           :: newvalue, allvalue
  REAL(dpp), ALLOCATABLE                    :: temp(:)

#include "epi_inc_resize_vector.F90"
  
  RETURN

#undef _Resize1Real
END SUBROUTINE Resize1Real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Int_ikn(vector, newlength, newvalue, allvalue)
  !! Changes size of a integer vector to newlength  
  IMPLICIT NONE
#define _Resize1Int_ikn
  CHARACTER(14), PARAMETER                  :: routine_name = "Resize1Int_ikn"
  INTEGER(ikn), ALLOCATABLE, INTENT(INOUT)  :: vector(:)
  INTEGER(ikn), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  INTEGER(ikn), ALLOCATABLE                 :: temp(:)

#include "epi_inc_resize_vector.F90"
  
  RETURN

#undef _Resize1Int_ikn
END SUBROUTINE Resize1Int_ikn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Int_iks(vector, newlength, newvalue, allvalue)
  !! Changes size of a integer vector to newlength  
  IMPLICIT NONE
#define _Resize1Int_iks
  CHARACTER(14), PARAMETER                  :: routine_name = "Resize1Int_iks"
  INTEGER(iks), ALLOCATABLE, INTENT(INOUT)  :: vector(:)
  INTEGER(iks), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  INTEGER(iks), ALLOCATABLE                 :: temp(:)

#include "epi_inc_resize_vector.F90"

  RETURN

#undef _Resize1Int_iks
END SUBROUTINE Resize1Int_iks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Int_ikb(vector, newlength, newvalue, allvalue)
!! Changes size of a integer vector to newlength  
  IMPLICIT NONE
#define _Resize1Int_ikb
  CHARACTER(14), PARAMETER                  :: routine_name = "Resize1Int_ikb"
  INTEGER(ikb), ALLOCATABLE, INTENT(INOUT)  :: vector(:)
  INTEGER(ikb), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  INTEGER(ikb), ALLOCATABLE                 :: temp(:)

#include "epi_inc_resize_vector.F90"

 RETURN

#undef _Resize1Int_ikb
END SUBROUTINE Resize1Int_ikb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! SUBROUTINE Resize1IntL(vector, newlength, newvalue, allvalue)
! !! Changes size of a integer vector to newlength  
!   IMPLICIT NONE
! #define _Resize1IntL
!   CHARACTER(11), PARAMETER :: routine_name = "Resize1IntL"
!   INTEGER(ikl), ALLOCATABLE, INTENT(INOUT)  :: vector(:)
!   INTEGER(ikl), INTENT(IN), OPTIONAL        :: newvalue, allvalue
!   INTEGER(ikl), ALLOCATABLE                 :: temp(:)
! 
! #include "epi_inc_resize_vector.F90"
! 
!  RETURN
! 
! #undef _Resize1IntL
! END SUBROUTINE Resize1IntL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Char(vector, newlength, newvalue, allvalue)
  !! Changes size of a character vector to newlength  
  IMPLICIT NONE
#define _Resize1Char
  CHARACTER(11), PARAMETER                  :: routine_name = "Resize1Char"
  CHARACTER(*), ALLOCATABLE, INTENT(INOUT)  :: vector(:)
  CHARACTER(*), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  CHARACTER(LEN(vector)), ALLOCATABLE       :: temp(:)

#include "epi_inc_resize_vector.F90"

  RETURN

#undef _Resize1Char
END SUBROUTINE Resize1Char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1Logi(vector, newlength, newvalue, allvalue)
  !! Changes size of a integer vector to newlength  
  IMPLICIT NONE
#define _Resize1Logi
  CHARACTER(11), PARAMETER                 :: routine_name = "Resize1Logi"
  LOGICAL, ALLOCATABLE, INTENT(INOUT)      :: vector(:)
  LOGICAL, INTENT(IN), OPTIONAL            :: newvalue, allvalue
  LOGICAL, ALLOCATABLE                     :: temp(:)

#include "epi_inc_resize_vector.F90"

  RETURN

#undef _Resize1Logi
END SUBROUTINE Resize1Logi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1TEMP(vector, newlength, newvalue, allvalue)
  !! Changes size of a character vector to newlength  
  IMPLICIT NONE
#define _Resize1TEMP
  CHARACTER(11), PARAMETER                  :: routine_name = "Resize1TEMP"
  TYPE(TMPTYPE), ALLOCATABLE, INTENT(INOUT) :: vector(:)
  TYPE(TMPTYPE), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  TYPE(TMPTYPE), ALLOCATABLE                :: temp(:)

#include "epi_inc_resize_vector.F90"

  RETURN

#undef _Resize1TEMP
END SUBROUTINE Resize1TEMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize1TEMPmin(vector, newlength, newvalue, allvalue)
  !! Changes size of a character vector to newlength  
  IMPLICIT NONE
#define _Resize1TEMPmin
  CHARACTER(14), PARAMETER                     :: routine_name = "Resize1TEMPmin"
  TYPE(TMPTYPEmin), ALLOCATABLE, INTENT(INOUT) :: vector(:)
  TYPE(TMPTYPEmin), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  TYPE(TMPTYPEmin), ALLOCATABLE                :: temp(:)

#include "epi_inc_resize_vector.F90"

  RETURN

#undef _Resize1TEMPmin
END SUBROUTINE Resize1TEMPmin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!                         TWO DIMENSIONAL OBJECTS                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize2Real(array, newdim1, newdim2, newvalue, allvalue)
  !! Changes size of a real array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize2Real
  CHARACTER(11), PARAMETER                 :: routine_name = "Resize2Real"
  REAL(dpp), ALLOCATABLE, INTENT(INOUT)    :: array(:,:)
  REAL(dpp), INTENT(IN), OPTIONAL          :: newvalue, allvalue
  REAL(dpp), ALLOCATABLE                   :: temp(:,:)
  
#include "epi_inc_resize_array.F90"

  RETURN

#undef _Resize2Real
END SUBROUTINE Resize2Real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize2Int_ikn(array, newdim1, newdim2, newvalue, allvalue)
  !! Changes size of a integer array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize2Int_ikn
  CHARACTER(14), PARAMETER                 :: routine_name = "Resize2Int_ikn"
  INTEGER(ikn), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
  INTEGER(ikn), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  INTEGER(ikn), ALLOCATABLE                :: temp(:,:)

#include "epi_inc_resize_array.F90"

  RETURN

#undef _Resize2Int_ikn
END SUBROUTINE Resize2Int_ikn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize2Int_iks(array, newdim1, newdim2, newvalue, allvalue)
  !! Changes size of a integer array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize2Int_iks
  CHARACTER(14), PARAMETER                 :: routine_name = "Resize2Int_iks"
  INTEGER(iks), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
  INTEGER(iks), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  INTEGER(iks), ALLOCATABLE                :: temp(:,:)
  
#include "epi_inc_resize_array.F90"

  RETURN

#undef _Resize2Int_iks
END SUBROUTINE Resize2Int_iks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize2Int_ikb(array, newdim1, newdim2, newvalue, allvalue)
!! Changes size of a integer array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize2Int_ikb
  CHARACTER(14), PARAMETER                 :: routine_name = "Resize2Int_ikb"
  INTEGER(ikb), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
  INTEGER(ikb), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  INTEGER(ikb), ALLOCATABLE                :: temp(:,:)

#include "epi_inc_resize_array.F90"

  RETURN

#undef _Resize2Int_ikb
END SUBROUTINE Resize2Int_ikb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! SUBROUTINE Resize2IntL(array, newdim1, newdim2, newvalue, allvalue)
! !! Changes size of a integer array to newdim2 x newdim2   
!   IMPLICIT NONE
! #define _Resize2IntL
!   CHARACTER(11), PARAMETER                 :: routine_name = "Resize2IntL"
!   INTEGER(ikl), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
!   INTEGER(ikl), INTENT(IN), OPTIONAL       :: newvalue, allvalue
!   INTEGER(ikl), ALLOCATABLE                :: temp(:,:)
! 
! #include "epi_inc_resize_array.F90"
! 
!   RETURN
! 
! #undef _Resize2IntL
! END SUBROUTINE Resize2IntL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize2Char(array, newdim1, newdim2, newvalue, allvalue)
  !! Changes size of a character array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize2Char
  CHARACTER(11), PARAMETER                 :: routine_name = "Resize2Char"
  CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
  CHARACTER(*), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  CHARACTER(LEN(array)), ALLOCATABLE       :: temp(:,:)
    
#include "epi_inc_resize_array.F90"

  RETURN

#undef _Resize2Char
END SUBROUTINE Resize2Char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!                       THREE DIMENSIONAL OBJECTS                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize3Real(array, newdim1, newdim2, newdim3, newvalue, &
                            allvalue)
  !! Changes size of a real array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize3Real
  CHARACTER(11), PARAMETER                :: routine_name = "Resize3Real"
  REAL(dpp), ALLOCATABLE, INTENT(INOUT)   :: array(:,:,:)
  REAL(dpp), INTENT(IN), OPTIONAL         :: newvalue, allvalue
  REAL(dpp), ALLOCATABLE                  :: temp(:,:,:)

#include "epi_inc_resize_array3.F90"

  RETURN

#undef _Resize3Real
END SUBROUTINE Resize3Real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize3Int_ikn(array, newdim1, newdim2, newdim3, newvalue, allvalue)
  !! Changes size of a integer array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize3Int_ikn
  CHARACTER(14), PARAMETER                  :: routine_name = "Resize3Int_ikn"
  INTEGER(ikn), ALLOCATABLE, INTENT(INOUT)  :: array(:,:,:)
  INTEGER(ikn), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  INTEGER(ikn), ALLOCATABLE                 :: temp(:,:,:)

#include "epi_inc_resize_array3.F90"

  RETURN

#undef _Resize3Int_ikn
END SUBROUTINE Resize3Int_ikn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize3Int_iks(array, newdim1, newdim2, newdim3, newvalue, allvalue)
  !! Changes size of a integer array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize3Int_iks
  CHARACTER(14), PARAMETER                  :: routine_name = "Resize3Int_iks"
  INTEGER(iks), ALLOCATABLE, INTENT(INOUT)  :: array(:,:,:)
  INTEGER(iks), INTENT(IN), OPTIONAL        :: newvalue, allvalue
  INTEGER(iks), ALLOCATABLE                 :: temp(:,:,:)

#include "epi_inc_resize_array3.F90"

  RETURN

#undef _Resize3Int_iks
END SUBROUTINE Resize3Int_iks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize3Int_ikb(array, newdim1, newdim2, newdim3, newvalue, allvalue)
 !! Changes size of a integer array to newdim2 x newdim2   
 IMPLICIT NONE
#define _Resize3Int_ikb
 CHARACTER(14), PARAMETER                 :: routine_name = "Resize3Int_ikb"
 INTEGER(ikb), ALLOCATABLE, INTENT(INOUT) :: array(:,:,:)
 INTEGER(ikb), INTENT(IN), OPTIONAL       :: newvalue, allvalue
 INTEGER(ikb), ALLOCATABLE                :: temp(:,:,:)

#include "epi_inc_resize_array3.F90"

 RETURN

#undef _Resize3Int_ikb
END SUBROUTINE Resize3Int_ikb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! SUBROUTINE Resize3IntL(array, newdim1, newdim2, newdim3, newvalue, allvalue)
!  !! Changes size of a integer array to newdim2 x newdim2   
!  IMPLICIT NONE
! #define _Resize3IntL
!  CHARACTER(11), PARAMETER                 :: routine_name = "Resize3IntL"
!  INTEGER(ikl), ALLOCATABLE, INTENT(INOUT) :: array(:,:,:)
!  INTEGER(ikl), INTENT(IN), OPTIONAL       :: newvalue, allvalue
!  INTEGER(ikl), ALLOCATABLE                :: temp(:,:,:)
! 
! #include "epi_inc_resize_array3.F90"
! 
!  RETURN
! 
! #undef _Resize3IntL
! END SUBROUTINE Resize3IntL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE Resize3Char(array, newdim1, newdim2, newdim3, newvalue, allvalue)
  !! Changes size of a character array to newdim2 x newdim2   
  IMPLICIT NONE
#define _Resize3Char
  CHARACTER(11), PARAMETER                 :: routine_name = "Resize3Char"
  CHARACTER(*), ALLOCATABLE, INTENT(INOUT) :: array(:,:,:)
  CHARACTER(*), INTENT(IN), OPTIONAL       :: newvalue, allvalue
  CHARACTER(LEN(array)), ALLOCATABLE       :: temp(:,:,:)

#include "epi_inc_resize_array3.F90"

  RETURN

#undef _Resize3Char
END SUBROUTINE Resize3Char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!**************************************************************************!!  
!!******************     TYPE CONVERSION ROUTINES     **********************!! 
!!**************************************************************************!!  

FUNCTION r2cp(rl) RESULT(char)
!! Converts integer number to a char
  IMPLICIT NONE
  INTEGER, PARAMETER             :: power = 11
  REAL(dpp), PARAMETER           :: eps = ten**(-power)
  REAL(dpp), INTENT(IN)          :: rl
  !CHARACTER(INT(LOG10(DBLE(MAX(1,ABS(INT(rl))))))+1-MAX(-1,MIN(INT(rl),0))) :: char
  CHARACTER(mstl)                :: char
  CHARACTER(mstl)                :: format1
  REAL(dpp)                      :: difference, rl1
  INTEGER                        :: digits1, digits2

  !! If rl is integer, just return integer
  IF(ABS(REAL(INT(rl, ikb), dpp) - rl) < eps) THEN
    char = i2cp(INT(rl, ikb))
  ELSE
    !! If true real then do some rounding to get rid of unnecessary digits
    digits1 = LEN_TRIM(i2cp(INT(rl)))
    digits2 = power
    difference = zero
    DO WHILE (difference<eps .AND. digits2>1)
      digits2 = digits2 - 1
      rl1 = round(rl, digits2-1)
      difference = ABS(rl-rl1)
    ENDDO
    !! And convert the real to a string
    WRITE(format1,'(A,I0,A,I0,A)') '(F',digits1+digits2+1,'.',digits2,')' 
    WRITE(char, format1) rl
  ENDIF
  
  RETURN
  
END FUNCTION r2cp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION real2char_sp(rl, iostat, c, s, mindigs) RESULT(char)
!! Converts single precision real number to a string
  IMPLICIT NONE
  REAL(4), INTENT(IN)            :: rl
  CHARACTER(mstl)                :: char
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  LOGICAL, INTENT(IN), OPTIONAL  :: c, s
  INTEGER, INTENT(IN), OPTIONAL  :: mindigs
  
  char = real2char_dp(REAL(rl, 8), iostat, c, s, mindigs)
  
  RETURN
  
END FUNCTION real2char_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION real2char_dp(rl, iostat, c, s, mindigs) RESULT(char)
!! Converts dpp real number to a string
  IMPLICIT NONE
  INTEGER, PARAMETER             :: power = 11 ! Keep this at most 11, otherwise the trimmed print doesn't work
  REAL(8), PARAMETER             :: eps = ten**(-power)
  REAL(8), INTENT(IN)            :: rl
  CHARACTER(mstl)                :: char
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  LOGICAL, INTENT(IN), OPTIONAL  :: c, s
  INTEGER, INTENT(IN), OPTIONAL  :: mindigs
  CHARACTER(mstl)                :: format1
  REAL(8)                        :: difference, rl0, rl1
  INTEGER                        :: ios, digits1, digits2, mindigs1, idec
  LOGICAL                        :: commas, spaces, negative

  !! Default ios is zero
  IF(PRESENT(iostat)) iostat = 0

  !! Minimum number of digits
  mindigs1 = -1
  IF(PRESENT(mindigs)) mindigs1 = mindigs

  ios = 0
  
  !! If rl is integer, just return integer
  IF(rl == zero) THEN
  
    char = "0"
  
  ELSEIF(LOG10(rl) < -power) THEN

    WRITE(char, '(ES11.4)', IOSTAT=ios) rl
    
  ELSEIF(ABS(REAL(INT(rl, ikb), dpp) - rl) < eps) THEN
    
    char = TRIM(i2cp(INT(rl, ikb)))
  
  !! If true real then do some rounding to get rid of unnecessary digits
  ELSE

    !! Store the sign separately
    negative = rl < zero
    rl0 = ABS(rl) 
  
    !! Figure out the number of digits
    digits1 = LEN_TRIM(i2cp(INT(rl0)))
    digits2 = power
    difference = zero
    DO WHILE (difference<eps .AND. digits2>1)
      digits2 = digits2 - 1
      rl1 = round(rl0, digits2 - 1)
      difference = ABS(rl0-rl1)
    ENDDO
    
    !! Do the conversion
    WRITE(format1,'(A,I0,A,I0,A)') '(F',digits1+digits2+1,'.',digits2,')' 
    WRITE(char, format1, IOSTAT=ios) rl0
    
    !! Restore the sign
    IF(negative) char = "-"//TRIM(char)
    
  ENDIF
  
  !! Error occured
  IF(ios /= 0) THEN
    IF(PRESENT(iostat)) THEN
      iostat = ios
      RETURN
    ENDIF 
    CALL PrntE("Argument must be a real (IOSTAT "//i2cp(ios)//").", Q=.TRUE.)
  ENDIF

  !! Add spaces and/or commas
  commas = add_commas
  spaces = add_spaces
  IF(PRESENT(c)) commas = c
  IF(PRESENT(s)) spaces = s
  CALL AddSpacesCommas(char, c, s)
  
  !! Append trailing zeros if mindigs1 positive
  IF(mindigs1 > 0) THEN
    
    !! Locate the decimal point
    idec = INDEX(char, ".")
  
    !! Determine how many digits after decimal point and attach zeros
    IF(idec==0) THEN
      CALL AttachString(char, "0", mindigs1, ".")
    ELSEIF(idec==LEN_TRIM(char)) THEN
      CALL AttachString(char, "0", mindigs1)
    ELSE
      mindigs1 = mindigs1 - LEN_TRIM(char(idec+1:)) 
      CALL AttachString(char, "0", mindigs1)
    ENDIF
    
  ENDIF
  
  RETURN
  
END FUNCTION real2char_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

PURE FUNCTION si2cp(in) RESULT(char)
!! PURE function that converts integer to a string
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: in
  CHARACTER(INT(LOG10(MAX(one,ABS(DBLE(in)))))+1-MAX(-1,MIN(INT(in),0))) :: char

  !! Convert integer into character
  WRITE(char, '(I0)') in
  
  RETURN
   
END FUNCTION si2cp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

PURE FUNCTION ni2cp(in) RESULT(char)
!! PURE function that converts integer to a string
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: in
  CHARACTER(INT(LOG10(DBLE(MAX(1,ABS(in)))),ikn)+1-MAX(-1,MIN(in,0))) :: char
  !CHARACTER :: char
  
  !! Convert integer into character
  WRITE(char, '(I0)') in
  
  RETURN
   
END FUNCTION ni2cp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

PURE FUNCTION bi2cp(in) RESULT(char)
!! PURE function that converts integer to a string
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: in
  CHARACTER(INT(LOG10(DBLE(MAX(one,ABS(DBLE(in))))))+1-INT(MAX(-one,MIN(DBLE(in),zero)))) :: char
  !CHARACTER(50) :: char

  !! Convert integer into character
  WRITE(char, '(I0)') in
  
  RETURN
   
END FUNCTION bi2cp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! PURE FUNCTION li2cp(in) RESULT(char)
! !! PURE function that converts integer to a string
!   IMPLICIT NONE
!   INTEGER(ikl), INTENT(IN) :: in
!   CHARACTER(INT(LOG10(DBLE(MAX(one,ABS(DBLE(in))))))+1-INT(MAX(-one,MIN(DBLE(in),zero)))) :: char
!   !CHARACTER(50) :: char
! 
!   !! Convert integer into character
!   WRITE(char, '(I0)') in
!   
!   RETURN
!    
! END FUNCTION li2cp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION si2c(in, iostat, c, s, length, ity) RESULT(char)
!! Converts integer number to a string
!! Length of the output string is determined either based on length of in,
!! based on the parameter length (by adding spaces at the beggining), or based
!! on integer type given in ity
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: in
  
#include "epi_inc_convert_i2c.F90"

  RETURN
   
END FUNCTION si2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION ni2c(in, iostat, c, s, length, ity) RESULT(char)
!! Converts integer number to a string
!! Length of the output string is determined either based on length of in,
!! based on the parameter length (by adding spaces at the beggining), or based
!! on integer type given in ity
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: in
  
#include "epi_inc_convert_i2c.F90"

  RETURN
   
END FUNCTION ni2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION li2c(in, iostat, c, s, length, ity) RESULT(char)
!! Converts integer number to a string
!! Length of the output string is determined either based on length of in,
!! based on the parameter length (by adding spaces at the beggining), or based
!! on integer type given in ity
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: in
  
#include "epi_inc_convert_i2c.F90"
  
  RETURN
   
END FUNCTION li2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION si2cN(in, iostat, c, s, length, ity) RESULT(char)
!! Multidimentional si2c
  IMPLICIT NONE
  INTEGER(iks), INTENT(IN) :: in(:)

#include "epi_inc_convert_i2c_array.F90"

  RETURN
   
END FUNCTION si2cN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION ni2cN(in, iostat, c, s, length, ity) RESULT(char)
!! Multidimentional ni2c
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: in(:)

#include "epi_inc_convert_i2c_array.F90"

  RETURN
   
END FUNCTION ni2cN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION li2cN(in, iostat, c, s, length, ity) RESULT(char)
!! Multidimentional li2c
  IMPLICIT NONE
  INTEGER(ikb), INTENT(IN) :: in(:)

#include "epi_inc_convert_i2c_array.F90"
  
  RETURN
   
END FUNCTION li2cN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

PURE FUNCTION c2rp(string) RESULT(number)
!! PURE function that converts a string to a double precision real
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)       :: string
  REAL(dpp)                      :: number

  !! If string empty, return
  IF(LEN_TRIM(string)==0) THEN
    number = zero
    RETURN
  ENDIF
  
  !! Do the conversion of string to real
  READ(string, *) number

  RETURN

END FUNCTION c2rp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION char2real(string, iostat) RESULT(number)
!! Converts number string to a double precision real number
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)       :: string
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  REAL(dpp)                      :: number
  INTEGER                        :: iostat1

  !! If string empty, return
  IF(LEN_TRIM(string)==0) THEN
    number = zero
    IF(PRESENT(iostat)) iostat = -9
    RETURN
  ENDIF
  
  !! Do the conversion of string to real
  READ(string, *, IOSTAT=iostat1) number

  !! Return iostat
  IF(PRESENT(iostat)) iostat = 0

  !! Check for error status
  IF(iostat1 /= 0) number = zero
  
  RETURN

END FUNCTION char2real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION char2real1(string, iostat) RESULT(number)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)       :: string
  REAL(dpp)                      :: number
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  INTEGER                        :: iostat1, exit_code
  CHARACTER(mstl)                :: string1

  10 CONTINUE

  !! Attempt conversion
  number = char2real(string, iostat1)
  
  !! Check for errors
  IF(iostat1/=0) THEN
    IF(PRESENT(iostat)) THEN
      iostat = iostat1
    ELSE
      CALL PrntE("Cannot convert string '"//TRIM(string)//"' to REAL.", &
                      log=silent, Q=silent)
      string1 = string
      exit_code = 0
      CALL AskForNewValue(string1, exit_code)
      IF(exit_code/=0) CALL Terminate(ask=.TRUE.)
      GOTO 10 
    ENDIF
  ENDIF

  RETURN

END FUNCTION char2real1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION char2realN(string, iostat) RESULT(number)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)       :: string(:)
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  REAL(dpp)                      :: number(SIZE(string))
  INTEGER                        :: i
  
  DO i=1,SIZE(string)
    number(i) = char2real(string(i), iostat)
  ENDDO

  RETURN

END FUNCTION char2realN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION nchar2int(string, iostat) RESULT(number)
!! Converts number string to a double precision real number
  IMPLICIT NONE
  INTEGER :: number

#include "epi_inc_convert_c2i.F90"
  
  RETURN

END FUNCTION nchar2int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION nchar2intN(string, iostat) RESULT(number)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: string(:)
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  INTEGER                        :: number(SIZE(string)), i
  
  DO i=1,SIZE(string)
    number(i) = nchar2int(string(i), iostat)
  ENDDO

  RETURN

END FUNCTION nchar2intN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION nchar2int8(string, iostat) RESULT(number)
!! Converts number string to a double precision real number
  IMPLICIT NONE
  INTEGER(8) :: number

#include "epi_inc_convert_c2i.F90"
  
  RETURN

END FUNCTION nchar2int8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION nchar2int8N(string, iostat) RESULT(number)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: string(:)
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  INTEGER(8)                     :: number(SIZE(string))
  INTEGER                        :: i
  
  DO i=1,SIZE(string)
    number(i) = nchar2int8(string(i), iostat)
  ENDDO

  RETURN

END FUNCTION nchar2int8N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION lchar2int(string, iostat) RESULT(number)
!! Converts number string to a double precision real number
  IMPLICIT NONE
  INTEGER(ikb)                   :: number

#include "epi_inc_convert_c2i.F90"
  
  RETURN

END FUNCTION lchar2int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION lchar2intN(string, iostat) RESULT(number)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)    :: string(:)
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  INTEGER(ikb)                   :: number(SIZE(string))
  INTEGER                        :: i
  
  DO i=1,SIZE(string)
    number(i) = lchar2int(string(i), iostat)
  ENDDO

  RETURN

END FUNCTION lchar2intN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION IsInteger(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: string 
  LOGICAL                  :: IsInteger 
  INTEGER                  :: i
  
  !! Initilize return value to false
  IsInteger = .FALSE.

  !! Look for invalid characters and if found any, exit with a false value
  DO i=1,LEN_TRIM(string)
    IF(i==1 .AND. (string(1:1)=="-" .OR. string(1:1)=="+")) CYCLE
    IF(string(i:i)<'0' .OR. string(i:i)>'9') RETURN
  END DO
  
  !! All characters were OK, the string is a (possibly negative) integer
  IsInteger = .TRUE.

  RETURN

END FUNCTION IsInteger

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION round_int(A)
  REAL(dpp), INTENT(IN) :: A
  INTEGER               :: round_int
  
  round_int = INT(A+half)
  
END FUNCTION round_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION round(A, n)
  REAL(dpp), INTENT(IN)         :: A
  INTEGER, OPTIONAL, INTENT(IN) :: n
  REAL(dpp)                     :: round, B
  
  B = INT(A)
  
  IF(PRESENT(n)) THEN
    IF(n>0) THEN
      round = B+REAL(INT(10**n*(A-B)+half), dpp)/10**n
    ELSE
      round = B+REAL(INT((A-B)+half), dpp)
    ENDIF
  ELSE
    round = B+REAL(INT((A-B)+half), dpp)
  ENDIF
  
END FUNCTION round

!!**************************************************************************!!  
!!******************       CHARACTER MANIPULATION       ********************!! 
!!**************************************************************************!!  

SUBROUTINE upcase1(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string 
  INTEGER                     :: i

  DO i=1,LEN(string)
    IF ( string(i:i) >= 'a' .AND. string(i:i) <= 'z' ) &
      string(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
  END DO
  
  RETURN

END SUBROUTINE upcase1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE upcaseN(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string(:) 
  INTEGER                     :: i

  DO i=1,SIZE(string)
    CALL upcase1(string(i))
  ENDDO
  
  RETURN

END SUBROUTINE upcaseN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION upcase1f(string) RESULT(char)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)   :: string 
  CHARACTER(LEN=LEN(string)) :: char 
  INTEGER :: i
  
  char = string

  DO i=1,LEN(string)
    IF ( string(i:i) >= 'a' .AND. string(i:i) <= 'z' ) &
      char(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
  ENDDO
  
  RETURN

END FUNCTION upcase1f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION upcaseNf(string) RESULT(char)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)   :: string(:) 
  CHARACTER(LEN=LEN(string)) :: char(SIZE(string)) 
  INTEGER                    :: i
  
  DO i=1,SIZE(string)
    char(i) = upcase1f(string(i))
  ENDDO
  
  RETURN

END FUNCTION upcaseNf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE lowcase1(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string
  INTEGER                     :: i

  DO i=1,LEN(string)
    IF ( string(i:i) >= 'A' .AND. string(i:i) <= 'Z' ) &
      string(i:i) = ACHAR ( IACHAR ( string(i:i) ) + 32 )
  ENDDO

  RETURN

END SUBROUTINE lowcase1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE lowcaseN(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string(:)
  INTEGER                     :: i

  DO i=1,LEN(string)
    CALL lowcase1(string(i))
  ENDDO

  RETURN

END SUBROUTINE lowcaseN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION lowcase1f(string) RESULT(char)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)   :: string
  CHARACTER(LEN=LEN(string)) :: char 
  INTEGER                    :: i

  char = string

  DO i=1,LEN(string)
    IF ( string(i:i) >= 'A' .AND. string(i:i) <= 'Z' ) &
      char(i:i) = ACHAR ( IACHAR ( string(i:i) ) + 32 )
  ENDDO

  RETURN

END FUNCTION lowcase1f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION lowcaseNf(string) RESULT(char)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)   :: string(:)
  CHARACTER(LEN=LEN(string)) :: char(SIZE(string)) 
  INTEGER                    :: i

  DO i=1,LEN(string)
    char(i) = lowcase1f(string(i))
  ENDDO

  RETURN

END FUNCTION lowcaseNf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE KeepNumbers(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string 
  INTEGER                     :: i, j
  CHARACTER(LEN(string))      :: string2 

  string2 = string
  string = ""
  j = 0
  DO i=1,LEN(string2)
    IF(string2(i:i) >= '0' .AND. string2(i:i) <= '9') THEN
      j = j+1
      string(j:j) = string2(i:i)
    ENDIF
  ENDDO
  
  RETURN

END SUBROUTINE KeepNumbers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE DelAllSubstr1(str, substr, lenstr)
!! Deletes all occurrences of substring 'substr' from string 'str' and
!! shifts characters left to fill holes.
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: str
  CHARACTER(*), INTENT(IN)      :: substr
  INTEGER, INTENT(IN), OPTIONAL :: lenstr
  INTEGER                       :: position, lensubstr, lenstr1
  
  IF(PRESENT(lenstr)) THEN
    lenstr1 = lenstr
  ELSE
    lenstr1 = LEN_TRIM(str)
  ENDIF
  lensubstr = LEN(substr)
  position = -1
  DO
     position=INDEX(str(1:lenstr1),substr)
     IF(position == 0) EXIT
     IF(position == 1) THEN
        str=str(lensubstr+1:)
     ELSE
        str=str(:position-1)//str(position+lensubstr:)
     ENDIF
     lenstr1 = lenstr1 - lensubstr
  ENDDO
  str = str(1:lenstr1)
  
  RETURN

END SUBROUTINE DelAllSubstr1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE DelAllSubstrN(str, substr, lenstr)
!! Deletes all occurrences of substring 'substr' from string 'str' and
!! shifts characters left to fill holes.
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: str(:)
  CHARACTER(*), INTENT(IN)      :: substr
  INTEGER, INTENT(IN), OPTIONAL :: lenstr
  INTEGER                       :: i
  
  DO i=1,SIZE(str)
    CALL DelAllSubstr1(str(i), substr, lenstr)
  ENDDO
  
  RETURN

END SUBROUTINE DelAllSubstrN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReplaceAllSubstr1(str, substr, newsubstr, lenstr)
!! Replaces all occurrences of substring 'substr' from string 'str' with
!! 'newsubstr', but only if 'substr' and 'newsubstr' have same length
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: str
  CHARACTER(*), INTENT(IN)      :: substr, newsubstr
  INTEGER, INTENT(IN), OPTIONAL :: lenstr
  INTEGER                       :: pos, lensubstr, lennewsubstr, lenstr1, k
  
  IF(LEN(substr)<1) THEN
    CALL PrntE("Cannot replace substring. Substring is empty.")
    RETURN
  ENDIF

  !! Return if no job needs to be done
  IF(INDEX(str, substr)==0) RETURN
  IF(substr==newsubstr .AND. LEN(substr)==LEN(newsubstr)) RETURN
  
  IF(PRESENT(lenstr)) THEN
    lenstr1 = lenstr
  ELSE
    lenstr1 = LEN_TRIM(str)
  ENDIF

  !! If newsubstr is empty then what is actually needed is removal
  IF(LEN(newsubstr)==0) THEN
    CALL DelAllSubstr(str, substr, lenstr1)
    RETURN
  ENDIF

  lensubstr = LEN(substr)
  lennewsubstr = LEN(newsubstr)
  IF(LEN_TRIM(substr)>0) lensubstr = LEN_TRIM(substr)
  IF(LEN_TRIM(newsubstr)>0) lennewsubstr = LEN_TRIM(newsubstr)

  k = 1
  DO
    pos = INDEX(str(k:lenstr1), substr)
    IF(pos == 0) EXIT
    IF(pos == 1) THEN
     IF(k==1) str = newsubstr(1:lennewsubstr)//str(k+lensubstr:) 
     IF(k>1)  str = str(:k-1)//newsubstr(1:lennewsubstr)//str(k+lensubstr:)
    ENDIF 
    IF(pos > 1 .AND. pos<lenstr1-lensubstr) &
     str = str(:k+pos-2)//newsubstr(1:lennewsubstr)//str(k-1+pos+lensubstr:) 
    IF(pos == lenstr1-lensubstr) str = str(:k+pos-2)//newsubstr(1:lennewsubstr)
    k = k-1+pos+lennewsubstr+1 
    lenstr1 = LEN_TRIM(str)
  ENDDO
  
  RETURN

END SUBROUTINE ReplaceAllSubstr1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReplaceAllSubstrN(str, substr, newsubstr, lenstr)
!! Replaces all occurrences of substring 'substr' from string 'str' with
!! 'newsubstr', but only if 'substr' and 'newsubstr' have same length
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)   :: str(:)
  CHARACTER(*), INTENT(IN)      :: substr, newsubstr
  INTEGER, INTENT(IN), OPTIONAL :: lenstr
  INTEGER                       :: i
  
  DO i=1,SIZE(str)
    CALL ReplaceAllSubstr1(str(i), substr, newsubstr, lenstr)
  ENDDO
  
  RETURN

END SUBROUTINE ReplaceAllSubstrN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE ReplaceAllSubstr2(str, substr, newsubstr)
!! Replaces all occurrences of substring 'substr' from string 'str' with
!! 'newsubstr', but only if 'substr' and 'newsubstr' have same length
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: str
  CHARACTER(*), INTENT(IN)    :: substr, newsubstr
  INTEGER                     :: position, lensubstr, lenstr
  
  IF(substr==newsubstr .OR. INDEX(str, substr)==0) RETURN

  IF(LEN(substr)/=LEN(newsubstr)) THEN
    CALL PrntE("Cannot replace substring. Strings have different size.")
    RETURN
  ENDIF
  IF(LEN(substr)<1) THEN
    CALL PrntE("Cannon replace substring. Substring is empty.")
    RETURN
  ENDIF

  lenstr = LEN_TRIM(str)
  lensubstr = LEN(substr)
  position = -1
  DO
     position = INDEX(str(1:lenstr),substr)
     IF(position == 0) EXIT
     IF(position > 0) THEN
       str(position:position+lensubstr-1) = newsubstr(1:lensubstr)
     ENDIF
  ENDDO
  
  RETURN

END SUBROUTINE ReplaceAllSubstr2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SplitString1(string, split, words, nwords, iostat)
!! A wrapper around SplitStringN (see below) which allows scalar input 'split'
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)               :: string, split
  CHARACTER(*), ALLOCATABLE, INTENT(OUT) :: words(:)
  INTEGER, INTENT(OUT), OPTIONAL         :: nwords, iostat
  CHARACTER(LEN(split))                  :: split1(1)
  
  IF(PRESENT(iostat)) iostat = 0
  split1 = split
  CALL SplitStringN(string, split1, words, nwords, iostat)
  
END SUBROUTINE SplitString1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SplitStringN(string, split, words, nwords, iostat)
!! Splits the strint in 'string' with the values in 'split' giving the 
!! separators. Returns the split output in words and its count in nwords 
  IMPLICIT NONE
  INTEGER, PARAMETER                     :: chunk = 50
  CHARACTER(*), INTENT(IN)               :: string
  CHARACTER(*), INTENT(IN)               :: split(:)
  CHARACTER(*), ALLOCATABLE, INTENT(OUT) :: words(:)
  INTEGER, INTENT(OUT), OPTIONAL         :: nwords, iostat
  INTEGER                                :: i, k, lenst, lensp, words_size, &
                                            firstsplit, prevsplitend, nwords1
  CHARACTER(LEN(string))                 :: string0
  
  IF(PRESENT(iostat)) iostat = 0
  
  !! Store the string locally so that changes are possible
  string0 = string
  
  !! Remember lengths
  lenst = LEN_TRIM(string0)
  lensp = LEN(split(1))
  
  !! If multiple splits given in 'split', replace 2nd, 3rd, etc. with 1st
  DO i=2,SIZE(split)
    CALL ReplaceAllSubstr(string0, split(i), split(1), lenst)
  ENDDO
  
  !! Check for too long split 
  !! THIS PROCEDURE NEEDS TO BE MODIFIED SO THAT THIS LIMITATION IS REMOVED!
  IF(LEN(split(1))>1) &
    CALL PrntE("String split(s) must be of length 1 (SplitStringN).", Q=.TRUE.)
                    
  !! Prepare words for first word
  words_size = 1 
  CALL ResizeVar(words, words_size, "")
  nwords1 = 0
  
  !! Check for lengths of string and split
  IF(lenst==0 .OR. lensp==0 .OR. lenst<lensp) THEN
    IF(PRESENT(iostat)) THEN
      iostat = -1
      IF(lenst==0)    iostat = -2
      IF(lensp==0)    iostat = -3
      IF(lenst<lensp) iostat = -4
    ENDIF 
    IF(PRESENT(nwords)) nwords = nwords1
    RETURN
  ENDIF
  
  !! Find first split and assign the first word
  firstsplit = INDEX(string0, split(1))
  IF(firstsplit==0) firstsplit = lenst+1
  IF(firstsplit>1) words(1) = string0(:firstsplit-1) 
  nwords1 = 1
  !! Find all other splits and assign the words
  prevsplitend = firstsplit+lensp-1
  k = prevsplitend
  DO WHILE (k+lensp<=lenst+1)
    IF(string0(k+1:MIN(k+lensp,lenst))/=split(1) .AND. k+lensp<=lenst) THEN
      k = k+1
      CYCLE
    ENDIF  
    nwords1 = nwords1 + 1
    !! If words is too small, enlarge it
    IF(nwords1 > words_size) THEN
      words_size = words_size + chunk
      CALL ResizeVar(words, words_size, "")
    ENDIF
    !! Assign the last word
    IF(k+lensp>=lenst) & 
      words(nwords1) = string0(prevsplitend+1:)
    !! Assign the next word
    IF(prevsplitend+1<=k .AND. k<=lenst) & 
      words(nwords1) = string0(prevsplitend+1:k)
    !! Remember where the current split ends
    prevsplitend = k+lensp
    k = prevsplitend  
  ENDDO
  
  !! Remove extra space for words
  CALL ResizeVar(words, nwords1)
  
  IF(PRESENT(nwords)) nwords = nwords1

  RETURN
    
END SUBROUTINE SplitStringN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION FindLastSubstr(string,substr) RESULT(pos)
  !! Returns the last occurence of substr within string
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: string, substr
  INTEGER :: pos, i
  
  IF(LEN_TRIM(string)>0 .AND. LEN_TRIM(substr)>0) THEN
    DO i=LEN_TRIM(string),1,-1
      IF(string(i:MIN(LEN_TRIM(string),i+LEN_TRIM(substr)-1)) == substr) THEN
        pos = i
        RETURN
      ENDIF
    ENDDO
  ENDIF 

  !! If advanced here then there is no substring of string matching substr  
  pos = 0
  
  RETURN
  
END FUNCTION FindLastSubstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION CountSubstr(string, substr) RESULT(count)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: string, substr
  CHARACTER(LEN(string))   :: str
  INTEGER                  :: count, pos

  count = 0
  IF(LEN_TRIM(string)==0 .OR. LEN_TRIM(substr)==0) RETURN

  !DO i=1,LEN_TRIM(string)
  !  IF(string(i:MIN(LEN_TRIM(string),i+LEN_TRIM(substr)-1)) == substr) THEN
  !    count = count + 1
  !  ENDIF
  !ENDDO

  str = string
  DO
    !! Look for substr
    pos = INDEX(str, substr)
    !! If none found exit
    IF(pos == 0) EXIT
    !! Otherwise increase counter
    count = count + 1
    !! Remove the part until the end of first substr from str
    IF(pos+1 > LEN(str)) EXIT
    str = str(pos+1:)
  ENDDO

  RETURN
  
END FUNCTION CountSubstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE AttachString(string, attach, ntimes, between)
!! Attaches 'attach' at the end of 'string' ntimes and if 'between' present
!! it is put once between 'string' and 'attach'
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT)        :: string
  CHARACTER(*), INTENT(IN)           :: attach
  INTEGER, INTENT(IN), OPTIONAL      :: ntimes
  CHARACTER(*), INTENT(IN), OPTIONAL :: between
  INTEGER                            :: ntimes1, k
  
  !! How many times 
  ntimes1 = 1
  IF(PRESENT(ntimes)) ntimes1 = ntimes
  
  !! Return if ntimes1 non-positive
  IF(ntimes1 <= 0) RETURN

  !! Attach the first 'attach' (and between)
  IF(ntimes1>0) THEN
    IF(PRESENT(between)) THEN
      string = TRIM(string) // between // TRIM(attach)
    ELSE
      string = TRIM(string) // TRIM(attach)
    ENDIF
  ENDIF

  !! Attach 'attach' the remaining ntimes1 - 1 times
  DO k = 2, ntimes1
    string = TRIM(string)//TRIM(attach)
  ENDDO
  
  RETURN
  
END SUBROUTINE AttachString
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE AddSpacesCommas(string, commas, spaces)
!! Ads spaces and/or commas to an integer part of a string (usually a number)
  IMPLICIT NONE
  INTEGER, PARAMETER            :: minlength = 4, groupsize = 3
  CHARACTER(*), INTENT(INOUT)   :: string
  LOGICAL, INTENT(IN), OPTIONAL :: commas, spaces
  LOGICAL                       :: add_commas1, add_spaces1, add
  CHARACTER(mltl)               :: temp
  INTEGER                       :: i, k, start

  !! Set variables based on optional input (spaces supersede commas)
  add_spaces1 = add_spaces
  add_commas1 = add_commas
  IF(PRESENT(commas)) add_commas1 = commas
  IF(PRESENT(spaces)) add_spaces1 = spaces
  add = .FALSE.
  IF(add_commas1 .OR. add_spaces1) add = .TRUE.

  !! If char too short, then don't add commas or spaces
  IF(LEN_TRIM(string) < minlength .OR. .NOT.add) RETURN

  !! Otherwise add commas and/or spaces
  temp = ""

  !! If string is a real number, modify only the integer part 
  start = INDEX(string, ".")
  IF(start == 0) THEN
    start = LEN_TRIM(string)
  ELSE
    temp = string(start:)
    start = start - 1
  ENDIF
  
  !! Do the modification
  i=0
  DO k=start,1,-1
    temp = string(k:k)//temp
    i=i+1
    !! Every third character if not the first add commas and spaces (or not)
    IF(i==groupsize .AND. k>1) THEN
      i=0
      IF(add_spaces1) temp = " "//TRIM(temp)
      IF(add_commas1) temp = ","//TRIM(temp)
    ENDIF 
  ENDDO
  
  !! Change the value of string
  string = TRIM(temp)
  
  RETURN

END SUBROUTINE AddSpacesCommas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION AC(str) RESULT(char)
!! Adds commas to str which contains an integer or a real number
!! The variable length char is commented because ifort doesn't want to take it
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: str
  !CHARACTER(MAX(LEN(str)+INT((LEN(ni2cp(ABS(INT(c2rp(str)))))-1)/3),1)) :: char
  CHARACTER(mstl) :: char
  
  !! Save str into char
  char = str
  
  !! Add spaces and/or commas
  CALL AddSpacesCommas(char, commas=.TRUE.)
  
  RETURN

END FUNCTION AC
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

FUNCTION AS(str) RESULT(char)
!! Adds commas to str which contains an integer or a real number
!! The variable length char is commented because ifort doesn't want to take it
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: str
  !CHARACTER(MAX(LEN(str)+INT((LEN(ni2cp(ABS(INT(c2rp(str)))))-1)/3),1)) :: char
  CHARACTER(mstl) :: char
  
  !! Save str into char
  char = str
  
  !! Add spaces and/or commas
  CALL AddSpacesCommas(char, spaces=.TRUE.)
  
  RETURN

END FUNCTION AS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE CheckSeparator(sep, nsep, line, findsep) 
!! Performs a check for the proper number of 'sep' in line and if the expected
!! number 'nsep' is not found, it tries using a different separator
  IMPLICIT NONE
  CHARACTER, INTENT(INOUT)      :: sep
  INTEGER, INTENT(INOUT)        :: nsep
  CHARACTER(*), INTENT(IN)      :: line
  LOGICAL, INTENT(IN), OPTIONAL :: findsep
  INTEGER                       :: ansep
  
  !! Check for proper count of sep in line
  ansep = CountSubstr(TRIM(line), sep)
  IF(ansep == nsep) RETURN
  
  !! Find new separator
  IF(.NOT.PRESENT(findsep)) GOTO 100
  IF(.NOT.findsep) GOTO 100
  
  !! Try tab
  ansep = CountSubstr(TRIM(line), tab)  
  IF(ansep == nsep) THEN
    sep = tab
    GOTO 100
  ENDIF
  
  !! Try space
  ansep = CountSubstr(TRIM(line), space)  
  IF(ansep == nsep) THEN
    sep = space
    GOTO 100
  ENDIF
  
  !! Unknown separator
  nsep = -1

  RETURN
  
  100 CONTINUE

  nsep = ansep 
  RETURN
      
END SUBROUTINE CheckSeparator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!**************************************************************************!!  
!!**************            PRINTING ROUTINES           ********************!! 
!!**************************************************************************!!  

SUBROUTINE PrintProgramHeader(unit, file, commented)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL      :: unit
  CHARACTER(*), INTENT(IN), OPTIONAL :: file
  LOGICAL, INTENT(IN), OPTIONAL      :: commented

  CHARACTER(60), SAVE                :: compiler_name, compiler_version, &
                                        compiler_build
#ifdef __DATE__
  CHARACTER(11), SAVE                :: compile_date = __DATE__
#else
  CHARACTER(11), SAVE                :: compile_date = ""
#endif

  CHARACTER(head_width)              :: top, bot, spc, prg, cl, cl1, cl2
  CHARACTER                          :: cc, ssid, ltop, lbot, rtop, rbot
  INTEGER                            :: i, ind

  cl1 = ""
  IF(LEN_TRIM(compile_date)>0) cl1 = "Compiled on "//TRIM(compile_date)
  cl2 = "Compiled under "//TRIM(sys_type)//" "//TRIM(sys_bit)//" "//TRIM(sys_envir)

#ifdef __INTEL_COMPILER
  compiler_name = "Intel Fortran Compiler"
  compiler_version = ""
  compiler_build = i2cp(__INTEL_COMPILER_BUILD_DATE)

  !! Check for single digit days and adjust by adding 0
  IF(compile_date(5:5)==" ") compile_date(5:5) = "0"
  !! Compiler version and build
  cl1 = "Compiled on "//TRIM(compile_date)//" with "//TRIM(compiler_name)//&
           " "//TRIM(compiler_version)//" (build "//TRIM(compiler_build)//")"
#else
#ifdef __GFORTRAN__
  compiler_name = "GNU Fortran Compiler"
  compiler_version = __VERSION__

  !! Check for single digit days and adjust by adding 0
  IF(compile_date(5:5)==" ") compile_date(5:5) = "0"
  !! Compiler version and build
  compiler_build = compiler_version(7:15)
  compiler_version = compiler_version(1:6)
  cl1 = "Compiled on "//TRIM(compile_date)//" with "//TRIM(compiler_name)//&
           " "//TRIM(compiler_version)
  IF(LEN_TRIM(compiler_build)>0) &
    cl1 = TRIM(cl1)//" (build "//TRIM(compiler_build)//")"
#endif
#endif

  ltop = "="        ! left top corner
  lbot = "="        ! left bottom corner 
  rtop = "="        ! right top corner 
  rbot = "="        ! right bottom corner
  ssid = "+"        ! sides

  cc = ""
  IF(PRESENT(commented)) cc = comment

  !! Print an empty line above the box
  CALL Prnt0(cc, unit=unit, file=file)

  !! Define the line with mostly empty space and the top bottom bottom lines
  spc = cc; spc(2:2) = ssid; spc(LEN(spc):) = ssid
  top = cc; top(2:2) = ltop; top(LEN(top):) = rtop
  bot = cc; bot(2:2) = lbot; bot(LEN(bot):) = rbot
  DO i=3,LEN(top)-1; top(i:i) = "="; bot(i:i) = "="; ENDDO

  !! Define the program full name line
  prg = program_fullname
  ind = MAX(INT((LEN(prg)-LEN_TRIM(prg))/2),2) + 1
  prg(ind:) = TRIM(prg)
  prg(:ind-1) = ""
  prg(1:1) = cc
  prg(2:2) = ssid
  prg(LEN(prg):) = ssid

  !! Print the header top line, i.e. boundary made out of '='
  CALL Prnt0(top, unit=unit, file=file)

  !! Print the "spacious" line
  CALL Prnt0(spc, unit=unit, file=file)

  CALL Prnt0(prg, unit=unit, file=file)

  !! Print the 2 compilation info lines
  DO i=1,2
    IF(i==1) cl = cl1
    IF(i==2) cl = cl2
    IF(LEN_TRIM(cl)>0) THEN
      ind = MAX(INT((LEN(cl)-LEN_TRIM(cl))/2),2) + 1
      cl(ind:) = TRIM(cl)
      cl(:ind-1) = ""
      cl(1:1) = cc
      cl(2:2) = ssid
      cl(LEN(cl):) = ssid
      CALL Prnt0(cl, unit=unit, file=file)
    ENDIF
  ENDDO
  
  !! Print the "spacious" line
  IF(LEN_TRIM(cl1)>0 .OR. LEN_TRIM(cl2)>0) CALL Prnt0(spc, unit=unit, file=file)
  
  !! Print the header bottom line, i.e. boundary made out of '='
  CALL Prnt0(bot, unit=unit, file=file)
  
  !! Print an empty line below the box
  CALL Prnt0(cc, unit=unit, file=file)

  RETURN

END SUBROUTINE PrintProgramHeader

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrntE(text, Q, skip1, skip2, nolead, lead, advance, trimit, &
                   log, screen, warning, premature, ask_to_continue, dellog)
!! Print error message
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)            :: text
  INTEGER, INTENT(IN), OPTIONAL       :: skip1, skip2, lead
  CHARACTER(*), INTENT(IN), OPTIONAL  :: advance
  LOGICAL, INTENT(IN), OPTIONAL       :: Q, nolead, log, screen, trimit, &
                                         warning, premature, ask_to_continue, &
                                         dellog
  INTEGER                             :: skip1a, skip2a
  CHARACTER                           :: response

  skip1a = 1
  skip2a = 1
  IF(PRESENT(skip1))     skip1a = skip1
  IF(PRESENT(skip2))     skip2a = skip2

  !! Print error message
  CALL Prnt(text="ERROR: "//text, Q=Q, skip1=skip1a, skip2=skip2a, &
                 nolead=nolead, lead=lead, advance=advance, trimit=trimit, &
                 log=log, screen=screen, warning=warning, error=.TRUE., &
                 premature=premature, dellog=dellog)
                 
  !! If Q was false or missing and ask_to_continue is present and true, ask user
  !! if the program should continue running
  IF(PRESENT(ask_to_continue)) THEN
    IF(ask_to_continue) THEN
      CALL AskUser(response, terminate=.TRUE.)
      IF(response/="N") CALL Terminate()
    ENDIF
  ENDIF

  RETURN 

END SUBROUTINE PrntE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PrntW(text, Q, skip1, skip2, nolead, lead, advance, trimit,&
                        log, screen, premature, dellog)
!! Print warning message
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)            :: text
  CHARACTER(*), INTENT(IN), OPTIONAL  :: advance
  INTEGER, INTENT(IN), OPTIONAL       :: skip1, skip2, lead
  LOGICAL, INTENT(IN), OPTIONAL       :: Q, nolead, log, screen, trimit, &
                                         premature, dellog

  CALL Prnt(text="WARNING: "//text, Q=Q, skip1=skip1, skip2=skip2, &
            nolead=nolead, lead=lead, advance=advance, trimit=trimit, &
            log=log, screen=screen, warning=.TRUE., premature=premature, &
            dellog=dellog)
  RETURN 

END SUBROUTINE PrntW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Prnt0(text, text1, text2, Q, skip1, skip2, lead, advance, &
                      trimit, log, screen, unit, file, warning, error, &
                      premature, wait, l, dots, dellog)
!! Print text without the lead text (' >')
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN), OPTIONAL  :: text, text1, text2, advance, file
  INTEGER, INTENT(IN), OPTIONAL       :: skip1, skip2, lead, unit, l
  LOGICAL, INTENT(IN), OPTIONAL       :: Q, log, screen, trimit, warning, &
                                         error, premature, wait, dots, dellog

  !! Call Prnt with nolead = .TRUE.
  CALL Prnt(text=text, text1=text1, text2=text2, Q=Q, skip1=skip1, &
            skip2=skip2, nolead=.TRUE., lead=lead, advance=advance, &
            trimit=trimit, log=log, screen=screen, unit=unit, file=file, &
            warning=warning, error=error, premature=premature, l=l, &
            wait=wait, dots=dots, dellog=dellog)

END SUBROUTINE Prnt0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PrntF(text, text1, text2, Q, skip1, skip2, lead, advance, &
                      trimit, log, screen, unit, file, warning, error, &
                      premature, wait, l, dots, dellog)
!! Print text with a leading star
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN), OPTIONAL  :: text, text1, text2, advance, file
  INTEGER, INTENT(IN), OPTIONAL       :: skip1, skip2, lead, unit, l
  LOGICAL, INTENT(IN), OPTIONAL       :: Q, log, screen, trimit, warning, &
                                         error, premature, wait, dots, dellog

  !! Call Prnt with nolead = .TRUE.
  CALL Prnt(text=" * "//text, text1=text1, text2=text2, Q=Q, skip1=skip1, &
            skip2=skip2, nolead=.TRUE., lead=lead, advance=advance, &
            trimit=trimit, log=log, screen=screen, unit=unit, file=file, &
            warning=warning, error=error, premature=premature, l=l, &
            wait=wait, dots=dots, dellog=dellog)

END SUBROUTINE PrntF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RECURSIVE SUBROUTINE Prnt(text, text1, text2, Q, skip1, skip2, nolead, lead, &
                     advance, trimit, log, screen, unit, file, warning, error, &
                     premature, wait, l, dots, dellog)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN), OPTIONAL  :: text, text1, text2, advance, file
  INTEGER, INTENT(IN), OPTIONAL       :: skip1, skip2, lead, unit, l
  LOGICAL, INTENT(IN), OPTIONAL       :: Q, nolead, log, screen, trimit, &
                                         warning, error, premature, wait, &
                                         dots, dellog
  INTEGER                             :: units(3), i, j, iscreen, ilog, & 
                                         iextra, skip1a, skip2a, lead1, &
                                         length0, length1, length2, &
                                         lempt, lempt1, lempt2
  LOGICAL                             :: printunit(3), Q1, nolead1, trimtext, &
                                         premature1, wait1, file_exists
  CHARACTER(mel)                      :: emptysp, text00, text0
  CHARACTER(3)                        :: advance1, adv
  CHARACTER(mfl)                      :: filename(3)
  CHARACTER                           :: stat(3)
  
  emptysp = ""
  filename = ""

  !! Decide whether dots instead of spaces are to be printed between text1, text
  !! and text2 (leave an empty space at the beginning of emptysp)
  IF(PRESENT(dots)) THEN
    IF(dots) CALL AttachString(emptysp, ".", LEN(emptysp), " ")
  ENDIF

  !! Store the units
  iscreen = 1
  ilog = 2
  iextra = 3
  
  units(iscreen) = usto
  units(ilog) = ulog
  units(iextra) = ufre
  IF(PRESENT(unit)) units(iextra) = unit

  !! Determine whether to PRINT on screen (1) and log file (2)
  printunit(iscreen) = .NOT.silent 
  printunit(ilog) = .NOT.no_log
  printunit(iextra) = .FALSE.
  filename(ilog) = log_file
  stat(ilog) = 'O'
  
  !! Change where printing happens based on user input
  IF(PRESENT(file)) THEN
    printunit(iscreen) = .FALSE.
    printunit(ilog) = .FALSE.
    printunit(iextra) = .TRUE.
    filename(iextra) = file
    stat(iextra) = 'U'
  ENDIF
  IF(PRESENT(screen) .AND. .NOT.silent) printunit(iscreen) = screen 
  IF(PRESENT(log) .AND. .NOT.no_log)    printunit(ilog) = log
  
   !! Create log file
  IF(printunit(ilog)) THEN
    IF(LEN_TRIM(log_file)==0) THEN
      CALL CreateLogFile(log_file)
      filename(ilog) = log_file
    ELSE
      INQUIRE(FILE=log_file, EXIST=file_exists)
      IF(.NOT.file_exists) CALL CreateLogFile(log_file)
    ENDIF
  ENDIF
  
  !! Assign initial values
  text00 = "" 
  Q1 = .FALSE.
  skip1a = 0
  skip2a = 0
  lead1 = 0
  nolead1 = .NOT..FALSE.
  trimtext = .TRUE.
  premature1 = .TRUE.
  wait1 = .FALSE.
  advance1 = 'YES'
  
  !! Assign values based on optional variables    
  IF(PRESENT(text))      text00 = text
  IF(PRESENT(Q))         Q1 = Q
  IF(PRESENT(skip1))     skip1a = skip1
  IF(PRESENT(skip2))     skip2a = skip2
  IF(PRESENT(lead))      lead1 = lead
  IF(PRESENT(nolead))    nolead1 = nolead
  IF(PRESENT(trimit))    trimtext = trimit
  IF(PRESENT(premature)) premature1 = premature
  IF(PRESENT(wait))      wait1 = wait
  IF(PRESENT(advance)) THEN
    IF(LEN_TRIM(advance)>0) THEN
      IF(upcasef(advance(1:1))=='N') advance1='NO'
    ENDIF
  ENDIF
  
  length0 = 0
  text0 = ""

  !! Increase warning counter and check for too many warnings
  IF(PRESENT(warning)) THEN
    IF(warning) THEN
      nwarnings = nwarnings + 1
      !! Store the warning
      CALL ResizeVar(warnings, nwarnings, TRIM(text00))
      IF(nwarnings > nwarnings_limit) printunit(iscreen) = .FALSE.
    ENDIF
  ENDIF
  
  !! Increase warning counter
  IF(PRESENT(error)) THEN
    IF(error) THEN
      nerrors = nerrors + 1
      CALL ResizeVar(errors, nerrors, TRIM(text00))
      IF(nerrors > nerrors_limit) printunit(iscreen) = .FALSE.
    ENDIF
  ENDIF

  !! If text1 or text2 give, no trimming and no lead text
  IF(PRESENT(text1) .OR. PRESENT(text2)) THEN
    nolead1 = .TRUE.
    trimtext = .FALSE.
  ENDIF

  !! Attach leadtext if nolead1 is .false.
  IF(.NOT.nolead1) THEN
    text0 = lt
    length0 = length0 + LEN(lt)
    IF(lead1==0 .AND. .NOT.PRESENT(lead)) lead1 = LEN(lt)
  ENDIF

  !! Finally assign text
  IF(PRESENT(text)) THEN
    !! Skip printing leading characters if nolead present and true 
    IF(nolead1) THEN
      text0 = text
    ELSE
      text0 = lt//text
    ENDIF

    !! Based on trimip1 attach text (only if text present)
    IF(trimtext) THEN 
      length0 = length0 + LEN_TRIM(text)
      IF(LEN(text)==LEN_TRIM(text)+1) length0 = length0 + 1
    ELSE
      length0 = length0 + LEN(text)
    ENDIF
  ENDIF

  !! Get lengths and positions
  length1 = 0
  length2 = 0
  IF(PRESENT(text1)) length1 = LEN_TRIM(text1)
  IF(PRESENT(text2)) length2 = LEN_TRIM(text2)
  !IF(PRESENT(text1)) length1 = LEN(text1)
  !IF(PRESENT(text2)) length2 = LEN(text2)

  lempt1 = 0
  lempt2 = 0
  IF(PRESENT(l)) THEN
    lempt = MAX(0,l - length1 - length0 - length2)
    lempt1 = MIN(LEN(emptysp), round_int(REAL(lempt, dpp)/2))
    lempt2 = MIN(LEN(emptysp), lempt - lempt1)
  ENDIF
  
!$OMP CRITICAL
  !! Do cycle over screen and log file
  DO i=1,SIZE(printunit)

    !! Skip current unit
    IF(.NOT.printunit(i)) CYCLE

    !! Open log file if currently writing to it
    !! If failure, don't announce and just quit
    IF(i==iscreen) THEN
      adv = advance1
    ELSE
      CALL OpenFile(UNIT=units(i), FILE=filename(i), ACTION='W', &
                    STATUS=stat(i), FORM='F', POSITION="A", ACCESS='S', &
                    announce=.FALSE., chck_miss=(i==ilog), islog=(i==ilog))
      adv = 'YES'
    ENDIF

    !! Skip skip1a lines
    DO j=1,skip1a; WRITE(units(i),'(A)'); ENDDO

    !! Print text1 if given
    IF(PRESENT(text1) .AND. length1>0) &
      CALL PrintLongText(text1(1:length1), unit=units(i), ADVANCE='NO', &
                         trimit=.FALSE., lead=lead1)

    !! Print some empty space
    IF(lempt1>0) &
      CALL PrintLongText(emptysp(1:lempt1), unit=units(i), ADVANCE='NO', &
                         trimit=.FALSE., lead=lead1)
                         
    !! PRINT THE MAIN TEXT
    IF(length0>0) &
      CALL PrintLongText(text0(1:length0), unit=units(i), ADVANCE='NO', &
                         trimit=.FALSE., lead=lead1)

    !! Print some empty space
    IF(lempt2>0) THEN
      IF(lempt2>1) &
        CALL PrintLongText(emptysp(2:lempt2), unit=units(i), ADVANCE='NO', &
                           trimit=.FALSE., lead=lead1)
      CALL PrintLongText(emptysp(1:1), unit=units(i), ADVANCE='NO', &
                         trimit=.FALSE., lead=lead1)
    ENDIF

    !! Print text2 if given
    IF(PRESENT(text2) .AND. length2>0) &
      CALL PrintLongText(text2(1:length2), unit=units(i), ADVANCE='NO', &
                         trimit=.FALSE., lead=lead1)

    WRITE(units(i),'(A)', ADVANCE=adv)
    
    !! Skip skip2a lines
    DO j=1,skip2a; WRITE(units(i),'(A)'); ENDDO
    
    !! Flush the output
    CALL FlushOutput(units(i))
    
    !! Close log file
    IF(i/=iscreen) CLOSE(units(i))

  ENDDO
!$OMP END CRITICAL

  IF(wait1) READ(*,'(A)')  
  IF(Q1) CALL Terminate(ask=.FALSE., premature=premature1, dellog=dellog)

  RETURN
  
END SUBROUTINE Prnt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PrintLongText(text, lead, flead, split, unit, width, advance, trimit)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)        :: text
  CHARACTER, INTENT(IN), OPTIONAL :: split, advance
  INTEGER, INTENT(IN), OPTIONAL   :: unit, lead, flead, width
  LOGICAL, INTENT(IN), OPTIONAL   :: trimit
  
  CHARACTER(1000)                 :: emptysp
  CHARACTER                       :: split1
  CHARACTER(3)                    :: advance1
  INTEGER                         :: lead1, flead1, unit1, width1, &
                                     k, lentext, lastbr, nextsplit,&
                                     splitcounter, lead2, spaceleft !, length
  LOGICAL                         :: trimtext, break

  !! Assign default values
  emptysp = ""
  lead1 = 0
  flead1 = 0
  unit1 = usto
  split1 = space
  trimtext = .TRUE.
  advance1 = 'YES'

  !! Read in optional arguments
  IF(PRESENT(lead))   lead1 = lead
  IF(PRESENT(flead))  flead1 = flead
  IF(PRESENT(unit))   unit1 = unit
  IF(PRESENT(split))  split1 = split
  IF(PRESENT(trimit)) trimtext = trimit
  IF(PRESENT(advance)) THEN
    IF(LEN_TRIM(advance)>0) THEN
      IF(upcasef(advance(1:1))=='N') advance1='NO'
    ENDIF
  ENDIF
  
  !! If width present, use it, otherwise determine width from constants
  IF(PRESENT(width)) THEN
    width1 = width
  ELSE
    IF(unit1 == ulog) THEN
      width1 = log_width
    ELSE
      width1 = stdo_width
    ENDIF
  ENDIF 
  
  !! Check for too long lead empty space
  IF(lead1>=stdo_width) lead1 = 0

  !! Get text length (with trimming or without) 
  IF(trimtext) THEN
    lentext = LEN_TRIM(text)
  ELSE
    lentext = LEN(text)
  ENDIF 
  lead2 = flead1
  splitcounter = 0
  lastbr = 0
  break = .FALSE.
  k=0
  
  !! Do the printing
  DO WHILE(k<lentext)
    
    k = k+1
    
    !! Print the next letter
    WRITE(unit1,'(A)', ADVANCE='NO') text(k:k)
    
    !! If entire text printed, exit
    IF(k==lentext) EXIT

    !! If the current character is split1, allow for breaking
    IF(text(k:k)==split1) THEN
      break = .TRUE.
      splitcounter = splitcounter + 1
    ENDIF

    !! If breaking not allowed and not the end of the line, 
    !! just cycle to print the next letter
    IF(.NOT.break .AND. k-lastbr < width1-lead1) CYCLE
    
    !! Check for exceeding of width1 if the next word is to be printed on the
    !! same line. If so, then print end of line and leading empty spaces
    nextsplit = INDEX(text(k+1:lentext),split1)
    IF(nextsplit==0) nextsplit = lentext - k + 1
    spaceleft = width1 - lead2 - (k - lastbr)
    !! If the next word is longer than one line, print until the end of line
    IF(nextsplit - 1 > width1 - lead1) THEN
      WRITE(unit1,'(A)') text(k+1:MIN(k+spaceleft,lentext))
      WRITE(unit1,'(A)', ADVANCE='NO') emptysp(1:lead1)
      lead2 = lead1
      lastbr = k + spaceleft
      k = lastbr
      IF(k < nextsplit) break = .FALSE.
      CYCLE
    ENDIF
    
    !! If the next word cannot fit on this line, start a new line
    !length = k - lastbr + nextsplit + 1
    spaceleft = width1 - lead2 - (k - lastbr)
    !IF(length > width1 - lead2) THEN
    IF(nextsplit - 1 > spaceleft) THEN
      WRITE(unit1,'(A)') ""
      WRITE(unit1,'(A)', ADVANCE='NO') emptysp(1:lead1)
      lead2 = lead1
      lastbr = k
      CYCLE
    ENDIF

  ENDDO
  
  !! Print end-of-line if advance1 is yes
  IF(advance1=='YES') WRITE(unit1,'(A)') ""
  
  RETURN


END SUBROUTINE PrintLongText

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE EraseTextL(n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER             :: k
  
  IF(n<=0) RETURN
  DO k=1,n
    WRITE(usto, '(A)', ADVANCE='NO') char(8)
    WRITE(usto, '(A)', ADVANCE='NO') " "
    WRITE(usto, '(A)', ADVANCE='NO') char(8)
  ENDDO
  RETURN

END SUBROUTINE EraseTextL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE EraseTextT(string)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: string
  
  CALL EraseTextL(LEN_TRIM(string))

  RETURN

END SUBROUTINE EraseTextT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrintCD1(text, istep, maxnsteps, percent, nerrors, testing, finish, &
                    show_pct, stext, cdstep, adv)
!! Print countdown in percentages
!! text ... used to store what was written on last iteration that needs to be
!!          deleted on the current iteration
!! istep ... contains "size of the update step" in percents
!! maxnsteps ... maximum number of progress refreshes
!! percent ... last value of the progress counter
!! nerrors ... if present the last progress refresh will contain error report
!! testing ... determines whether progress line formulated in terms of tests (if true)
!! finish ... determines whether then the last progress refresh will not show 
!!            text "Finished!" instead of "100%" (if true)
!! stext ... starting progress line text 
!! NOTE: Countdown can be initialized by stext being present and having non-"" 
!! value or by the value in 'percent' being negative
  IMPLICIT NONE
#if defined __INTEL_COMPILER
  INTEGER(ikb), PARAMETER               :: df_cd_jump_pct = 5, &
                                           df_cd_jump_cnt = 1
#else
  INTEGER(ikb), PARAMETER               :: df_cd_jump_pct = 1, &
                                           df_cd_jump_cnt = 1
#endif
  INTEGER, PARAMETER                    :: limit = 10**7
  CHARACTER(*), INTENT(OUT)             :: text
  INTEGER(ikb), INTENT(IN)              :: istep, maxnsteps
  INTEGER(ikb), INTENT(INOUT)           :: percent
  INTEGER(ikb), INTENT(IN), OPTIONAL    :: nerrors
  LOGICAL, INTENT(IN), OPTIONAL         :: testing, finish, show_pct, adv
  CHARACTER(*), INTENT(INOUT), OPTIONAL :: stext
  INTEGER(ikb), INTENT(IN), OPTIONAL    :: cdstep
  INTEGER(ikb)                          :: prcnt, D, i, cd_jump1, df_cd_jump
  CHARACTER(3)                          :: advnc
  LOGICAL                               :: show_pct1, adv1

  !! Quit if no countdown should printed
  IF(.NOT.countdown .OR. istep > maxnsteps) RETURN
  !IF(.NOT.countdown) RETURN
  
  !! Resolve whether percentage or counter will be shown
  show_pct1 = .TRUE.
  adv1 = .TRUE.
  IF(PRESENT(show_pct)) show_pct1 = show_pct
  IF(PRESENT(adv)) adv1 = adv
  
  !! Resolve how big the step between successive updates will be
  IF(show_pct1) THEN
    df_cd_jump = df_cd_jump_pct 
  ELSE
    df_cd_jump = df_cd_jump_cnt
  ENDIF
  cd_jump1 = df_cd_jump 
  IF(PRESENT(cdstep)) cd_jump1 = cdstep
  
  !! Make sure cd_jump1 has sensefull value 
  IF(cd_jump1<=0 .OR. (show_pct1 .AND. cd_jump1>100)) &
    cd_jump1 = df_cd_jump

  !! Save how far along the progress is 
  i = istep
  
  !! If the start of countdown is indicated by stext, initialize counters and print stext
  IF(PRESENT(stext)) THEN
    IF(LEN_TRIM(stext)>0) THEN
      CALL Prnt(TRIM(stext)//" ", ADVANCE='NO')
      stext = ""
      text = "" 
      percent = -1
      i = 1
    ENDIF    
  ENDIF

  !! Another way to indicate start of countdown
  IF(percent == -1) THEN
    IF(PRESENT(stext)) stext = ""
    text = "" 
    i = 1
  ENDIF

  !! Print initial testing progress announcement
  IF(PRESENT(testing)) THEN
    IF(testing .AND. percent<0) &
      CALL Prnt("Running test ", log=.FALSE., ADVANCE='NO')
  ENDIF
  
  !! Compute progress
  IF(show_pct1) THEN
    !! Decide whether decimal will be printed or not
    IF(maxnsteps <= limit) D = 100
    IF(maxnsteps > limit)  D = 1000  
    !! Calculate current progress
    prcnt = MIN(MAX(0, INT(DBLE(i)/MAX(1,INT(maxnsteps))*INT(D))), INT(D))
  ELSE
    D = maxnsteps
    prcnt = i
  ENDIF
  
  !! If no change to report, just return
  IF(percent>=0 .AND. prcnt<D .AND. prcnt-percent<cd_jump1) RETURN

  !! Erase the previous text
  CALL EraseText(LEN_TRIM(text))

  advnc = 'NO'
  text = TRIM(i2c(prcnt))
  
  !! Progress appeared, select text to be printed
  IF(.NOT.show_pct1) THEN
    text = TRIM(text)//" (out of "//TRIM(i2c(maxnsteps))//")"
  ELSEIF(maxnsteps<=limit) THEN
    text = TRIM(text)//" %"
  ELSE
    IF(prcnt==0)                    text = "0.0 %"  
    IF(prcnt>0 .AND. prcnt<10)      text = "0."//text(1:1)//" %"  
    IF(prcnt>=10 .AND. prcnt<100)   text = text(1:1)//"."//text(2:2)//" %"  
    IF(prcnt>=100 .AND. prcnt<1000) text = text(1:2)//"."//text(3:3)//" %"
    IF(prcnt == 1000)               text = "100.0 %"
  ENDIF  

  !! Countdown finished and if 'finished' present, announcement is printed
  IF(prcnt >= D) THEN
    IF(adv1) advnc = 'YES'
    IF(PRESENT(finish)) THEN
      IF(finish) text = "Finished!"
    ENDIF
  ENDIF

  !! Print the announcement
  CALL Prnt0(TRIM(text), ADVANCE=advnc, log=.FALSE.)
  
  !! Return progress counter
  percent = prcnt
  
  !! If the countdown is of tests and all test finished, print announcement
  if(.false.) then
    IF(prcnt == D) THEN
      IF(.NOT.PRESENT(testing)) GOTO 10 
      IF(.NOT.testing)          GOTO 10
      IF(.NOT.PRESENT(finish))  GOTO 10 
      IF(.NOT.finish)           GOTO 10 
      IF(.NOT.PRESENT(nerrors)) THEN
        CALL Prnt("Testing finished!")
        GOTO 10
      ENDIF
      IF(nerrors==0) THEN
        CALL Prnt("Testing finished successfully without errors.")
      ELSE
        CALL Prnt("Testing finished with errors.")
      ENDIF
    ENDIF
  endif
  
  10 CONTINUE
  
  !! Flush the output
  !CALL FlushOutput(usto)

  RETURN
  
END SUBROUTINE PrintCD1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE PrintCD2(text, istep, maxnsteps, percent, nerrors, testing, finish, &
                    show_pct, stext, cdstep, adv)
!! Alias for PrintCD1 which takes regular integer arguments istep and
!! maxnsteps
  IMPLICIT NONE
  CHARACTER(*), INTENT(OUT)             :: text
  INTEGER, INTENT(IN)                   :: istep, maxnsteps
  INTEGER(ikb), INTENT(INOUT)           :: percent
  INTEGER(ikb), INTENT(IN), OPTIONAL    :: nerrors
  LOGICAL, INTENT(IN), OPTIONAL         :: testing, finish, show_pct, adv
  CHARACTER(*), INTENT(INOUT), OPTIONAL :: stext
  INTEGER(ikb), INTENT(IN), OPTIONAL    :: cdstep

  CALL PrintCD1(text, INT(istep, ikb), INT(maxnsteps, ikb), percent, nerrors, &
                testing, finish, show_pct, stext, cdstep, adv)
                           
  RETURN
END SUBROUTINE PrintCD2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

SUBROUTINE FlushOutput(unit)
  IMPLICIT NONE
  !EXTERNAL :: FLUSH
  INTEGER, INTENT(IN), OPTIONAL :: unit
  INTEGER                       :: unit1

  unit1 = usto
  IF(PRESENT(unit)) unit1 = unit

#ifdef __INTEL_COMPILER
  IF(unit1==usto) THEN
    CLOSE(unit1)
    OPEN(UNIT=unit1)
    RETURN
  ENDIF
#endif

  CALL FLUSH(unit1)

  RETURN

END SUBROUTINE FlushOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE WriteVarIntoFile(filename, s, as)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN)           :: filename
  CHARACTER(*), INTENT(IN), OPTIONAL :: s
  CHARACTER(*), INTENT(IN), OPTIONAL :: as(:)
  !INTEGER(iks), INTENT(IN), OPTIONAL  :: as_ik(:)
  !INTEGER, INTENT(IN), OPTIONAL      :: as_i(:)
  !INTEGER(ikb), INTENT(IN), OPTIONAL :: as_ikb(:)
  !REAL(dpp), INTENT(IN), OPTIONAL    :: as_r(:)
  LOGICAL                            :: file_exists
  CHARACTER(10)                      :: position
  INTEGER                            :: k
  
  IF(.NOT.PRESENT(s) .AND. .NOT.PRESENT(as)) RETURN

  !! Check existence of the given file
  INQUIRE(FILE=filename, EXIST=file_exists)
  IF(file_exists) THEN
    position = 'APPEND'
  ELSE
    position = 'REWIND'
  ENDIF
  
  !! If file doesn't exist create new, if file exists append it
  OPEN(UNIT=ufre, FILE=filename, ACTION='WRITE', FORM='FORMATTED', &
       POSITION=position)

  !! Write given string(s) into the file
  IF(PRESENT(s)) WRITE(ufre,'(A)', ADVANCE='NO') TRIM(s)
  IF(PRESENT(as)) THEN
    DO k=1,SIZE(as)
      IF(k>1 .OR. PRESENT(s)) WRITE(ufre,'(A)', ADVANCE='NO') tab 
      WRITE(ufre,'(A)', ADVANCE='NO') TRIM(as(k))
    ENDDO
  ENDIF
  
  !! Print end of line
  WRITE(ufre,'(A)') ""

  !! Close the file
  CLOSE(ufre)

  RETURN
 
END SUBROUTINE WriteVarIntoFile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GetEffectiveNRow(x) RESULT(n)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: x(:,:)
  INTEGER :: n, i
  
  n = SIZE(x,1)
  DO i=SIZE(x,1),1,-1
    IF(ANY(LEN_TRIM(x(n,:))>0)) EXIT
    n = n-1
  ENDDO

  RETURN

END FUNCTION GetEffectiveNRow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GetEffectiveNCol(x) RESULT(n)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: x(:,:)
  INTEGER :: n, i
  
  n = SIZE(x,2)
  DO i=SIZE(x,2),1,-1
    IF(ANY(LEN_TRIM(x(:,n))>0)) EXIT
    n = n-1
  ENDDO

  RETURN

END FUNCTION GetEffectiveNCol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GetNumberFormat(x)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: x
  CHARACTER(LEN(x))        :: GetNumberFormat
  INTEGER                  :: before, after, dec, length
  
  length = LEN_TRIM(x)
  IF(length==0) THEN
    GetNumberFormat="(F0.0)"
    RETURN
  ENDIF

  dec = INDEX(x,".")
  IF(dec==0) THEN
    before = length
    after = 1 
  ELSEIF(dec==length) THEN
    before = length-1
    after = 1
  ELSE
    before = LEN_TRIM(x(1:dec-1))
    after = LEN_TRIM(x(dec+1:length))
  ENDIF
  
  GetNumberFormat = "(F"//TRIM(i2cp(before+after+1))//"."//TRIM(i2cp(after))//")"
  
  RETURN

END FUNCTION GetNumberFormat

!!**************************************************************************!!  
!!**************         MISCELLANEOUS ROUTINES         ********************!! 
!!**************************************************************************!!  

SUBROUTINE SetValReal1(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: val
  REAL(dpp), INTENT(OUT)           :: x1
  REAL(dpp), INTENT(OUT), OPTIONAL :: x2, x3, x4, x5, x6, x7, x8, x9, x10, &
      x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, & 
      x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, &
      x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValReal1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValRealN(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: val
  REAL(dpp), INTENT(OUT)           :: x1(:)
  REAL(dpp), INTENT(OUT), OPTIONAL :: x2(:), x3(:), x4(:), x5(:), x6(:), &
      x7(:), x8(:), x9(:), x10(:), x11(:), x12(:), x13(:), x14(:), x15(:), & 
      x16(:), x17(:), x18(:), x19(:), x20(:), x21(:), x22(:), x23(:), x24(:), &
      x25(:), x26(:), x27(:), x28(:), x29(:), x30(:), x31(:), x32(:), x33(:), &
      x34(:), x35(:), x36(:), x37(:), x38(:), x39(:), x40(:), x41(:), x42(:), &
      x43(:), x44(:), x45(:), x46(:), x47(:), x48(:), x49(:), x50(:)
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValRealN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValRealNM(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  REAL(dpp), INTENT(IN)            :: val
  REAL(dpp), INTENT(OUT)           :: x1(:,:)
  REAL(dpp), INTENT(OUT), OPTIONAL :: x2(:,:), x3(:,:), x4(:,:), x5(:,:), &
      x6(:,:), x7(:,:), x8(:,:), x9(:,:), x10(:,:), x11(:,:), x12(:,:), & 
      x13(:,:), x14(:,:), x15(:,:), x16(:,:), x17(:,:), x18(:,:), x19(:,:), &
      x20(:,:), x21(:,:), x22(:,:), x23(:,:), x24(:,:), x25(:,:), x26(:,:), &
      x27(:,:), x28(:,:), x29(:,:), x30(:,:), x31(:,:), x32(:,:), x33(:,:), &
      x34(:,:), x35(:,:), x36(:,:), x37(:,:), x38(:,:), x39(:,:), x40(:,:), &
      x41(:,:), x42(:,:), x43(:,:), x44(:,:), x45(:,:), x46(:,:), x47(:,:), &
      x48(:,:), x49(:,:), x50(:,:)
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValRealNM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValInt1(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: val
  INTEGER, INTENT(OUT)           :: x1
  INTEGER, INTENT(OUT), OPTIONAL :: x2, x3, x4, x5, x6, x7, x8, x9, x10, &
      x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, & 
      x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, &
      x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValInt1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValIntN(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  INTEGER, INTENT(IN)            :: val
  INTEGER, INTENT(OUT)           :: x1(:)
  INTEGER, INTENT(OUT), OPTIONAL :: x2(:), x3(:), x4(:), x5(:), x6(:), &
      x7(:), x8(:), x9(:), x10(:), x11(:), x12(:), x13(:), x14(:), x15(:), & 
      x16(:), x17(:), x18(:), x19(:), x20(:), x21(:), x22(:), x23(:), x24(:), &
      x25(:), x26(:), x27(:), x28(:), x29(:), x30(:), x31(:), x32(:), x33(:), &
      x34(:), x35(:), x36(:), x37(:), x38(:), x39(:), x40(:), x41(:), x42(:), &
      x43(:), x44(:), x45(:), x46(:), x47(:), x48(:), x49(:), x50(:)
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValIntN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValBool1(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  LOGICAL, INTENT(IN)            :: val
  LOGICAL, INTENT(OUT)           :: x1
  LOGICAL, INTENT(OUT), OPTIONAL :: x2, x3, x4, x5, x6, x7, x8, x9, x10, &
      x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, & 
      x25, x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, &
      x39, x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValBool1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SetValBoolN(val, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, &
      x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, & 
      x26, x27, x28, x29, x30, x31, x32, x33, x34, x35, x36, x37, x38, x39, &
      x40, x41, x42, x43, x44, x45, x46, x47, x48, x49, x50)
!! Subroutine to assign val into x1, x2, etc. Up to 50 variables can be modified
!! this way. All except x1 can be missing.
  IMPLICIT NONE
  LOGICAL, INTENT(IN)            :: val
  LOGICAL, INTENT(OUT)           :: x1(:)
  LOGICAL, INTENT(OUT), OPTIONAL :: x2(:), x3(:), x4(:), x5(:), x6(:), &
      x7(:), x8(:), x9(:), x10(:), x11(:), x12(:), x13(:), x14(:), x15(:), & 
      x16(:), x17(:), x18(:), x19(:), x20(:), x21(:), x22(:), x23(:), x24(:), &
      x25(:), x26(:), x27(:), x28(:), x29(:), x30(:), x31(:), x32(:), x33(:), &
      x34(:), x35(:), x36(:), x37(:), x38(:), x39(:), x40(:), x41(:), x42(:), &
      x43(:), x44(:), x45(:), x46(:), x47(:), x48(:), x49(:), x50(:)
                      
#include "epi_inc_assignvalue.F90"
  
  RETURN  

END SUBROUTINE SetValBoolN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AskForNewValue(string, exit_code)
  IMPLICIT NONE
  CHARACTER(*), INTENT(INOUT) :: string
  INTEGER, INTENT(OUT)        :: exit_code
  CHARACTER(100)              :: response
  
  exit_code = -1

  11 CONTINUE 
  CALL Prnt("Do you wish to change the value of string '"//&
                    TRIM(string)//"'? (Y or N) ", ADVANCE='NO', log=.FALSE.)
  READ(*,'(A)') response
  IF(LEN_TRIM(response)==0) GOTO 11
  IF(upcasef(response(1:1))/="Y" .AND. upcasef(response(1:1))/="N") GOTO 11
  IF(upcasef(response(1:1))=="N") RETURN

  CALL Prnt("Enter new value: ", ADVANCE='NO', log=.FALSE.)
  READ(*,'(A)') response
  IF(LEN_TRIM(response)==0) GOTO 11
  string = TRIM(response)
  exit_code = 0 
  
  RETURN
  
END SUBROUTINE AskForNewValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AskUser(response, question, retry, cont, newfile, confirm_delete, &
                   overwrite, overwriteall, terminate, allow_no_response, skip1)
  IMPLICIT NONE
  CHARACTER, INTENT(OUT), OPTIONAL   :: response
  CHARACTER(*), INTENT(IN), OPTIONAL :: question
  LOGICAL, INTENT(IN), OPTIONAL      :: retry, cont, newfile, confirm_delete, &
                                        overwrite, terminate, overwriteall, &
                                        allow_no_response
  INTEGER, INTENT(IN), OPTIONAL      :: skip1
  CHARACTER(mmtl)                    :: question1
  CHARACTER                          :: response1, answers(5), def_answer
  INTEGER                            :: qc, i
  LOGICAL                            :: allow_no_response1, valid_response
  
  allow_no_response1 = .FALSE.
  IF(PRESENT(allow_no_response)) allow_no_response1 = allow_no_response
  
  !! LIST OF QUESTIONS

  def_answer = ""
  
  IF(PRESENT(terminate)) THEN
    qc = 2
    question1 = "> Terminate program? (Y or N)"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "N"
  ENDIF

  IF(PRESENT(overwriteall)) THEN
    qc = 3
    question1 = "> Overwrite? (Y(es), N(o), A(ll))"
    answers(1:qc) = (/"Y","N","A"/)
    def_answer = "Y"
  ENDIF

  IF(PRESENT(overwrite)) THEN
    qc = 2
    question1 = "> Overwrite? (Y(es), N(o))"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "Y"
  ENDIF

  IF(PRESENT(confirm_delete)) THEN
    qc = 2
    question1 = "> Are you sure you want to delete the file? (Y(es), N(o))"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "N"
  ENDIF

  IF(PRESENT(retry)) THEN
    qc = 2
    question1 = "> Retry? (Y or N)"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "Y"
  ENDIF
  
  IF(PRESENT(cont)) THEN
    qc = 2
    question1 = "> Are you sure you want to continue? (Y or N)"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "Y"
  ENDIF

  IF(PRESENT(newfile)) THEN
    qc = 2
    question1 = "> Would you like to enter a different filename? (Y or N)"
    answers(1:qc) = (/"Y","N"/)
    def_answer = "N"
  ENDIF
  
  IF(PRESENT(question)) question1 = question

  10 CONTINUE

  !! Leave skip1 lines empty
  IF(PRESENT(skip1)) THEN
    DO i=1,skip1; WRITE(usto,'(A)'); ENDDO
  ENDIF
  
  !! Ask the question and read the answer
  WRITE (usto,'(A)', ADVANCE='NO') TRIM(question1)//" " 
  
  IF(LEN_TRIM(def_answer)>0) &
    WRITE (usto,'(A)', ADVANCE='NO') "(Default answer is '"//TRIM(def_answer)//"') "
   
  READ(*,'(A)') response1
  
  !! If empty response then use the default answer
  IF(LEN_TRIM(response1)==0) response1 = def_answer

  !! Make the response upper case
  CALL upcase(response1)

  !! Evaluate validity of the response
  valid_response = ANY(response1==answers(1:qc))
  
  !! If no answer is not allowed and the answer does not match 
  IF(.NOT.allow_no_response1 .AND. .NOT.valid_response) THEN
    CALL Prnt("Invalid response.", log=.FALSE.)
    GOTO 10
  ENDIF

  IF(PRESENT(response)) response = response1
  RETURN
  
END SUBROUTINE AskUser

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Terminate(announce, premature, ask, exit_code, log, dellog)
!! Terminates the run of the program possibly with some announcements
!! Variables: ask ... Ask user for termination confirmation
!!            announce ... Announce termination
!!            warn ... Warn user about possible premature termination
!!            log ... Whether warnings should be writen into log file
  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL  :: ask, announce, premature, log, dellog
  INTEGER, INTENT(OUT), OPTIONAL :: exit_code
  CHARACTER                      :: response
  LOGICAL                        :: announce1, premature1, ask1, log1
  INTEGER                        :: ios
  CHARACTER(mltl)                :: text

  !! Exit code if not terminated
  IF(PRESENT(exit_code)) exit_code = 1

  !! Set default values for announcements
  ask1       = .TRUE. 
  announce1  = .TRUE. 
  premature1 = .TRUE.
  log1       = .TRUE.
  
  !! Set values to variables
  IF(PRESENT(ask))       ask1 = ask
  IF(PRESENT(announce))  announce1 = announce
  IF(PRESENT(premature)) premature1 = premature
  IF(PRESENT(log))       log1 = log

  !! Ask whether user wants to quit
  IF(ask1) THEN
    CALL AskUser(response, terminate=.TRUE.) 
    IF(response=="N") GOTO 20
  ENDIF

  !! Announce quiting
  IF(announce1) THEN
    IF(premature1) &
      CALL PrntW("Premature termination probably occurred! Inspect"//& 
                        " the log and output files for more details.", &
                        skip1=1, log=log1)
    CALL Prnt("Program terminated.", skip2=1, log=log1)
  ENDIF

  !! Delete log file
  IF(PRESENT(dellog)) THEN
    IF(dellog) THEN
      CALL DeleteFile(log_file, ios=ios)
      IF(ios==0) THEN
        text = "Log file ["//TRIM(log_file)//"] deleted."
      ELSE
        text = "Log file ["//TRIM(log_file)//"] COULD NOT be deleted."
      ENDIF
      CALL Prnt(text, skip2=1, log=.FALSE.)
    ENDIF
  ENDIF

  !!! TERMINATE !!!
  STOP
  !!! TERMINATE !!!
  
  20 CONTINUE
  
  !! If reached here, then no termination will happen
  IF(PRESENT(exit_code)) exit_code = 0
  
  RETURN

END SUBROUTINE Terminate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

END MODULE EPI_UTILS
