!! This is intentended to be included in epi_utilities.F90 into every
!! subroutine that converts strings to integers (of selected type)

  CHARACTER(*), INTENT(INOUT)    :: string
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  INTEGER                        :: iostat1, exit_code

  10 CONTINUE

  !! If string empty, return
  IF(LEN_TRIM(string)==0) THEN
    number = 0
    IF(PRESENT(iostat)) iostat = -9
    RETURN
  ENDIF
  
  !! Do the conversion of string to integer
  READ(string, *, IOSTAT=iostat1) number
  
  !! Check for error status
  IF(iostat1 /= 0) THEN
    number = 0
    IF(PRESENT(iostat)) THEN
      iostat = iostat1
      RETURN
    ENDIF
    CALL PrntE("Cannot convert '"//TRIM(string)//"' to integer.", log=silent, Q=silent)
    CALL AskForNewValue(string, exit_code)
    IF(exit_code/=0) STOP
    IF(exit_code==0) GOTO 10
  ENDIF
  
  IF(PRESENT(iostat)) iostat = 0
