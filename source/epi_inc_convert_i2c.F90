!! This is intentended to be included in epi_utilities.F90 into every
!! subroutine that converts integers to characters
!! INT(LOG10(DBLE(MAX(1,ABS(in)))))+1 is the exact length of string can contain
!! 'in' (without +- sign), INT(INT(LOG10(DBLE(MAX(1,ABS(in)))))/3) is the number
!! of commas that will be added 

  !! Input variables
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  LOGICAL, INTENT(IN), OPTIONAL  :: c, s
  INTEGER, INTENT(IN), OPTIONAL  :: length, ity
  CHARACTER(INT(LOG10(DBLE(MAX(one,ABS(DBLE(in)))))) + 1 - INT(MAX(-one,MIN(DBLE(in),zero))) + &
    IND(add_commas,add_spaces)*INT(INT(LOG10(MAX(one,ABS(DBLE(in)))))/three)) :: char
  CHARACTER(mstl)                :: temp
  INTEGER                        :: iostatus, k, length1

  !! Initialize
  length1 = 0

  !! Determine length of output based on given integer type
  IF(PRESENT(ity)) THEN
    IF(ity==iks) length1 = iks_len  
    IF(ity==ikn) length1 = ikn_len  
    IF(ity==ikb) length1 = ikb_len  
  ENDIF
  
  !! If length present, ity is irrelevant
  IF(PRESENT(length)) length1 = length
   
  !! Convert integer into character
  WRITE(char, '(I0)', IOSTAT=iostatus) in
  
  !! Check for errors in coversion
  IF(iostatus /= 0) THEN
    IF(PRESENT(iostat)) THEN
      iostat = iostatus
      char = ""
      RETURN
    ENDIF
    PRINT *,"" 
    PRINT *," ERROR: Problem with argument", in, "!" 
    CALL PrntE("Argument must be an integer.", Q=.TRUE.)
  ENDIF

  !! Add spaces and/or commas to the string at every 3rd character from back
  CALL AddSpacesCommas(char, c, s)

  !! Set length of output by adding trailing spaces at the beggining of output
  IF(length1>0) THEN
    k = MIN(LEN(temp), MIN(length1,LEN(char)) - LEN_TRIM(char))
    IF(k>0) THEN
      temp = ""
      char = temp(1:k)//char
    ENDIF
  ENDIF
  

    
