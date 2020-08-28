!! This is intentended to be included in epi_utilities.F90 into every
!! subroutine that resizes vectors

  INTEGER, INTENT(IN) :: newlength
  INTEGER(ikb)        :: memory_required
  INTEGER             :: length, oldlength, unit

  unit = 1000
#ifdef _Resize1Real
  unit = dpp
#elif defined _Resize1Int_ikn
  unit = ikn
#elif defined _Resize1Int_iks
  unit = iks
#elif defined _Resize1Int_ikb
  unit = ikb
#elif defined _Resize1Logi
  unit = 1
#elif defined _Resize1Char
  unit = LEN(vector)
#elif defined _Resize1TEMP
  unit = 1000
#elif defined _Resize1TEMPmin
  unit = 200
#endif

  !! Check for non-positive length
  IF(newlength<=0) &
    CALL Prnt("ERROR: Cannot resize vector to nonpositive length"//&
                   " ("//routine_name//").", Q=.TRUE.)
                   
  !! Check for too large vector
  memory_required = INT(newlength,ikb)*unit
  IF(memory_required > mem_warn_limit) &
    CALL PrntW("The size of vector to be allocated by "//routine_name//&
               " is very large ("//TRIM(i2c(memory_required,c=.TRUE.))//&
               " bytes). This may cause memory allocation problems on"//&
               " some machines.", skip1=1)

  !! Resize vector to newlength
  IF(ALLOCATED(vector)) THEN
    oldlength = SIZE(vector)
    !! If the old and new sizes are the same, don't do anything
    IF(newlength==oldlength) RETURN
    !! Otherwise change the size
    length = MIN(oldlength,newlength)
    ALLOCATE(temp(oldlength))
    temp = vector
    DEALLOCATE(vector)
    ALLOCATE(vector(newlength))
    IF(PRESENT(newvalue)) vector = newvalue
    vector(1:length) = temp(1:length)
    DEALLOCATE(temp)
  ELSE
    oldlength = 0
    ALLOCATE(vector(newlength))
    IF(PRESENT(newvalue)) vector = newvalue
  ENDIF
  
  !! Assign allvalue to all elements of the vector (this should be modified to
  !! be more efficient, i.e. assign allvalue sooner without using temp)
  IF(PRESENT(allvalue)) vector = allvalue

