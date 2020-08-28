!! This is intentended to be included in epi_utilities.F90 into every
!! subroutine that resizes arrays

  INTEGER, INTENT(IN) :: newdim1, newdim2
  INTEGER(ikb)        :: memory_required
  INTEGER             :: unit, dim1, dim2, olddim1, olddim2

  unit = 1000
#ifdef _Resize2Real
  unit = dpp
#elif defined _Resize2Int_ikn
  unit = ikn
#elif defined _Resize2Int_iks
  unit = iks
#elif defined _Resize2Int_ikb
  unit = ikb
#elif defined _Resize2Logi
  unit = 1
#elif defined _Resize2Char
  unit = LEN(array)
#elif defined _Resize2TEMP
  unit = 100
#endif

  !! Check for non-positive dimensions
  IF(newdim1 <= 0) &
    CALL PrntE("Dimension 1 must be positive ("//routine_name//").", Q=.TRUE.)
  IF(newdim2 <= 0) &
    CALL PrntE("Dimension 2 must be positive ("//routine_name//").", Q=.TRUE.)

  !! Check for too large array
  memory_required = INT(newdim1,ikb)*INT(newdim2,ikb)*unit
  IF(memory_required > mem_warn_limit) &
    CALL PrntW("The size of array to be allocated by "//routine_name//&
               " is very large ("//TRIM(i2c(memory_required,c=.TRUE.))//&
               " bytes). This may cause memory allocation problems on"//&
               " some machines.", skip1=1)

  !! Resize array to newlength and assign old values
  IF(ALLOCATED(array)) THEN
    !! Remember the old dimensions
    olddim1 = SIZE(array,1) 
    olddim2 = SIZE(array,2) 

    !! If the old and new sizes are the same, don't do anything
    IF(newdim1==olddim1 .AND. newdim2==olddim2) RETURN

    !! Otherwise change the size
    dim1 = MIN(olddim1, newdim1)
    dim2 = MIN(olddim2, newdim2)

    ALLOCATE(temp(olddim1, olddim2))
    temp = array
    DEALLOCATE(array)
    ALLOCATE(array(newdim1,newdim2))
    IF(PRESENT(newvalue)) array = newvalue 
    array(1:dim1,1:dim2) = temp(1:dim1,1:dim2)
    DEALLOCATE(temp)
  !! Resize array to newlength only, because there are no old values
  ELSE
    olddim1 = 0
    olddim2 = 0
    ALLOCATE(array(newdim1,newdim2))
    IF(PRESENT(newvalue)) array = newvalue 
  ENDIF

  !! Assign allvalue to all elements of the array (this should be modified to
  !! be more efficient, i.e. assign allvalue sooner without using temp)
  IF(PRESENT(allvalue)) array = allvalue

