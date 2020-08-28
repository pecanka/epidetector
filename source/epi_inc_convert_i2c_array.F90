!! This is intentended to be included in epi_utilities.F90 into every
!! subroutine that converts integers to characters

  !! Input variables
  INTEGER, INTENT(OUT), OPTIONAL :: iostat
  LOGICAL, INTENT(IN), OPTIONAL  :: c, s
  INTEGER, INTENT(IN), OPTIONAL  :: length, ity
  CHARACTER(mttl)                :: char(SIZE(in))
  INTEGER                        :: k

  DO k=1,SIZE(in)
    char(k) = i2c(in(k), iostat, c, s, length, ity)
  ENDDO  

    
