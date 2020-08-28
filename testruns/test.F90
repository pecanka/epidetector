MODULE PARAMETERS
!! This file is intended for inclusion whenever the below variables and 
!! constants are needed. The idea is to have them defined in one place only
!! for convenience of changing their values.

!! Default accessibility is PUBLIC  
  PUBLIC

  !! PROGRAM NAME AND VERSION
  CHARACTER(11), PARAMETER    :: program_name = "EpiDetector"
  CHARACTER(4), PARAMETER     :: program_version = "0.5"
  CHARACTER(50), PARAMETER    :: contact_email = "j.pecanka@vu.nl"
  
  !! PROGRAM FULL NAME
  INTEGER, PARAMETER          :: fnl = LEN(program_name)+LEN(program_version)+1
  CHARACTER(fnl), PARAMETER   :: program_fullname = program_name//" "//&
                                                    program_version
  
  !! Define a constant which identifies the system (taken from 
  !! http://sourceforge.net/p/predef/wiki/OperatingSystems/ or
  !! http://nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system
  !! List of predefined macros can be obtained by calling: gcc -dM -E - < /dev/null
  CHARACTER(20), PARAMETER    :: &
!#ifdef __gnu_linux__
                                 system_type = "GNU Linux", &
                                 system_bit = "32/64bit", &
                                 system_envir = ""
!#endif

#if defined(__MINGW64__)
                                 system_type = "Windows", &
                                 system_bit = "64bit", &
                                 system_envir = "with MinGW"
#elif defined(__MINGW32__)
                                 system_type = "Windows", &
                                 system_bit = "32bit", &
                                 system_envir = "with MinGW"
#elif defined(__CYGWIN64__)
                                 system_type = "Windows", &
                                 system_bit = "64bit", &
                                 system_envir = "with Cygwin"
#elif defined(__CYGWIN32__)
                                 system_type = "Windows", &
                                 system_bit = "32bit", &
                                 system_envir = "with Cygwin"
#elif defined(__CYGWIN__)
                                 system_type = "Windows", &
                                 system_bit = "32bit", &
                                 system_envir = "with Cygwin"
#elif defined(_WIN64)
                                 system_type = "Windows", &
                                 system_bit = "64bit", &
                                 system_envir = ""
#elif defined(_WIN32)
                                 system_type = "Windows", &
                                 system_bit = "32bit", &
                                 system_envir = ""
#else                                                   
                                 system_type = "unknown system", &
                                 system_bit = "", &
                                 system_envir = ""
#endif


END MODULE PARAMETERS