@rem g95 -cpp -O2 -ftrace=full -ftrace=frame epi_logistic3.F90 epi_utilities.F90 epi_math.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testingproc.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector

@rem gfortran -cpp -O2 -fopenmp -lgomp -fbacktrace -fbounds-check epi_logistic3.F90 epi_utilities.F90 epi_par_zig_mod.F90 epi_math.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testingproc.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector

@REM gfortran -cpp -O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ftrapv -ffpe-trap=invalid,zero,overflow ^
@REM   epi_cdflib_biomath_constants_mod.f90 ^
@REM   epi_cdflib_biomath_strings_mod.f90 epi_cdflib_biomath_sort_mod.f90 ^
@REM   epi_cdflib_biomath_interface_mod.f90 epi_cdflib_biomath_mathlib_mod.f90 ^
@REM   epi_cdflib_zero_finder.f90 epi_cdflib_cdf_aux_mod.f90 ^
@REM   epi_cdflib_cdf_normal_mod.f90 epi_cdflib_cdf_gamma_mod.f90 ^
@REM   epi_cdflib_cdf_chisq_mod.f90 epi_cdflib_cdf_nc_chisq_mod.f90 ^
@REM   epi_fmin.f epi_logistic3_mod.F90 epi_parameters_mod.F90 epi_utilities_mod.F90 ^
@REM   epi_par_zig_mod.F90 epi_math_mod.F90 epi_simulation_mod.F90 epi_initialize_mod.F90 ^
@REM   epi_analysis_mod.F90 epi_procedures_mod.F90 epi_dataoutputinput_mod.F90 ^
@REM   epi_testingproc_mod.F90 epi_testing_mod.F90 EpiDetector.F90 -o EpiDetector

rem gfortran -cpp -O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ftrapv -ffpe-trap=invalid,zero,overflow ^
rem gfortran -cpp -O0 -fopenmp -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ftrapv -ffpe-trap=invalid,zero,overflow ^
rem gfortran -cpp -O2 -fopenmp ^
rem gfortran -cpp -O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ftrapv -ffpe-trap=invalid,zero,overflow ^
gfortran -cpp -O2 -fopenmp ^
  epi_fmin.f epi_logistic3_mod.F90 epi_parameters_mod.F90 epi_utilities_mod.F90 ^
  epi_par_zig_mod.F90 epi_math_mod.F90 epi_simulation_mod.F90 epi_testingproc_mod.F90 ^
  epi_initialize_mod.F90 epi_analysis_mod.F90 epi_procedures_mod.F90 epi_dataoutputinput_mod.F90 ^
  epi_testing_mod.F90 EpiDetector.F90 -o EpiDetector

del *.mod

@rem gfortran -cpp -O2 -fopenmp -lgomp epi_fmin.f epi_logistic3.F90 epi_utilities.F90 epi_par_zig_mod.F90 epi_math.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testingproc.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -static -O2 -fopenmp -lgomp -std=f95 -Wextra -Wall -pedantic utilities.F90 epi_fmin.f epi_logistic3.F90 epi_initialize.F90 epi_procedures.F90 epi_dataoutputinput.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -static -O2 -fopenmp -lgomp -pedantic -fbacktrace -fbounds-check epi_fmin.f epi_logistic3.F90 epi_utilities.F90 epi_initialize.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -O2 -fopenmp -lgomp -static -fbounds-check epi_fmin.f epi_logistic3.F90 epi_utilities.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -O2 epi_fmin.f epi_logistic3.F90 epi_utilities.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -O2 epi_fmin.f epi_logistic3.F90 epi_utilities.F90 epi_simulation.F90 epi_initialize.F90 epi_analysis.F90 epi_procedures.F90 epi_dataoutputinput.F90 epi_testing.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -static -O2 -pedantic epi_utilities.F90 epi_fmin.f epi_logistic3.F90 epi_initialize.F90 epi_procedures.F90 epi_testing.F90 epi_dataoutputinput.F90 EpiDetector.F90 -o EpiDetector
@rem gfortran -cpp -O2 utilities.F90 epi_fmin.f epi_logistic3.F90 epi_initialize.F90 epi_procedures.F90 epi_dataoutputinput.F90 EpiDetector.F90 -o EpiDetector
