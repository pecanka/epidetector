cd ..
cd source
mkoctfile -I. --link-stand-alone epi_fmin.f epi_logistic3_mod.F90 epi_eispack.F90 epi_params_mod.F90 epi_utils_mod.F90 epi_par_zig_mod.F90 epi_math_mod.F90 epi_simul_mod.F90 epi_testprocs_mod.F90 epi_init_mod.F90 epi_analysis_mod.F90 epi_procs_mod.F90 epi_dataio_mod.F90 epi_test_mod.F90 EpiDetector.F90 -o EpiDetector.exe -lgfortran
delete *.mod
cd ..
cd compile
