
### ANALYSIS ###
if false; then
echo "Running ..."
./EpiDetector.exe --pretest 3 --maintest 0 --delta 0.8 --level1 0.0001 --outreport1 0.0001 --nthreads 1 \
  --bfile PelvicOrganProlapse_CaCo_zCall-mind-0.05-geno-0.05-maf-0.005-hwe-0.0001-sexed-dupl-het-missphen-0.001-misshap-1e-6-mds 
fi

### DEBUGGING ###
# --maxntests 42554 --seed 43207 --keep-temp  
if false; then
./EpiDetector.exe --maxntest 1234567 --pretest-same-chr --keep-temp --pretest 4 --delta 0.8 \
--center-group 2 --level1 0.005 --simulate-input --nloci 10 --ncontrols 500 --ncases 1000 \
--report-all --reject-all --prev 1 --out-debug --group-var --group-cov --out-minimalistic
fi
#./EpiDetector.exe --seed 43207 --pretest 3 --maintest 0 --delta 0.9 --level1 0.9 --outreport1 1 --tmp-to-out-only \
#  --overwrite --tmp PelvicOrganProlapse_CaCo_zCall-mind-0.05-geno-0.05-maf-0.005-hwe-0.0001-sexed-dupl-het-missphen-0.001-misshap-1e-6-mds_aug082014161323834_pid=9400_amodel=A_delta=0.9_cntrgrp=0.tmp
#./plink-1.9 --epistasis --bfile PelvicOrganProlapse_CaCo_zCall-mind-0.05-geno-0.05-maf-0.005-hwe-0.0001-sexed-dupl-het-missphen-0.001-misshap-1e-6-mds
#--pretest 3 --maintest 0 --delta 0.9 --level1 0.0001 --outreport1 0.0001 --nthreads 2 --bfile PelvicOrganProlapse_CaCo_zCall-mind-0.05-geno-0.05-maf-0.005-hwe-0.0001-sexed-dupl-het-missphen-0.001-misshap-1e-6-mds 
#gprof EpiDetector > gprof.out

./EpiDetector.exe --simulateinput --test 1 --pretest 4 --amodel A \
  --smodel A --fixedmaf --maf1 0.35 --maf2 0.35 --maf3 0.35 --delta-min 0.1 \
  --delta-max 0.1 --delta-step 0.2 --delta-simple-frac 1 --OR-min 1 --OR-max 1 \
  --OR-step 1 --nrepeats max --nloci 10 --ncases 1000 --ncontrols 500 \
  --nsamples 1 --prevalence 0.05 --OR1 1 --OR2 1 --LD 0 --mincellcount 0 \
  --cellcorrection 1 --reportdisjoint --reportscore --reportall --alwaysscore \
  --rejectall --out-minimalistic --ndigit-pval 1 --hide-errors --nocountdown \
  --center-group 2 --group-var --group-cov --max-runtime-s 7200 --HWE
  
Rscript get_means.R

