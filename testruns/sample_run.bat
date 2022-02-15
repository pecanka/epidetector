copy ..\source\EpiDetector.exe .

@rem EpiDetector

@rem EpiDetector --bfile PDcontrols_APOE-GAB2 --report-all --out test --overwrite

@rem EpiDetector --bfile PDcontrols_APOE-GAB2 --reportall

@rem EpiDetector --simulateinput --test 0 --pretest 3 --delta 1 ^
@rem   --OR 1 --amodel A --smodel A --maf 0.33 --fixed-maf --nloci 500 --ncases 500 --ncontrols 1000 ^
@rem   --nsamples 1 --cellcorrection 1 --reportdisjoint --reportscore --reportall --rejectall ^
@rem   --mincellcount 0 --out-maf --HWE --group-var --group-cov --center-group 1 --prevalence 0.01 ^
@rem   --rseed 43232 --nolog --nocountdown --nooutput

@rem EpiDetector --simulateinput --test 1 --pretest 3 --amodel O --smodel A --fixedmaf --maf1 auto --maf2 0.25 \
@rem  --delta-min 0.3 --delta-max 0.3 --delta-step 0.4 --OR-min 2 --OR-max 2 --OR-step 1 --nrepeats max \
@rem  --ncases 1500 --ncontrols 2500 --nloci 50 --nsamples 1 --prevalence 0.01 --OR1 1 --OR2 1 \
@rem  --mincellcount 0 --cellcorrection 1 --reportdisjoint --reportscore --reportall --alwaysscore \
@rem  --rejectall --ndigits-pval 3 --HWE --nthreads 1 --out-minimalistic --seed 47897

EpiDetector --simulateinput --test 1 --pretest 3 --delta 0.5 --nloci 10 --ncases 500 --ncontrols 1000 --nsamples 10 --OR 1.5 --reportdisjoint --reportscore --reportall --rejectall --out-maf --out test_simul --out-loc
