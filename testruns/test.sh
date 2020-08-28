#./EpiDetector --simulateinput --smodel A --fixedmaf --maf1 0.35 --maf2 0.35 --maf3 0.35 \
  --OR 1.2 --nloci 100 --ncases 772  --ncontrols 2024 --nsamples 1 \
  --prevalence -3.3 --OR1 1.0 --OR2 1.0 --LD 0 \
  --test 1 --pretest 4 --amodel A --level1 0.5 --delta 0.7 --nthreads 1 \
  --tempcycle 400000 --non-group-var --group-cov --center-by-self \
  --report-adjusted --report-disjoint --report-score --always-score--keep-all-chr \
  --keep-all-bp --seed 444444 --out-minimal --out-maf --out-loc --out-errcode-single
  
./EpiDetector --bfile data_PT=4_sm=A_am=A_LD=0_OR=3-3_OR1=1_OR2=1_prev=0.05_d=0.05-0.05_a1=0.00000001-0.95_M1=10_nca=1000_nco=3000_nL1=0_nL2=0_pid26116_72378_OR=3 \
  --test 1 --pretest 4 --amodel A --level1 0.5 --delta 0.7 --nthreads 1 \
  --tempcycle 400000 --non-group-var --group-cov --center-by-self \
  --report-adjusted --report-disjoint --report-score --always-score--keep-all-chr \
  --keep-all-bp --seed 444444 --out-minimal --out-maf --out-loc --out-errcode-single
