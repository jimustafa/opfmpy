#!/bin/bash

PW='pw.x'

ln -sf ../psp/si_lda_v1.uspp.F.UPF Si.UPF
ln -sf ../psp/o_lda_v1.2.uspp.F.UPF O.UPF

$PW < scf.in > scf.out
$PW < bands.in > bands.out

cd opfm_1-16
./run.sh
cd ../

cd opfm_9-16
./run.sh
cd ../
