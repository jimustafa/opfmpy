#!/bin/bash

PW='pw.x'

ln -sf ../psp/si_lda_v1.uspp.F.UPF Si.UPF

$PW < scf.in > scf.out
$PW < bands.in > bands.out

cd opfm
./run.sh
cd ../

cd opfm_projwfc
./run.sh
cd ../
