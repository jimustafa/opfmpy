#!/bin/bash

PW2WANNIER90='pw2wannier90.x'
WANNIER90='wannier90.x'
OPFM='opfm.py'

$OPFM create-cfg --force

$OPFM setup

$WANNIER90 -pp opfm
$PW2WANNIER90 < pw2opfm.in > pw2opfm.out

$OPFM optimize

$OPFM opfm2w90

$WANNIER90
