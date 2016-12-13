#!/bin/bash

PW2WANNIER90='pw2wannier90.x'
PROJWFC='projwfc.x'
WANNIER90='wannier90.x'
OPFM='opfm.py'

$OPFM create-cfg --force

$OPFM setup

$WANNIER90 -pp opfm
$PW2WANNIER90 < pw2opfm.in > pw2opfm.out

$PROJWFC < projwfc.in > projwfc.out
$OPFM projwfc2amn

$OPFM optimize

$OPFM opfm2w90

$WANNIER90
