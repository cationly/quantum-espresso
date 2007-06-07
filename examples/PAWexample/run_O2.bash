#!/bin/bash

PWROOT=$HOME/QE/PAW/espresso/

PW=$PWROOT/bin/pw.x
SCRATCH=/scratch/degironc/o2/
echo $PW, $SCRATCH


fileout=O2_dist.dat

DISTLIST="0.8 1.0 1.2 1.4"
DISTLIST="0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 "
PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW3.PAW PAW4.PAW" 
#ECUTWFCLIST="100 70 50 35 25 18 13"


echo -n "#dist" > $fileout
for pseudo in $PSEUDOLIST; do
echo -n " $pseudo " >> $fileout
done
echo >> $fileout

ecutwfc=70
ecutrho=400

for dist in $DISTLIST; do
name=O2_$dist
echo -n " $dist " >> $fileout

for pseudo in $PSEUDOLIST; do
echo Now: $ecutwfc $pseudo $dist

pwout=${name}_${ecutwfc}_${pseudo}.out
cat <<EOF | $PW > $pwout
${name}
Oxygen-hydrogen
 &control
    calculation = 'scf'
    prefix='o2',
    pseudo_dir = './PSEUDO/'
    outdir='$SCRATCH/'
    disk_io='low'
 /
 &system
    ibrav=  1, celldm(1) =10.0, nat=  2, ntyp= 1,
    ecutwfc= $ecutwfc
    ecutrho= $ecutrho
    occupations = 'smearing'
    smearing = 'gauss'
    degauss = 0.01
 /
 &electrons
    electron_maxstep=20,
    mixing_mode = 'plain'
    mixing_beta = 0.7
    mixing_ndim = 1
    conv_thr =  1.0d-6
 /
ATOMIC_SPECIES
 O   1.000  $pseudo
 H   1.000  H.pbe-rrkjus.UPF
ATOMIC_POSITIONS {bohr}
O  -$dist   -$dist   -$dist
O  +$dist   +$dist   +$dist
K_POINTS automatic
 1 1 1 0 0 0
EOF

energy=$(awk '/\!/ {print $5}' $pwout)
echo -n " $energy " >> $fileout

done
echo >> $fileout
done

cat <<EOF | gnuplot
set term post
set out 'E_${fileout}.ps'
set title "${fileout}"
set xlabel "ecutwfc (Ryd)  (no dual)"
set ylabel "Total energy (Ryd)"
set da s lp
plot 'E_${fileout}.dat' u 1:2 t "NC1", 'E_${fileout}.dat' u 1:3 t "NC2", 'E_${fileout}.dat' u 1:4 t "US1", 'E_${fileout}.dat' u 1:5 t "US2", 'E_${fileout}.dat' u 1:6 t "PAW1", 'E_${fileout}.dat' u 1:7 t "PAW3", 'E_${fileout}.dat' u 1:8 t "PAW4"
EOF
